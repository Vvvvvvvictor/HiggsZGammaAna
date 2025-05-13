#!/usr/bin/env python
import numpy as np
import pandas as pd
import uproot
import pickle as pkl
import json
import os
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from xgboost import XGBClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import QuantileTransformer
from tqdm import tqdm

from pdb import set_trace

# General settings
USE_STORED_MODEL = True
TRANSFORM = True  # Use transformed scores with suffix '_t'
REWEIGHT = True   # Apply reweighting to background for better modeling
MINBKG = 2       # Minimum background events required in a bin
MASS_WINDOW = (120, 130)  # Signal mass window for significance calculation
SIDEBAND = ((100, 120), (130, 180))  # Sideband regions for background estimation

"""
Function to calculate significance considering weights.
Inputs:
  - signal: sum of weighted signal events
  - background: sum of weighted background events
  - sig_err: error on signal estimation (optional)
  - bkg_err: error on background estimation (optional)
Returns:
  - significance value and its error (if errors are provided)
"""
def calculate_significance(signal, background, sig_err=None, bkg_err=None):
    # Avoid division by zero
    if background <= 0 or signal <= 0:
        return 0, 0
    
    ntot = signal + background
    # Calculate significance using likelihood ratio formula
    significance = 2 * ((ntot * np.log(ntot / background)) - signal)
    
    # If errors are not provided, return only significance
    if sig_err is None or bkg_err is None:
        return significance, 0
    
    # Calculate error propagation on significance
    error = (np.log(ntot/background)*sig_err)**2 + ((np.log(ntot/background) - (signal/background))*bkg_err)**2
    return significance, error

"""
Function to calculate bin weights within given boundaries.
Inputs:
  - data: histogram data array
  - bl: lower bin boundary
  - br: upper bin boundary
  - data_err: error on data (optional)
Returns:
  - sum of weights and its error (if error is provided)
"""
def calculate_weights(data, bl, br, data_err=None):
    weight = data[bl:br].sum()
    if data_err is None:
        return weight, 0
    weight_err = np.sqrt(data_err[bl:br].sum())
    return weight, weight_err

"""
Function to optimize thresholds for maximum significance.
This is the core algorithm for finding optimal bin boundaries.
Inputs:
  - sig: signal histogram data
  - sig_err: error on signal data (optional)
  - bkg: background histogram data
  - bkg_err: error on background data (optional)
  - bl: lower bin boundary
  - br: upper bin boundary
  - nbin: number of bins to divide into
  - minN: minimum number of background events required in a bin
  - early_stop: stop optimization if no improvement after this many iterations
  - floatB: whether to allow floating bin boundaries
  - pbar: whether to show progress bar
  - num: numerator histogram for reweighting (optional)
  - den: denominator histogram for reweighting (optional)
Returns:
  - list of optimized bin boundaries, total significance, and error on significance
"""
def fit(sig, sig_err, bkg, bkg_err, bl, br, nbin, minN=MINBKG, early_stop=-1, floatB=False, pbar=False, num=1, den=1):
    # Base case: single bin
    if nbin == 1:
        if floatB: 
            return [], 0, 0

        # Calculate signal and background weights within boundaries
        signal, signal_err = calculate_weights(sig, bl, br, sig_err)
        background, background_err = calculate_weights(bkg, bl, br, bkg_err)
        
        # Apply data-driven reweighting if requested
        if REWEIGHT:
            rw_num, _ = calculate_weights(num, bl, br)
            rw_den, _ = calculate_weights(den, bl, br)
            # Skip if insufficient statistics for reweighting
            if rw_num < 10 or rw_den < 10: 
                return -1, -1, -1
            background = background * rw_num / rw_den
            background_err = background_err * (rw_num / rw_den)
        
        # Skip if not enough background events
        if background < minN:
            return -1, -1, -1
            
        # Calculate significance for this bin
        significance, error = calculate_significance(signal, background, signal_err, background_err)
        return [bl], significance, error

    # Recursive case: multiple bins
    elif nbin > 1:
        # Binary splitting for efficient bin optimization
        L = int(np.ceil(np.log2(nbin)))
        N2 = 2**(L-1)  # Number of bins in second half
        N1 = nbin - N2  # Number of bins in first half
        
        b_opt, s_opt, e_opt, stop = -1, -1, -1, 0
        
        # Try each possible splitting point
        for b in (tqdm(np.arange(bl+1, br)) if pbar else np.arange(bl+1, br)):
            # Optimize first half
            b1, s1, e1 = fit(sig, sig_err, bkg, bkg_err, bl, b, N1, minN, floatB=floatB, num=num, den=den)
            if b1 == -1:  # Skip if invalid
                continue
                
            # Optimize second half
            b2, s2, e2 = fit(sig, sig_err, bkg, bkg_err, b, br, N2, minN, num=num, den=den)
            if b2 == -1:  # Break if invalid (can't improve further)
                break
                
            # Update if better significance found
            if s1 + s2 > s_opt:
                s_opt = s1 + s2
                b_opt = sorted(list(set(b1 + [b] + b2)))
                e_opt = (e1*s1) + (e2*s2)
            else:
                # Early stopping mechanism
                stop += 1
                if early_stop > 0 and stop == early_stop:
                    break

        return b_opt, s_opt, e_opt

"""
Function to optimize the number of bins for maximum significance.
This function tests different bin counts and returns the configuration that gives the best significance.
Inputs:
  - sig_hist: signal histogram data
  - sig_hist_err: error on signal histogram
  - bkg_hist: background histogram data
  - bkg_hist_err: error on background histogram 
  - edges: bin edges for histograms
  - min_bins: minimum number of bins to try (default 2)
  - max_bins: maximum number of bins to try (default 7)
  - **kwargs: additional arguments to pass to fit function
Returns:
  - best_n_bins: optimal number of bins
  - best_boundaries: list of bin boundaries (indices)
  - best_thresholds: list of bin boundaries (actual values)
  - best_significance: significance value for optimal configuration
  - best_significance_error: error on significance for optimal configuration
"""
def optimize_bin_count(
    sig_hist, sig_hist_err, bkg_hist, bkg_hist_err, edges, min_bins=2, max_bins=7, **kwargs
):
    best_n_bins = 1
    best_significance = 0
    best_significance_error = 0
    best_boundaries = []
    best_thresholds = []
    
    # Try different numbers of bins from min_bins to max_bins
    for n_bins in range(min_bins, max_bins + 1):
        print(f"Testing {n_bins} bins...")
        boundaries, significance, significance_error = fit(
            sig_hist, sig_hist_err, bkg_hist, bkg_hist_err,
            0, len(sig_hist), n_bins, pbar=True, **kwargs
        )
        
        if boundaries == -1:
            print(f"Invalid result for {n_bins} bins. Skipping.")
            continue
            
        # Add the upper boundary
        boundaries.append(len(sig_hist))
        thresholds = [edges[idx] for idx in boundaries]
        
        # Convert significance to more familiar form
        significance = np.sqrt(significance)
        significance_error = np.sqrt(significance_error) if significance_error > 0 else 0
        
        # Check if new binning is significantly better
        if significance - best_significance > significance_error:
            best_significance = significance
            best_significance_error = significance_error
            best_n_bins = n_bins
            best_boundaries = boundaries
            best_thresholds = thresholds
            print(f"New best: {n_bins} bins with significance {significance:.4f} ± {significance_error:.4f}")
        else:
            # Stop early if we don't improve for more than 1 bin past the best
            if n_bins > best_n_bins + 1:
                print(f"No significant improvement after {n_bins} bins. Stopping.")
                break
    
    return best_n_bins, best_boundaries, best_thresholds, best_significance, best_significance_error

"""
Main function to load data and perform the optimization.
This function handles the full workflow from loading data to finding optimal categorization.
"""
def main():
    # Directory where the input ROOT files are stored
    filepath = "/eos/user/j/jiehan/root/outputs/test/all_jet"
    
    # Define background, signal, and data samples
    bkg_names = ["ZGToLLG", "DYJetsToLL", "EWKZ2J"]
    sig_names = ["ggH_M125", "VBF_M125"]
    data_names = ["Data"]
    
    # File path for output figures
    if not os.path.exists("figs"):
        os.makedirs("figs")
    
    print("Loading data...")
    
    # Load background data
    bkg_data = pd.DataFrame()
    for bkg in bkg_names:
        try:
            # Load data and apply basic selection
            tree = uproot.open(f"{filepath}/{bkg}.root")["all_jet"]
            df = tree.arrays(["bdt_score", "bdt_score_t", "vbf_score", "vbf_score_t",
                             "bdt_2j_score", "bdt_2j_score_t", "weight", "H_mass", 
                             "event", "n_jets"], library="pd")
            df = df.query("H_mass > 95 & H_mass < 180")
            bkg_data = pd.concat([bkg_data, df], ignore_index=True)
        except Exception as e:
            print(f"Error loading {bkg}: {e}")
    
    bkg_data["label"] = 0  # Mark as background
    
    # Load signal data
    sig_data = pd.DataFrame()
    ggH_data = pd.DataFrame()  # Separate ggH for VBF purity calculation
    VBF_data = pd.DataFrame()  # Separate VBF for VBF purity calculation
    
    for sig in sig_names:
        try:
            sig_file = f"{filepath}/{sig}.root"
            if os.path.exists(sig_file):
                tree = uproot.open(sig_file)["all_jet"]
                df = tree.arrays(["bdt_score", "bdt_score_t", "vbf_score", "vbf_score_t",
                                 "bdt_2j_score", "bdt_2j_score_t", "weight", "H_mass", 
                                 "event", "n_jets"], library="pd")
                df = df.query("H_mass > 95 & H_mass < 180")
                
                # Store separately for VBF purity calculation
                if sig == "ggH_M125":
                    ggH_data = df.copy()
                elif sig == "VBF_M125":
                    VBF_data = df.copy()
                
                sig_data = pd.concat([sig_data, df], ignore_index=True)
        except Exception as e:
            print(f"Error loading {sig}: {e}")
    
    sig_data["label"] = 1  # Mark as signal
    
    # Load data (for data-driven background estimation)
    data_data = pd.DataFrame()
    for data in data_names:
        try:
            data_file = f"{filepath}/{data}.root"
            if os.path.exists(data_file):
                tree = uproot.open(data_file)["all_jet"]
                df = tree.arrays(["bdt_score", "bdt_score_t", "vbf_score", "vbf_score_t",
                                 "bdt_2j_score", "bdt_2j_score_t", "weight", "H_mass", 
                                 "event", "n_jets"], library="pd")
                # Use sideband regions for background estimation
                df = df.query("(H_mass > 95 & H_mass < 120) | (H_mass > 130 & H_mass < 180)")
                data_data = pd.concat([data_data, df], ignore_index=True)
        except Exception as e:
            print(f"Error loading {data}: {e}")
    
    print("Data loaded successfully")
    
    # Extract samples with n_jets >= 2
    bkg_2jets = bkg_data.query("n_jets >= 2")
    sig_2jets = sig_data.query("n_jets >= 2")
    ggH_2jets = ggH_data.query("n_jets >= 2") if not ggH_data.empty else pd.DataFrame()
    VBF_2jets = VBF_data.query("n_jets >= 2") if not VBF_data.empty else pd.DataFrame()
    
    # Define the score columns to use
    score_col = "bdt_score_t" if TRANSFORM else "bdt_score"
    vbf_score_col = "vbf_score_t" if TRANSFORM else "vbf_score"
    score_2j_col = "bdt_2j_score_t" if TRANSFORM else "bdt_2j_score"
    
    print("Finding VBF-enriched region...")
    
    # Calculate VBF purity in bins of vbf_score
    n_bins = 100  # Number of bins for vbf_score
    vbf_edges = np.linspace(0, 1, n_bins+1)
    vbf_centers = (vbf_edges[:-1] + vbf_edges[1:]) / 2
    
    # Arrays to store purity values
    vbf_purity = np.zeros(n_bins)
    vbf_yield = np.zeros(n_bins)
    ggH_yield = np.zeros(n_bins)
    
    # Calculate VBF purity in each bin
    for i in range(n_bins):
        # bin_low, bin_high = vbf_edges[i], vbf_edges[i+1]
        bin_low, bin_high = vbf_edges[i], vbf_edges[-1]
        
        # Sum of weights for VBF in this bin
        if not VBF_2jets.empty:
            vbf_sum = VBF_2jets[(VBF_2jets[vbf_score_col] >= bin_low) & 
                             (VBF_2jets[vbf_score_col] < bin_high)]["weight"].sum()
            vbf_yield[i] = vbf_sum
        else:
            vbf_sum = 0
            
        # Sum of weights for ggH in this bin
        if not ggH_2jets.empty:
            ggH_sum = ggH_2jets[(ggH_2jets[vbf_score_col] >= bin_low) & 
                              (ggH_2jets[vbf_score_col] < bin_high)]["weight"].sum()
            ggH_yield[i] = ggH_sum
        else:
            ggH_sum = 0
            
        # Calculate purity
        total = vbf_sum + ggH_sum
        vbf_purity[i] = vbf_sum / total if total > 0 else 0

    # Find the vbf_score threshold for 70% VBF purity
    vbf_purity_threshold = 0.7
    vbf_purity_index = 0
    vbf_score_threshold = 0
    
    for i in range(n_bins):
        if vbf_purity[i] >= vbf_purity_threshold:
            vbf_score_threshold = vbf_edges[i]
            vbf_purity_index = i
            break
    
    # If no bin exceeds the threshold, use the bin with maximum purity
    if vbf_score_threshold == 0:
        max_purity_idx = np.argmax(vbf_purity)
        vbf_score_threshold = vbf_edges[max_purity_idx]
        print(f"Warning: No bin exceeds {vbf_purity_threshold*100:.1f}% VBF purity. Using threshold with max purity: {vbf_purity[max_purity_idx]:.4f}")
    
    print(f"VBF-enriched region threshold: vbf_score{'_t' if TRANSFORM else ''} >= {vbf_score_threshold:.4f}")
    print(f"Expected VBF purity in this region: {vbf_purity[vbf_purity_index]:.4f} (VBF / (VBF + ggH))")
    
    # Plot VBF purity vs vbf_score
    plt.figure(figsize=(10, 6))
    plt.plot(vbf_centers, vbf_purity, 'b-', label='VBF Purity')
    plt.axhline(y=vbf_purity_threshold, color='r', linestyle='--', label=f'Threshold ({vbf_purity_threshold:.1f})')
    plt.axvline(x=vbf_score_threshold, color='g', linestyle='--', label=f'vbf_score threshold ({vbf_score_threshold:.4f})')
    plt.xlabel('vbf_score' + ('_t' if TRANSFORM else ''))
    plt.ylabel('VBF Purity (VBF / (VBF+ggH))')
    plt.title('VBF Purity vs vbf_score')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('figs/vbf_purity.png')
    
    # Plot event yields
    plt.figure(figsize=(10, 6))
    plt.bar(vbf_centers, vbf_yield, width=vbf_edges[1]-vbf_edges[0], alpha=0.5, label='VBF')
    plt.bar(vbf_centers, ggH_yield, width=vbf_edges[1]-vbf_edges[0], bottom=vbf_yield, alpha=0.5, label='ggH')
    plt.axvline(x=vbf_score_threshold, color='r', linestyle='--', label=f'Threshold ({vbf_score_threshold:.4f})')
    plt.xlabel('vbf_score' + ('_t' if TRANSFORM else ''))
    plt.ylabel('Event Yield')
    plt.title('Signal Event Yields vs vbf_score')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('figs/vbf_yield.png')
    
    # ----- VBF-enriched region analysis -----
    print("\nAnalyzing VBF-enriched region...")
    
    # Select events in VBF-enriched region (n_jets >= 2 and high vbf_score)
    vbf_region_sig = sig_data.query(f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold}")
    vbf_region_bkg = bkg_data.query(f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold}")
    
    # Calculate optimal bin boundary for VBF-enriched region using bdt_2j_score
    print("Finding optimal number of bins and boundaries for VBF-enriched region...")
    
    # Create histograms for bdt_2j_score in the VBF-enriched region
    n_scan = 100  # Number of bins for scanning
    
    # Signal histogram
    sig_hist, sig_edges = np.histogram(
        vbf_region_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_2j_col], 
        bins=n_scan, range=(0, 1), 
        weights=vbf_region_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    sig_hist_err = np.histogram(
        vbf_region_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_2j_col], 
        bins=n_scan, range=(0, 1), 
        weights=vbf_region_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]**2
    )[0]**0.5
    
    # Background histogram
    bkg_hist, bkg_edges = np.histogram(
        vbf_region_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_2j_col], 
        bins=n_scan, range=(0, 1), 
        weights=vbf_region_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    bkg_hist_err = np.histogram(
        vbf_region_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_2j_col], 
        bins=n_scan, range=(0, 1), 
        weights=vbf_region_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]**2
    )[0]**0.5
    
    # For reweighting
    if REWEIGHT:
        sideband_query = " | ".join([f"(H_mass > {sb[0]} & H_mass < {sb[1]})" for sb in SIDEBAND])
        sideband_hist = np.histogram(
            vbf_region_bkg.query(sideband_query)[score_2j_col], 
            bins=n_scan, range=(0, 1), 
            weights=vbf_region_bkg.query(sideband_query)["weight"]
        )[0]
        
        data_hist = np.histogram(
            data_data.query(f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold}")[score_2j_col], 
            bins=n_scan, range=(0, 1), 
            weights=data_data.query(f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold}")["weight"]
        )[0]

    # Find the optimal number of bins for VBF-enriched region
    vbf_best_n_bins, vbf_region_boundaries, vbf_region_thresholds, vbf_region_significance, _ = optimize_bin_count(
        sig_hist, sig_hist_err, bkg_hist, bkg_hist_err, sig_edges,
        min_bins=2, max_bins=7, 
        minN=MINBKG, 
        num=data_hist if REWEIGHT else 1,
        den=sideband_hist if REWEIGHT else 1
    )

    print(f"\nOptimal number of bins for VBF-enriched region: {vbf_best_n_bins}")
    print(f"VBF-enriched region boundaries: {vbf_region_thresholds}")
    print(f"VBF-enriched region significance: {vbf_region_significance:.4f}")
    
    # ----- Remaining events analysis -----
    print("\nAnalyzing remaining events...")
    
    # Select remaining events (n_jets < 2 OR (n_jets >= 2 AND vbf_score < threshold))
    remaining_sig = sig_data.query(f"n_jets < 2 or (n_jets >= 2 and {vbf_score_col} < {vbf_score_threshold})")
    remaining_bkg = bkg_data.query(f"n_jets < 2 or (n_jets >= 2 and {vbf_score_col} < {vbf_score_threshold})")
    remaining_data = data_data.query(f"n_jets < 2 or (n_jets >= 2 and {vbf_score_col} < {vbf_score_threshold})")
    
    # Optimize binning for remaining events
    print("Finding optimal number of bins and boundaries for remaining events...")
    
    # Create histograms for bdt_score in the remaining events region
    remaining_sig_hist, remaining_sig_edges = np.histogram(
        remaining_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_col], 
        bins=n_scan, range=(0, 1), 
        weights=remaining_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    remaining_sig_hist_err = np.histogram(
        remaining_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_col], 
        bins=n_scan, range=(0, 1), 
        weights=remaining_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]**2
    )[0]**0.5
    
    remaining_bkg_hist, remaining_bkg_edges = np.histogram(
        remaining_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_col], 
        bins=n_scan, range=(0, 1), 
        weights=remaining_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    remaining_bkg_hist_err = np.histogram(
        remaining_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_col], 
        bins=n_scan, range=(0, 1), 
        weights=remaining_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]**2
    )[0]**0.5
    
    # For data-driven background estimation
    if REWEIGHT:
        sideband_query = " | ".join([f"(H_mass > {sb[0]} & H_mass < {sb[1]})" for sb in SIDEBAND])
        remaining_sideband_hist = np.histogram(
            remaining_bkg.query(sideband_query)[score_col], 
            bins=n_scan, range=(0, 1), 
            weights=remaining_bkg.query(sideband_query)["weight"]
        )[0]
        
        remaining_data_hist = np.histogram(
            remaining_data[score_col], 
            bins=n_scan, range=(0, 1), 
            weights=remaining_data["weight"]
        )[0]
    
    # Find the optimal number of bins for remaining events
    best_n_bins, best_boundaries, best_thresholds, best_significance, best_significance_error = optimize_bin_count(
        remaining_sig_hist, remaining_sig_hist_err, 
        remaining_bkg_hist, remaining_bkg_hist_err, 
        remaining_sig_edges,
        min_bins=2, max_bins=7,
        minN=MINBKG,
        num=remaining_data_hist if REWEIGHT else 1,
        den=remaining_sideband_hist if REWEIGHT else 1
    )
    
    print(f"\nOptimal number of bins for remaining events: {best_n_bins}")
    print(f"Optimal thresholds for remaining events: {best_thresholds}")
    print(f"Remaining events significance: {best_significance:.4f}")
    
    # ----- Calculate significances for each category -----
    print("\nCalculating significance for each category...")
    
    # First, calculate significances for VBF-enriched region categories
    vbf_region_cats = []
    
    for i in range(len(vbf_region_thresholds) - 1):
        cat_name = f"VBF_Cat{i}"
        low, high = vbf_region_thresholds[i], vbf_region_thresholds[i+1]
        
        # Select events in this category
        cat_sig = vbf_region_sig.query(f"{score_2j_col} >= {low} and {score_2j_col} < {high} and H_mass > {MASS_WINDOW[0]} and H_mass < {MASS_WINDOW[1]}")
        cat_bkg = vbf_region_bkg.query(f"{score_2j_col} >= {low} and {score_2j_col} < {high} and H_mass > {MASS_WINDOW[0]} and H_mass < {MASS_WINDOW[1]}")
        
        # Calculate VBF events in this category
        cat_vbf = VBF_data.query(f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold} and {score_2j_col} >= {low} and {score_2j_col} < {high} and H_mass > {MASS_WINDOW[0]} and H_mass < {MASS_WINDOW[1]}")
        vbf_yield = cat_vbf["weight"].sum() if not cat_vbf.empty else 0
        
        # Calculate yields and significance
        sig_yield = cat_sig["weight"].sum()
        sig_err = np.sqrt((cat_sig["weight"]**2).sum())
        bkg_yield = cat_bkg["weight"].sum()
        bkg_err = np.sqrt((cat_bkg["weight"]**2).sum())
        
        # Apply data-driven reweighting if enabled
        if REWEIGHT:
            sb_query = " | ".join([f"(H_mass > {sb[0]} & H_mass < {sb[1]})" for sb in SIDEBAND])
            sb_bkg = vbf_region_bkg.query(f"{score_2j_col} >= {low} and {score_2j_col} < {high} and ({sb_query})")
            sb_data = data_data.query(f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold} and {score_2j_col} >= {low} and {score_2j_col} < {high}")
            
            if not sb_bkg.empty and not sb_data.empty:
                rw_factor = sb_data["weight"].sum() / sb_bkg["weight"].sum() if sb_bkg["weight"].sum() > 0 else 1
                bkg_yield *= rw_factor
                bkg_err *= rw_factor
        
        # Calculate significance
        significance, sig_error = calculate_significance(sig_yield, bkg_yield, sig_err, bkg_err)
        
        vbf_region_cats.append({
            'name': cat_name,
            'low_threshold': low,
            'high_threshold': high,
            'signal': sig_yield,
            'background': bkg_yield,
            'vbf_yield': vbf_yield,  # Add VBF yield
            'vbf_fraction': vbf_yield / sig_yield if sig_yield > 0 else 0,  # VBF purity in this category
            'significance': np.sqrt(significance),
            'S/B': sig_yield/bkg_yield if bkg_yield > 0 else 0
        })
        
        print(f"{cat_name}: S={sig_yield:.2f}±{sig_err:.2f}, B={bkg_yield:.2f}±{bkg_err:.2f}, VBF={vbf_yield:.2f} ({vbf_yield/sig_yield*100:.1f}%), S/B={sig_yield/bkg_yield if bkg_yield > 0 else 0:.3f}, Significance={np.sqrt(significance):.3f}")
    
    # Then calculate significances for the remaining events categories
    remaining_cats = []
    
    for i in range(len(best_thresholds) - 1):
        cat_name = f"Remaining_Cat{i}"
        low, high = best_thresholds[i], best_thresholds[i+1]
        
        # Select events in this category
        cat_sig = remaining_sig.query(f"{score_col} >= {low} and {score_col} < {high} and H_mass > {MASS_WINDOW[0]} and H_mass < {MASS_WINDOW[1]}")
        cat_bkg = remaining_bkg.query(f"{score_col} >= {low} and {score_col} < {high} and H_mass > {MASS_WINDOW[0]} and H_mass < {MASS_WINDOW[1]}")
        
        # Calculate VBF events in this category
        cat_vbf = VBF_data.query(f"(n_jets < 2 or (n_jets >= 2 and {vbf_score_col} < {vbf_score_threshold})) and {score_col} >= {low} and {score_col} < {high} and H_mass > {MASS_WINDOW[0]} and H_mass < {MASS_WINDOW[1]}")
        vbf_yield = cat_vbf["weight"].sum() if not cat_vbf.empty else 0
        
        # Calculate yields and significance
        sig_yield = cat_sig["weight"].sum()
        sig_err = np.sqrt((cat_sig["weight"]**2).sum())
        bkg_yield = cat_bkg["weight"].sum()
        bkg_err = np.sqrt((cat_bkg["weight"]**2).sum())
        
        # Apply data-driven reweighting if enabled
        if REWEIGHT:
            sb_query = " | ".join([f"(H_mass > {sb[0]} & H_mass < {sb[1]})" for sb in SIDEBAND])
            sb_bkg = remaining_bkg.query(f"{score_col} >= {low} and {score_col} < {high} and ({sb_query})")
            sb_data = remaining_data.query(f"{score_col} >= {low} and {score_col} < {high}")
            
            if not sb_bkg.empty and not sb_data.empty:
                rw_factor = sb_data["weight"].sum() / sb_bkg["weight"].sum() if sb_bkg["weight"].sum() > 0 else 1
                bkg_yield *= rw_factor
                bkg_err *= rw_factor
        
        # Calculate significance
        significance, sig_error = calculate_significance(sig_yield, bkg_yield, sig_err, bkg_err)
        
        remaining_cats.append({
            'name': cat_name,
            'low_threshold': low,
            'high_threshold': high,
            'signal': sig_yield,
            'background': bkg_yield,
            'vbf_yield': vbf_yield,  # Add VBF yield
            'vbf_fraction': vbf_yield / sig_yield if sig_yield > 0 else 0,  # VBF purity in this category
            'significance': np.sqrt(significance),
            'S/B': sig_yield/bkg_yield if bkg_yield > 0 else 0
        })
        
        print(f"{cat_name}: S={sig_yield:.2f}±{sig_err:.2f}, B={bkg_yield:.2f}±{bkg_err:.2f}, VBF={vbf_yield:.2f} ({vbf_yield/sig_yield*100:.1f}%), S/B={sig_yield/bkg_yield if bkg_yield > 0 else 0:.3f}, Significance={np.sqrt(significance):.3f}")
    
    # ----- Calculate total significance -----
    total_sig = sum(cat['signal'] for cat in vbf_region_cats + remaining_cats)
    total_bkg = sum(cat['background'] for cat in vbf_region_cats + remaining_cats)
    total_vbf = sum(cat['vbf_yield'] for cat in vbf_region_cats + remaining_cats)
    total_significance, _ = calculate_significance(total_sig, total_bkg)
    sum_quadrature = np.sqrt(sum(cat['significance']**2 for cat in vbf_region_cats + remaining_cats))
    
    print(f"\nTotal signal: {total_sig:.2f}")
    print(f"Total background: {total_bkg:.2f}")
    print(f"Total VBF: {total_vbf:.2f} ({total_vbf/total_sig*100:.1f}%)")
    print(f"Total S/B: {total_sig/total_bkg if total_bkg > 0 else 0:.3f}")
    print(f"Total significance (direct): {np.sqrt(total_significance):.4f}")
    print(f"Total significance (quadrature sum): {sum_quadrature:.4f}")
    
    # ----- Save configuration to JSON -----
    def convert_to_python_float(obj):
        """Convert numpy float types to standard Python float for JSON serialization"""
        if isinstance(obj, dict):
            return {k: convert_to_python_float(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_python_float(i) for i in obj]
        elif isinstance(obj, np.number):
            return float(obj)
        else:
            return obj

    config = {
        "vbf_enriched_region": {
            "score_column": vbf_score_col,
            "threshold": vbf_score_threshold,
            "binning": {
                "score_column": score_2j_col,
                "thresholds": vbf_region_thresholds
            },
            "categories": vbf_region_cats
        },
        "remaining_region": {
            "score_column": score_col,
            "thresholds": best_thresholds,
            "categories": remaining_cats
        },
        "total_significance": sum_quadrature
    }
    
    # Convert all numpy float values to Python float before saving to JSON
    config = convert_to_python_float(config)
    
    with open("figs/optimal_categorization.json", "w") as f:
        json.dump(config, f, indent=2)
    
    # ----- Plot the score distributions -----
    plt.figure(figsize=(12, 6))
    
    # VBF-enriched region
    ax1 = plt.subplot(1, 2, 1)
    sig_hist, bins = np.histogram(
        vbf_region_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_2j_col], 
        bins=50, range=(0, 1), 
        weights=vbf_region_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    sig_hist = sig_hist / np.sum(sig_hist) if np.sum(sig_hist) > 0 else sig_hist
    
    bkg_hist, _ = np.histogram(
        vbf_region_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_2j_col], 
        bins=bins, range=(0, 1), 
        weights=vbf_region_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    bkg_hist = bkg_hist / np.sum(bkg_hist) if np.sum(bkg_hist) > 0 else bkg_hist
    
    centers = (bins[:-1] + bins[1:]) / 2
    width = bins[1] - bins[0]
    
    plt.bar(centers, sig_hist, width=width, alpha=0.5, label="Signal", color="blue")
    plt.bar(centers, bkg_hist, width=width, alpha=0.5, label="Background", color="red")
    
    for threshold in vbf_region_thresholds:
        plt.axvline(x=threshold, color='green', linestyle='--')
    
    plt.xlabel(score_2j_col)
    plt.ylabel("Normalized Events")
    plt.title("VBF-enriched Region")
    plt.legend()
    plt.grid(alpha=0.3)
    
    # Remaining region
    ax2 = plt.subplot(1, 2, 2)
    sig_hist, bins = np.histogram(
        remaining_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_col], 
        bins=50, range=(0, 1), 
        weights=remaining_sig.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    sig_hist = sig_hist / np.sum(sig_hist) if np.sum(sig_hist) > 0 else sig_hist
    
    bkg_hist, _ = np.histogram(
        remaining_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")[score_col], 
        bins=bins, range=(0, 1), 
        weights=remaining_bkg.query(f"H_mass > {MASS_WINDOW[0]} & H_mass < {MASS_WINDOW[1]}")["weight"]
    )
    bkg_hist = bkg_hist / np.sum(bkg_hist) if np.sum(bkg_hist) > 0 else bkg_hist
    
    centers = (bins[:-1] + bins[1:]) / 2
    width = bins[1] - bins[0]
    
    plt.bar(centers, sig_hist, width=width, alpha=0.5, label="Signal", color="blue")
    plt.bar(centers, bkg_hist, width=width, alpha=0.5, label="Background", color="red")
    
    for threshold in best_thresholds:
        plt.axvline(x=threshold, color='green', linestyle='--')
    
    plt.xlabel(score_col)
    plt.ylabel("Normalized Events")
    plt.title("Remaining Events Region")
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figs/score_distributions.png')
    
    print("\nAnalysis completed. Results saved to 'figs/optimal_categorization.json'")
    print("Score distributions plot saved to 'figs/score_distributions.png'")

    # ----- Plot Higgs mass distributions for each category -----
    print("\nPlotting H_mass distributions for each category...")
    
    # Create directory for mass distribution plots if it doesn't exist
    if not os.path.exists("figs/mass_distributions"):
        os.makedirs("figs/mass_distributions")
    
    # Define color dictionary for different background processes
    color_dict = {
        r"Z$+\gamma$": "#3f90da", 
        "Z+Fake Photon": "#ffa90e", 
        r"VBSZ+$\gamma$": "#92dadd", 
        r"t$\bar{t}$": "#e76300", 
        r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01",
        "multiboson": "#832db6", 
        r"t$\bar{t}$+X": "#94a4a2"
    }
    
    # Define background labels mapping
    background_labels = {
        "ZGToLLG": r"Z$+\gamma$", 
        "DYJetsToLL": "Z+Fake Photon", 
        "EWKZ2J": r"VBSZ+$\gamma$"
    }
    
    # Define mass range and binning for histograms
    mass_bins = 85
    mass_range = (95, 180)
    
    # Function to plot H_mass distribution for a category
    def plot_category_mass_distribution(category_name, category_query, output_filename):
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Prepare data for plotting
        bkg_hists = []
        bkg_labels = []
        bkg_colors = []
        
        # Process background samples
        for bkg_name in bkg_names:
            try:
                # Load background samples directly from the original file
                bkg_file = f"{filepath}/{bkg_name}.root"
                if os.path.exists(bkg_file):
                    tree = uproot.open(bkg_file)["all_jet"]
                    # Get the required variables
                    cols = ["bdt_score", "bdt_score_t", "vbf_score", "vbf_score_t",
                            "bdt_2j_score", "bdt_2j_score_t", "weight", "H_mass", 
                            "event", "n_jets"]
                    bkg_df = tree.arrays(cols, library="pd")
                    
                    # Apply category query conditions
                    query_parts = category_query.replace("and", "&")
                    # Handle score column names
                    score_col_actual = score_col if TRANSFORM else score_col.replace("_t", "")
                    vbf_score_col_actual = vbf_score_col if TRANSFORM else vbf_score_col.replace("_t", "")
                    score_2j_col_actual = score_2j_col if TRANSFORM else score_2j_col.replace("_t", "")
                    
                    # Replace column names in the query
                    query_parts = query_parts.replace(score_col, score_col_actual)
                    query_parts = query_parts.replace(vbf_score_col, vbf_score_col_actual)
                    query_parts = query_parts.replace(score_2j_col, score_2j_col_actual)
                    
                    # Apply mass window
                    final_query = f"({query_parts}) & (H_mass >= {mass_range[0]}) & (H_mass <= {mass_range[1]})"
                    bkg_sample = bkg_df.query(final_query)
                    
                    if len(bkg_sample) > 0:
                        # Create histogram
                        hist, bin_edges = np.histogram(
                            bkg_sample["H_mass"], 
                            bins=mass_bins, 
                            range=mass_range, 
                            weights=bkg_sample["weight"]
                        )
                        bkg_hists.append(hist)
                        bkg_labels.append(background_labels.get(bkg_name, bkg_name))
                        bkg_colors.append(color_dict.get(background_labels.get(bkg_name, bkg_name), "#999999"))
            except Exception as e:
                print(f"Error processing background {bkg_name}: {e}")
        
        # Create signal histogram
        sig_hist = np.zeros(mass_bins)
        sig_edges = None
        vbf_hist = None  # Separate VBF signal for comparison
        
        # Process signal samples directly from files
        for sig_name in sig_names:
            try:
                # Load signal samples directly from the original file
                sig_file = f"{filepath}/{sig_name}.root"
                if os.path.exists(sig_file):
                    tree = uproot.open(sig_file)["all_jet"]
                    # Get the required variables
                    cols = ["bdt_score", "bdt_score_t", "vbf_score", "vbf_score_t",
                            "bdt_2j_score", "bdt_2j_score_t", "weight", "H_mass", 
                            "event", "n_jets"]
                    sig_df = tree.arrays(cols, library="pd")
                    
                    # Apply category query conditions
                    query_parts = category_query.replace("and", "&")
                    # Handle score column names
                    score_col_actual = score_col if TRANSFORM else score_col.replace("_t", "")
                    vbf_score_col_actual = vbf_score_col if TRANSFORM else vbf_score_col.replace("_t", "")
                    score_2j_col_actual = score_2j_col if TRANSFORM else score_2j_col.replace("_t", "")
                    
                    # Replace column names in the query
                    query_parts = query_parts.replace(score_col, score_col_actual)
                    query_parts = query_parts.replace(vbf_score_col, vbf_score_col_actual)
                    query_parts = query_parts.replace(score_2j_col, score_2j_col_actual)
                    
                    # Apply mass window
                    final_query = f"({query_parts}) & (H_mass >= {mass_range[0]}) & (H_mass <= {mass_range[1]})"
                    sig_sample = sig_df.query(final_query)
                    
                    if len(sig_sample) > 0:
                        # Create histogram
                        if sig_edges is None:
                            hist, sig_edges = np.histogram(
                                sig_sample["H_mass"], 
                                bins=mass_bins, 
                                range=mass_range, 
                                weights=sig_sample["weight"]
                            )
                            sig_hist = hist
                        else:
                            hist = np.histogram(
                                sig_sample["H_mass"], 
                                bins=sig_edges, 
                                weights=sig_sample["weight"]
                            )[0]
                            sig_hist += hist
                        
                        # If this is VBF signal, store separately for purity calculation
                        if sig_name == "VBF_M125":
                            vbf_hist = hist
            except Exception as e:
                print(f"Error processing signal {sig_name}: {e}")
        
        # Create data histogram
        data_hist = np.zeros(mass_bins)
        
        # Process data samples directly from files
        for data_name in data_names:
            try:
                # Load data samples directly from the original file
                data_file = f"{filepath}/{data_name}.root"
                if os.path.exists(data_file):
                    tree = uproot.open(data_file)["all_jet"]
                    # Get the required variables
                    cols = ["bdt_score", "bdt_score_t", "vbf_score", "vbf_score_t",
                            "bdt_2j_score", "bdt_2j_score_t", "weight", "H_mass", 
                            "event", "n_jets"]
                    data_df = tree.arrays(cols, library="pd")
                    
                    # Apply category query conditions
                    query_parts = category_query.replace("and", "&")
                    # Handle score column names
                    score_col_actual = score_col if TRANSFORM else score_col.replace("_t", "")
                    vbf_score_col_actual = vbf_score_col if TRANSFORM else vbf_score_col.replace("_t", "")
                    score_2j_col_actual = score_2j_col if TRANSFORM else score_2j_col.replace("_t", "")
                    
                    # Replace column names in the query
                    query_parts = query_parts.replace(score_col, score_col_actual)
                    query_parts = query_parts.replace(vbf_score_col, vbf_score_col_actual)
                    query_parts = query_parts.replace(score_2j_col, score_2j_col_actual)
                    
                    # Apply mass window
                    final_query = f"({query_parts}) & (H_mass >= {mass_range[0]}) & (H_mass <= {mass_range[1]}) & (H_mass < 120 or H_mass > 130)"
                    data_sample = data_df.query(final_query)
                    
                    if len(data_sample) > 0:
                        # Create histogram using signal bins for consistency
                        if sig_edges is not None:
                            hist = np.histogram(
                                data_sample["H_mass"], 
                                bins=sig_edges, 
                                weights=data_sample["weight"]
                            )[0]
                            data_hist += hist
            except Exception as e:
                print(f"Error processing data {data_name}: {e}")
        
        # Calculate scale factor for background normalization (optional)
        if np.sum(data_hist) > 0 and np.sum(np.sum(bkg_hists, axis=0)) > 0:
            sf = np.sum(data_hist) / (np.sum(np.sum(bkg_hists, axis=0)) - np.sum(bkg_hists, axis=0)[25:35].sum())
            bkg_hists = [bkg * sf for bkg in bkg_hists]
            print(f"Applied scale factor: {sf:.4f} to background in {category_name}")
            
        # Plot stacked background histograms
        if bkg_hists and sig_edges is not None:
            hep.histplot(bkg_hists, sig_edges, stack=True, label=bkg_labels, histtype="fill", color=bkg_colors, ax=ax)
            
            # Plot data points with error bars, masking signal region
            bin_centers = (sig_edges[:-1] + sig_edges[1:]) / 2
            data_hist_masked = np.where((bin_centers > 120) & (bin_centers < 130), 0, data_hist)
            ax.errorbar(bin_centers, data_hist_masked, yerr=np.sqrt(data_hist_masked), fmt="o", label="Data", color="black")
            
            # Plot signal lines
            if vbf_hist is not None:
                ax.plot(bin_centers, vbf_hist, label="VBF", color="yellow", linestyle="--", linewidth=2)
            ax.plot(bin_centers, sig_hist, label="Signal", color="red", linewidth=3)

            vbf_purity = np.sum(vbf_hist) / np.sum(sig_hist) if np.sum(sig_hist) > 0 else 0
            ax.text(0.05, 0.95, f"VBF Purity: {vbf_purity:.2f}", transform=ax.transAxes, fontsize=20, verticalalignment='top')
            
            # Set plot properties
            ax.set_xlim(mass_range)
            max_height = max(np.max(np.sum(bkg_hists, axis=0) if bkg_hists else 0), np.max(data_hist if len(data_hist) > 0 else 0))
            ax.set_ylim(0, max_height * 1.2)
            
            # Add labels and title
            ax.set_xlabel("Higgs Mass [GeV]")
            ax.set_ylabel("Events")
            ax.title.set_text(f"Events in {category_name}")
            ax.legend(loc="upper right")
            
            # Save figure
            plt.tight_layout()
            plt.savefig(output_filename)
            plt.close()
            print(f"Saved mass distribution plot for {category_name} to {output_filename}")
    
    # Plot VBF-enriched categories
    for i, cat in enumerate(vbf_region_cats):
        cat_name = cat['name']
        low = cat['low_threshold']
        high = cat['high_threshold']
        
        # Define query for this category
        category_query = f"n_jets >= 2 and {vbf_score_col} >= {vbf_score_threshold} and {score_2j_col} >= {low} and {score_2j_col} < {high}"
        
        # Plot and save
        output_file = f"figs/mass_distributions/{cat_name}_hmass.png"
        plot_category_mass_distribution(cat_name, category_query, output_file)
    
    # Plot remaining categories
    for i, cat in enumerate(remaining_cats):
        cat_name = cat['name']
        low = cat['low_threshold']
        high = cat['high_threshold']
        
        # Define query for this category
        category_query = f"(n_jets < 2 or (n_jets >= 2 and {vbf_score_col} < {vbf_score_threshold})) and {score_col} >= {low} and {score_col} < {high}"
        
        # Plot and save
        output_file = f"figs/mass_distributions/{cat_name}_hmass.png"
        plot_category_mass_distribution(cat_name, category_query, output_file)
    
    # Plot combined categories
    output_file = f"figs/mass_distributions/all_categories_hmass.png"
    plot_category_mass_distribution("All Categories", "H_mass > 0", output_file)  # Simple condition that applies to all events
    
    print("\nFinished plotting H_mass distributions for all categories")

if __name__ == "__main__":
    main()
