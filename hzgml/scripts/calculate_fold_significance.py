#!/usr/bin/env python
"""Calculate significance for each fold and val/test dataset based on existing boundary values"""

import os
import json
import argparse
import numpy as np
from ROOT import TFile, TH1F, TH1, gROOT
from categorizer import calc_sig, categorizer
import matplotlib.pyplot as plt
import copy

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Calculate significance for each fold and dataset")
    parser.add_argument('-i', '--input', default='test/significances/0_0_two_jet_1D_4.json', help='JSON file with boundaries')
    parser.add_argument('-r', '--region', default='two_jet', 
                        choices=['all_jet', 'two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH'],
                        help='Region to process')
    parser.add_argument('-d', '--datapath', default='/eos/home-j/jiehan/root/outputs/', help='Path to root files')
    parser.add_argument('-v', '--variable', default='bdt', choices=['bdt', 'NN'], help='MVA variable')
    parser.add_argument('-t', '--transform', type=bool, default=True, help='Use transformed scores')
    parser.add_argument('-p', '--prime', type=int, default=314159, help='Prime number for fold division')
    parser.add_argument('-o', '--output', default=None, help='Output JSON file path')
    parser.add_argument('-e', '--estimate', default='fullSimrw', 
                        choices=['fullSim', 'fullSimrw', 'data_sid'], help='Significance estimation method')
    parser.add_argument('--skip', action='store_true', default=False, help='Skip the hadd part')
    parser.add_argument('-n', '--nfold', type=int, default=4, help='Number of folds')
    parser.add_argument('-s', '--shield', type=int, default=-1, help='Shield variable index')
    parser.add_argument('-a', '--add', type=int, default=-1, help='Add variable index')
    
    return parser.parse_args()

def calculate_significance_for_fold(data_dir, region, variable, boundaries, transform, fold, nfold, prime, estimate):
    """Calculate significance for a specific fold"""
    # Open ROOT files based on estimation method
    f_sig = TFile(f'{data_dir}/{region}/sig.root')
    t_sig = f_sig.Get(region)
    
    f_bkgmc = t_bkgmc = f_data_sid = t_data_sid = None
    if estimate in ["fullSim", "fullSimrw"]:
        f_bkgmc = TFile(f'{data_dir}/{region}/bkgmc.root')
        t_bkgmc = f_bkgmc.Get(region)
    if estimate in ["fullSimrw", "data_sid"]:
        f_data_sid = TFile(f'{data_dir}/{region}/Data.root')
        t_data_sid = f_data_sid.Get(region)
    
    # Setup for analysis
    nscan = 100  # Fixed bin count
    fold_condition = "1" if fold == -1 else f"((event%{prime})%{nfold}=={fold})"
    signal_window = f"weight*((H_mass>=120&&H_mass<=130)&&{fold_condition})"
    sideband = f"weight*((H_mass>=100&&H_mass<=180)&&!(H_mass>=120&&H_mass<=130)&&{fold_condition})"
    score_var = f"{variable}_score{'_t' if transform else ''}"
    
    # Fill histograms
    h_sig = TH1F('h_sig', 'h_sig', nscan, 0, 1.)
    t_sig.Draw(f"{score_var}>>h_sig", signal_window)
    
    # Create background histograms
    h_bkgmc_cen = h_bkgmc_sid = h_data_sid = None
    if estimate in ["fullSim", "fullSimrw"]:
        h_bkgmc_cen = TH1F('h_bkgmc_cen', 'h_bkgmc_cen', nscan, 0, 1.)
        t_bkgmc.Draw(f"{score_var}>>h_bkgmc_cen", signal_window)
    
    if estimate == "fullSimrw":
        h_bkgmc_sid = TH1F('h_bkgmc_sid', 'h_bkgmc_sid', nscan, 0, 1.)
        t_bkgmc.Draw(f"{score_var}>>h_bkgmc_sid", sideband)
    
    if estimate in ["fullSimrw", "data_sid"]:
        h_data_sid = TH1F('h_data_sid', 'h_data_sid', nscan, 0, 1.)
        t_data_sid.Draw(f"{score_var}>>h_data_sid", sideband)
        if estimate == "data_sid":
            h_data_sid.Scale(0.20)
    
    # Create categorizer objects
    if estimate == "data_sid":
        cgz = categorizer(h_sig, h_data_sid)
    elif estimate == "fullSimrw":
        cgz = categorizer(h_sig, h_bkgmc_cen, h_bkg_rw_num=h_data_sid, h_bkg_rw_den=h_bkgmc_sid)
    else:  # fullSim
        cgz = categorizer(h_sig, h_bkgmc_cen)
    
    # Store original background
    orig_bkg_values = [cgz.h_bkg.GetBinContent(i) for i in range(1, nscan+1)]
    orig_bkg_errors = [cgz.h_bkg.GetBinError(i) for i in range(1, nscan+1)]
    cgz_orig = copy.deepcopy(cgz)
    
    # Setup bin boundaries
    boundaries_with_end = boundaries + [1.0] if boundaries[-1] != 1.0 else boundaries
    bin_boundaries = [int(b * nscan) + 1 for b in boundaries_with_end]
    
    # Helper function to calculate bin integrals and handle reweighting
    def get_bin_data(categorizer_obj, bin_boundaries, do_reweight=False):
        bin_sigs, bin_bkgs, bin_sig_errs, bin_bkg_errs = [], [], [], []
        
        for i in range(len(bin_boundaries) - 1):
            bl, br = bin_boundaries[i], bin_boundaries[i+1] - 1
            nsig, dsig = categorizer_obj.h_sig.IntegralAndError(bl, br)
            nbkg, dbkg = categorizer_obj.h_bkg.IntegralAndError(bl, br)
            
            # Apply reweighting if needed
            if do_reweight and estimate == "fullSimrw" and nbkg != 0:
                rw_num = cgz.h_bkg_rw_num.IntegralAndError(bl, br)[0]
                rw_den = cgz.h_bkg_rw_den.IntegralAndError(bl, br)[0]
                if rw_den > 0:
                    nbkg *= rw_num / rw_den
                    dbkg *= rw_num / rw_den
            
            bin_sigs.append(nsig)
            bin_bkgs.append(nbkg)
            bin_sig_errs.append(dsig)
            bin_bkg_errs.append(dbkg)
        
        return bin_sigs, bin_bkgs, bin_sig_errs, bin_bkg_errs
    
    # Calculate for original background
    bin_sigs_orig, bin_bkgs_orig, bin_sig_errs_orig, bin_bkg_errs_orig = get_bin_data(cgz_orig, bin_boundaries, True)
    
    # Calculate significance for original background
    zs_orig = calc_sig(np.array(bin_sigs_orig), np.array(bin_bkgs_orig), 
                      np.array(bin_sig_errs_orig), np.array(bin_bkg_errs_orig))
    z_orig = np.sqrt((zs_orig[0]**2).sum())
    u_orig = np.sqrt((zs_orig[0]**2 * zs_orig[1]**2).sum()) / z_orig if z_orig > 0 else 0
    
    # Smooth background for standard calculation
    cgz.smooth_sim(1, SorB='B')
    
    # Ensure valid bin boundaries
    if bin_boundaries[0] != 1:
        bin_boundaries = [1] + bin_boundaries[1:]
    if bin_boundaries[-1] != nscan + 1:
        bin_boundaries[-1] = nscan + 1
    
    # Calculate for smoothed background
    bin_sigs, bin_bkgs, bin_sig_errs, bin_bkg_errs = get_bin_data(cgz, bin_boundaries, True)
    
    # Calculate significance for smoothed background
    zs = calc_sig(np.array(bin_sigs), np.array(bin_bkgs), 
                 np.array(bin_sig_errs), np.array(bin_bkg_errs))
    z = np.sqrt((zs[0]**2).sum())
    u = np.sqrt((zs[0]**2 * zs[1]**2).sum()) / z if z > 0 else 0
    
    # Get histogram data for plotting
    sig_bins, sig_values, sig_errors = [], [], []
    bkg_bins, bkg_values, bkg_errors, bkg_smooth_values = [], [], [], []
    
    for i in range(1, nscan+1):
        bin_center = float(i-0.5)/nscan
        sig_bins.append(bin_center)
        sig_values.append(cgz.h_sig.GetBinContent(i))
        sig_errors.append(cgz.h_sig.GetBinError(i))
        bkg_bins.append(bin_center)
        bkg_values.append(orig_bkg_values[i-1])
        bkg_errors.append(orig_bkg_errors[i-1])
        bkg_smooth_values.append(cgz.h_bkg.GetBinContent(i))
    
    # Clean up resources
    for f in [f_sig, f_bkgmc, f_data_sid]:
        if f: f.Close()
    
    return (z, u, bin_sigs, bin_bkgs, bin_sig_errs, bin_bkg_errs, zs[0], 
            sig_bins, sig_values, sig_errors, bkg_bins, bkg_values, bkg_errors, bkg_smooth_values,
            z_orig, u_orig, zs_orig[0])

def plot_fold_distribution(dataset, fold, boundaries, sig_bins, sig_values, sig_errors, 
                          bkg_bins, bkg_values, bkg_errors, bkg_smooth_values, significance, output_dir, 
                          bin_zs=None, bin_z_errs=None, significance_err=None, significance_orig=None, 
                          significance_orig_err=None, bin_zs_orig=None, bin_zs_orig_errs=None):
    """Plot signal and background distributions for a specific fold"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    plt.figure(figsize=(12, 10))
    plt.rcParams.update({'font.size': 18})
    
    # Normalize distributions
    sig_sum = sum(sig_values)
    bkg_sum = sum(bkg_values)
    bkg_smooth_sum = sum(bkg_smooth_values) if bkg_smooth_values else 0
    
    # Normalize values for plotting
    sig_values_norm = [v/sig_sum for v in sig_values] if sig_sum > 0 else sig_values
    sig_errors_norm = [e/sig_sum for e in sig_errors] if sig_sum > 0 else sig_errors
    bkg_values_norm = [v/bkg_sum for v in bkg_values] if bkg_sum > 0 else bkg_values
    bkg_errors_norm = [e/bkg_sum for e in bkg_errors] if bkg_sum > 0 else bkg_errors
    bkg_smooth_norm = [v/bkg_smooth_sum for v in bkg_smooth_values] if bkg_smooth_sum > 0 else bkg_smooth_values
    
    # Plot data
    plt.errorbar(sig_bins, sig_values_norm, yerr=sig_errors_norm, fmt='o', markersize=6,
                 color='blue', ecolor='blue', label='Signal', capsize=4, linewidth=3)
    plt.errorbar(bkg_bins, bkg_values_norm, yerr=bkg_errors_norm, fmt='o', markersize=6,
                 color='red', ecolor='red', label='Background', capsize=4, linewidth=3)
    plt.plot(bkg_bins, bkg_smooth_norm, color='darkred', linewidth=4, label='Background (smoothed)')
    
    # Set plot attributes
    plt.xlim(0, 1)
    plt.ylim(bottom=0)
    boundaries_with_end = boundaries + [1.0] if boundaries[-1] != 1.0 else boundaries
    
    # Add category significance values
    if bin_zs is not None:
        for i in range(len(boundaries_with_end) - 1):
            if i >= len(bin_zs):
                continue
                
            center = (boundaries_with_end[i] + boundaries_with_end[i+1]) / 2
            cat_sig = bin_zs[i]
            
            # Format text with error if available
            if fold == -1 and bin_z_errs is not None and i < len(bin_z_errs):
                text = f'significance: {cat_sig:.3f}±{bin_z_errs[i]:.3f}'
            else:
                text = f'significance: {cat_sig:.3f}'
                
            # Add original significance info
            if bin_zs_orig is not None and i < len(bin_zs_orig):
                if fold == -1 and bin_zs_orig_errs is not None and i < len(bin_zs_orig_errs):
                    text += f'\nw/o smooth: {bin_zs_orig[i]:.3f}±{bin_zs_orig_errs[i]:.3f}'
                else:
                    text += f'\nw/o smooth: {bin_zs_orig[i]:.3f}'
            
            # Alternate text positions for readability
            y_pos = plt.ylim()[1] * (0.625 if i % 2 == 1 else 0.7)
            plt.text(center, y_pos, text, horizontalalignment='center', fontsize=16,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='lightgray'))
    
    # Add boundary lines
    for boundary in boundaries:
        plt.axvline(x=boundary, color='red', linestyle='--', alpha=0.7, linewidth=3)
    
    # Set labels and title
    plt.xlabel('BDT Score', fontsize=20)
    plt.ylabel('Normalized Events', fontsize=20)
    title = f"{dataset.capitalize()} - Fold {fold}" if fold >= 0 else f"{dataset.capitalize()} - All Events"
    plt.title(title, fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(loc='upper right', fontsize=18)
    
    # Add total significance text
    if significance_orig is not None:
        if fold == -1 and significance_err is not None:
            sig_text = f"Significance: {significance:.3f}±{significance_err:.3f}\nWithout smoothing: {significance_orig:.3f}±{significance_orig_err:.3f}"
        else:
            sig_text = f"Significance: {significance:.3f}\nWithout smoothing: {significance_orig:.3f}"
    else:
        sig_text = f"Significance: {significance:.3f}" + (f"±{significance_err:.3f}" if fold == -1 and significance_err else "")
    
    plt.figtext(0.95, 0.75, sig_text, fontsize=20, horizontalalignment='right',
               bbox=dict(facecolor='white', alpha=0.8, edgecolor='lightgray'))
    
    # Save plot
    fold_str = str(fold) if fold >= 0 else "all"
    output_file = os.path.join(output_dir, f"{dataset}_fold_{fold_str}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Plot saved: {output_file}")

def plot_comparison(test_data, val_data, fold, boundaries, significance_test, significance_val, output_dir, 
                    delta_test=None, delta_val=None, test_bin_zs=None, val_bin_zs=None,
                    test_bin_z_errs=None, val_bin_z_errs=None, test_bin_zs_orig=None, val_bin_zs_orig=None,
                    test_bin_zs_orig_errs=None, val_bin_zs_orig_errs=None, significance_test_orig=None, 
                    significance_val_orig=None, delta_test_orig=None, delta_val_orig=None):
    """Plot comparison of test and validation distributions"""
    # Unpack data
    test_sig_bins, test_sig_values, test_bkg_bins, test_bkg_smooth_values = test_data
    val_sig_bins, val_sig_values, val_bkg_bins, val_bkg_smooth_values = val_data
    
    # Setup plot
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(14, 10))
    plt.rcParams.update({'font.size': 18})
    
    # Normalize distributions
    test_sig_sum = sum(test_sig_values)
    test_bkg_sum = sum(test_bkg_smooth_values)
    val_sig_sum = sum(val_sig_values)
    val_bkg_sum = sum(val_bkg_smooth_values)
    
    # Create normalized values for plotting
    test_sig_norm = [v/test_sig_sum for v in test_sig_values] if test_sig_sum > 0 else test_sig_values
    val_sig_norm = [v/val_sig_sum for v in val_sig_values] if val_sig_sum > 0 else val_sig_values
    test_bkg_norm = [v/test_bkg_sum for v in test_bkg_smooth_values] if test_bkg_sum > 0 else test_bkg_smooth_values
    val_bkg_norm = [v/val_bkg_sum for v in val_bkg_smooth_values] if val_bkg_sum > 0 else val_bkg_smooth_values
    
    # Colors for validation dataset
    val_signal_color = 'deepskyblue'
    val_bkg_color = 'magenta'
    
    # Plot signals and backgrounds
    plt.plot(test_sig_bins, test_sig_norm, 'o-', markersize=6, linewidth=3, 
             color='blue', label='Test Signal')
    plt.plot(val_sig_bins, val_sig_norm, 's-', markersize=6, linewidth=3,
             color=val_signal_color, label='Val Signal')
    plt.plot(test_bkg_bins, test_bkg_norm, linewidth=4, 
             color='red', label='Test Background(smoothed)')
    plt.plot(val_bkg_bins, val_bkg_norm, linewidth=4,
             color=val_bkg_color, label='Val Background(smoothed)')
    
    # Add significance text
    is_overall = fold == -1 and delta_test is not None and delta_val is not None
    sig_text = ""
    
    if is_overall:
        sig_text += f'Test Sig: {significance_test:.3f}±{delta_test:.3f}'
        if significance_test_orig is not None:
            sig_text += f', w/o smooth: {significance_test_orig:.3f}'
            if delta_test_orig is not None:
                sig_text += f'±{delta_test_orig:.3f}'
        sig_text += '\n'
        sig_text += f'Val Sig: {significance_val:.3f}±{delta_val:.3f}'
        if significance_val_orig is not None:
            sig_text += f', w/o smooth: {significance_val_orig:.3f}'
            if delta_val_orig is not None:
                sig_text += f'±{delta_val_orig:.3f}'
    else:
        sig_text += f'Test Sig: {significance_test:.3f}'
        if significance_test_orig is not None:
            sig_text += f', w/o smooth: {significance_test_orig:.3f}'
        sig_text += '\n'
        sig_text += f'Val Sig: {significance_val:.3f}'
        if significance_val_orig is not None:
            sig_text += f', w/o smooth: {significance_val_orig:.3f}'
    
    plt.figtext(0.12, 0.85, sig_text, fontsize=18,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='lightgray'))
    
    # Set plot attributes
    plt.xlim(0, 1)
    plt.ylim(bottom=0)
    boundaries_with_end = boundaries + [1.0] if boundaries[-1] != 1.0 else boundaries
    
    # Add category significance values
    if test_bin_zs is not None:
        add_category_texts(plt, boundaries_with_end, test_bin_zs, test_bin_z_errs, fold, 
                          0.75, 'Test', 'blue', test_bin_zs_orig, test_bin_zs_orig_errs)
    
    if val_bin_zs is not None:
        add_category_texts(plt, boundaries_with_end, val_bin_zs, val_bin_z_errs, fold, 
                          0.60, 'Val', val_signal_color, val_bin_zs_orig, val_bin_zs_orig_errs)
    
    # Add boundary lines
    for boundary in boundaries:
        plt.axvline(x=boundary, color='red', linestyle='--', alpha=0.7, linewidth=3)
    
    # Set labels and save plot
    plt.xlabel('BDT Score', fontsize=20)
    plt.ylabel('Normalized Events', fontsize=20)
    plt.title(f"Test vs Val - Fold {fold}" if fold >= 0 else "Test vs Val - All Events", fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(loc='upper right', fontsize=18)
    
    fold_str = str(fold) if fold >= 0 else "all"
    output_file = os.path.join(output_dir, f"comparison_fold_{fold_str}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Comparison plot saved: {output_file}")

def add_category_texts(plt_obj, boundaries, bin_zs, bin_z_errs, fold, y_factor, prefix, color, 
                      bin_zs_orig=None, bin_zs_orig_errs=None):
    """Add category significance texts to plots"""
    for i in range(len(boundaries) - 1):
        if i >= len(bin_zs):
            continue
        
        center = (boundaries[i] + boundaries[i+1]) / 2
        cat_sig = bin_zs[i]
        
        # Create significance text with error if available
        if fold == -1 and bin_z_errs is not None and i < len(bin_z_errs):
            text = f'{prefix} significance: {cat_sig:.3f}±{bin_z_errs[i]:.3f}'
        else:
            text = f'{prefix} significance: {cat_sig:.3f}'
        
        # Add original significance info when available
        if fold == -1 and bin_zs_orig is not None and i < len(bin_zs_orig):
            if bin_zs_orig_errs is not None and i < len(bin_zs_orig_errs):
                text += f'\nw/o smooth: {bin_zs_orig[i]:.3f}±{bin_zs_orig_errs[i]:.3f}'
            else:
                text += f'\nw/o smooth: {bin_zs_orig[i]:.3f}'
        
        # Alternate text positions for readability
        y_pos = plt_obj.ylim()[1] * (y_factor - 0.075 if i % 2 == 1 else y_factor)
        plt_obj.text(center, y_pos, text, horizontalalignment='center', fontsize=14, 
                    color=color, bbox=dict(facecolor='white', alpha=0.7, edgecolor='lightgray'))

def prepare_data_files(data_dir, region, skip=False):
    """Prepare required data files"""
    if skip:
        return
    
    sigs = ['ggH_M125', 'VBF_M125'] #'WminusH_M125', 'WplusH_M125', 'ZH_M125', 'ttH_M125'
    bkgs = ["ZGToLLG", "DYJetsToLL", "EWKZ2J"] #"ZG2JToG2L2J", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "WWG", "WZG", "ZZG", "ttZJets", "ttWJets"
    
    # Merge signal files
    siglist = ''
    for sig in sigs:
        if os.path.isfile('%s/%s/%s.root' % (data_dir, region, sig)):
            siglist += ' %s/%s/%s.root' % (data_dir, region, sig)
    if siglist:
        os.system(f"hadd -f {data_dir}/{region}/sig.root{siglist}")
    
    # Merge background files
    bkglist = ''
    for bkg in bkgs:
        if os.path.isfile('%s/%s/%s.root' % (data_dir, region, bkg)):
            bkglist += ' %s/%s/%s.root' % (data_dir, region, bkg)
    if bkglist:
        os.system(f"hadd -f {data_dir}/{region}/bkgmc.root{bkglist}")

def read_config(datapath, input_file):
    """Read the JSON configuration file"""
    try:
        with open(f"{datapath}{input_file}", 'r') as f:
            content = f.read()
            # Find the last closing brace of the JSON object for error tolerance
            last_brace_pos = content.rstrip().rfind('}')
            if last_brace_pos > 0:
                # Only parse the valid JSON part
                valid_json = content[:last_brace_pos+1]
                config = json.loads(valid_json)
            else:
                print(f"Error: Could not find a valid JSON object in {datapath}{input_file}")
                return None, None
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return None, None
    
    # Extract boundaries correctly
    if isinstance(config['boundaries_values'], list) and len(config['boundaries_values']) > 0:
        if isinstance(config['boundaries_values'][0], list):
            boundaries_values = config['boundaries_values'][0]
        else:
            boundaries_values = config['boundaries_values']
    else:
        print("Error: Invalid boundaries_values in config")
        return None, None
    
    print(f"Using boundaries: {boundaries_values}")
    return config, boundaries_values

def calculate_bin_errors(fold_bin_zs, fold_bin_zs_orig):
    """Calculate standard errors for bin significances across folds"""
    bin_z_errs = []
    bin_z_orig_errs = []
    
    if not fold_bin_zs:
        return bin_z_errs, bin_z_orig_errs
    
    # Make sure all folds have the same number of bins
    min_bins = min(len(bz) for bz in fold_bin_zs)
    min_bins_orig = min(len(bz) for bz in fold_bin_zs_orig)
    
    # Calculate standard deviation for each bin in smoothed background
    for bin_idx in range(min_bins):
        bin_values = [fold_bz[bin_idx] for fold_bz in fold_bin_zs if bin_idx < len(fold_bz)]
        if bin_values:
            std_dev = np.std(bin_values)
            std_err = std_dev / np.sqrt(len(bin_values))
            bin_z_errs.append(2*std_err)
    
    # Calculate standard deviation for each bin in original background
    for bin_idx in range(min_bins_orig):
        bin_values = [fold_bz[bin_idx] for fold_bz in fold_bin_zs_orig if bin_idx < len(fold_bz)]
        if bin_values:
            std_dev = np.std(bin_values)
            std_err = std_dev
            bin_z_orig_errs.append(2*std_err)
    
    return bin_z_errs, bin_z_orig_errs

def process_dataset(dataset, data_dir, region, variable, boundaries_values, transform, nfold, prime, estimate, plots_dir, skip=False):
    """Process a single dataset (test or validation)"""
    prepare_data_files(data_dir, region, skip)
    
    fold_results = []
    fold_significances = []
    fold_bin_zs = []  # Store bin significances for each fold
    fold_significances_orig = []  # Store original (non-smoothed) significances
    fold_bin_zs_orig = []  # Store original bin significances for each fold
    distribution_data = {}
    
    # Calculate significance for each fold
    for fold in range(nfold):
        print(f"Calculating significance for {dataset} dataset, fold {fold}...")
        z, u, bin_sigs, bin_bkgs, bin_sig_errs, bin_bkg_errs, bin_zs, sig_bins, sig_values, sig_errors, bkg_bins, bkg_values, bkg_errors, bkg_smooth_values, z_orig, u_orig, bin_zs_orig = calculate_significance_for_fold(
            data_dir, region, variable, 
            boundaries_values, transform, fold, nfold, prime, estimate
        )
        
        # Store data for comparison plots
        distribution_data[fold] = {
            'significance': z,
            'significance_orig': z_orig,
            'sig_bins': sig_bins,
            'sig_values': sig_values,
            'sig_errors': sig_errors,
            'bkg_bins': bkg_bins,
            'bkg_values': bkg_values,
            'bkg_errors': bkg_errors,
            'bkg_smooth_values': bkg_smooth_values,
            'bin_zs': bin_zs,
            'bin_zs_orig': bin_zs_orig
        }
        
        # Plot and save distribution
        plot_fold_distribution(dataset, fold, boundaries_values, sig_bins, sig_values, sig_errors, 
                              bkg_bins, bkg_values, bkg_errors, bkg_smooth_values, z, plots_dir, 
                              bin_zs, significance_orig=z_orig, bin_zs_orig=bin_zs_orig)
        
        fold_results.append({
            'significance': float(z),
            'significance_orig': float(z_orig),
            'delta': float(u),
            'delta_orig': float(u_orig),
            'bin_sigs': [float(s) for s in bin_sigs],
            'bin_bkgs': [float(b) for b in bin_bkgs],
            'bin_sig_errs': [float(se) for se in bin_sig_errs],
            'bin_bkg_errs': [float(be) for be in bin_bkg_errs],
            'bin_zs': [float(bz) for bz in bin_zs],
            'bin_zs_orig': [float(bz) for bz in bin_zs_orig]
        })
        fold_significances.append(z)
        fold_significances_orig.append(z_orig)
        fold_bin_zs.append(bin_zs)
        fold_bin_zs_orig.append(bin_zs_orig)
        
        print(f"Fold {fold} significance: {z:.3f} ± {u:.3f} (without smoothing: {z_orig:.3f} ± {u_orig:.3f})")
    
    # Calculate total significance (all events)
    print(f"Calculating total significance for {dataset} dataset (all folds combined)...")
    total_z, total_u, total_bin_sigs, total_bin_bkgs, total_bin_sig_errs, total_bin_bkg_errs, total_bin_zs, sig_bins, sig_values, sig_errors, bkg_bins, bkg_values, bkg_errors, bkg_smooth_values, total_z_orig, total_u_orig, total_bin_zs_orig = calculate_significance_for_fold(
        data_dir, region, variable, 
        boundaries_values, transform, -1, nfold, prime, estimate
    )
    
    # Calculate standard errors for bin significances
    bin_z_errs, bin_z_orig_errs = calculate_bin_errors(fold_bin_zs, fold_bin_zs_orig)
    
    # Calculate standard deviation for original significances
    significance_orig_err = np.std(fold_significances_orig) if fold_significances_orig else None
    
    # Store overall distribution data
    distribution_data[-1] = {
        'significance': total_z,
        'significance_orig': total_z_orig,
        'significance_orig_err': significance_orig_err,
        'sig_bins': sig_bins,
        'sig_values': sig_values,
        'sig_errors': sig_errors,
        'bkg_bins': bkg_bins,
        'bkg_values': bkg_values,
        'bkg_errors': bkg_errors,
        'bkg_smooth_values': bkg_smooth_values,
        'bin_zs': total_bin_zs,
        'bin_z_errs': bin_z_errs,
        'bin_zs_orig': total_bin_zs_orig,
        'bin_zs_orig_errs': bin_z_orig_errs
    }
    
    # Calculate average significance across all folds
    avg_significance = sum(fold_significances) / len(fold_significances)
    std_dev_significance = np.std(fold_significances)
    significance_err = std_dev_significance
    
    print(f"{dataset} dataset average significance: {avg_significance:.3f} ± {std_dev_significance:.3f}")
    print(f"{dataset} dataset total significance: {total_z:.3f} ± {significance_err:.3f} (without smoothing: {total_z_orig:.3f} ± {significance_orig_err:.3f})")
    
    # Create dataset results dictionary
    dataset_results = {
        'avg_significance': float(avg_significance),
        'total_significance': float(total_z),
        'total_significance_orig': float(total_z_orig),
        'total_deviation': float(std_dev_significance)/np.sqrt(len(fold_significances)),
        'total_deviation_orig': float(significance_orig_err) if significance_orig_err else 0,
        'total_delta': float(total_u),
        'total_delta_orig': float(total_u_orig),
        'total_bin_sigs': [float(s) for s in total_bin_sigs],
        'total_bin_bkgs': [float(b) for b in total_bin_bkgs],
        'total_bin_zs': [float(bz) for bz in total_bin_zs],
        'total_bin_zs_orig': [float(bz) for bz in total_bin_zs_orig],
        'bin_zs': [float(bz) for bz in total_bin_zs],
        'bin_z_errs': [float(err) for err in bin_z_errs] if bin_z_errs else [],
        'bin_zs_orig': [float(bz) for bz in total_bin_zs_orig],
        'bin_zs_orig_errs': [float(err) for err in bin_z_orig_errs] if bin_z_orig_errs else [],
        'total_z': float(total_z),
        'significance_err': float(significance_err) if significance_err is not None else 0.0,
        'total_z_orig': float(total_z_orig),
        'significance_orig_err': float(significance_orig_err) if significance_orig_err is not None else 0.0
    }
    
    return fold_results, dataset_results, distribution_data

def calculate_propagated_error(bin_zs, bin_z_errs):
    """Calculate propagated error using error propagation formula"""
    # Check if inputs are None or empty
    if bin_z_errs is None or bin_zs is None:
        return 0.0
    
    # Convert to numpy arrays if they aren't already
    bin_zs_array = np.array(bin_zs)
    bin_z_errs_array = np.array(bin_z_errs)
    
    # Check if arrays are empty
    if bin_zs_array.size == 0 or bin_z_errs_array.size == 0:
        return 0.0
    
    try:
        # Calculate propagated error using proper numpy operations
        numerator = np.sum((bin_z_errs_array * bin_zs_array)**2)
        denominator = np.sum(bin_zs_array**2)
        
        if denominator <= 0:
            return 0.0
            
        delta = np.sqrt(numerator / denominator)
        return float(delta)  # Ensure we return a Python float, not numpy type
    except Exception as e:
        print(f"Error in calculating propagated error: {e}")
        return 0.0

def create_comparison_plots(distribution_data, boundaries_values, plots_dir, nfold):
    """Create comparison plots between test and validation datasets"""
    print("Creating comparison plots for test and validation datasets...")
    
    # Create fold comparison plots
    for fold in range(nfold):
        if fold in distribution_data['test'] and fold in distribution_data['val']:
            test_data = (
                distribution_data['test'][fold]['sig_bins'],
                distribution_data['test'][fold]['sig_values'],
                distribution_data['test'][fold]['bkg_bins'],
                distribution_data['test'][fold]['bkg_smooth_values']
            )
            val_data = (
                distribution_data['val'][fold]['sig_bins'],
                distribution_data['val'][fold]['sig_values'],
                distribution_data['val'][fold]['bkg_bins'],
                distribution_data['val'][fold]['bkg_smooth_values']
            )
            plot_comparison(
                test_data, val_data, fold, boundaries_values,
                distribution_data['test'][fold]['significance'],
                distribution_data['val'][fold]['significance'],
                plots_dir,
                test_bin_zs=distribution_data['test'][fold]['bin_zs'],
                val_bin_zs=distribution_data['val'][fold]['bin_zs'],
                significance_test_orig=distribution_data['test'][fold].get('significance_orig'),
                significance_val_orig=distribution_data['val'][fold].get('significance_orig')
            )

def calculate_errors(distribution_data, results):
    """Calculate propagated errors for the test and validation datasets"""
    if -1 not in distribution_data['test'] or -1 not in distribution_data['val']:
        return {}
    
    # Prepare error dictionary
    error_dict = {}
    
    # Calculate errors for test dataset
    test_bin_zs = distribution_data['test'][-1].get('bin_zs', None)
    test_bin_z_errs = distribution_data['test'][-1].get('bin_z_errs', None)
    test_bin_zs_orig = distribution_data['test'][-1].get('bin_zs_orig', None)
    test_bin_zs_orig_errs = distribution_data['test'][-1].get('bin_zs_orig_errs', None)
    
    # Calculate for validation dataset
    val_bin_zs = distribution_data['val'][-1].get('bin_zs', None)
    val_bin_z_errs = distribution_data['val'][-1].get('bin_z_errs', None)
    val_bin_zs_orig = distribution_data['val'][-1].get('bin_zs_orig', None)
    val_bin_zs_orig_errs = distribution_data['val'][-1].get('bin_zs_orig_errs', None)
    
    # Calculate propagated errors using error propagation formula
    test_delta = calculate_propagated_error(test_bin_zs, test_bin_z_errs) or results['datasets']['test']['total_deviation']
    val_delta = calculate_propagated_error(val_bin_zs, val_bin_z_errs) or results['datasets']['val']['total_deviation']
    
    # Calculate errors for original (unsmoothed) results
    delta_test_orig = calculate_propagated_error(test_bin_zs_orig, test_bin_zs_orig_errs) or distribution_data['test'][-1].get('significance_orig_err', 0)
    delta_val_orig = calculate_propagated_error(val_bin_zs_orig, val_bin_zs_orig_errs) or distribution_data['val'][-1].get('significance_orig_err', 0)
    
    # Store in dictionary
    error_dict = {
        'test_delta': test_delta,
        'val_delta': val_delta,
        'delta_test_orig': delta_test_orig,
        'delta_val_orig': delta_val_orig
    }
    
    return error_dict

def main():
    """Main function"""
    args = parse_args()
    
    # Configure ROOT
    gROOT.SetBatch(True)
    TH1.SetDefaultSumw2(1)
    
    # Read configuration file
    config, boundaries_values = read_config(args.datapath, args.input)
    if not config or not boundaries_values:
        return
    
    # Initialize results dictionary
    results = {
        'boundaries': config['boundaries'],
        'boundaries_values': boundaries_values,
        'original_significance': config.get('significance', 0),
        'original_delta': config.get('Delta_significance', 0),
        'variable': args.variable,
        'transform': args.transform,
        'prime': args.prime,
        'nfold': args.nfold,
        'folds': {},
        'datasets': {}
    }
    
    # Create plots directory
    plots_dir = f"{args.datapath}/plots/fold_distributions/{args.region}"
    
    # Dictionary to store distribution data
    distribution_data = {'test': {}, 'val': {}}
    
    # Process each dataset
    for dataset in ['val', 'test']:
        data_dir = f"{args.datapath}/{dataset}"
        fold_results, dataset_results, dataset_dist_data = process_dataset(
            dataset, data_dir, args.region, args.variable, boundaries_values,
            args.transform, args.nfold, args.prime, args.estimate, plots_dir, args.skip
        )
        
        # Store results
        results['folds'][dataset] = fold_results
        results['datasets'][dataset] = dataset_results
        distribution_data[dataset] = dataset_dist_data
    
    # Create comparison plots
    create_comparison_plots(distribution_data, boundaries_values, plots_dir, args.nfold)
    
    # Calculate propagated errors
    error_dict = calculate_errors(distribution_data, results)
    
    # If we have overall (-1) data for both test and val, create overall comparison
    if -1 in distribution_data['test'] and -1 in distribution_data['val']:
        # Update error values in results
        for dataset in ['test', 'val']:
            results['datasets'][dataset].update({
                'total_delta': float(error_dict[f'{dataset}_delta']),
                'total_delta_orig': float(error_dict[f'delta_{dataset}_orig'])
            })
        
        # Create plots for each dataset with errors
        for dataset in ['val', 'test']:
            plot_fold_distribution(
                dataset, -1, boundaries_values, 
                distribution_data[dataset][-1]['sig_bins'],
                distribution_data[dataset][-1]['sig_values'],
                distribution_data[dataset][-1]['sig_errors'],
                distribution_data[dataset][-1]['bkg_bins'],
                distribution_data[dataset][-1]['bkg_values'],
                distribution_data[dataset][-1]['bkg_errors'],
                distribution_data[dataset][-1]['bkg_smooth_values'],
                results['datasets'][dataset]['total_z'], plots_dir, 
                results['datasets'][dataset]['bin_zs'],
                results['datasets'][dataset].get('bin_z_errs', []),
                results['datasets'][dataset]['total_delta'],
                results['datasets'][dataset]['total_z_orig'], 
                results['datasets'][dataset]['total_delta_orig'],
                results['datasets'][dataset]['bin_zs_orig'],
                results['datasets'][dataset].get('bin_zs_orig_errs', [])
            )
        
        # Create overall comparison plot
        test_data = (
            distribution_data['test'][-1]['sig_bins'],
            distribution_data['test'][-1]['sig_values'],
            distribution_data['test'][-1]['bkg_bins'],
            distribution_data['test'][-1]['bkg_smooth_values']
        )
        val_data = (
            distribution_data['val'][-1]['sig_bins'],
            distribution_data['val'][-1]['sig_values'],
            distribution_data['val'][-1]['bkg_bins'],
            distribution_data['val'][-1]['bkg_smooth_values']
        )
        
        plot_comparison(
            test_data, val_data, -1, boundaries_values,
            distribution_data['test'][-1]['significance'],
            distribution_data['val'][-1]['significance'],
            plots_dir,
            delta_test=error_dict['test_delta'],
            delta_val=error_dict['val_delta'],
            test_bin_zs=distribution_data['test'][-1].get('bin_zs', None),
            val_bin_zs=distribution_data['val'][-1].get('bin_zs', None),
            test_bin_z_errs=distribution_data['test'][-1].get('bin_z_errs', None),
            val_bin_z_errs=distribution_data['val'][-1].get('bin_z_errs', None),
            test_bin_zs_orig=distribution_data['test'][-1].get('bin_zs_orig', None),
            val_bin_zs_orig=distribution_data['val'][-1].get('bin_zs_orig', None),
            test_bin_zs_orig_errs=distribution_data['test'][-1].get('bin_zs_orig_errs', None),
            val_bin_zs_orig_errs=distribution_data['val'][-1].get('bin_zs_orig_errs', None),
            significance_test_orig=distribution_data['test'][-1].get('significance_orig'),
            significance_val_orig=distribution_data['val'][-1].get('significance_orig'),
            delta_test_orig=error_dict['delta_test_orig'],
            delta_val_orig=error_dict['delta_val_orig']
        )
    
    # Save results to JSON file
    if args.output:
        output_path = args.output
    else:
        # Match the directory structure from categorization_1D.py
        test_dir = f"{args.datapath}/test"
        if not os.path.isdir(f'{test_dir}/significances/{args.region}'):
            print(f'INFO: Creating output folder: "{test_dir}/significances/{args.region}"')
            os.makedirs(f'{test_dir}/significances/{args.region}')
        output_path = f"{test_dir}/significances/{args.shield+1}_{args.add+1}_{args.region}_1D_fold_significance.json"
    
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)
    
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    main()
