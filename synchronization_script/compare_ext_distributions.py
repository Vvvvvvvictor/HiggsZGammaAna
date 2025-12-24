#!/usr/bin/env python3
"""
Compare distributions between DYJetsToLL and DYJetsToLL_ext files
across different years.
"""

import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

from pdb import set_trace as bp

plt.style.use(hep.style.CMS)

# Configuration
TREE = "zero_to_one_jet"  # Tree name for both folders
VAR = "H_mass"  # Variable to compare
WEIGHT = "weight_corr"  # Weight variable (will check if exists)
PATH1 = "/eos/home-j/jiehan/root/skimmed_ntuples_rui_new/DYJetsToLL/"
PATH2 = "/eos/home-j/jiehan/root/skimmed_ntuples_rui_new/DYJetsToLL_ext/"

# Years mapping based on file listing
years_mapping = {
    "2016preVFP": "2016preVFP",
    "2016postVFP": "2016postVFP", 
    "2017": "2017",
    "2018": "2018",
    "2022preEE": "2022preEE",
    "2022postEE": "2022postEE",
    "2023postBPix": "2023postBPix"
}

# Color scheme for different sample types
color_dict = {
    "DYJetsToLL": "#3f90da",
    "DYJetsToLL_ext": "#ffa90e"
}

# Variable configuration
var_config = {
    "range": (95, 180),
    "bins": 85,
    "title": r"$m_{\ell\ell\gamma}$"
}

def read_root_file(filepath, tree_name, variables):
    """Read ROOT file and return pandas DataFrame"""
    try:
        with uproot.open(f"{filepath}:{tree_name}") as tree:
            data = tree.arrays(variables, library="pd")
            data = data.query(f"H_mass>95 & H_mass<180")
        return data
    except Exception as e:
        print(f"Error reading {filepath}:{tree_name} - {e}")
        return None

def get_histogram_data(path, year, tree_name):
    """Get histogram data for a specific path and year"""
    filename = f"{year}.root"
    filepath = os.path.join(path, filename)
    
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None, None
    
    print(f"Reading {filepath}:{tree_name}")
    
    variables = [VAR]
    if VAR != "H_mass":
        variables.append("H_mass")
    
    # Check if weight exists
    try:
        with uproot.open(f"{filepath}:{tree_name}") as tree:
            if WEIGHT in tree.keys():
                variables.append(WEIGHT)
                use_weight = True
                print(f"Using weight: {WEIGHT}")
            else:
                use_weight = False
                print(f"Warning: Weight '{WEIGHT}' not found in {filename}, using unweighted")
    except Exception as e:
        print(f"Error checking weight in {filepath}: {e}")
        use_weight = False
    
    data = read_root_file(filepath, tree_name, variables)
    if data is None:
        return None, None
    
    # Filter valid data
    data = data.query(f"({VAR} > -900) & ({VAR} < 900)")
    
    # Get weights
    weights = data[WEIGHT].values if use_weight else np.ones(len(data))
    
    return data[VAR].values, weights

def create_combined_plot(years):
    """Create a single comparison plot combining all specified years"""
    print(f"\n=== Processing all years: {years} ===")
    
    all_dy_data = []
    all_dy_weights = []
    all_dy_ext_data = []
    all_dy_ext_weights = []

    for year in years:
        print(f"--- Reading data for year {year} ---")
        # Read data from both folders
        dy_data, dy_weights = get_histogram_data(PATH1, year, TREE)
        dy_ext_data, dy_ext_weights = get_histogram_data(PATH2, year, TREE)
        
        if dy_data is not None:
            all_dy_data.append(dy_data)
            all_dy_weights.append(dy_weights)
        
        if dy_ext_data is not None:
            all_dy_ext_data.append(dy_ext_data)
            all_dy_ext_weights.append(dy_ext_weights)

    # Concatenate data from all years
    if not all_dy_data and not all_dy_ext_data:
        print("No data found for any year, skipping plot generation.")
        return

    final_dy_data = np.concatenate(all_dy_data) if all_dy_data else None
    final_dy_weights = np.concatenate(all_dy_weights) if all_dy_weights else None
    final_dy_ext_data = np.concatenate(all_dy_ext_data) if all_dy_ext_data else None
    final_dy_ext_weights = np.concatenate(all_dy_ext_weights) if all_dy_ext_weights else None

    # Create histogram bins
    bins = np.linspace(var_config["range"][0], var_config["range"][1], var_config["bins"] + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_width = bins[1] - bins[0]
    
    # Create histograms
    dy_hist = None
    dy_hist_err = None
    dy_ext_hist = None
    dy_ext_hist_err = None
    
    if final_dy_data is not None:
        dy_hist, _ = np.histogram(final_dy_data, bins=bins, weights=final_dy_weights)
        dy_hist_err = np.sqrt(np.abs(dy_hist))
    
    if final_dy_ext_data is not None:
        dy_ext_hist, _ = np.histogram(final_dy_ext_data, bins=bins, weights=final_dy_ext_weights)
        dy_ext_hist_err = np.sqrt(np.abs(dy_ext_hist))
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), dpi=200, 
                                   gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
    
    # Main plot
    ax1.tick_params(axis='x', which='both', bottom=True, top=True, 
                    labelbottom=False, left=True, right=True, labelright=False)
    
    # Plot DYJetsToLL
    if dy_hist is not None:
        ax1.errorbar(bin_centers, dy_hist, yerr=dy_hist_err, 
                     xerr=bin_width/2, fmt="o", 
                     label=f"DYJetsToLL [N={np.sum(dy_hist):.1f}]", 
                     color=color_dict["DYJetsToLL"], markersize=5, linewidth=2)
    
    # Plot DYJetsToLL_ext
    if dy_ext_hist is not None:
        ax1.errorbar(bin_centers, dy_ext_hist, yerr=dy_ext_hist_err, 
                     xerr=bin_width/2, fmt="s", 
                     label=f"DYJetsToLL_ext [N={np.sum(dy_ext_hist):.1f}]", 
                     color=color_dict["DYJetsToLL_ext"], markersize=5, linewidth=2)
    
    # Formatting
    y_label = f"Events/{bin_width:.1f}"
    ax1.set_ylabel(y_label, fontsize=14)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(var_config["range"][0], var_config["range"][1])
    
    # Set y-axis limit
    max_y = 0
    if dy_hist is not None:
        max_y = max(max_y, dy_hist.max())
    if dy_ext_hist is not None:
        max_y = max(max_y, dy_ext_hist.max())
    ax1.set_ylim(0, 1.4 * max_y)
    
    # Add year annotation
    ax1.annotate(f"All Years Combined", xy=(0.02, 0.98), xycoords='axes fraction', 
                 fontsize=14, ha="left", va="top",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Ratio plot
    ax2.tick_params(axis='x', which='both', bottom=True, top=True, 
                    labelbottom=True, left=True, right=True, labelright=False)
    
    if dy_hist is not None and dy_ext_hist is not None:
        # Calculate ratio: DYJetsToLL / DYJetsToLL_ext
        dy_ext_hist_safe = np.where(dy_ext_hist == 0, 1e-8, dy_ext_hist)  # Avoid division by zero
        
        ratio = dy_hist / dy_ext_hist_safe
        ratio_err = dy_hist_err / dy_ext_hist_safe
        
        # Plot ratio
        mask = (ratio > 0) & (ratio < 5)  # Reasonable ratio range
        ax2.errorbar(bin_centers[mask], ratio[mask], yerr=ratio_err[mask], 
                     xerr=(bin_width/2), fmt="o", color="black", 
                     markersize=5, linewidth=2)
        
        ax2.axhline(1, color="red", linestyle="--", alpha=0.7)
        ax2.set_ylim(0, 3)
        ax2.set_ylabel("DY/DY_ext", fontsize=14)
    else:
        ax2.text(0.5, 0.5, "No ratio plot\n(missing data)", 
                 transform=ax2.transAxes, ha="center", va="center",
                 fontsize=12, alpha=0.7)
        ax2.set_ylim(0, 1)
    
    ax2.set_xlim(var_config["range"][0], var_config["range"][1])
    ax2.set_xlabel(var_config["title"], fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # Save plot
    output_dir = Path("plots/dy_comparison")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"H_mass_comparison_all_years.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    print(f"Saved combined plot: {output_file}")
    
    # Print summary
    print(f"Summary for all years combined:")
    if dy_hist is not None:
        print(f"  DYJetsToLL events: {np.sum(dy_hist):.1f}")
    if dy_ext_hist is not None:
        print(f"  DYJetsToLL_ext events: {np.sum(dy_ext_hist):.1f}")

def main():
    """Main function to process all years"""
    print("Starting DYJetsToLL vs DYJetsToLL_ext comparison analysis...")
    print(f"Input directory 1: {PATH1}")
    print(f"Input directory 2: {PATH2}")
    print(f"Variable: {VAR}")
    print(f"Bins: {var_config['bins']}, Range: {var_config['range']}")
    
    # Get available years from the years_mapping
    available_years = list(years_mapping.keys())
    
    print(f"\nProcessing years: {available_years}")
    
    # Process all years together
    try:
        create_combined_plot(available_years)
    except Exception as e:
        print(f"Error processing combined years: {e}")
    
    print("\nAnalysis completed!")
    print("Check the 'plots/dy_comparison/' directory for output plots.")

if __name__ == "__main__":
    main()