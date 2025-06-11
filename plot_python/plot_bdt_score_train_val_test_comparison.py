import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import os

plt.style.use(hep.style.CMS) # <--- 步骤 1: 尝试注释掉这一行

# --- Configuration ---
# Path to the ROOT files after BDT application
BASE_PATH = "/eos/home-j/jiehan/root/outputs/" 

# Samples
SIG_SAMPLES = ["ggH_M125", "VBF_M125"] 
BKG_SAMPLES = ["ZGToLLG", "DYJetsToLL", "EWKZ2J"] 
DATA_SAMPLES = ["Data"] # Assuming Data.root exists, or will be handled if missing

# TTree name in the ROOT files
TREE_NAME = "all_jet" # This was in the prompt's path, assuming it's the TTree name

# BDT score column name
BDT_SCORE_COL = "bdt_score_t" # Or "bdt_score" if not transformed

# Output directory for plots (relative to script location)
OUTPUT_DIR = "figs/bdt_score_train_val_test_comparison"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")

# Mass selections
SIG_MASS_QUERY = "H_mass > 120 & H_mass < 130"
SIDEBAND_QUERY = "(H_mass > 100 & H_mass < 120) | (H_mass > 130 & H_mass < 180)"

# Columns to load from ROOT files
COLUMNS = [BDT_SCORE_COL, "H_mass", "weight", "event"]

def get_data_for_set(sample_names, base_path, tree_name, columns_to_load, set_type, mass_query):
    """
    Loads data for specified samples and applies mass query.
    - For 'train' set_type, reads from 'train1' and 'train2' subdirectories.
    - For 'validation' set_type, reads from 'val' subdirectory.
    - For 'test' set_type, reads from 'test' subdirectory.
    Assumes data in these directories are correctly pre-split and weighted.
    """
    full_df_list = []
    
    # Determine the actual subdirectories to read from based on set_type
    # These are the first level subdirectories under BASE_PATH 
    # (e.g., /eos/home-j/jiehan/root/outputs/train1/)
    sub_dirs_for_set_type = []
    if set_type == "train":
        sub_dirs_for_set_type = ["train1", "train2"]
    elif set_type == "validation":
        sub_dirs_for_set_type = ["val"]
    elif set_type == "test":
        sub_dirs_for_set_type = ["test"] 

    for sub_dir_name in sub_dirs_for_set_type:
        for sample_name in sample_names:
            # Path construction: BASE_PATH / sub_dir_name / tree_name / sample_name.root
            file_path = os.path.join(base_path, sub_dir_name, tree_name, f"{sample_name}.root")
            try:
                tree = uproot.open(file_path)[tree_name]
                df = tree.arrays(columns_to_load, library="pd") # Load all requested columns directly

                if 'weight' not in df.columns: # Simplified weight handling
                     df['weight'] = 1.0

                if mass_query:
                    df = df.query(mass_query).copy() 
                
                full_df_list.append(df)
            except Exception as e:
                print(f"Error loading or processing {sample_name} from {file_path}: {e}") # Keep basic error logging

    if not full_df_list:
        return pd.DataFrame() # Return empty DataFrame if nothing loaded
    
    combined_df = pd.concat(full_df_list, ignore_index=True)
    
    return combined_df

def plot_bdt_comparison(data_dict, title_prefix, filename, bdt_col, bins=50, range_bdt=(0,1), use_log_y=True):
    """
    Plots BDT score distributions for train, validation, and test sets.
    data_dict = {'Train': df_train, 'Validation': df_val, 'Test': df_test}
    use_log_y: If True, use logarithmic scale for y-axis.
    """
    plt.figure(figsize=(10, 10))
    
    for set_name, df in data_dict.items():
        weights = df['weight'] 
        bdt_scores = df[bdt_col] 
        
        current_sum_w = np.sum(weights)
        
        # Calculate histogram density
        hist_density, bins_edges = np.histogram(bdt_scores, bins=bins, range=range_bdt, weights=weights, density=True)
        
        weights_squared = df['weight']**2 
        sum_w2_hist, _ = np.histogram(bdt_scores, bins=bins, range=range_bdt, weights=weights_squared, density=False)
        
        errors_counts = np.sqrt(sum_w2_hist)
        bin_widths = np.diff(bins_edges)
        
        # Handle potential division by zero if current_sum_w is zero or bin_widths contains zero.
        # A robust solution would add checks here, but per request, keeping it minimal.
        # If current_sum_w is 0, errors_normalized_density will be inf or nan.
        errors_normalized_density = np.zeros_like(errors_counts, dtype=float)
        if current_sum_w != 0: # Basic check to prevent division by zero
            errors_normalized_density = errors_counts / current_sum_w
        
        bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2
        
        plt.errorbar(bin_centers, hist_density, yerr=errors_normalized_density, label=set_name, marker='o', linestyle='-', markersize=3, capsize=2, elinewidth=1)

    plt.xlabel(f"BDT Score ({bdt_col})")
    plt.ylabel("Normalized Events")
    plt.title(title_prefix)
    
    # Call plt.legend() to generate/update the legend
    plt.legend(title="Dataset Split")
    if use_log_y:
        plt.yscale('log')  # Use logarithmic scale for y-axis if specified
        
    plt.grid(True, alpha=0.3)
    # plt.ylim(bottom=0) # Removed: Not ideal with log scale. Log scale handles lower bound.
    
    full_plot_path = os.path.join(OUTPUT_DIR, filename)
    print(full_plot_path)
    plt.savefig(full_plot_path)
    # plt.close()
    print(f"Saved plot: {full_plot_path}")

def main():
    # 1. Signal MC in mass window [120, 130]
    print("Processing Signal MC...")
    sig_train_df = get_data_for_set(SIG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "train", SIG_MASS_QUERY)
    sig_val_df   = get_data_for_set(SIG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "validation", SIG_MASS_QUERY)
    sig_test_df  = get_data_for_set(SIG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "test", SIG_MASS_QUERY)
    
    plot_bdt_comparison(
        {'Train': sig_train_df, 'Validation': sig_val_df, 'Test': sig_test_df},
        f"Signal MC (H mass region) BDT Score",
        "signal_mc_bdt_comparison.png",
        BDT_SCORE_COL,
        use_log_y=False  # Do not use log scale for signal
    )

    # 2. Background MC in sidebands [100-120] U [130-180]
    print("\nProcessing Background MC...")
    bkg_train_df = get_data_for_set(BKG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "train", SIDEBAND_QUERY)
    bkg_val_df   = get_data_for_set(BKG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "validation", SIDEBAND_QUERY)
    bkg_test_df  = get_data_for_set(BKG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "test", SIDEBAND_QUERY)

    plot_bdt_comparison(
        {'Train': bkg_train_df, 'Validation': bkg_val_df, 'Test': bkg_test_df},
        f"Background MC (Sidebands) BDT Score",
        "background_mc_bdt_comparison.png",
        BDT_SCORE_COL
        # use_log_y defaults to True, so no need to specify here
    )

    # 3. Data in sidebands [100-120] U [130-180]
    print("\nProcessing Data...")
    data_train_df = get_data_for_set(DATA_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "train", SIDEBAND_QUERY)
    data_val_df   = get_data_for_set(DATA_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "validation", SIDEBAND_QUERY)
    data_test_df  = get_data_for_set(DATA_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "test", SIDEBAND_QUERY)

    plot_bdt_comparison(
        {'Train': data_train_df, 'Validation': data_val_df, 'Test': data_test_df},
        f"Data (Sidebands) BDT Score",
        "data_sidebands_bdt_comparison.png",
        BDT_SCORE_COL
        # use_log_y defaults to True
    )
     # 4. Background MC in mass window [120, 130]
    print("\nProcessing Background MC in mass window...")
    bkg_train_df_mass = get_data_for_set(BKG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "train", SIG_MASS_QUERY)
    bkg_val_df_mass   = get_data_for_set(BKG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "validation", SIG_MASS_QUERY)
    bkg_test_df_mass  = get_data_for_set(BKG_SAMPLES, BASE_PATH, TREE_NAME, COLUMNS, "test", SIG_MASS_QUERY)
    plot_bdt_comparison(
        {'Train': bkg_train_df_mass, 'Validation': bkg_val_df_mass, 'Test': bkg_test_df_mass},
        f"Background MC (H mass region) BDT Score",
        "background_mc_mass_bdt_comparison.png",
        BDT_SCORE_COL
        # use_log_y defaults to True
    )

    print("\nAll plots generated.")

if __name__ == "__main__":
    main()
