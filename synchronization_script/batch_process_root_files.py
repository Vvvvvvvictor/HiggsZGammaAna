#!/usr/bin/env python3
import os
import glob
import re
import argparse
import uproot
import numpy as np
import pandas as pd
import awkward as ak
from tqdm import tqdm

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Batch process VBF input ROOT files')
    parser.add_argument('--input-dir', type=str, default='/eos/user/j/jiehan/root/VBF_input/',
                        help='Directory of input files')
    parser.add_argument('--output-dir', type=str, 
                        default='/eos/user/j/jiehan/root/skimmed_ntuples_rui/',
                        help='Directory to save output files')
    parser.add_argument('--dry-run', action='store_true',
                        help='If set, only list files without processing them')
    return parser.parse_args()

def main():
    """Main function"""
    args = parse_arguments()
    
    # Process name mapping
    proc_mapping = {
        "DY0": "DYJetsToLL",
        "EWK": "EWKZ2J",
        "SM1": "ZGToLLG",
        "GGF": "ggH_M125",
        "VBF": "VBF_M125"
    }
    
    # Map tree names based on process
    tree_mapping = {
        "DY0": "TreeB",
        "EWK": "TreeB",
        "SM1": "TreeB",
        "GGF": "TreeS",
        "VBF": "TreeS"
    }
    
    # Branch renaming mapping
    branch_mapping = {
        "llphoton_m": "H_mass",
        "ll_m": "Z_mass",
        "photon_pt_ratio": "gamma_relpt",
        "dijet_m": "mass_jj",
        "dijet_deta": "delta_eta_jj",
        "dijet_dphi": "delta_phi_jj",
        "llphoton_dijet_dphi": "delta_phi_zgjj",
        "photon_rapidity": "gamma_eta",
        "photon_mva": "gamma_mvaID",
        "llphoton_pTt": "H_ptt",
        "photon_pt": "gamma_pt",
        "j1_pt": "jet_1_pt",
        "j2_pt": "jet_2_pt",
        "j1_eta": "jet_1_eta",
        "j2_eta": "jet_2_eta",
        "photon_jet1_dr": "jet1G_deltaR",
        "photon_jet2_dr": "jet2G_deltaR",
        "max_dR": "l1g_deltaR",
        "min_dR": "l2g_deltaR",
        "costheta": "lep_cos_theta",
        "phi": "lep_phi",
        "photon_zeppenfeld": "photon_zeppenfeld",
        "llphoton_dijet_balance": "pt_balance",
        "cosTheta": "Z_cos_theta",
        "l1_rapidity": "Z_lead_lepton_eta",
        "l2_rapidity": "Z_sublead_lepton_eta",
        "pt_mass": "H_relpt",
        "photon_res": "gamma_ptRelErr",
        "njet": "n_jets"
    }
    
    # Find all matching files
    input_files = glob.glob(os.path.join(args.input_dir, "*_*_pinnacles_fixed.root"))
    print(f"Found {len(input_files)} files")
    
    # Use regex to match filename pattern with year captured as either "2016APV" or a 4-digit year
    pattern = re.compile(r'([^_]+)_((?:2016APV)|[0-9]{4})_pinnacles_fixed\.root')
    
    for file_path in tqdm(input_files):
        file_name = os.path.basename(file_path)
        match = pattern.match(file_name)
        
        if not match:
            print(f"Warning: Filename {file_name} does not match expected format, skipping")
            continue
            
        proc, year_raw = match.groups()
        # Map year: "2016" -> "2016postVFP", "2016APV" -> "2016preVFP"
        if year_raw == "2016":
            year = "2016postVFP"
        elif year_raw == "2016APV":
            year = "2016preVFP"
        else:
            year = year_raw
        
        # Check if the process is in our mapping
        if proc not in proc_mapping:
            print(f"Warning: Process {proc} is not in the mapping list, skipping")
            continue
            
        # Get the new process name and tree name
        new_proc = proc_mapping[proc]
        tree_name = tree_mapping[proc]
        
        # Create output directory
        output_dir = os.path.join(args.output_dir, new_proc)
        os.makedirs(output_dir, exist_ok=True)
        
        # Output file path using the transformed year
        output_path = os.path.join(output_dir, f"{year}.root")
        
        if args.dry_run:
            print(f"Will process: {file_path} -> {output_path} (Tree: {tree_name} -> two_jet)")
            continue
            
        try:
            # Open input file
            with uproot.open(file_path) as f:
                if tree_name not in f:
                    print(f"Error: Tree {tree_name} not found in file {file_name}")
                    continue
                    
                # Read data
                data = f[tree_name].arrays()
                
                # Create a new dict to store renamed branches
                new_data = {}
                
                # Rename branches if needed
                for old_name in data.fields:
                    if old_name in branch_mapping:
                        new_name = branch_mapping[old_name]
                        new_data[new_name] = data[old_name]
                    else:
                        # Retain branches that do not require renaming
                        new_data[old_name] = data[old_name]
                
                # Convert data to an Array (or DataFrame)
                ak_array = ak.Array(new_data)
                
                # Write to new file with tree named "two_jet"
                with uproot.recreate(output_path) as out:
                    out["two_jet"] = ak_array
                
                print(f"Successfully processed: {file_path} -> {output_path}")
                
        except Exception as e:
            print(f"Error processing file {file_name}: {str(e)}")

if __name__ == "__main__":
    main()
