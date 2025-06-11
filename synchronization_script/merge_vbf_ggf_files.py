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
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Merge VBF and GGF input ROOT files')
    parser.add_argument('--vbf-input-dir', type=str, default='/eos/project/h/htozg-dy-privatemc/rzou/bdt/VBF_input/',
                        help='Directory of VBF input files')
    parser.add_argument('--ggf-input-dir', type=str, default='/eos/project/h/htozg-dy-privatemc/rzou/bdt/GGF_input/',
                        help='Directory of GGF input files')
    parser.add_argument('--output-dir', type=str, 
                        default='/eos/user/j/jiehan/root/merged_ntuples_vbf_ggf/',
                        help='Directory to save output files')
    parser.add_argument('--dry-run', action='store_true',
                        help='If set, only list files without processing them')
    parser.add_argument('--recover', action='store_true',
                        help='Attempt to recover corrupted files')
    parser.add_argument('--years', type=str, nargs='+', default=["2016", "2016APV", "2017", "2018"],
                        help='List of years to process')
    return parser.parse_args()

def find_matching_files(vbf_dir, ggf_dir, years):
    """Find matching files in VBF and GGF directories based on process and year"""
    # Dictionary to store file pairs (vbf, ggf)
    file_pairs = {}
    
    # Pattern to match VBF files
    vbf_pattern = re.compile(r'([^_]+)_((?:2016APV)|(?:2022EE)|(?:2023BPix)|[0-9]{4})_pinnacles_fixed\.root')
    
    # Pattern to match GGF files
    ggf_pattern = re.compile(r'([^_]+)_((?:2016APV)|(?:2022EE)|(?:2023BPix)|[0-9]{4})_pinnacles_ggf_fixed\.root')
    
    # Find all VBF files
    vbf_files = glob.glob(os.path.join(vbf_dir, "*_*_pinnacles_fixed.root"))
    logging.info(f"Found {len(vbf_files)} VBF files")
    
    # Find all GGF files
    ggf_files = glob.glob(os.path.join(ggf_dir, "*_*_pinnacles_ggf_fixed.root"))
    logging.info(f"Found {len(ggf_files)} GGF files")
    
    # Process VBF files
    for vbf_file in vbf_files:
        vbf_name = os.path.basename(vbf_file)
        match = vbf_pattern.match(vbf_name)
        
        if match:
            proc, year = match.groups()
            
            # Skip if year not in specified years
            if year not in years:
                continue
                
            key = (proc, year)
            if key not in file_pairs:
                file_pairs[key] = {"vbf": None, "ggf": None}
                
            file_pairs[key]["vbf"] = vbf_file
    
    # Process GGF files
    for ggf_file in ggf_files:
        ggf_name = os.path.basename(ggf_file)
        match = ggf_pattern.match(ggf_name)
        
        if match:
            proc, year = match.groups()
            
            # Skip if year not in specified years
            if year not in years:
                continue
                
            key = (proc, year)
            if key not in file_pairs:
                file_pairs[key] = {"vbf": None, "ggf": None}
                
            file_pairs[key]["ggf"] = ggf_file
    
    return file_pairs

def main():
    """Main function"""
    args = parse_arguments()
    
    # Process name mapping
    proc_mapping = {
        "DY0": "DYJetsToLL",
        "EWK": "EWKZ2J",
        "SM1": "ZGToLLG",
        "SM2": "ZGToLLG2",
        "data": "Data",
        "GGF": "ggH_M125",
        "VBF": "VBF_M125",
        "WH": "WH_M125",
        "ZH": "ZH_M125",
        "ttH": "ttH_M125"
    }
    
    # Map tree names based on process
    tree_mapping = {
        "DY0": "TreeB",
        "EWK": "TreeB",
        "SM1": "TreeB",
        "SM2": "TreeB",
        "data": "TreeB",
        "GGF": "TreeS",
        "VBF": "TreeS",
        "WH": "TreeS",
        "ZH": "TreeS",
        "ttH": "TreeS"
    }
    
    # Year mapping for output file naming
    year_mapping = {
        "2016": "2016postVFP",
        "2016APV": "2016preVFP",
        "2022": "2022preEE",
        "2022EE": "2022postEE", 
        "2023": "2023preBPix",
        "2023BPix": "2023postBPix"
    }
    
    # Branch renaming mapping
    branch_mapping = {
        "llphoton_refit_pTt": "H_ptt",
        "pt_mass_refit": "H_relpt",
        "llphoton_refit_m": "H_mass",
        "ll_refit_m": "Z_mass",
        "llphoton_refit_dijet_dphi": "delta_phi_zgjj",
        "llphoton_refit_costheta": "lep_cos_theta",
        "llphoton_refit_psi": "lep_phi",
        "llphoton_refit_dijet_balance": "pt_balance",
        "llphoton_refit_cosTheta": "Z_cos_theta",

        "photon_pt_ratio": "gamma_relpt",
        "dijet_m": "mass_jj",
        "dijet_deta": "delta_eta_jj",
        "dijet_dphi": "delta_phi_jj",
        "photon_rapidity": "gamma_eta",
        "photon_mva": "gamma_mvaID",
        "photon_pt": "gamma_pt",
        "j1_pt": "jet_1_pt",
        "j2_pt": "jet_2_pt",
        "j1_eta": "jet_1_eta",
        "j2_eta": "jet_2_eta",
        "photon_jet1_dr": "jet1G_deltaR",
        "photon_jet2_dr": "jet2G_deltaR",
        "max_dR": "l1g_deltaR",
        "min_dR": "l2g_deltaR",
        "photon_zeppenfeld": "photon_zeppenfeld",
        "l1_rapidity": "Z_lead_lepton_eta",
        "l2_rapidity": "Z_sublead_lepton_eta",
        "photon_res": "gamma_ptRelErr",
        "njet": "n_jets"
    }
    
    # Statistics for processing results
    stats = {"success": 0, "skipped": 0, "recovered": 0, "failed": 0}
    
    # Find matching files from VBF and GGF directories
    file_pairs = find_matching_files(args.vbf_input_dir, args.ggf_input_dir, args.years)
    logging.info(f"Found {len(file_pairs)} matching file pairs")
    
    # Process each pair of files
    for (proc, year), files in tqdm(file_pairs.items()):
        vbf_file = files["vbf"]
        ggf_file = files["ggf"]
        
        # Skip if either file is missing
        if vbf_file is None or ggf_file is None:
            logging.warning(f"Missing file for {proc}_{year}: VBF={vbf_file}, GGF={ggf_file}")
            stats["skipped"] += 1
            continue
        
        # Check if the process is in our mapping
        if proc not in proc_mapping:
            logging.warning(f"Process {proc} is not in the mapping list, skipping")
            stats["skipped"] += 1
            continue
        
        # Get the new process name and tree name
        new_proc = proc_mapping[proc]
        tree_name = tree_mapping[proc]
        
        # Map year using the year mapping dictionary
        out_year = year_mapping.get(year, year)  # Fall back to year if not in mapping
        
        # Create output directory
        output_dir = os.path.join(args.output_dir, new_proc)
        os.makedirs(output_dir, exist_ok=True)
        
        # Output file path
        output_path = os.path.join(output_dir, f"{out_year}.root")
        
        if args.dry_run:
            logging.info(f"Will merge: VBF={vbf_file}, GGF={ggf_file} -> {output_path}")
            continue
        
        try:
            # Options for file recovery if needed
            options = {}
            if args.recover:
                options = {"handler": lambda x: logging.warning(f"Recovered from error: {x}")}
            
            # Load VBF data
            vbf_data = None
            with uproot.open(vbf_file, **options) as f:
                if tree_name not in f:
                    logging.error(f"Tree {tree_name} not found in VBF file {vbf_file}")
                    stats["failed"] += 1
                    continue
                
                vbf_data = f[tree_name].arrays()
                logging.info(f"VBF data loaded: {len(vbf_data)} events")
            
            # Load GGF data
            ggf_data = None
            with uproot.open(ggf_file, **options) as f:
                if tree_name not in f:
                    logging.error(f"Tree {tree_name} not found in GGF file {ggf_file}")
                    stats["failed"] += 1
                    continue
                
                ggf_data = f[tree_name].arrays()
                logging.info(f"GGF data loaded: {len(ggf_data)} events")
            
            # Add a source flag to each dataset
            vbf_data_with_flag = vbf_data.copy()
            vbf_data_with_flag["source"] = ak.full_like(ak.zeros_like(vbf_data["njet"]), 0)  # 0 for VBF
            
            ggf_data_with_flag = ggf_data.copy()
            ggf_data_with_flag["source"] = ak.full_like(ak.zeros_like(ggf_data["njet"]), 1)  # 1 for GGF
            
            # Merge the datasets
            # First check if all fields match between datasets
            vbf_fields = set(vbf_data.fields)
            ggf_fields = set(ggf_data.fields)
            
            if vbf_fields != ggf_fields:
                logging.warning(f"Fields mismatch between VBF and GGF files for {proc}_{year}")
                logging.warning(f"VBF only: {vbf_fields - ggf_fields}")
                logging.warning(f"GGF only: {ggf_fields - vbf_fields}")
                
                # Use intersection of fields
                common_fields = vbf_fields.intersection(ggf_fields)
                common_fields.add("source")  # Add our flag field
                
                # Filter the datasets to only include common fields
                vbf_filtered = {field: vbf_data_with_flag[field] for field in common_fields if field in vbf_data_with_flag}
                ggf_filtered = {field: ggf_data_with_flag[field] for field in common_fields if field in ggf_data_with_flag}
                
                vbf_data_with_flag = ak.Array(vbf_filtered)
                ggf_data_with_flag = ak.Array(ggf_filtered)
            
            # Concatenate the arrays
            merged_data = ak.concatenate([vbf_data_with_flag, ggf_data_with_flag])
            logging.info(f"Merged data: {len(merged_data)} events")
            
            # Create a new dict for renamed branches
            new_data = {}
            
            # Rename branches if needed
            for old_name in merged_data.fields:
                if old_name in branch_mapping:
                    new_name = branch_mapping[old_name]
                    new_data[new_name] = merged_data[old_name]
                else:
                    # Retain branches that do not require renaming
                    new_data[old_name] = merged_data[old_name]
            
            # Convert to awkward array
            final_array = ak.Array(new_data)
            print(f"Final array event count before VBF selection: {len(final_array)}")
            
            # Add VBF channel selection if 'n_jets' field exists
            if "n_jets" in new_data:
                print(f"Final array event count after njets selection: {len(final_array[final_array['n_jets'] >= 2])}")
                final_array = final_array[(final_array["n_jets"] >= 2) & (final_array["nbdfm"] == 0)]
                print(f"Final array event count after VBF selection: {len(final_array)}")
            
            # Write to new file with tree named "two_jet"
            with uproot.recreate(output_path) as out:
                out["two_jet"] = final_array
            
            logging.info(f"Successfully processed and merged: {output_path}")
            stats["success"] += 1
            
        except Exception as e:
            logging.error(f"Error processing files for {proc}_{year}: {str(e)}")
            stats["failed"] += 1
    
    # Output processing statistics
    logging.info(f"Processing completed: Success {stats['success']}, Recovered {stats['recovered']}, Skipped {stats['skipped']}, Failed {stats['failed']}")

if __name__ == "__main__":
    main()
