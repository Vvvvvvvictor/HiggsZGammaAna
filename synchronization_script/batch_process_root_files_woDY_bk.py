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
    parser = argparse.ArgumentParser(description='Batch process VBF input ROOT files')
    parser.add_argument('--input-dir', type=str, default='/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_input_redwood/',
                        help='Directory of input files')
    parser.add_argument('--output-dir', type=str, 
                        default='/eos/user/j/jiehan/root/skimmed_ntuples_rui_new/',
                        help='Directory to save output files')
    parser.add_argument('--dry-run', action='store_true',
                        help='If set, only list files without processing them')
    parser.add_argument('--recover', action='store_true',
                        help='Attempt to recover corrupted files')
    parser.add_argument('--years', type=str, nargs='+', default=["2016", "2016APV", "2017", "2018", "2022", "2022EE", "2023", "2023BPix"],
                        help='List of years to process (default: "2016", "2016APV", "2017", "2018", "2022", "2022EE", "2023", "2023BPix")')
    return parser.parse_args()

def main():
    """Main function"""
    args = parse_arguments()
    
    # Process name mapping
    proc_mapping = {
        "DY0": "DYJetsToLL",
        "EWK": "EWKZ2J",
        "SM1": "ZGToLLG1",
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

    # Define year to integer mapping
    year_to_int_mapping = {
        "2016preVFP": 0,
        "2016postVFP": 1,
        "2017": 2,
        "2018": 3,
        "2022preEE": 4,
        "2022postEE": 5,
        "2023preBPix": 6,
        "2023postBPix": 7
    }
    
    # Branch renaming mapping
    branch_mapping = {
        "llphoton_pTt": "H_ptt",
        "pt_mass": "H_relpt",
        "llphoton_m": "H_mass",
        "ll_m": "Z_mass",
        "llphoton_dijet_dphi": "delta_phi_zgjj",
        "costheta": "lep_cos_theta",
        "phi": "lep_phi",
        "llphoton_dijet_balance": "pt_balance",
        "cosTheta": "Z_cos_theta",

        "j1_m": "jet_1_mass",
        "j1_eta": "jet_1_eta",

        # "llphoton_refit_pTt": "H_ptt",
        # "pt_mass_refit": "H_relpt",
        # "llphoton_refit_m": "H_mass",
        # "ll_refit_m": "Z_mass",
        # "llphoton_refit_dijet_dphi": "delta_phi_zgjj",
        # "llphoton_refit_pTt": "H_ptt",
        # "llphoton_refit_costheta": "lep_cos_theta",
        # "llphoton_refit_psi": "lep_phi",
        # "llphoton_refit_dijet_balance": "pt_balance",
        # "llphoton_refit_cosTheta": "Z_cos_theta",

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
    
    # Find all matching files
    input_files = glob.glob(os.path.join(args.input_dir, "*_*_redwood_v1_*.root"))
    logging.info(f"Found {len(input_files)} files")
    
    # 修改正则表达式以匹配更多年份格式，包括2022EE和2023BPix等特殊年份
    pattern = re.compile(r'([^_]+)_((?:2016APV)|(?:2022EE)|(?:2023BPix)|[0-9]{4})_redwood_v1_(ggf|vbf)\.root')
    
    # Statistics for processing results
    stats = {"success": 0, "skipped": 0, "recovered": 0, "failed": 0}
    
    # Get the list of years to process
    years_to_process = args.years
    logging.info(f"Processing only these years: {years_to_process}")

    # Group files by process and year
    files_to_process = {}
    
    for file_path in tqdm(input_files):
        file_name = os.path.basename(file_path)
        match = pattern.match(file_name)
        
        if not match:
            # logging.warning(f"Filename {file_name} does not match expected format, skipping")
            stats["skipped"] += 1
            continue
            
        proc, year_raw, ext = match.groups()
        
        # Skip if the year is not in the list of years to process
        if year_raw not in years_to_process:
            # logging.info(f"Skipping {file_name} as year {year_raw} is not in the specified years list")
            stats["skipped"] += 1
            continue
            
        # Map year using the year mapping dictionary
        year = year_mapping.get(year_raw, year_raw)  # Fall back to year_raw if not in mapping
        
        # Check if the process is in our mapping
        if proc not in proc_mapping:
            # logging.warning(f"Process {proc} is not in the mapping list, skipping")
            stats["skipped"] += 1
            continue

        if proc not in files_to_process:
            files_to_process[proc] = {}
        if year not in files_to_process[proc]:
            files_to_process[proc][year] = []
        
        files_to_process[proc][year].append(file_path)

    for proc, year_files in files_to_process.items():
        for year, file_list in year_files.items():
            # Get the new process name and tree name
            new_proc = proc_mapping[proc]
            tree_name = tree_mapping[proc]
            
            # Create output directory
            output_dir = os.path.join(args.output_dir, new_proc)
            os.makedirs(output_dir, exist_ok=True)
            
            # Output file path using the transformed year
            output_path = os.path.join(output_dir, f"{year}.root")
            
            if args.dry_run:
                logging.info(f"Will process and merge for {new_proc}/{year}: {file_list} -> {output_path}")
                continue

            all_data = {"two_jet": [], "zero_to_one_jet": []}
            
            try:
                for file_path in file_list:
                    file_name = os.path.basename(file_path)
                    match = pattern.match(file_name)
                    _, _, ext = match.groups()

                    # Try to open file with recovery options
                    options = {}
                    if args.recover:
                        options = {"handler": lambda x: logging.warning(f"Recovered from error: {x}")}
                    
                    with uproot.open(file_path, **options) as f:
                        if tree_name not in f:
                            logging.error(f"Tree {tree_name} not found in file {file_name}")
                            continue
                        
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
                        
                        # Get the integer value for the current year
                        year_int_value = year_to_int_mapping.get(year, -1) # Default to -1 if year not in map
                        new_data["year"] = np.full(len(data), year_int_value, dtype=np.int32)

                        # Convert data to an Array
                        ak_array = ak.Array(new_data)

                        ak_array = ak_array[ak_array["nel"] + ak_array["nmu"] == 2]

                        # Ensure zg_categorizationBitMap is integer for bitwise operations
                        if "zg_categorizationBitMap" in ak_array.fields:
                            ak_array["zg_categorizationBitMap"] = ak.values_astype(ak_array["zg_categorizationBitMap"], np.int64)

                        if ext == "vbf":
                            ak_array_vbf = ak_array[(ak_array["n_jets"] >= 2) & (ak_array["nbdfm"] == 0) & (ak_array["zg_categorizationBitMap"] & 64 == 64)]
                            all_data["two_jet"].append(ak_array_vbf)
                        elif ext == "ggf":
                            ak_array_ggf = ak_array[(ak_array["n_jets"] < 2) & (ak_array["met"] < 90) & (ak_array["zg_categorizationBitMap"] & 128 == 128)]
                            all_data["zero_to_one_jet"].append(ak_array_ggf)

                with uproot.recreate(output_path) as out:
                    # Merge and write two_jet
                    if all_data["two_jet"]:
                        merged_two_jet = ak.concatenate(all_data["two_jet"])
                        out["two_jet"] = merged_two_jet
                        print(f"{new_proc} {year} two_jet event count: {len(merged_two_jet)}")
                    
                    # Merge and write zero_to_one_jet
                    if all_data["zero_to_one_jet"]:
                        merged_zero_to_one_jet = ak.concatenate(all_data["zero_to_one_jet"])
                        out["zero_to_one_jet"] = merged_zero_to_one_jet
                        print(f"{new_proc} {year} zero_to_one_jet event count: {len(merged_zero_to_one_jet)}")

                    # Merge and write all_jet
                    all_jet_arrays = []
                    if all_data["two_jet"]:
                        all_jet_arrays.extend(all_data["two_jet"])
                    if all_data["zero_to_one_jet"]:
                        all_jet_arrays.extend(all_data["zero_to_one_jet"])
                    
                    if all_jet_arrays:
                        # Align fields before concatenating
                        if len(all_jet_arrays) > 1:
                            common_fields = set(all_jet_arrays[0].fields)
                            for arr in all_jet_arrays[1:]:
                                common_fields &= set(arr.fields)
                            
                            aligned_arrays = []
                            sorted_common_fields = sorted(list(common_fields))
                            for arr in all_jet_arrays:
                                field_dict = {field: arr[field] for field in sorted_common_fields}
                                aligned_arrays.append(ak.Array(field_dict))
                            merged_all_jet = ak.concatenate(aligned_arrays)
                        else:
                            merged_all_jet = all_jet_arrays[0]

                        out["all_jet"] = merged_all_jet
                        print(f"{new_proc} {year} all_jet event count: {len(merged_all_jet)}")

                logging.info(f"Successfully processed and merged for {new_proc}/{year} -> {output_path}")
                stats["success"] += 1

            except Exception as e:
                logging.error(f"Unexpected error processing {new_proc}/{year}: {str(e)}")
                stats["failed"] += 1
    
    # Output processing statistics
    logging.info(f"Processing completed: Success {stats['success']}, Recovered {stats['recovered']}, Skipped {stats['skipped']}, Failed {stats['failed']}")

if __name__ == "__main__":
    main()

# rm ZGToLLG/2022preEE.root; hadd ZGToLLG/2022preEE.root ZGToLLG1/2022preEE.root ZGToLLG2/2022preEE.root; rm ZGToLLG/2022postEE.root; hadd ZGToLLG/2022postEE.root ZGToLLG1/2022postEE.root ZGToLLG2/2022postEE.root