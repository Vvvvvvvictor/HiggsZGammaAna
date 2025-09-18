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
from pdb import set_trace

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Batch process DY input ROOT files')
    parser.add_argument('--input-dir', type=str, default='/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_input_redwood/',
                        help='Directory of input files')
    parser.add_argument('--output-dir', type=str, 
                        default='/eos/user/j/jiehan/root/skimmed_ntuples_rui_new/',
                        help='Directory to save output files')
    parser.add_argument('--dry-run', action='store_true',
                        help='If set, only list files without processing them')
    parser.add_argument('--years', type=str, nargs='+', default=["2022", "2022EE", "2023", "2023BPix"],
                        help='List of years to process (default: "2016", "2016APV", "2017", "2018", "2022", "2022EE", "2023", "2023BPix")')
    return parser.parse_args()

def main():
    """Main function"""
    args = parse_arguments()
    
    # Process name mapping - only DY0 related
    proc_mapping = {
        "DY0": "DYJetsToLL",
        "DYfilter": "DYJetsToLL",
        "DYmix": "DYJetsToL"
    }
    
    # Map tree names based on process - only DY related
    tree_mapping = {
        "DY0": "TreeB",
        "DYfilter": "TreeB", 
        "DYmix": "TreeB"
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
        # "pt_mass": "H_relpt",
        # "llphoton_m": "H_mass",
        # "ll_m": "Z_mass",
        # "llphoton_dijet_dphi": "delta_phi_zgjj",
        # "llphoton_pTt": "H_ptt",
        # "costheta": "lep_cos_theta",
        # "phi": "lep_phi",
        # "llphoton_dijet_balance": "pt_balance",
        # "cosTheta": "Z_cos_theta",

        "j1_m": "jet_1_mass",
        "j1_eta": "jet_1_eta",

        "llphoton_refit_pTt": "H_ptt",
        "pt_mass_refit": "H_relpt",
        "llphoton_refit_m": "H_mass",
        "ll_refit_m": "Z_mass",
        "llphoton_refit_dijet_dphi": "delta_phi_zgjj",
        "llphoton_refit_pTt": "H_ptt",
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
    
    # Find all matching files - including special DY files
    input_files = glob.glob(os.path.join(args.input_dir, "*_*_redwood_v1_*.root"))
    logging.info(f"Found {len(input_files)} files")
    
    # 修改正则表达式以匹配DY相关文件和其他格式
    # 匹配格式如: DYfilter_2016_redwood_v1_vbf.root
    pattern = re.compile(r'(DY0|DYfilter|DYmix|[^_]+)_((?:2016APV)|(?:2022EE)|(?:2023BPix)|[0-9]{4})_redwood_v1_(?:(?:all|ggf|vbf)_)?(?:fixedjet_)?(ggf_ext|vbf_ext)\.root')
    
    # Statistics for processing results
    stats = {"success": 0, "skipped": 0, "recovered": 0, "failed": 0}
    
    # Get the list of years to process
    years_to_process = args.years
    logging.info(f"Processing only these years: {years_to_process}")
    
    # Group DY files by year and extension for merging
    dy_files_by_year = {}
    
    for file_path in input_files:
        file_name = os.path.basename(file_path)
        match = pattern.match(file_name)
        
        if not match:
            logging.warning(f"Filename {file_name} does not match expected format, skipping")
            stats["skipped"] += 1
            continue
            
        proc, year_raw, ext = match.groups()
        
        # Skip if the year is not in the list of years to process
        if year_raw not in years_to_process:
            logging.info(f"Skipping {file_name} as year {year_raw} is not in the specified years list")
            stats["skipped"] += 1
            continue
        
        # Group DY files for merging - only process files with ggf_ext or vbf_ext suffix
        if proc in ["DY0", "DYfilter", "DYmix"]:
            year = year_mapping.get(year_raw, year_raw)
            if year not in dy_files_by_year:
                dy_files_by_year[year] = {"ggf_ext": [], "vbf_ext": []}
            
            if ext == "ggf_ext":
                dy_files_by_year[year]["ggf_ext"].append(file_path)
            elif ext == "vbf_ext":
                dy_files_by_year[year]["vbf_ext"].append(file_path)
            # Skip files without special suffix
        else:
            logging.info(f"Skipping non-DY file: {file_name}")
            stats["skipped"] += 1
    
    # Process DY files (merge them)
    for year, file_groups in dy_files_by_year.items():
        output_dir = os.path.join(args.output_dir, "DYJetsToLL_ext")
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{year}.root")
        
        if args.dry_run:
            logging.info(f"Will merge DY files for year {year}:")
            for ext_type, files in file_groups.items():
                if files:
                    logging.info(f"  {ext_type}: {files}")
            continue
        
        try:
            all_data = {"two_jet": [], "zero_to_one_jet": []}
            
            logging.info(f"Processing DY files for year {year}")
            # Process each group of files
            for ext_type, files in file_groups.items():
                logging.info(f"Processing {ext_type} files: {len(files)} files")
                for file_path in files:
                    with uproot.open(file_path) as f:
                        if "TreeB" not in f:
                            logging.error(f"TreeB not found in file {file_path}")
                            continue
                        
                        data = f["TreeB"].arrays()
                        # if "DYmix" in file_path:
                        #     set_trace()
                        
                        # Create a new dict to store renamed branches
                        new_data = {}
                        
                        # Rename branches if needed
                        for old_name in data.fields:
                            if old_name in branch_mapping:
                                new_name = branch_mapping[old_name]
                                new_data[new_name] = data[old_name]
                            else:
                                new_data[old_name] = data[old_name]
                        
                        # Get the integer value for the current year
                        year_int_value = year_to_int_mapping.get(year, -1)
                        new_data["year"] = np.full(len(data), year_int_value, dtype=np.int32)
                        
                        ak_array = ak.Array(new_data)
                        
                        # Check if nel and nmu fields exist
                        if "nel" in ak_array.fields and "nmu" in ak_array.fields:
                            try:
                                nel = ak_array["nel"]
                                nmu = ak_array["nmu"]
                                
                                # Ensure nel and nmu are not None and have proper values
                                if nel is None or nmu is None:
                                    logging.error(f"nel or nmu is None in {file_path}")
                                    continue
                                
                                # Convert any potential None values to 0
                                nel = ak.fill_none(nel, 0) 
                                nmu = ak.fill_none(nmu, 0)
                                
                                # Ensure they are integers
                                nel = ak.values_astype(nel, np.int32)
                                nmu = ak.values_astype(nmu, np.int32)
                                
                                # Apply lepton selection
                                lepton_sum = nel + nmu
                                ak_array = ak_array[lepton_sum == 2]
                                logging.info(f"Processed {file_path}: {len(ak_array)} events after lepton selection")
                                
                            except Exception as e:
                                logging.error(f"Error in lepton selection for {file_path}: {e}")
                                logging.error(f"nel type: {type(nel) if 'nel' in locals() else 'undefined'}")
                                logging.error(f"nmu type: {type(nmu) if 'nmu' in locals() else 'undefined'}")
                                continue
                        else:
                            logging.warning(f"nel or nmu fields not found in {file_path}, skipping lepton selection")
                        
                        # Apply selections based on file type
                        if ext_type == "ggf_ext":
                            # Only zero_to_one_jet selection
                            # Check both original and renamed field names
                            njet_field = None
                            if "n_jets" in ak_array.fields:
                                njet_field = "n_jets"
                            elif "njet" in ak_array.fields:
                                njet_field = "njet"
                            
                            if njet_field:
                                n_jets = ak_array[njet_field]
                                n_jets = ak.where(ak.is_none(n_jets), 0, n_jets)
                                ak_array_filtered = ak_array[n_jets < 2]
                            else:
                                logging.warning(f"Neither n_jets nor njet field found in {file_path}")
                                ak_array_filtered = ak_array
                            
                            if "met" in ak_array_filtered.fields and "zg_categorizationBitMap" in ak_array_filtered.fields:
                                zg_categorization = ak_array_filtered["zg_categorizationBitMap"]
                                zg_categorization = ak.where(ak.is_none(zg_categorization), 0, zg_categorization)
                                # ak_array_filtered = ak_array_filtered[(zg_categorization == 128)]
                                met = ak_array_filtered["met"]
                                met = ak.where(ak.is_none(met), 0, met)
                                ak_array_filtered = ak_array_filtered[met < 90]
                            else:
                                logging.warning(f"met or zg_categorizationBitMap field not found in {file_path}")

                            all_data["zero_to_one_jet"].append(ak_array_filtered)
                        elif ext_type == "vbf_ext":
                            # Only two_jet selection
                            # Check both original and renamed field names
                            njet_field = None
                            if "n_jets" in ak_array.fields:
                                njet_field = "n_jets"
                            elif "njet" in ak_array.fields:
                                njet_field = "njet"
                                
                            if njet_field and "nbdfm" in ak_array.fields and "zg_categorizationBitMap" in ak_array.fields:
                                n_jets = ak_array[njet_field]
                                nbdfm = ak_array["nbdfm"]
                                zg_categorization = ak_array["zg_categorizationBitMap"]
                                n_jets = ak.where(ak.is_none(n_jets), 0, n_jets)
                                nbdfm = ak.where(ak.is_none(nbdfm), 0, nbdfm)
                                # zg_categorization = ak.where(ak.is_none(zg_categorization), 0, zg_categorization)
                                ak_array_filtered = ak_array[(n_jets >= 2) & (nbdfm == 0)]
                            else:
                                logging.warning(f"Required fields not found in {file_path}: njet_field={njet_field}, nbdfm={'nbdfm' in ak_array.fields}")
                                ak_array_filtered = ak_array
                            all_data["two_jet"].append(ak_array_filtered)
            
            # Merge and write output
            # Define expected fields to ensure consistency
            expected_fields = ['weight_corr', 'llphoton_m', 'llphoton_eta', 'llphoton_phi', 'llphoton_pt', 'part', 'll_m', 'gamma_mvaID', 'pt_mass', 'H_relpt', 'gamma_ptRelErr', 'gamma_eta', 'Z_lead_lepton_eta', 'Z_sublead_lepton_eta', 'cosTheta', 'costheta', 'phi', 'l2g_deltaR', 'l1g_deltaR', 'weight', 'weight_dm', 'dm', 'mllg', 'gamma_pt', 'photon_eta', 'photon_phi', 'photon_reliso', 'jet1G_deltaR', 'jet2G_deltaR', 'photon_jet_mindr', 'photon_jet_maxdr', 'photon_zeppenfeld', 'gamma_relpt', 'll_lepid', 'll_refit_pt', 'Z_mass', 'll_refit_l1_pt', 'll_refit_l2_pt', 'll_refit_eta', 'll_refit_phi', 'llphoton_refit_pt', 'llphoton_refit_eta', 'llphoton_refit_phi', 'H_mass', 'llphoton_refit_dphi', 'llphoton_refit_deta', 'llphoton_refit_dr', 'llphoton_refit_l1_masserr', 'llphoton_refit_l2_masserr', 'llphoton_refit_ph_masserr', 'lep_cos_theta', 'Z_cos_theta', 'lep_phi', 'delta_phi_zgjj', 'pt_balance', 'llphoton_refit_dijet_dr', 'H_ptt', 'llphoton_refit_pTt_an_hig019014', 'llphoton_dijet_dphi', 'llphoton_dijet_balance', 'llphoton_dijet_dr', 'llphoton_hmiss_photon_dphi', 'llphoton_pTt', 'llphoton_pTt_an_hig019014', 'event', 'n_jets', 'nbdfm', 'jet_1_pt', 'jet_1_eta', 'j1_phi', 'jet_1_mass', 'j1_isgood', 'jet_2_pt', 'jet_2_eta', 'j2_phi', 'j2_isgood', 'j2_m', 'j3_pt', 'j3_eta', 'j3_phi', 'j3_isgood', 'j3_m', 'dijet_pt', 'dijet_eta', 'dijet_phi', 'mass_jj', 'dijet_dr', 'delta_phi_jj', 'delta_eta_jj', 'j1_qgl', 'j1_puid_disc', 'j1_ne_emef', 'j2_qgl', 'j2_puid_disc', 'j2_ne_emef', 'j3_qgl', 'j3_puid_disc', 'j3_zeppenfeld', 'zg_cutBitMap', 'met', 'isrun3', 'nel', 'nmu', 'year']
            with uproot.recreate(output_path) as out:
                for tree_name, arrays in all_data.items():
                    # if "zero_to_one_jet" in tree_name:
                    #     continue
                    try:
                        logging.info(f"Merging {len(arrays)} arrays for {tree_name}")
                        
                        # Ensure all arrays have the same fields by keeping only common fields
                        if len(arrays) > 1:
                            # Find intersection of all fields
                            common_fields = set(arrays[0].fields)
                            for arr in arrays[0:]:
                                common_fields &= set(arr.fields)
                            
                            logging.info(f"Common fields: {len(common_fields)} out of expected {len(expected_fields)}")
                            
                            # Create aligned arrays with common fields in consistent order
                            aligned_arrays = []
                            sorted_common_fields = sorted(common_fields)
                            
                            for i, arr in enumerate(arrays):
                                try:
                                    # Create new array with only common fields
                                    field_dict = {}
                                    for field in sorted_common_fields:
                                        field_dict[field] = arr[field]
                                    aligned_arrays.append(ak.Array(field_dict))
                                except Exception as e:
                                    logging.error(f"Error aligning array {i} for {tree_name}: {e}")
                                    continue
                            
                            if aligned_arrays:
                                merged_array = ak.concatenate(aligned_arrays)
                                logging.info(f"DY {year} {tree_name}: {len(merged_array)} events")
                                out[tree_name] = merged_array
                            else:
                                logging.error(f"No arrays could be aligned for {tree_name}")
                                out[tree_name] = ak.Array({})
                        else:
                            # Single array, use directly
                            merged_array = arrays[0]
                            logging.info(f"DY {year} {tree_name}: {len(merged_array)} events")
                            out[tree_name] = merged_array
                            
                    except Exception as e:
                        logging.error(f"Error merging arrays for {tree_name}: {e}")
                        import traceback
                        traceback.print_exc()
                        # Create empty array as fallback
                        out[tree_name] = ak.Array({})
                
                # # Also create all_jet tree by combining two_jet and zero_to_one_jet  
                # try:
                #     if all_data["two_jet"] or all_data["zero_to_one_jet"]:
                #         combined_arrays = []
                #         if all_data["two_jet"]:
                #             combined_arrays.extend(all_data["two_jet"])
                #         if all_data["zero_to_one_jet"]:
                #             combined_arrays.extend(all_data["zero_to_one_jet"])
                        
                #         # Use the same alignment logic
                #         common_fields = set(combined_arrays[0].fields)
                #         for arr in combined_arrays[1:]:
                #             common_fields &= set(arr.fields)
                        
                #         aligned_arrays = []
                #         sorted_common_fields = sorted(common_fields)
                        
                #         for arr in combined_arrays:
                #             field_dict = {field: arr[field] for field in sorted_common_fields}
                #             aligned_arrays.append(ak.Array(field_dict))
                        
                #         merged_all = ak.concatenate(aligned_arrays)
                #         logging.info(f"DY {year} all_jet: {len(merged_all)} events")
                #         out["all_jet"] = merged_all
                # except Exception as e:
                #     logging.error(f"Error creating all_jet tree: {e}")
                #     import traceback
                #     traceback.print_exc()
            
            logging.info(f"Successfully merged DY files for year {year} -> {output_path}")
            stats["success"] += 1
            
        except Exception as e:
            logging.error(f"Error merging DY files for year {year}: {str(e)}")
            stats["failed"] += 1
    
    # Output processing statistics
    logging.info(f"Processing completed: Success {stats['success']}, Recovered {stats['recovered']}, Skipped {stats['skipped']}, Failed {stats['failed']}")

if __name__ == "__main__":
    main()