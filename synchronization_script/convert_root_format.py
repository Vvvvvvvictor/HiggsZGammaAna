import os
import json
import uproot
import pandas as pd
import numpy as np
import awkward as ak
from argparse import ArgumentParser
import logging

from pdb import set_trace

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=logging.DEBUG)

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Convert ROOT files to specific format')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-j/jiehan/root/VBF_output/', help='Path to the input folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/home-j/jiehan/root', help='Path to the output folder')
    parser.add_argument('-t', '--type', default='data', choices=['data', 'signal'], help='Type of samples to process (data or signal)')
    return parser.parse_args()

def map_year(year):
    """Map year names to their corresponding format in the input files."""
    year_mapping = {
        '2016postVFP': '2016',
        '2016preVFP': '2016APV',
        '2022preEE': '2022',
        '2022postEE': '2022EE',
        '2023preBPix': '2023',
        '2023postBPix': '2023BPix'
    }
    return year_mapping.get(year, year)

def map_sample_name(sample):
    """Map sample names to their corresponding format in the input files."""
    sample_mapping = {
        'ggH_M125': 'GGF',
        'VBF_M125': 'VBF',
        'WH_M125': 'WH',
        'ZH_M125': 'ZH',
        'ttH_M125': 'ttH',
        'Data': 'data'
    }
    return sample_mapping.get(sample, sample)

class DataConverter:
    def __init__(self):
        """Initialize the Data Converter with category boundaries."""
        self.cat_boundaries = [-float('inf'), 0.083, 0.402, 0.518, float('inf')]
        self.bdt_branch = 'BDT_score_2j'
        self.lepton_branch = 'll_lepid'
        
    def split_by_categories(self, data):
        """Split the data based on category boundaries using BDT_score_2j."""
        categories = []
        lepton_cats = {}
        
        # First split by BDT score
        for i in range(len(self.cat_boundaries) - 1):
            lower = self.cat_boundaries[i]
            upper = self.cat_boundaries[i+1]
            
            # Handle -inf and inf boundaries
            if lower == -float('inf'):
                if upper == float('inf'):
                    cat_data = data
                else:
                    cat_data = data[data[self.bdt_branch] < upper]
            elif upper == float('inf'):
                cat_data = data[data[self.bdt_branch] >= lower]
            else:
                cat_data = data[(data[self.bdt_branch] >= lower) & (data[self.bdt_branch] < upper)]
            
            # Further split by lepton ID
            ele_data = cat_data[cat_data[self.lepton_branch] == 11]
            mu_data = cat_data[cat_data[self.lepton_branch] == 13]
            
            lepton_cats[f"{i}_ele"] = ele_data
            lepton_cats[f"{i}_mu"] = mu_data
            
            # Keep the original category split for backward compatibility
            categories.append(cat_data)
        
        return categories, lepton_cats

def process_data_files(converter, output_folder, input_folder):
    """Process data files and write the results to output ROOT files."""
    procductions = ['Data']
    years = ['2016preVFP', '2016postVFP', '2017', '2018', '2022preEE', '2022postEE', '2023preBPix', '2023postBPix']
    region = 'VBF'  # Using VBF as default region

    # Create all necessary directories in one go
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all process-year combinations and write output
    for proc, year in [(p, y) for p in procductions for y in years]:
        if not os.path.exists(f"{output_folder}/fitting_bkg/{proc}_{year}"): 
            os.makedirs(f"{output_folder}/fitting_bkg/{proc}_{year}", exist_ok=True)
        
        with uproot.recreate(f"{output_folder}/fitting_bkg/{proc}_{year}/output_{proc}_{year}.root") as outfile:
            logging.info(f"Output file: {output_folder}/fitting_bkg/{proc}_{year}/output_{proc}_{year}.root")
            logging.info(f"Processing {proc} {year}")
            
            # Map the year to the corresponding format in input files
            input_year = map_year(year)
            input_sample = map_sample_name(proc)
            input_path = f'{input_folder}/{input_sample}_{input_year}_output.root'
            
            with uproot.open(input_path) as infile:
                # Read data using awkward array to avoid issues during pandas conversion
                ak_data = infile["outtree"].arrays(library='ak')
                try:
                    # Try converting to pandas DataFrame using the correct ak.to_pandas() function
                    data = ak.to_pandas(ak_data)
                except Exception as e:
                    logging.warning(f"Error converting awkward array to pandas directly: {e}")
                    # If direct conversion fails, try manual conversion
                    columns = ak_data.fields
                    data = pd.DataFrame({col: ak_data[col].to_numpy() for col in columns})
                
                # Rename llphoton_refit_m to CMS_hgg_mass
                data = data.rename(columns={"llphoton_refit_m": "CMS_hgg_mass"})
                
                # Add new dZ branch filled with 0.0
                data['dZ'] = 0.0
                data = data.query('(CMS_hgg_mass > 95) & (CMS_hgg_mass < 170)')
                
                # Split data by BDT categories and lepton types
                data_splits, lepton_splits = converter.split_by_categories(data)
                
                # Write categories to output file - both original split and lepton split
                for i, split_data in enumerate(data_splits):
                    if region == 'VBF':
                        if i == 0:
                            split_data = split_data.query('CMS_hgg_mass > 100')
                        elif i == 1 or i == 2 or i == 3:
                            split_data = split_data.query('CMS_hgg_mass > 95')
                    logging.info(f"Number of events in {region}{i}: {len(split_data)}")
                    outfile[f'DiphotonTree/Data_13TeV_{region}{i}'] = split_data
    
    # Create combined output for all years
    if not os.path.exists(f"{output_folder}/fitting_bkg/Data"): 
        os.makedirs(f"{output_folder}/fitting_bkg/Data", exist_ok=True)
    
    # Combine Run2 data (2016-2018)
    logging.info("Merging Run2 data (2016-2018)")
    os.system(f"hadd -f {output_folder}/fitting_bkg/Data/output_Data_Run2.root "
              f"{output_folder}/fitting_bkg/Data_2016preVFP/output_*.root "
              f"{output_folder}/fitting_bkg/Data_2016postVFP/output_*.root "
              f"{output_folder}/fitting_bkg/Data_2017/output_*.root "
              f"{output_folder}/fitting_bkg/Data_2018/output_*.root")
    
    # Combine Run3 data (2022-2023)
    logging.info("Merging Run3 data (2022-2023)")
    os.system(f"hadd -f {output_folder}/fitting_bkg/Data/output_Data_Run3.root "
              f"{output_folder}/fitting_bkg/Data_2022preEE/output_*.root "
              f"{output_folder}/fitting_bkg/Data_2022postEE/output_*.root "
              f"{output_folder}/fitting_bkg/Data_2023preBPix/output_*.root "
              f"{output_folder}/fitting_bkg/Data_2023postBPix/output_*.root")

    # Combine all data
    logging.info("Merging all data")
    os.system(f"hadd -f {output_folder}/fitting_bkg/Data/output_Data_all.root "
              f"{output_folder}/fitting_bkg/Data/output_Data_Run2.root "
              f"{output_folder}/fitting_bkg/Data/output_Data_Run3.root")

def process_signal_files(converter, output_folder, input_folder):
    """Process signal files and write the results to output ROOT files."""
    region = 'VBF'  # Using VBF as default region
    syst_variations = [
        "nominal"
        # , "FNUF_up", "FNUF_down", "Material_up", "Material_down",
        # "Scale_up", "Scale_down", "Smearing_up", "Smearing_down", "JER_up", 
        # "JER_down", "JES_up", "JES_down", "MET_JES_up", "MET_JES_down",
        # "MET_Unclustered_up", "MET_Unclustered_down", "Muon_pt_up", "Muon_pt_down"
    ]
    procductions = ['ggH_M125', 'VBF_M125', 'WH_M125', 'ZH_M125', 'ttH_M125']
    years = ['2016preVFP', '2016postVFP', '2017', '2018', '2022preEE', '2022postEE', '2023preBPix', '2023postBPix']

    # Create all necessary directories in one go
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all process-year combinations and write output
    for proc, year in [(p, y) for p in procductions for y in years]:
        proc_name = proc.replace('M125', '125')
        if not os.path.exists(f"{output_folder}/fitting_signal/{proc}_{year}/output_{proc}.root"): 
            os.makedirs(f"{output_folder}/fitting_signal/{proc}_{year}", exist_ok=True)
        
        with uproot.recreate(f"{output_folder}/fitting_signal/{proc}_{year}/output_{proc}.root") as outfile:
            logging.info(f"Output file: {output_folder}/fitting_signal/{proc}_{year}/output_{proc}.root")
            
            for syst in syst_variations:
                logging.info(f"Processing {proc} {year} {syst}")
                
                # Map the year to the corresponding format in input files
                input_year = map_year(year)
                input_sample = map_sample_name(proc)
                
                if syst == 'nominal':
                    input_path = f'{input_folder}{input_sample}_{input_year}_output.root'
                    syst_name = ''
                else:
                    input_path = f'{input_folder}{proc}_{syst}/{input_year}.root'
                    syst_parts = syst.split('_')
                    syst_type = syst_parts[0] + ''.join(part.capitalize() for part in syst_parts[1:-1])
                    syst_uod = syst_parts[-1].capitalize()
                    syst_name = f'_{syst_type}{syst_uod}01sigma'
                
                try:
                    with uproot.open(input_path) as infile:
                        # Use the same method to read signal files
                        ak_data = infile["outtree"].arrays(library='ak')
                        try:
                            # Use the correct ak.to_pandas() function
                            data = ak.to_pandas(ak_data)
                        except Exception as e:
                            logging.warning(f"Error converting awkward array to pandas directly: {e}")
                            columns = ak_data.fields
                            data = pd.DataFrame({col: ak_data[col].to_numpy() for col in columns})
                        
                        data.columns = [col.replace('up', 'Up').replace('down', 'Down') for col in data.columns]
                        
                        # Rename llphoton_refit_m to CMS_hgg_mass
                        data = data.rename(columns={"llphoton_refit_m": "CMS_hgg_mass", "weight_corr": "weight"})
                        
                        # Add new dZ branch filled with 0.0
                        data['dZ'] = 0.0
                        
                        # Split data by BDT categories and lepton types
                        data_splits, lepton_splits = converter.split_by_categories(data)
                        
                        # Write categories to output file
                        for i, split_data in enumerate(data_splits):
                            logging.info(f"Number of events in {region}{i}{syst_name}: {len(split_data)}")
                            tree_name = f'{proc_name}_13TeV_{region}{i}{syst_name}'
                            outfile[f'DiphotonTree/{tree_name}'] = split_data
                        
                        # Write lepton-specific categories
                        for cat_name, split_data in lepton_splits.items():
                            i, flav = cat_name.split('_')
                            logging.info(f"Number of events in {region}{i}{flav}{syst_name}: {len(split_data)}")
                            tree_name = f'{proc_name}_{flav}_13TeV_{region}{i}{syst_name}'
                            outfile[f'DiphotonTree/{tree_name}'] = split_data
                except Exception as e:
                    logging.error(f"Error processing {input_path}: {e}")
                    continue

if __name__ == "__main__":
    args = get_args()
    converter = DataConverter()
    
    if args.type.lower() == 'data':
        process_data_files(converter, args.outputFolder, args.inputFolder)
    elif args.type.lower() == 'signal':
        process_signal_files(converter, args.outputFolder, args.inputFolder)
    else:
        logging.error(f"Unknown process type: {args.type}")
