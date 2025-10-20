#!/usr/bin/env python3
"""
Apply BDT models to all skimmed ntuples produced by generate_condor_sub.py.
This script processes signal and background samples, applying models and
organizing output according to apply_bdt_*.py conventions.
"""

import os
import json
import uproot
import pandas as pd
import numpy as np
import xgboost as xgb
from argparse import ArgumentParser
import logging

# Set up logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=logging.INFO)

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Apply BDT models to all skimmed ntuples')
    parser.add_argument('-c', '--config', default='model_config.json', 
                        help='Path to the model config file')
    parser.add_argument('-i', '--inputFolder', 
                        default='/eos/home-j/jiehan/root/skimmed_ntuples/', 
                        help='Path to the input folder with skimmed ntuples')
    parser.add_argument('-m', '--modelFolder', 
                        default='/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/', 
                        help='Path to the model folder')
    parser.add_argument('-o', '--outputFolder', 
                        default='/eos/home-j/jiehan/root/fitting_signal', 
                        help='Path to the output folder')
    parser.add_argument('--process-type', choices=['signal', 'background'], 
                        default='signal', help='Type of process to apply models to')
    parser.add_argument('--region', choices=['zero_to_one_jet', 'two_jet'], 
                        default='zero_to_one_jet', help='Analysis region')
    parser.add_argument('--max-events', type=int, default=None,
                        help='Maximum number of events to process (for testing)')
    return parser.parse_args()

class ModelApplicator:
    def __init__(self, region, args):
        """Initialize the Model Applicator."""
        # Map region names
        self.region_mapping = {
            'zero_to_one_jet': 'ggf',
            'two_jet': 'vbf'
        }
        # BDT score thresholds for classification
        self.bdt_thresholds = {
            'ggf': [0.94, 0.83, 0.57],  # ggF categories
            'vbf': [0.91, 0.81, 0.48]   # VBF categories
        }
        self.region = region
        self.args = args
        self.mapped_region = self.region_mapping.get(region, region)
        self.config = self._load_config()
        self.models = self._load_models()
        self.variable_mapping = self.config.get('variable_mapping', {})
        self.model_variables = self.config['model_variables'].get(self.mapped_region, [])
        
        logging.info(f"Initialized model applicator for region: {region} -> {self.mapped_region}")

    def _load_config(self):
        """Load the model configuration."""
        # Try different possible config paths
        possible_paths = [
            self.args.config,
            os.path.join(os.path.dirname(__file__), self.args.config),
            os.path.join(os.path.dirname(__file__), 'data', 'model_config.json'),
            'model_config.json'
        ]
        
        config_path = None
        for path in possible_paths:
            if os.path.exists(path):
                config_path = path
                break
        
        if config_path is None:
            raise FileNotFoundError(f"Config file not found. Tried: {possible_paths}")
            
        logging.info(f"Loading model config from {config_path}")
        with open(config_path, 'r') as f:
            return json.load(f)

    def _load_models(self):
        """Load BDT models."""
        models = []
        
        # Determine model type based on region mapping
        model_type = self.mapped_region  # 'ggf' or 'vbf'
        
        for i in range(4):  # 4-fold cross-validation
            model_path = f'{self.args.modelFolder}/model_xgb_{i}_{model_type}_redwood_v1_ext_val.json'
            if os.path.exists(model_path):
                logging.info(f"Loading model: {model_path}")
                model = xgb.Booster()
                model.load_model(model_path)
                models.append(model)
            else:
                logging.warning(f"Model not found: {model_path}")
        
        if not models:
            raise FileNotFoundError(f"No models found for region {self.mapped_region}")
        
        return models

    def map_variables(self, data):
        """Map variable names to model naming convention."""
        mapped_data = data.copy()
        
        # Apply variable mapping
        for model_var, internal_var in self.variable_mapping.items():
            if internal_var in data.columns and model_var in self.model_variables:
                mapped_data[model_var] = data[internal_var]
        
        # Check for missing variables and fill with zeros
        missing_vars = [var for var in self.model_variables if var not in mapped_data.columns]
        
        if missing_vars:
            logging.warning(f"Missing variables: {missing_vars}, filling with 0.0")
            for var in missing_vars:
                mapped_data[var] = 0.0
        
        return mapped_data

    def classify_bdt_category(self, bdt_scores):
        """Classify events into BDT categories based on score thresholds."""
        thresholds = self.bdt_thresholds[self.mapped_region]
        categories = np.zeros(len(bdt_scores), dtype=int)
        
        for cat_idx, threshold in enumerate(thresholds):
            categories[bdt_scores >= threshold] = cat_idx
        
        # The last category (lowest score) gets index 3
        categories[bdt_scores < thresholds[-1]] = 3
        
        return categories

    def apply_model(self, data):
        """Apply BDT models and return data with scores and categories."""
        if data.empty:
            logging.warning("Empty dataset provided")
            return data
        
        # Map variables to model format
        mapped_data = self.map_variables(data)
        
        # Prepare data for model
        features = mapped_data[self.model_variables]
        
        # Apply models using cross-validation approach
        scores = []
        
        for i in range(len(self.models)):
            # Get subset for this fold
            if 'event' in data.columns:
                subset_mask = data['event'] % len(self.models) == i
            else:
                subset_mask = np.arange(len(data)) % len(self.models) == i
            
            subset_features = features[subset_mask]
            
            if len(subset_features) > 0:
                # Apply model
                fold_scores = self.models[i].predict(xgb.DMatrix(subset_features))
                
                # Store scores with their original indices
                scores.extend([(idx, score) for idx, score in zip(data.index[subset_mask], fold_scores)])
        
        # Sort by original index and extract scores
        scores.sort(key=lambda x: x[0])
        
        # Add scores to the original data
        result_data = data.copy()
        result_data['bdt_score'] = [score for _, score in scores]
        
        # Add BDT category classification
        result_data['bdt_category'] = self.classify_bdt_category(result_data['bdt_score'].values)
        
        logging.info(f"Applied model to {len(result_data)} events")
        return result_data

def get_systs(year):
    """Get systematic uncertainties for a given year."""
    if any(run2_year in year for run2_year in ["2016", "2017", "2018"]):
        return ["Photon_scale", "Photon_smear", "Electron_scale", "Electron_smear", "JER", "JES", "MET_JES", "MET_Unclustered", "Muon_pt_smear"]
    return ["Photon_scale", "Photon_smear", "Electron_scale", "Electron_smear", "JER", "JES", "MET_JES", "MET_Unclustered", "Muon_pt_scale", "Muon_pt_smear"]

def get_sample_output_name(sample_name):
    """Convert sample names for output according to special naming rules."""
    # Handle theory variation samples: *_up/*_down -> *_Tune_up/*_Tune_down
    if sample_name.endswith('_up') or sample_name.endswith('_down'):
        base_name = sample_name.rsplit('_', 1)[0]
        direction = sample_name.rsplit('_', 1)[1]
        return f"{base_name}_Tune_{direction}"
    
    # Handle renormalization scale samples: *_mu -> *_Mmu
    if sample_name.endswith('_mu'):
        base_name = sample_name.rsplit('_', 1)[0]
        return f"{base_name}_Mmu"
    
    return sample_name

def filter_columns(df):
    """Filter columns of the dataframe based on a predefined list."""
    columns_to_keep = []
    for col in df.columns:
        if (col in ["CMS_hgg_mass", "weight", "dZ", "bdt_score", "bdt_category"] or 
            col.endswith("Up") or col.endswith("Down") or col.endswith("central")):
            columns_to_keep.append(col)
            
    return df[columns_to_keep]

def process_input_files(applicator, input_folder, output_folder, process_type):
    """Process input files and apply models."""
    
    # Define years
    years = ['2016preVFP', '2016postVFP', '2017', '2018', '2022preEE', '2022postEE', '2023preBPix', '2023postBPix']
    
    # Create output directory
    os.makedirs(output_folder, exist_ok=True)
    
    if process_type == 'signal':
        # Signal samples - consistent with generate_condor_sub.py
        signal_samples = ['ggH_M125', 'VBF_M125', 'WplusH_M125', 'WminusH_M125', 'ZH_M125', 'ttH_M125']
        nominal_only_samples = [
            'ggH_M120', 'VBF_M120', 'WplusH_M120', 'WminusH_M120', 'ZH_M120', 'ttH_M120',
            'ggH_M130', 'VBF_M130', 'WplusH_M130', 'WminusH_M130', 'ZH_M130', 'ttH_M130',
            'ggH_up', 'VBF_up', 'WplusH_up', 'WminusH_up', 'ZH_up', 'ttH_up',
            'ggH_down', 'VBF_down', 'WplusH_down', 'WminusH_down', 'ZH_down', 'ttH_down',
            'ggH_mu', 'VBF_mu', 'WplusH_mu', 'WminusH_mu', 'ZH_mu', 'ttH_mu'
        ]
        
        # Process nominal signal samples
        for sample in signal_samples + nominal_only_samples:
            for year in years:
                input_path = f"{input_folder}/{sample}/{year}.root"
                
                if not os.path.exists(input_path):
                    logging.warning(f"Input file not found: {input_path}")
                    continue
                
                # Apply naming rules for output
                output_sample_name = get_sample_output_name(sample)
                output_dir = f"{output_folder}/{output_sample_name}_{year}"
                output_path = f"{output_dir}/output_{output_sample_name}.root"
                
                logging.info(f"Processing {sample} {year} -> {output_sample_name}")
                _process_single_file(applicator, input_path, output_path, output_sample_name, year, 'signal')
        
        # Process systematic variations (only for M125 signal samples)
        for sample in signal_samples:
            for year in years:
                systs = get_systs(year)
                for syst in systs:
                    for uod in ["up", "down"]:
                        syst_sample = f"{sample}_{syst}_{uod}"
                        input_path = f"{input_folder}/{syst_sample}/{year}.root"
                        
                        if not os.path.exists(input_path):
                            logging.warning(f"Input file not found: {input_path}")
                            continue
                        
                        output_sample_name = sample.replace('M125', '125')
                        output_dir = f"{output_folder}/{output_sample_name}_{year}"
                        output_path = f"{output_dir}/output_{output_sample_name}.root"
                        
                        logging.info(f"Processing {syst_sample} {year}")
                        _process_single_file(applicator, input_path, output_path, output_sample_name, year, 'signal', syst, uod)
    
    else:  # Background/Data
        processes = ['Data']
        for proc in processes:
            for year in years:
                input_path = f"{input_folder}/{proc}/{year}.root"
                
                if not os.path.exists(input_path):
                    logging.warning(f"Input file not found: {input_path}")
                    continue
                
                output_dir = f"{output_folder}/{proc}_{year}"
                output_path = f"{output_dir}/output_{proc}_{year}.root"
                
                logging.info(f"Processing {proc} {year}")
                _process_single_file(applicator, input_path, output_path, proc, year, 'background')

def _process_single_file(applicator, input_path, output_path, proc_name, year, process_type, syst=None, direction=None):
    """Process a single ROOT file with model application."""
    # Create output directory
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Region mapping for output tree names
    region_map = {
        'zero_to_one_jet': 'ggH', 
        'two_jet': 'VBF'
    }
    
    try:
        with uproot.open(input_path) as infile:
            # Handle existing files
            if os.path.exists(output_path):
                backup_path = output_path + ".backup"
                if os.path.exists(backup_path):
                    os.remove(backup_path)
                os.rename(output_path, backup_path)
            
            with uproot.recreate(output_path) as outfile:
                # Process all trees in the input file
                for key in infile.keys():
                    tree_name = key.split(';')[0]
                    
                    if tree_name not in region_map:
                        continue
                        
                    logging.info(f"Processing tree: {tree_name}")
                    
                    # Load data
                    data = infile[tree_name].arrays(library='pd', 
                                                  entry_stop=applicator.args.max_events)
                    
                    if len(data) == 0:
                        logging.warning(f"Empty tree: {tree_name}")
                        continue
                    
                    # Apply model
                    data_with_score = applicator.apply_model(data)
                    
                    # Prepare data according to process type
                    if process_type == 'signal':
                        _process_signal_data(outfile, data_with_score, proc_name, year, tree_name, region_map, syst, direction)
                    else:
                        _process_background_data(outfile, data_with_score, proc_name, year, tree_name, region_map)
        
        logging.info(f"Successfully processed: {output_path}")
        
    except Exception as e:
        logging.error(f"Error processing {input_path}: {e}")
        # Clean up partial files
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except:
                pass
        raise

def _process_signal_data(outfile, data, proc_name, year, tree_name, region_map, syst=None, direction=None):
    """Process signal data and write to output file."""
    # Add required columns
    data.columns = [col.replace('up', 'Up').replace('down', 'Down') for col in data.columns]
    if 'dZ' not in data.columns:
        data['dZ'] = np.zeros(len(data))
    
    # Process for both electron and muon channels
    for lep in ['ele', 'mu']:
        # Apply systematic variation naming if applicable
        syst_name = ''
        if syst and direction:
            syst_parts = syst.split('_')
            syst_type = syst_parts[0] + ''.join(part.capitalize() for part in syst_parts[1:])
            syst_uod = direction.capitalize()
            syst_name = f'_{syst_type}{syst_uod}01sigma'
        
        # Rename columns and filter
        processed_data = data.rename(columns={"H_mass": "CMS_hgg_mass"})
        processed_data = filter_columns(processed_data)
        
        # Create tree name
        output_tree_name = f'{proc_name}_{lep}_13TeV_{region_map[tree_name]}{syst_name}'
        outfile[f'DiphotonTree/{output_tree_name}'] = processed_data
        
        logging.info(f"Written {len(processed_data)} events to {output_tree_name}")

def _process_background_data(outfile, data, proc_name, year, tree_name, region_map):
    """Process background data and write to output file."""
    # Rename columns and filter
    processed_data = data.rename(columns={"H_mass": "CMS_hgg_mass"})
    processed_data = filter_columns(processed_data)
    
    # Create tree name following apply_bdt_bkg.py convention
    output_tree_name = f'Data_13TeV_{region_map[tree_name]}'
    outfile[f'DiphotonTree/{output_tree_name}'] = processed_data
    
    logging.info(f"Written {len(processed_data)} events to {output_tree_name}")

def main():
    args = get_args()
    
    logging.info("Starting model application")
    logging.info(f"Input folder: {args.inputFolder}")
    logging.info(f"Model folder: {args.modelFolder}")
    logging.info(f"Output folder: {args.outputFolder}")
    logging.info(f"Process type: {args.process_type}")
    logging.info(f"Region: {args.region}")
    
    try:
        # Initialize applicator
        applicator = ModelApplicator(args.region, args)
        
        # Process files
        process_input_files(applicator, args.inputFolder, args.outputFolder, args.process_type)
        
        logging.info("Model application completed successfully!")
        
    except Exception as e:
        logging.error(f"Error during model application: {e}")
        raise

if __name__ == "__main__":
    main()