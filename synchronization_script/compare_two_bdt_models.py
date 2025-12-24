#!/usr/bin/env python3
"""
Script to apply two different BDT models and compare their scores
- Your model (from hzgml/models)
- External model (from /eos/project/h/htozg-dy-privatemc/rzou/bdt/XGB_scores/)
"""

import os
import sys
import json
import glob
import numpy as np
import pandas as pd
import uproot
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from tqdm import tqdm
import logging
import gc
import warnings
import pickle
import joblib
warnings.filterwarnings('ignore')

from pdb import set_trace

# Try to import WeightedQuantileTransformer
try:
    # Add hzgml scripts to path
    hzgml_scripts_path = '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts'
    if hzgml_scripts_path not in sys.path:
        sys.path.insert(0, hzgml_scripts_path)
    
    from weighted_quantile_transformer import WeightedQuantileTransformer
    logging.info("Successfully imported WeightedQuantileTransformer")
except ImportError as e:
    logging.warning(f"Could not import WeightedQuantileTransformer: {e}")
    WeightedQuantileTransformer = None

# Set up logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def get_args():
    """Get command line arguments."""
    parser = ArgumentParser(description='Compare two BDT models')
    parser.add_argument('-i', '--input-dir', type=str, 
                        default='/eos/user/j/jiehan/root/skimmed_ntuples_rui_new/',
                        help='Directory containing input ROOT files')
    parser.add_argument('-o', '--output-dir', type=str, 
                        default='/eos/user/j/jiehan/root/outputs/bdt_comparison/',
                        help='Output directory for results')
    parser.add_argument('-r', '--region', type=str, choices=['two_jet', 'zero_to_one_jet', 'all_jet'],
                        default='two_jet', help='Region to process')
    parser.add_argument('--your-model-dir', type=str, 
                        default='/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/models',
                        help='Directory containing your trained models')
    parser.add_argument('--external-model-dir', type=str,
                        default='/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/',
                        help='Directory containing external models')
    parser.add_argument('--config-dir', type=str,
                        default='/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data',
                        help='Directory containing configuration files')
    parser.add_argument('--sample', type=str, nargs='+', 
                        default=['DYJetsToLL'], #'VBF_M125', 'ggH_M125', 'ZGToLLG', 'EWKZ2J'
                        help='Samples to process')
    parser.add_argument('--year', type=str, default='2016preVFP', 
                        help='Year to process')
    
    return parser.parse_args()

class BDTModelComparator:
    """Class to compare two different BDT models."""
    
    def __init__(self, args):
        self.args = args
        self.region = args.region
        self.input_dir = args.input_dir
        self.output_dir = args.output_dir
        self.your_model_dir = args.your_model_dir
        self.external_model_dir = args.external_model_dir
        self.config_dir = args.config_dir
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Load configurations
        self.load_configs()
        self.load_external_config()
        
        # Your model variables
        self.your_model_vars = []
        self.your_models = {}
        self.your_score_transformers = {}  # For transforming BDT output scores (not input features)
        
        # External model variables (need to be mapped)
        self.external_model_vars = []
        self.external_models = {}
        
        # Random index for fold splitting (like in apply_bdt.py)
        self.randomIndex = 'event'
        
        # Preselections (like in apply_bdt.py)
        self.preselections = []

        # Load models
        self.load_your_models()
        self.load_external_models()
    
    def get_external_variable_mapping(self):
        """Map external model variables to your variable names."""
        return self.external_config['variable_mapping']
    
    def load_configs(self):
        """Load training and apply configurations."""
        train_config_path = os.path.join(self.config_dir, 'training_config_BDT.json')
        apply_config_path = os.path.join(self.config_dir, 'apply_config_BDT.json')
        
        with open(train_config_path, 'r') as f:
            self.train_config = json.load(f)
            
        with open(apply_config_path, 'r') as f:
            self.apply_config = json.load(f)
            
        # Load preselections and randomIndex from apply config
        if 'common' in self.apply_config:
            config = self.apply_config['common']
            if 'preselections' in config:
                self.preselections = config['preselections']
            if 'randomIndex' in config:
                self.randomIndex = config['randomIndex']
                
        # Load region-specific preselections
        if self.region in self.apply_config:
            config = self.apply_config[self.region]
            if '+preselections' in config:
                self.preselections += config['+preselections']
    
    def load_external_config(self):
        """Load external model configuration."""
        external_config_path = os.path.join(os.path.dirname(__file__), 'external_model_config.json')
        
        with open(external_config_path, 'r') as f:
            self.external_config = json.load(f)
    
    def load_your_models(self):
        """Load your trained models and transformers."""
        try:
            # Get variables for your model
            common_vars = self.train_config['common']['train_variables']
            region_vars = self.train_config[self.region].get('+train_variables', [])
            self.your_model_vars = common_vars + region_vars
            
            logging.info(f"Your model variables: {self.your_model_vars}")
            
            # Try different model file patterns (prefer H5 over JSON due to format issues)
            model_patterns = [
                f'BDT_{self.region}_{{fold}}.h5',
                f'optuna_{self.region}/BDT_region_{self.region}_fold{{fold}}.json',
                f'optuna_{self.region}_*/BDT_region_{self.region}_fold{{fold}}.json',
                f'BDT_{self.region}_{{fold}}.json',
                f'{self.region}_{{fold}}.json',
                f'{self.region}_{{fold}}.h5'
            ]
            
            # Find the correct pattern
            model_dir_found = None
            for pattern in model_patterns:
                test_path = os.path.join(self.your_model_dir, pattern.format(fold=0))
                
                # Handle glob patterns for optuna directories
                if '*' in pattern:
                    matches = glob.glob(test_path)
                    if matches:
                        model_dir_found = pattern
                        break
                elif os.path.exists(test_path):
                    model_dir_found = pattern
                    break
            
            if not model_dir_found:
                logging.error(f"No model files found for region {self.region}")
                return
            
            logging.info(f"Using model pattern: {model_dir_found}")
            
            # Load your models (XGBoost format - json or h5)
            for fold in range(4):
                model_path = os.path.join(self.your_model_dir, model_dir_found.format(fold=fold))
                
                # Handle glob patterns
                if '*' in model_path:
                    matches = glob.glob(model_path)
                    if matches:
                        model_path = matches[0]  # Take the first match
                
                if os.path.exists(model_path):
                    try:
                        # Load XGBoost model
                        model = xgb.Booster()
                        model.load_model(model_path)
                        self.your_models[fold] = model
                        logging.info(f"Loaded your model fold {fold} from {model_path}")
                    except Exception as e:
                        logging.error(f"Failed to load model {model_path}: {e}")
                
                # Load score transformer (WeightedQuantileTransformer for BDT output) if exists
                transformer_patterns = [
                    f'BDT_tsf_{self.region}_{fold}.pkl',
                    f'optuna_{self.region}/BDT_tsf_{self.region}_{fold}.pkl'
                ]
                
                for transformer_pattern in transformer_patterns:
                    transformer_path = os.path.join(self.your_model_dir, transformer_pattern)
                    if os.path.exists(transformer_path):
                        try:
                            with open(transformer_path, 'rb') as f:
                                # This is the WeightedQuantileTransformer for transforming BDT output scores
                                self.your_score_transformers[fold] = pickle.load(f)
                            logging.info(f"Loaded output score transformer fold {fold}")
                            break
                        except Exception as e:
                            logging.error(f"Failed to load output score transformer {transformer_path}: {e}")
            
        except Exception as e:
            logging.error(f"Error loading your models: {e}")
    
    def load_external_models(self):
        """Load external models."""
        try:
            # Determine model type based on region
            if self.region == 'two_jet':
                model_type = 'vbf'
            else:
                model_type = 'ggf'
            
            # Get external model variables from config
            self.external_model_vars = self.external_config['external_model_variables'][model_type]
            logging.info(f"External model variables ({model_type}): {self.external_model_vars}")
            
            # Load external models (4 folds)
            for fold in range(4):
                model_path = os.path.join(self.external_model_dir, f'model_xgb_{fold}_{model_type}_redwood_v1_ext_val.json')
                if os.path.exists(model_path):
                    model = xgb.Booster()
                    model.load_model(model_path)
                    self.external_models[fold] = model
                    logging.info(f"Loaded external model fold {fold} ({model_type})")
            
        except Exception as e:
            logging.error(f"Error loading external models: {e}")
    
    def map_variables(self, data, target_vars, is_external=False):
        """Map variables in data to target variable names."""
        mapped_data = {}
        if is_external:
            # Map from external variable names to your data column names  
            var_mapping = self.get_external_variable_mapping()
            for ext_var in target_vars:
                if ext_var in var_mapping:
                    your_var = var_mapping[ext_var]
                    if your_var in data.columns:
                        mapped_data[ext_var] = data[your_var]
                    else:
                        logging.warning(f"Variable {your_var} not found in data for external model variable {ext_var}")
                elif ext_var in data.columns:
                    # Direct mapping if variable exists in data
                    mapped_data[ext_var] = data[ext_var]
                else:
                    logging.warning(f"External model variable {ext_var} not found")
        else:
            # For your model, use variables directly
            for var in target_vars:
                if var in data.columns:
                    mapped_data[var] = data[var]
                else:
                    logging.warning(f"Your model variable {var} not found in data")
        
        return pd.DataFrame(mapped_data, columns=target_vars)
    
    def preselect(self, data):
        """Apply preselections to data (similar to apply_bdt.py)."""
        for p in self.preselections:
            try:
                data = data.query(p)
            except Exception as e:
                logging.warning(f"Failed to apply preselection '{p}': {e}")
        return data
    
    def apply_models(self, sample, year):
        """Apply both models to a sample and save to separate ROOT file."""
        input_file = os.path.join(self.input_dir, sample, f'{year}.root')
        if not os.path.exists(input_file):
            logging.warning(f"Input file not found: {input_file}")
            return None
        
        # Output file for this sample
        output_file = os.path.join(self.output_dir, f'{sample}_{year}_{self.region}_bdt_comparison.root')
        
        logging.info(f"Processing {sample} {year}")
        
        try:
            # All variables to be saved, including all variables used by BDT models
            all_vars = set(self.your_model_vars) | set(self.external_model_vars)
            # Also save scores, H_mass, weight and randomIndex
            extra_vars = {'H_mass', 'weight_corr', self.randomIndex}
            all_vars = list(all_vars | extra_vars)
            
            # Store all data chunks
            all_data_chunks = []

            # Read and process in chunks
            for data in uproot.iterate(f"{input_file}:{self.region}", library="pd", step_size=100000):
                if len(data) == 0:
                    continue
                
                # Handle missing weight_corr/H_mass
                if 'weight_corr' not in data.columns:
                    data['weight_corr'] = 1.0
                if 'H_mass' not in data.columns:
                    data['H_mass'] = np.nan
                
                # Preselect data (similar to apply_bdt.py)
                data = self.preselect(data)
                
                # Process each fold separately (following apply_bdt.py pattern)
                for fold in range(4):
                    # Split data by fold using the same logic as apply_bdt.py
                    data_fold = data[data[self.randomIndex] % 314159 % 4 == fold]
                    if data_fold.shape[0] == 0: 
                        continue
                    
                    data_output = data_fold.copy()
                    # Calculate your model scores for this fold
                    if fold in self.your_models and self.your_model_vars:
                        your_data = self.map_variables(data_fold, self.your_model_vars, is_external=False)
                        if not your_data.empty:
                            try:
                                dmatrix = xgb.DMatrix(your_data[self.your_model_vars])
                                your_scores = self.your_models[fold].predict(dmatrix)
                                
                                # Apply score transformer if available
                                if fold in self.your_score_transformers and len(your_scores) > 0:
                                    try:
                                        your_scores_t = self.your_score_transformers[fold].transform(your_scores.reshape(-1,1)).reshape(-1)
                                        data_output['your_bdt_score_t'] = your_scores_t
                                    except Exception as e:
                                        logging.warning(f"Score transformer failed for fold {fold}: {e}")
                                        data_output['your_bdt_score_t'] = your_scores
                                else:
                                    data_output['your_bdt_score_t'] = your_scores
                                
                                data_output['your_bdt_score'] = your_scores
                            except Exception as e:
                                logging.error(f"Error predicting with your model fold {fold}: {e}")
                                data_output['your_bdt_score'] = np.full(len(data_fold), np.nan)
                                data_output['your_bdt_score_t'] = np.full(len(data_fold), np.nan)
                    else:
                        data_output['your_bdt_score'] = np.full(len(data_fold), np.nan)
                        data_output['your_bdt_score_t'] = np.full(len(data_fold), np.nan)
                    
                    # Calculate external model scores for this fold
                    if fold in self.external_models and self.external_model_vars:
                        external_data = self.map_variables(data_fold, self.external_model_vars, is_external=True)
                        if not external_data.empty:
                            try:
                                # # for column_name in self.external_model_vars:
                                # #     external_data.loc[:, column_name] = 1.0
                                # external_data = external_data[:1]
                                # event1 = [1.0, -0.668896, 0.111204, -0.506836, -1.183105, 0.186646, -1.970355, 0.447510, 0.080442, 2.203866, 2.207445, 0.264801, -999.0, -999.0, -999.0, 0.270044, -999.0, -999.0, -999.0, 0]
                                # for i, column_name in enumerate(self.external_model_vars):
                                #     # external_data_temp = external_data.copy(deep=True)
                                #     # if column_name in external_data.columns:
                                #     external_data.loc[:, column_name] = event1[i]
                                #     # print("fold: ", fold, " dataframe: ", external_data_temp[self.external_model_vars])
                                #     # for col in self.external_model_vars:
                                #     #     print(f"fold: {fold}, col: {col}, data: ", external_data_temp[col].values)
                                # dmatrix = xgb.DMatrix(external_data[self.external_model_vars])
                                # preds = self.external_models[fold].predict(dmatrix)
                                # print("fold: ", fold, " preds: ", preds)
                                # break


                                # Ensure column order matches external_model_vars
                                dmatrix = xgb.DMatrix(external_data[self.external_model_vars])
                                external_scores = self.external_models[fold].predict(dmatrix)
                                data_output['external_bdt_score'] = external_scores
                            except Exception as e:
                                logging.error(f"Error predicting with external model fold {fold}: {e}")
                                data_output['external_bdt_score'] = np.full(len(data_fold), np.nan)
                    else:
                        data_output['external_bdt_score'] = np.full(len(data_fold), np.nan)
                    
                    # Assemble content to save (keep only required variables)
                    save_dict = {}
                    var_mapping = self.get_external_variable_mapping()
                    for v in all_vars:
                        if v in data_output.columns:
                            save_dict[v] = data_output[v].values
                        elif v in var_mapping:
                            save_dict[v] = data_output[var_mapping[v]].values
                        else:
                            save_dict[v] = np.full(len(data_output), np.nan)
                    
                    # Add BDT scores
                    save_dict['your_bdt_score'] = data_output['your_bdt_score'].values
                    save_dict['your_bdt_score_t'] = data_output['your_bdt_score_t'].values
                    save_dict['external_bdt_score'] = data_output['external_bdt_score'].values
                    save_dict[self.randomIndex] = data_output[self.randomIndex].values
                    
                    if len(save_dict[list(save_dict.keys())[0]]) > 0:
                        all_data_chunks.append(save_dict)
                    
                    del data_output, save_dict
                    gc.collect()
                
                del data, data_fold
                gc.collect()
            
            # Merge all chunks and save to a single ROOT file
            if all_data_chunks:
                # Merge all chunks
                combined_dict = {}
                for key in all_data_chunks[0].keys():
                    combined_dict[key] = np.concatenate([chunk[key] for chunk in all_data_chunks])
                
                # Save to ROOT file
                with uproot.recreate(output_file) as fout:
                    fout['tree'] = combined_dict
                
                logging.info(f"Saved {sample} {year} to {output_file} with {len(combined_dict[list(combined_dict.keys())[0]])} events")
                
                del all_data_chunks, combined_dict
                gc.collect()
                
                return output_file
            else:
                logging.warning(f"No data processed for {sample} {year}")
                return None
                
        except Exception as e:
            logging.error(f"Error processing {sample} {year}: {e}")
            return None
    
    def plot_2d_comparison(self, all_results):
        """Create 2D plots comparing the two BDT scores with smooth curves and quantile bands."""
        from scipy import interpolate
        from scipy.ndimage import gaussian_filter1d
        
        plt.style.use('default')
        
        # Create plots for each sample
        for sample in set([r['sample'] for r in all_results if r is not None]):
            sample_results = [r for r in all_results if r is not None and r['sample'] == sample]
            
            if not sample_results:
                continue
            
            # Combine all years for this sample
            your_scores = np.concatenate([r['your_scores'] for r in sample_results if len(r['your_scores']) > 0])
            external_scores = np.concatenate([r['external_scores'] for r in sample_results if len(r['external_scores']) > 0])
            weights = np.concatenate([r['weight'] for r in sample_results])
            
            if len(your_scores) == 0 or len(external_scores) == 0:
                logging.warning(f"No scores available for {sample}")
                continue
            
            # Create figure
            plt.figure(figsize=(10, 8))
            
            # Define x bins for mean calculation
            n_bins = 50
            x_min, x_max = np.min(your_scores), np.max(your_scores)
            x_bins = np.linspace(x_min, x_max, n_bins + 1)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            
            # Calculate mean values and quantiles for each bin
            x_means = []
            y_means = []
            y_quantiles = {0.023: [], 0.159: [], 0.841: [], 0.977: []}
            
            for i in range(len(x_bins) - 1):
                # Find points in this bin
                mask = (your_scores >= x_bins[i]) & (your_scores < x_bins[i + 1])
                
                if np.sum(mask) < 5:  # Skip bins with too few points
                    continue
                
                x_bin_data = your_scores[mask]
                y_bin_data = external_scores[mask]
                w_bin_data = weights[mask]
                
                # Calculate weighted means
                x_mean = np.average(x_bin_data, weights=w_bin_data)
                y_mean = np.average(y_bin_data, weights=w_bin_data)
                
                x_means.append(x_mean)
                y_means.append(y_mean)
                
                # Calculate weighted quantiles
                # Sort by y values
                sorted_indices = np.argsort(y_bin_data)
                sorted_y = y_bin_data[sorted_indices]
                sorted_w = w_bin_data[sorted_indices]
                
                # Calculate cumulative weights
                cumulative_weights = np.cumsum(sorted_w)
                total_weight = cumulative_weights[-1]
                
                # Find quantiles
                for q in y_quantiles.keys():
                    target_weight = q * total_weight
                    idx = np.searchsorted(cumulative_weights, target_weight)
                    if idx < len(sorted_y):
                        y_quantiles[q].append(sorted_y[idx])
                    else:
                        y_quantiles[q].append(sorted_y[-1])
            
            # Convert to numpy arrays
            x_means = np.array(x_means)
            y_means = np.array(y_means)
            
            if len(x_means) < 3:
                logging.warning(f"Not enough points for smooth curves in {sample}")
                continue
            
            # Sort by x_means for smooth interpolation
            sort_idx = np.argsort(x_means)
            x_means_sorted = x_means[sort_idx]
            y_means_sorted = y_means[sort_idx]
            
            # Apply smoothing to quantiles
            for q in y_quantiles.keys():
                y_quantiles[q] = np.array(y_quantiles[q])[sort_idx]
                # Apply Gaussian smoothing
                y_quantiles[q] = gaussian_filter1d(y_quantiles[q], sigma=1.0)
            
            # Apply smoothing to means
            y_means_smooth = gaussian_filter1d(y_means_sorted, sigma=1.0)
            
            # Fill between quantiles
            # 2σ region (yellow): between 0.023 and 0.977 quantiles
            plt.fill_between(x_means_sorted, y_quantiles[0.023], y_quantiles[0.977], 
                           color='yellow', alpha=0.3, label='2σ')
            
            # 1σ region (green): between 0.159 and 0.841 quantiles  
            plt.fill_between(x_means_sorted, y_quantiles[0.159], y_quantiles[0.841], 
                           color='green', alpha=0.4, label='1σ')
            
            # Plot quantile curves
            plt.plot(x_means_sorted, y_quantiles[0.023], 'orange', linewidth=1, alpha=0.8)
            plt.plot(x_means_sorted, y_quantiles[0.159], 'darkgreen', linewidth=1, alpha=0.8)
            plt.plot(x_means_sorted, y_quantiles[0.841], 'darkgreen', linewidth=1, alpha=0.8)
            plt.plot(x_means_sorted, y_quantiles[0.977], 'orange', linewidth=1, alpha=0.8)
            
            # Plot mean curve
            plt.plot(x_means_sorted, y_means_smooth, 'blue', linewidth=2, label='Mean')
            
            plt.xlabel('Your BDT Score')
            plt.ylabel('External BDT Score')
            plt.title(f'BDT Score Comparison with Quantile Bands - {sample} ({self.region})')
            
            # Add correlation coefficient
            correlation = np.corrcoef(your_scores, external_scores)[0, 1]
            plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                    transform=plt.gca().transAxes, bbox=dict(boxstyle='round', facecolor='white'))
            
            # Add legend
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # Save plot
            output_path = os.path.join(self.output_dir, f'bdt_comparison_2d_{sample}_{self.region}.png')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"Saved 2D comparison plot: {output_path}")
            
            # Create 2D histogram heatmap with rainbow colormap
            plt.figure(figsize=(10, 8))
            
            # Create 2D histogram
            n_bins = 50
            hist, xbins, ybins = np.histogram2d(your_scores, external_scores, bins=n_bins, weights=weights)
            
            # Create meshgrid for plotting
            X, Y = np.meshgrid(xbins[:-1], ybins[:-1])
            
            # Plot 2D histogram as heatmap with rainbow colormap
            im = plt.pcolormesh(X, Y, hist.T, cmap='rainbow', shading='auto')
            
            # Add colorbar
            cbar = plt.colorbar(im)
            cbar.set_label('Weighted Event Count', fontsize=12)
            
            plt.xlabel('Your BDT Score')
            plt.ylabel('External BDT Score')
            plt.title(f'BDT Score 2D Histogram - {sample} ({self.region})')
            plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                    transform=plt.gca().transAxes, bbox=dict(boxstyle='round', facecolor='white'))
            
            # Add diagonal line
            min_score = min(np.min(your_scores), np.min(external_scores))
            max_score = max(np.max(your_scores), np.max(external_scores))
            plt.plot([min_score, max_score], [min_score, max_score], 'k--', alpha=0.8, linewidth=2, label='y=x')
            plt.legend()
            
            output_path = os.path.join(self.output_dir, f'bdt_comparison_heatmap_{sample}_{self.region}.png')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"Saved 2D heatmap plot: {output_path}")
    
    def load_results_for_plotting(self, output_files):
        """Load results from individual ROOT files for plotting."""
        all_results = []
        
        for output_file in output_files:
            if not output_file or not os.path.exists(output_file):
                continue
                
            try:
                # Extract sample name from filename
                basename = os.path.basename(output_file)
                # Format: SAMPLE_YEAR_REGION_bdt_comparison.root
                sample_name = basename.split('_')[0]
                year = basename.split('_')[1]
                
                with uproot.open(output_file) as f:
                    data = f['tree'].arrays(library='pd')
                    
                    if len(data) == 0:
                        continue
                        
                    result = {
                        'sample': sample_name,
                        'year': year,
                        'n_events': len(data),
                        'your_scores': data['your_bdt_score'].values,
                        'your_scores_t': data['your_bdt_score_t'].values if 'your_bdt_score_t' in data.columns else data['your_bdt_score'].values,
                        'external_scores': data['external_bdt_score'].values,
                        'H_mass': data['H_mass'].values if 'H_mass' in data.columns else None,
                        'weight': data['weight_corr'].values if 'weight_corr' in data.columns else np.ones(len(data))
                    }
                    all_results.append(result)
                    
            except Exception as e:
                logging.warning(f"Failed to load results from {output_file}: {e}")
                
        return all_results
    
    def run(self):
        """Run the comparison."""
        output_files = []
        
        for sample in self.args.sample:
            output_file = self.apply_models(sample, self.args.year)
            if output_file:
                output_files.append(output_file)
        
        if output_files:
            # 读取数据用于画图
            try:
                all_results = self.load_results_for_plotting(output_files)
                if all_results:
                    self.plot_2d_comparison(all_results)
            except Exception as e:
                logging.warning(f"Failed to load data for plotting: {e}")
        else:
            logging.warning("No results to process")

def main():
    args = get_args()
    
    logging.info("Starting BDT model comparison")
    logging.info(f"Region: {args.region}")
    logging.info(f"Samples: {args.sample}")
    logging.info(f"Year: {args.year}")
    
    comparator = BDTModelComparator(args)
    comparator.run()
    
    logging.info("Comparison completed")

if __name__ == '__main__':
    main()
