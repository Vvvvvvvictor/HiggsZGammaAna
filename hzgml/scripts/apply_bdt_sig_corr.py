import os
import json
import uproot
import pandas as pd
import numpy as np
import xgboost as xgb
import pickle
from argparse import ArgumentParser
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import logging

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=logging.DEBUG)

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Apply BDT signal correction to data')
    parser.add_argument('-c', '--config', default='data/training_config_BDT.json', help='Path to the training config file')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-j/jiehan/root/skimmed_ntuples', help='Path to the input folder')
    parser.add_argument('-m', '--modelFolder', default='models', help='Path to the model folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/home-j/jiehan/root/fitting_signal', help='Path to the output folder')
    parser.add_argument('-s', '--catSplitFolder', default='/eos/home-j/jiehan/root/outputs/test/significances', help='Path to the category split folder')
    return parser.parse_args()

class BDTApplicator:
    def __init__(self, region, args):
        """Initialize the BDT Applicator with models and transformers."""
        self.region = region
        self.args = args
        self.models = self._load_models()
        self.transformers = self._load_transformers()
        self.train_variables, self.preselection = self._load_config()
        # logging.info(f"preselection: {self.preselection}")
        self.cat_boundaries = self._load_category_boundaries()
        self.random_index = 'event'

    def _load_models(self):
        """Load all BDT models for the specified region."""
        logging.info(f"Loading BDT models for {self.region}")
        return [xgb.Booster(model_file=f'{self.args.modelFolder}/BDT_{self.region}_{i}.h5') for i in range(4)]

    def _load_transformers(self):
        """Load all transformers for the specified region."""
        logging.info(f"Loading transformers for {self.region}")
        return [
            pickle.load(open(f'{self.args.modelFolder}/BDT_tsf_{self.region}_{i}.pkl', 'rb')) 
            for i in range(4)
        ]

    def _load_config(self):
        """Load the training configuration."""
        logging.info(f"Loading config from {self.args.config}")
        with open(self.args.config, 'r') as f:
            config = json.load(f)
        common_vars = config['common']['train_variables']
        region_vars = config[self.region].get('+train_variables', [])
        preselection = config['common']['preselections'] + config[self.region].get('+preselections', [])
        return common_vars + region_vars, preselection

    def _load_category_boundaries(self):
        """Load category boundaries from a file."""
        with open(f'{self.args.catSplitFolder}/bin_boundaries_1D_{self.region}.txt', 'r') as f:
            return list(map(float, f.readline().strip().split()))

    def preselect(self, data):
        """Apply preselection filters to the input data."""
        return data.query(' and '.join(self.preselection))

    def apply_bdt(self, data):
        """Apply BDT models and return split data by categories."""
        # data = self.preselect(data)
        out_data = pd.concat(
            [self._process_subset(data, i) for i in range(4)], ignore_index=True
        )
        return self._split_by_categories(out_data)

    def _process_subset(self, data, i):
        """Process a subset of data with a specific BDT model and transformer."""
        subset = data[data[self.random_index] % 4 == i].copy()
        scores = self.models[i].predict(xgb.DMatrix(subset[self.train_variables]))
        # Use .loc to safely assign new columns
        subset.loc[:, 'bdt_score'] = scores
        subset.loc[:, 'bdt_score_t'] = self.transformers[i].transform(scores.reshape(-1, 1)).flatten()
        return subset

    def _split_by_categories(self, data):
        """Split the data based on category boundaries."""
        return [
            data.query(f'{self.cat_boundaries[i]} <= bdt_score_t < {self.cat_boundaries[i+1]}')
            for i in range(len(self.cat_boundaries) - 1)
        ]

def filter_columns(df):
    """Filter columns of the dataframe based on a predefined list."""
    columns_to_keep = []
    for col in df.columns:
        if col in ["CMS_hgg_mass", "weight", "dZ"] or col.endswith("Up") or col.endswith("Down") or col.endswith("central"):
            columns_to_keep.append(col)
    
    # Ensure essential columns are not accidentally dropped if they don't match the pattern
    for essential_col in ["CMS_hgg_mass", "weight", "dZ"]:
        if essential_col in df.columns and essential_col not in columns_to_keep:
            columns_to_keep.append(essential_col)
            
    return df[columns_to_keep]

def process_files(applicators, output_folder, input_folder):
    """Process input files and write the results to output ROOT files."""
    region_map = {
        'zero_to_one_jet': 'ggH', 'two_jet': 'VBF',
        # 'VH': 'VHlep', 'ZH': 'ZHinv', 'ttH_had': 'ttHh', 'ttH_lep': 'ttHl'
        # 'two_jet': 'VBF'
        # 'all_jet': 'Incl'
    }
    base_systs = ["nominal"] + [f"{syst}_{direction}" for syst in ["Photon_scale", "Photon_smear", "Electron_scale", "Electron_smear", "JER", "JES", "MET_JES", "MET_Unclustered", "Muon_pt_smear"] for direction in ["up", "down"]]
    run3_systs = [f"Muon_pt_scale_{direction}" for direction in ["up", "down"]]
    
    procductions = ['ggH_M125', 'VBF_M125', 'WplusH_M125', 'WminusH_M125', 'ZH_M125', 'ttH_M125'] # 'ggH_M125', 'VBF_M125', 'WplusH_M125', 'WminusH_M125', 'ZH_M125', 'ttH_M125'
    years = ['2016preVFP', '2016postVFP'] #'2016preVFP', '2016postVFP', '2017', '2018', '2022preEE', '2022postEE', '2023preBPix', '2023postBPix'

    # Create all necessary directories in one go
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all process-year combinations and write output
    for proc, year in [(p, y) for p in procductions for y in years]:
        systs_for_year = base_systs[:]
        if '2022' in year or '2023' in year:
            systs_for_year.extend(run3_systs)
        
        syst_variations = ["nominal"] + systs_for_year

        proc_name = proc.replace('M125', '125')
        if not os.path.exists(f"{output_folder}/{proc}_{year}/output_{proc}.root"): os.makedirs(f"{output_folder}/{proc}_{year}", exist_ok=True)
        with uproot.recreate(f"{output_folder}/{proc}_{year}/output_{proc}.root") as outfile:
            logging.info(f"Output file: {output_folder}/{proc}_{year}/output_{proc}.root")
            for syst in syst_variations:
                logging.info(f"Processing {proc} {year} {syst}")
                if syst == 'nominal':
                    input_path = f'{input_folder}/{proc}/{year}.root'
                    syst_name = ''
                else:
                    input_path = f'{input_folder}/{proc}_{syst}/{year}.root'
                    syst_parts = syst.split('_')
                    syst_type = syst_parts[0] + ''.join(part.capitalize() for part in syst_parts[1:-1])
                    syst_uod = syst_parts[-1].capitalize()
                    syst_name = f'_{syst_type}{syst_uod}01sigma'
                with uproot.open(input_path) as infile:
                    for channel in region_map.keys():
                        all_data = infile[channel].arrays(library='pd')
                        all_data.columns = [col.replace('up', 'Up').replace('down', 'Down') for col in all_data.columns]
                        if 'dZ' not in all_data.columns:
                            all_data['dZ'] = np.zeros(len(all_data))
                        for lep in ['ele', 'mu']:
                            # data = all_data.query('nel==2' if lep == 'ele' else 'nmu==2')
                            data = all_data
                            data.loc[:, 'llphoton_hmiss_photon_dphi'] = data['photon_mht_deltaphi']
                            if channel in applicators:
                                applicator = applicators[channel]
                                data_splits = applicator.apply_bdt(data)
                                for i, split_data in enumerate(data_splits):
                                    logging.info(f"Number of events in {region_map[channel]}{i}{syst_name}: {len(split_data)}")
                                    split_data = split_data.rename(columns={"H_mass": "CMS_hgg_mass"})
                                    split_data = filter_columns(split_data)
                                    tree_name = f'{proc_name}_{lep}_13TeV_{region_map[channel]}{i}{syst_name}'
                                    outfile[f'DiphotonTree/{tree_name}'] = split_data
                            else:
                                logging.info(f"Number of events in {region_map[channel]}{syst_name}: {len(data)}")
                                data = data.rename(columns={"H_mass": "CMS_hgg_mass"})
                                data = filter_columns(data)
                                tree_name = f'{proc_name}_{lep}_13TeV_{region_map[channel]}{syst_name}'
                                outfile[f'DiphotonTree/{tree_name}'] = data

if __name__ == "__main__":
    args = get_args()
    applicators = {
        'zero_to_one_jet': BDTApplicator('zero_to_one_jet', args),
        'two_jet': BDTApplicator('two_jet', args),
        # 'all_jet': BDTApplicator('all_jet', args)
    }
    process_files(applicators, args.outputFolder, args.inputFolder)
