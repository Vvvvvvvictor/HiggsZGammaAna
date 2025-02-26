import os
import json
import uproot
import pandas as pd
import xgboost as xgb
import pickle
from argparse import ArgumentParser
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import logging

from pdb import set_trace

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=logging.DEBUG)

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Apply BDT signal correction to data')
    parser.add_argument('-c', '--config', default='data/training_config_BDT.json', help='Path to the training config file')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-j/jiehan/root/skimmed_ntuples_run2/', help='Path to the input folder')
    parser.add_argument('-m', '--modelFolder', default='models', help='Path to the model folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/home-j/jiehan/root/fitting_bkg', help='Path to the output folder')
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
        logging.info(f"preselection: {self.preselection}")
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

def process_files(applicators, output_folder, input_folder):
    """Process input files and write the results to output ROOT files."""
    region_map = {
        'zero_to_one_jet': 'ggH', 'two_jet': 'VBF',
        'VH': 'VHlep', 'ZH': 'ZHinv', 'ttH_had': 'ttHh', 'ttH_lep': 'ttHl'
        # 'two_jet': 'VBF'
    }
    syst_variations = [
        "nominal", "FNUF_up", "FNUF_down", "Material_up", "Material_down",
        "Scale_up", "Scale_down", "Smearing_up", "Smearing_down", "JER_up", 
        "JER_down", "JES_up", "JES_down", "MET_JES_up", "MET_JES_down",
        "MET_Unclustered_up", "MET_Unclustered_down", "Muon_pt_up", "Muon_pt_down"
    ]
    procductions = ['Data']
    years = ['2016preVFP', '2016postVFP', '2017', '2018']

    # Create all necessary directories in one go
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all process-year combinations and write output
    for proc, year in [(p, y) for p in procductions for y in years]:
        if not os.path.exists(f"{output_folder}/{proc}_{year}/output_{proc}_{year}.root"): os.makedirs(f"{output_folder}/{proc}_{year}", exist_ok=True)
        with uproot.recreate(f"{output_folder}/{proc}_{year}/output_{proc}_{year}.root") as outfile:
            logging.info(f"Output file: {output_folder}/{proc}_{year}/output_{proc}_{year}.root")
            logging.info(f"Processing {proc} {year}")
            input_path = f'{input_folder}/{proc}/{year}.root'
            with uproot.open(input_path) as infile:
                for channel in region_map.keys():
                    data = infile[channel].arrays(library='pd')
                    if channel in applicators:
                        applicator = applicators[channel]
                        data_splits = applicator.apply_bdt(data)
                        for i, split_data in enumerate(data_splits):
                            logging.info(f"Number of events in {region_map[channel]}{i}: {len(split_data)}")
                            split_data = split_data.rename(columns={"H_mass": "CMS_hgg_mass"})
                            outfile[f'DiphotonTree/Data_13TeV_{region_map[channel]}{i}'] = split_data
                    else:
                        logging.info(f"Number of events in {region_map[channel]}: {len(data)}")
                        data = data.rename(columns={"H_mass": "CMS_hgg_mass"})
                        outfile[f'DiphotonTree/Data_13TeV_{region_map[channel]}'] = data
    if not os.path.exists(f"{output_folder}/Data"): os.makedirs(f"{output_folder}/Data", exist_ok=True)
    os.system(f"hadd -f {output_folder}/Data/output_Data_Run2.root {output_folder}/Data_*/output_{proc}_*.root")

if __name__ == "__main__":
    args = get_args()
    applicators = {
        'zero_to_one_jet': BDTApplicator('zero_to_one_jet', args),
        'two_jet': BDTApplicator('two_jet', args)
    }
    process_files(applicators, args.outputFolder, args.inputFolder)
