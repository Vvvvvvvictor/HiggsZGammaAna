import xgboost as xgb
import numpy as np
import pandas as pd
import uproot
import pickle
import json
import os
from scipy.stats import gaussian_kde
from argparse import ArgumentParser

from pdb import set_trace
import gc

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', nargs=2, default=['/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/training_config_BDT.json', '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/apply_config_BDT.json'], help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', default='/eos/home-j/jiehan/data_for_norm_float_v1', help='directory of training inputs')
    parser.add_argument('-m', '--modelFolder', action='store', default='/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/models', help='directory of BDT models')
    parser.add_argument('-o', '--outputFolder', action='store', default='/eos/home-j/jiehan/data_for_norm_float_v1/output_WP80', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH', 'all_jet'], default='zero_to_one_jet', help='Region to process')
    parser.add_argument('-cat', '--category', action='store', nargs='+', help='apply only for specific categories')
    parser.add_argument('-s', '--sample', action='store', nargs='+', help='apply only for specific samples')

    return parser.parse_args()

class ApplyXGBHandler(object):
    """
    Class for applying XGBoost model to data.
    """
    
    def __init__(self, region=''):
        """
        Initialize the class.
        """
        args = getArgs()
        self.args = args
        self.configPath = args.config
        self.inputFolder = args.inputFolder
        self.modelFolder = args.modelFolder
        self.outputFolder = args.outputFolder
        
        self.region = args.region if args.region else "inclusive"
        self.readTrainConfig(self.configPath[0])
        self.readApplyConfig(self.configPath[1])
        self.model = self.loadModel()
        
    def readTrainConfig(self, path):
        """
        Read training configuration.
        """
        self.trainVariables = ["Z_cos_theta", "lep_cos_theta", "H_ptt", "Z_lead_lepton_eta", "Z_sublead_lepton_eta", "gamma_eta", "lep_phi", "gamma_mvaID", "gamma_ptRelErr", "l1g_deltaR", "l2g_deltaR"]
        self.outputVariables = ["bdt_score", "H_mass", "weight"] + self.trainVariables
        self.preselections = "H_mass > 100 & H_mass < 180 & H_mass + Z_mass > 185"
        
        # select control region or signal region
        self.postselctions = "Z_mass > 80 & Z_mass < 100 & ((regions == 2) | (regions == 1))"
        # self.postselctions = "Z_mass > 80 & Z_mass < 100 & regions == 0"
        
    def readApplyConfig(self, path):
        """
        Read application configuration.
        """
        if self.args.sample:
            self.applySamples = list(self.args.sample)
        else:
            # self.applySamples = ["Data", "ZGToLLG", "ZG2JToG2L2J", "DYJets_01J"]
            self.applySamples = ["Data", "ZGToLLG", "ZG2JToG2L2J", "DYJets_01J"]
    
    def loadModel(self):
        """
        Load the trained model.
        """
        modelPath = f"{self.modelFolder}/BDT_{self.region}.pkl"
        model = pickle.load(open(modelPath, 'rb'))
        return model
    
    def applyModel(self):
        """
        Apply the model to data.
        """
        for sample in self.applySamples:
            print(f"Applying model to {sample}...")
            
            # set the output file
            with uproot.recreate(f"{self.outputFolder}/{self.region}/{sample}_fake.root") as output_file:
            # with uproot.recreate(f"{self.outputFolder}/{self.region}/{sample}.root") as output_file:
                data_out = pd.DataFrame()
                for f in os.listdir(f"{self.inputFolder}/{sample}"):
                    year = f.split('/')[-1].split('_')[-1].split('.')[0]
                    
                    print(f"Processing {f}...")
                    sample_data = uproot.open(f"{self.inputFolder}/{sample}/{f}")[self.region].arrays(library="pd")
                    sample_data = sample_data.query(self.preselections)
                    
                    # only control region need to change mvaID
                    sample_data = self.changeMvaID(sample_data)
                    # sample_data["gamma_mvaID"] = sample_data["gamma_mvaID"].apply(lambda x: 1)
                    
                    sample_data = sample_data.query(self.postselctions)
                    sample_data["bdt_score"] = self.model.predict_proba(sample_data[self.trainVariables].to_numpy())[:, 1]
                    # sample_data = sample_data[self.outputVariables]
                    data_out = pd.concat([data_out, sample_data], ignore_index=True, sort=False)
                
                    # # only control region need to check new mvaID distribution
                    # if not os.path.isdir(f'figures/fake_{self.region}'): 
                    #     os.makedirs(f'figures/fake_{self.region}')
                    # plt.figure()
                    # plt.hist(data_out['gamma_mvaID'], bins=100, histtype='step', label='gamma_mvaID')
                    # plt.legend()
                    # plt.savefig(f'figures/fake_{self.region}/{sample}_{year}_mimic_gamma_mvaID.png')
                    
                output_file[self.region] = data_out
                del data_out, sample_data
                gc.collect()
            
    def changeMvaID(self, data):
        """
        Change the mvaID to a value pass WP80
        """
        # 核密度估计并重新采样函数
        def kde_resample(data, condition, gamma_eta_selection):
            # gamma_mvaID_values = data[condition]['gamma_mvaID'].values
            gamma_mvaID_values = pd.read_pickle(f"/eos/user/j/jiehan/data_for_norm_float_v1/gamma_mvaID_01J_{condition}.pkl")['gamma_mvaID'].values
            
            kde = gaussian_kde(gamma_mvaID_values)
            gamma_mvaID_range = np.linspace(min(gamma_mvaID_values), max(gamma_mvaID_values), 1000)
            prob_dist = kde(gamma_mvaID_range)
            return np.random.choice(gamma_mvaID_range, size=data.query(gamma_eta_selection).shape[0], p=prob_dist/prob_dist.sum())

        # 重新采样并应用
        data.loc[(data['gamma_eta'] < 1.5) & (data['gamma_eta'] > -1.5), 'gamma_mvaID'] = kde_resample(data, "B", "(gamma_eta < 1.5) & (gamma_eta > -1.5)")
        data.loc[(data['gamma_eta'] > 1.5) | (data['gamma_eta'] < -1.5), 'gamma_mvaID'] = kde_resample(data, "E", "(gamma_eta > 1.5) | (gamma_eta < -1.5)")
        
        return data
        
            
if __name__ == "__main__":
    handler = ApplyXGBHandler()
    handler.applyModel()
    
        