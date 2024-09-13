import uproot
import numpy as np
import os
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
import gc

import copy
from argparse import ArgumentParser
import json
import pandas as pd
# from root_pandas import *
import pickle
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import xgboost as xgb
from tqdm import tqdm
import logging
from pdb import set_trace
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError + 1
pd.options.mode.chained_assignment = None

from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', nargs=2, default=['/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/training_config_BDT.json', '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/apply_config_BDT.json'], help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', default='/eos/home-j/jiehan/data_for_norm_float_v1', help='directory of training inputs')
    parser.add_argument('-m', '--modelFolder', action='store', default='/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/models', help='directory of BDT models')
    parser.add_argument('-o', '--outputFolder', action='store', default='/eos/home-j/jiehan/data_for_norm_float_v1/output', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH', 'all_jet'], default='zero_to_one_jet', help='Region to process')
    parser.add_argument('-cat', '--category', action='store', nargs='+', help='apply only for specific categories')

    return parser.parse_args()

class ApplyXGBHandler(object):
    "Class for applying XGBoost"

    def __init__(self, configPath, region=''):

        print('===============================')
        print('  ApplyXGBHandler initialized')
        print('===============================')

        args=getArgs()
        self._region = region
        self._inputFolder = ''
        self._inputTree = region if region else 'inclusive'
        self._modelFolder = ''
        self._outputFolder = ''
        self._chunksize = 500000
        self._category = []
        self._branches = []
        self._outbranches = []

        self.m_models = {}
        self.m_tsfs = {}

        self.train_variables = {}
        self.randomIndex = 'eventNumber'

        self.models = {}
        self.observables = []
        self.preselections = []

        self.readApplyConfig(configPath[1])
        self.readTrainConfig(configPath[0])
        self.arrangeBranches()
        self.arrangePreselections()

    def readApplyConfig(self, configPath):
        """Read configuration file formated in json to extract information to fill TemplateMaker variables."""
        try:
            member_variables = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("_") and not attr.startswith('m_')]

            stream = open(configPath, 'r')
            configs = json.loads(stream.read())

            # read from the common settings
            config = configs["common"]
            for member in config.keys():
                if member in member_variables:
                    setattr(self, member, config[member])

            # read from the region specific settings
            if self._region:
                config = configs[self._region]
                for member in config.keys():
                    if member in member_variables:
                        setattr(self, member, config[member])
                if '+preselections' in config.keys():
                    self.preselections += config['+preselections']
                if '+observables' in config.keys():
                    self.observables += config['+observables']

        except Exception as e:
            logging.error("Error reading apply configuration '{config}'".format(config=configPath))
            logging.error(e)

    def readTrainConfig(self, configPath):

        try:
            stream = open(configPath, 'r')
            configs = json.loads(stream.read())
   
            config = configs["common"]

            if 'randomIndex' in config.keys(): self.randomIndex = config['randomIndex']
 
            if self.models:
                for model in self.models:
    
                    # read from the common settings
                    config = configs["common"]
                    if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables'][:]
    
                    # read from the region specific settings
                    if model in configs.keys():
                        config = configs[model]
                        if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables'][:]
                        if '+train_variables' in config.keys(): self.train_variables[model] += config['+train_variables']

                        print("\n\n")
                        print(self.train_variables[model])
                        print(len(self.train_variables[model]))
                        print("\n\n")

        except Exception as e:
            logging.error("Error reading training configuration '{config}'".format(config=configPath))
            logging.error(e)

    def arrangeBranches(self):

        self._branches = set()
        for model in self.models:
            self._branches = self._branches | set(self.train_variables[model])

        self._branches = self._branches | set([self.randomIndex]) | set([p.split()[0] for p in self.preselections]) | set(self.observables)
        self._branches = list(self._branches)

        for model in self.models:
            self.train_variables[model] = [x.replace('noexpand:', '') for x in self.train_variables[model]]
        self.preselections = [x.replace('noexpand:', '') for x in self.preselections]
        self.randomIndex = self.randomIndex.replace('noexpand:', '')

        self._outbranches = [branch for branch in self._branches if 'noexpand' not in branch]

    def arrangePreselections(self):

        if self.preselections:
            self.preselections = ['data.' + p for p in self.preselections]

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setModelFolder(self, modelFolder):
        self._modelFolder = modelFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def preselect(self, data):

        for p in self.preselections:
            data = data[eval(p)]

        return data

    def loadModels(self):

        if self.models:
            for model in self.models:
                print('XGB INFO: Loading BDT model: ', model)
                self.m_models[model] = []
                for i in range(4):
                    bst = xgb.Booster()
                    bst.load_model('%s/BDT_%s_%d.h5'%(self._modelFolder, model, i))
                    # bst = pickle.load(open('%s/BDT_%s_%d.pkl'%(self._modelFolder, model, i), "rb" ), encoding = 'latin1' ) # FIXME: need to change depending on region
                    self.m_models[model].append(bst)
                    del bst

    def loadTransformer(self):
        
        if self.models:
            for model in self.models:
                print('XGB INFO: Loading score transformer for model: ', model)
                self.m_tsfs[model] = []
                for i in range(4):
                    tsf = pickle.load(open('%s/BDT_tsf_%s_%d.pkl'%(self._modelFolder, model, i), "rb" ), encoding = 'latin1' )
                    self.m_tsfs[model].append(tsf)

    def applyBDT(self, category, scale=1):
        outputbraches = copy.deepcopy(self._outbranches)
        branches = copy.deepcopy(self._branches)
        branches += ["gamma_mvaID_WPL", 'regions', "gamma_eta", "Z_mass", "H_mass"]
        outputbraches += []
        
        outputContainer = self._outputFolder + '/' + self._region
        output_path = outputContainer + '/%s_fake.root' % category
        if not os.path.isdir(outputContainer): os.makedirs(outputContainer)
        if os.path.isfile(output_path): os.remove(output_path)

        f_list = []
        cat_folder = self._inputFolder + '/' + category
        for f in os.listdir(cat_folder):
            if f.endswith('.root'): f_list.append(cat_folder + '/' + f)

        print('-------------------------------------------------')
        for f in f_list: print('XGB INFO: Including sample: ', f)

        #TODO put this to the config
        # for data in tqdm(read_root(sorted(f_list), key=self._inputTree, columns=branches, chunksize=self._chunksize), bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}', desc='XGB INFO: Applying BDTs to %s samples' % category):
        with uproot.recreate(output_path) as output_file:
            out_data = pd.DataFrame()
            for filename in tqdm(sorted(f_list), desc='XGB INFO: Applying BDTs to %s samples' % category, bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
                print(filename)
                file = uproot.open(filename)
                print(self._inputTree)
                print(branches)
                data = file[self._inputTree].arrays(library='pd')
                data = data[(data["H_mass"]>100) & (data["H_mass"]<180) & (data["Z_mass"] + data["H_mass"] > 185)]
                
                output_selection = ((data['regions'] == 2) | (data['regions'] == 1))
                
                if not os.path.isdir('figures/fake'): 
                    os.makedirs('figures/fake')
                # get the distribution of 'gamma_mvaID' with 'region==0', and mimic the distribution to 'region==2'
                
                # gamma_mvaID_dist_B = data[(data['regions'] == 0) & (data['gamma_eta'] < 1.5)]['gamma_mvaID'].values
                # gamma_mvaID_dist_E = data[(data['regions'] == 0) & (data['gamma_eta'] > 1.5)]['gamma_mvaID'].values
                # prob_dist_B, gamma_mvaID_dist_B = np.histogram(gamma_mvaID_dist_B, bins=100, density=True)
                # prob_dist_E, gamma_mvaID_dist_E = np.histogram(gamma_mvaID_dist_E, bins=100, density=True)
                # gamma_mvaID_dist_B, gamma_mvaID_dist_E = (gamma_mvaID_dist_B[1:] + gamma_mvaID_dist_B[:-1])/2, (gamma_mvaID_dist_E[1:] + gamma_mvaID_dist_E[:-1])/2
                
                
                # Load the data
                gamma_mvaID_dist_B = pd.read_pickle("/eos/user/j/jiehan/data_for_norm_float_v1/gamma_mvaID_2J_B.pkl")['gamma_mvaID']
                gamma_mvaID_dist_E = pd.read_pickle("/eos/user/j/jiehan/data_for_norm_float_v1/gamma_mvaID_2J_E.pkl")['gamma_mvaID']

                # KDE for each distribution
                kde_B = gaussian_kde(gamma_mvaID_dist_B, bw_method='scott')
                kde_E = gaussian_kde(gamma_mvaID_dist_E, bw_method='scott')

                # Define the ranges for sampling
                min_B, max_B = gamma_mvaID_dist_B.min(), gamma_mvaID_dist_B.max()
                min_E, max_E = gamma_mvaID_dist_E.min(), gamma_mvaID_dist_E.max()

                # Generate probability distributions using KDE
                gamma_mvaID_x_B = np.linspace(0.42, 1, 101)[1:] - (1-0.42)/200
                gamma_mvaID_x_E = np.linspace(0.14, 1, 101)[1:] - (1-0.14)/200

                prob_dist_B = kde_B(gamma_mvaID_x_B)
                prob_dist_E = kde_E(gamma_mvaID_x_E)

                # Plot the original and smoothed distributions
                plt.figure()
                plt.clf()
                
                plt.hist(gamma_mvaID_dist_B, bins=100, histtype='step', label='B_original', density=True)
                plt.hist(gamma_mvaID_dist_E, bins=100, histtype='step', label='E_original', density=True)

                plt.plot(gamma_mvaID_x_B, prob_dist_B, label='B_smooth')
                plt.plot(gamma_mvaID_x_E, prob_dist_E, label='E_smooth')

                plt.legend()
                plt.xlabel('gamma_mvaID', fontsize=28)
                plt.ylabel('a.u.', fontsize=28)
                plt.title("gamma_mvaID distribution", fontsize=36)
                plt.savefig('figures/fake/smooth_gamma_mvaID.png')

                plt.clf()

                # Plot histogram of the original data
                plt.hist(data[(data['regions'] == 2)]['gamma_mvaID'], bins=100, histtype='step', label='Original')

                # Sampling from KDE
                def sample_from_kde(kde, size, min_val, max_val):
                    sampled_values = np.empty(size)
                    samples = kde.resample(20*size).flatten()
                    mask = (samples >= min_val) & (samples <= max_val)
                    if np.sum(mask) >= size:
                        sampled_values = np.random.choice(samples[mask], size=size, replace=False)

                    return sampled_values

                # Sample values for each category
                data_filtered_B = data[output_selection & (data['gamma_eta'] < 1.5)]
                data_filtered_E = data[output_selection & (data['gamma_eta'] > 1.5)]

                sampled_values_B = sample_from_kde(kde_B, data_filtered_B.shape[0], min_B, max_B)
                sampled_values_E = sample_from_kde(kde_E, data_filtered_E.shape[0], min_E, max_E)
                print(sampled_values_B)

                # Assign sampled values to the DataFrame
                data.loc[output_selection & (data['gamma_eta'] < 1.5), 'gamma_mvaID'] = sampled_values_B
                data.loc[output_selection & (data['gamma_eta'] > 1.5), 'gamma_mvaID'] = sampled_values_E

                # Plot histogram of the transformed data
                plt.hist(data[(data['regions'] == 2)]['gamma_mvaID'], bins=100, histtype='step', label='Transformed')
                plt.legend(fontsize=28)
                plt.xlabel('gamma_mvaID', fontsize=28)
                plt.ylabel('a.u.', fontsize=28)
                plt.title("Transformed gamma_mvaID", fontsize=36)
                plt.yscale('log')
                name = filename.split('/')[-1].split('.')[0]
                plt.savefig(f'figures/fake/{name}_mimic_gamma_mvaID.png')

                # Further filter data
                data = data[output_selection & (data["Z_mass"] > 80) & (data["Z_mass"] < 100)]

                for i in range(4):

                    data_s = data[data[self.randomIndex]%4 == i]
                    data_o = data_s
                    # data_o = data_s[outputbraches]

                    for model in self.train_variables.keys():
                        x_Events = data_s[self.train_variables[model]]
                        dEvents = xgb.DMatrix(x_Events)
                        scores = self.m_models[model][i].predict(dEvents)
                        if len(scores) > 0:
                            scores_t = self.m_tsfs[model][i].transform(scores.reshape(-1,1)).reshape(-1)
                        else:
                            scores_t = scores
                    
                        xgb_basename = self.models[model]
                        data_o[xgb_basename] = scores
                        data_o[xgb_basename+'_t'] = scores_t

                    out_data = pd.concat([out_data, data_o], ignore_index=True, sort=False)

                # out_data.to_root(output_path, key='test', mode='a', index=False)
            print(out_data)
            output_file[self._inputTree] = out_data

            # for tree in "inclusive", "zero_jet", "one_jet", "two_jet", "zero_to_one_jet", "VH_ttH":
            #     print(tree)
            #     out_data = pd.DataFrame()
            #     if tree == self._inputTree: continue
            #     for filename in tqdm(sorted(f_list), desc='XGB INFO: Applying BDTs to %s samples' % category, bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
            #         print(filename)
            #         file = uproot.open(filename)
            #         data = file[tree].arrays(library='pd')
            #         out_data = pd.concat([out_data, data], ignore_index=True, sort=False)
            #         print(data)
            #     output_file[tree] = out_data
            #     print(out_data)

            del out_data, data_s, data_o, data
            gc.collect()


def main():

    args=getArgs()
    
    configPath = args.config
    xgb = ApplyXGBHandler(configPath, args.region)

    xgb.setInputFolder(args.inputFolder)
    xgb.setModelFolder(args.modelFolder)
    xgb.setOutputFolder(args.outputFolder)

    xgb.loadModels()
    xgb.loadTransformer() # FIXME: need to change depending on region

    # with open('/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/inputs_config.json') as f:
    #     config = json.load(f)
    # sample_list = config['sample_list']
    sample_list = ["Data", "ZGToLLG", "ZG2JToG2L2J", "DYJets_2J"] # "Data", "ZGToLLG", "ZG2JToG2L2J", "DYJets_2J"

    for category in sample_list:
        if args.category and category not in args.category: continue
        xgb.applyBDT(category)

    return

if __name__ == '__main__':
    main()
 
