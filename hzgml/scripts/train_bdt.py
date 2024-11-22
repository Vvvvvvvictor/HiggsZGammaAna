#!/usr/bin/env python
import os
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd
import uproot
import pickle
from sklearn.metrics import roc_curve, auc, confusion_matrix, roc_auc_score
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import xgboost as xgb
from tabulate import tabulate
import matplotlib.pyplot as plt
from tqdm import tqdm
from pdb import set_trace
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError + 1

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', default='data/training_config_BDT.json', help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', help='directory of training inputs')
    parser.add_argument('-o', '--outputFolder', action='store', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH', 'VBF', 'all_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-f', '--fold', action='store', type=int, nargs='+', choices=[0, 1, 2, 3], default=[0, 1, 2, 3], help='specify the fold for training')
    parser.add_argument('-p', '--params', action='store', type=dict, default=None, help='json string.') #type=json.loads
    parser.add_argument('--hyperparams_path', action='store', default=None, help='path of hyperparameters json') 
    parser.add_argument('--save', action='store_true', help='Save model weights to HDF5 file')
    parser.add_argument('--corr', action='store_true', default=False, help='Plot corelation between each training variables')
    parser.add_argument('--importance', action='store_true', default=True, help='Plot importance of variables, parameter "gain" is recommanded')
    parser.add_argument('--roc', action='store_true', default=True, help='Plot ROC')
    parser.add_argument('--optuna', action='store_true', default=False, help='Run hyperparameter tuning using optuna')
    parser.add_argument('--optuna_metric', action='store', default='auc', choices=['eval_auc', 'sqrt_eval_auc_minus_train_auc', 'eval_auc_minus_train_auc', 'eval_auc_over_train_auc'], help='Optuna metric to optimize')
    parser.add_argument('--n-calls', action='store', type=int, default=36, help='Steps of hyperparameter tuning using optuna')
    parser.add_argument('--continue-optuna', action='store', type=int, default=0, help='Continue tuning hyperparameters using optuna')

    parser.add_argument('-s', '--shield', action='store', type=int, default=-1, help='Which variables needs to be shielded')
    parser.add_argument('-a', '--add', action='store', type=int, default=-1, help='Which variables needs to be added')

    return parser.parse_args()

class XGBoostHandler(object):
    "Class for running XGBoost"

    def __init__(self, configPath, region=''):

        print("""
========================================================================
|| **     *     ** ###### $$        ######    ****    #       # %%%%% ||
||  **   ***   **  ##     $$       ##       ***  ***  ##     ## %     ||
||   ** ** ** **   ###### $$      ##       **      ** # #   # # %%%%% ||
||    ***   ***    ##     $$       ##       ***  ***  #  # #  # %     || 
||     *     *     ###### $$$$$$$   ######    ****    #   #   # %%%%% || 
||                                                                    ||
||    $$$$$$$$$$    ####         XX      XX     GGGGGG     BBBBB      ||
||        $$      ###  ###         XX  XX     GGG          B    BB    ||
||        $$     ##      ##          XX      GG     GGGG   BBBBB      ||
||        $$      ###  ###         XX  XX     GGG     GG   B    BB    ||
||        $$        ####         XX      XX     GGGGGG     BBBBB      ||
========================================================================
              """)

        args=getArgs()
        self.continue_optuna = args.continue_optuna
        self.optuna_metric = args.optuna_metric
        self.n_calls = args.n_calls
        self._shield = args.shield
        self._add = args.add

        self._region = region

        self._inputFolder = '/eos/home-j/jiehan/root/skimmed_ntuples_rui'
        self._outputFolder = 'models'
        self._chunksize = 500000
        self._branches = []
        self._sig_branches = []
        self._mc_branches = []
        self._data_branches = []

        self.m_data_sig = pd.DataFrame()
        self.m_data_bkg = pd.DataFrame()
        self.m_train_wt = {}
        self.m_val_wt = {}
        self.m_test_wt = {}
        self.m_y_train = {}
        self.m_y_val = {}
        self.m_y_test = {}
        self.m_dTrain = {}
        self.m_dVal = {}
        self.m_dTest = {}
        self.m_dTest_sig = {}
        self.m_dTest_bkg = {}
        self.m_bst = {}
        self.m_score_train = {}
        self.m_score_val = {}
        self.m_score_test = {}
        self.m_score_test_sig = {}
        self.m_score_test_bkg = {}
        self.m_tsf = {}

        self.inputTree = 'inclusive'
        self.train_signal = []
        self.train_data_background = []
        self.train_mc_background = []
        self.train_dd_background = []
        self.train_variables = []
        self.preselections = [] 
        self.mc_preselections = []
        self.data_preselections = []
        self.signal_preselections = []
        self.background_preselections = []
        self.randomIndex = 'eventNumber'
        self.weight = 'weight'
        self.params = [{'eval_metric': ['auc', 'logloss']}]
        self.early_stopping_rounds = 10
        self.numRound = 10000
        self.SF = -1
        self.readConfig(configPath)
        self.checkConfig()

    def readConfig(self, configPath):
        """Read configuration file formated in json to extract information to fill TemplateMaker variables."""
        try:
            member_variables = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("_") and not attr.startswith('m_')]

            stream = open(configPath, 'r')
            configs = json.loads(stream.read())

            # read from the common settings
            config = configs["common"]
            # if self._add >= 0:
            #     config["train_variables"].append(config["+train_variables"][self._add])

            for member in config.keys():
                if member in member_variables:
                    setattr(self, member, config[member])


            # read from the region specific settings
            if self._region:
                config = configs[self._region]
                if self._add >= 0:
                    config["+train_variables"].append(config["test_variables"][self._add])
                for member in config.keys():
                    if member in member_variables:
                        setattr(self, member, config[member])
                if '+train_mc_background' in config.keys():
                    self.train_mc_background += config['+train_mc_background']
                if '+train_signal' in config.keys():
                    self.train_signal += config['+train_signal']
                if '+train_variables' in config.keys():
                    self.train_variables += config['+train_variables']
                if '+preselections' in config.keys():
                    self.preselections += config['+preselections']
                if '+signal_preselections' in config.keys():
                    self.signal_preselections += config['+signal_preselections']
                if '+background_preselections' in config.keys():
                    self.background_preselections += config['+background_preselections']

            if self._shield >= 0:
                self.train_variables.pop(self._shield)

            print("\n\n")
            print(self.train_variables)
            print(len(self.train_variables))
            print("\n\n")

            self._branches = list( set(self.train_variables) | set([p.split()[0].replace("(", "") for p in self.preselections]) | set([p.split()[0].replace("(", "") for p in self.background_preselections]) | set([self.randomIndex, self.weight]))

            self._mc_branches = list( set(self.train_variables) | set([p.split()[0].replace("(", "") for p in self.preselections]) | set([p.split()[0].replace("(", "") for p in self.mc_preselections]) | set([p.split()[0].replace("(", "") for p in self.signal_preselections]) | set([p.split()[0].replace("(", "") for p in self.background_preselections]) | set([self.randomIndex, self.weight]))

            self._data_branches = list( set(self.train_variables) | set([p.split()[0].replace("(", "") for p in self.preselections]) | set([p.split()[0].replace("(", "") for p in self.data_preselections]) | set([p.split()[0].replace("(", "") for p in self.signal_preselections]) | set([p.split()[0].replace("(", "") for p in self.background_preselections]) | set([self.randomIndex, self.weight]))

            self._sig_branches = list( set(self.train_variables) | set([p.split()[0].replace("(", "") for p in self.preselections]) | set([p.split()[0].replace("(", "") for p in self.mc_preselections]) | set([p.split()[0].replace("(", "") for p in self.signal_preselections]) | set([self.randomIndex, self.weight]))



            self.train_variables = [x.replace('noexpand:', '') for x in self.train_variables]
            self.preselections = [x.replace('noexpand:', '') for x in self.preselections]
            self.mc_preselections = [x.replace('noexpand:', '') for x in self.mc_preselections]
            self.data_preselections = [x.replace('noexpand:', '') for x in self.data_preselections]
            self.signal_preselections = [x.replace('noexpand:', '') for x in self.signal_preselections]
            self.background_preselections = [x.replace('noexpand:', '') for x in self.background_preselections]
            self.randomIndex = self.randomIndex.replace('noexpand:', '')
            self.weight = self.weight.replace('noexpand:', '')

            # if self.preselections:
            #     self.preselections = ['data.' + p for p in self.preselections]
            # if self.signal_preselections:
            #     self.signal_preselections = ['data.' + p for p in self.signal_preselections]
            # if self.mc_preselections:
            #     self.mc_preselections = ['data.' + p for p in self.mc_preselections]
            # if self.data_preselections:
            #     self.data_preselections = ['data.' + p for p in self.data_preselections]
            # if self.background_preselections:
            #     self.background_preselections = ['data.' + p for p in self.background_preselections]

        except Exception as e:
            logging.error("Error reading configuration '{config}'".format(config=configPath))
            logging.error(e)

    def checkConfig(self):
        if not self.train_signal: print('ERROR: no training signal!!')
        if not self.train_data_background: print('ERROR: no data side band training background!!')
        if not self.train_mc_background: print('ERROR: no MC training background!!')
        if not self.train_dd_background: print('ERROR: no data-driven training background!!')
        if not self.train_variables: print('ERROR: no training variables!!')

    def setParams(self, params, fold=-1):

        print('XGB INFO: setting hyperparameters...')
        if fold == -1:
            self.params = [{'eval_metric': ['auc', 'logloss']}]
            fold = 0
        else:
            self.params[fold] = {'eval_metric': ['auc', 'logloss']}
        for key in params:
            self.params[fold][key] = params[key]

    def set_early_stopping_rounds(self, rounds):
        self.early_stopping_rounds = rounds

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def preselect(self, data, sample=''):
        if sample == 'signal':
            for p in self.signal_preselections:
                data = data.query(p)
            for p in self.mc_preselections:
                data = data.query(p)
        elif sample == 'background':
            for p in self.background_preselections:
                data = data.query(p)
        elif sample == 'data':
            for p in self.background_preselections:
                data = data.query(p)
            for p in self.data_preselections:
                data = data.query(p)
        elif sample == 'mc_background':
            for p in self.background_preselections:
                data = data.query(p)
            for p in self.mc_preselections:
                data = data.query(p)
        elif sample == 'DYJetsToLL':
            # data = data.query('n_iso_photons == 0')
            for p in self.background_preselections:
                data = data.query(p)
            for p in self.mc_preselections:
                data = data.query(p)
        for p in self.preselections:
            data = data.query(p)

        return data

    def readData(self):

        sig_list, bkg_mc_list, bkg_dd_list, bkg_data_list = [], [], [], []
        for sig_cat in self.train_signal:
            sig_cat_folder = self._inputFolder + '/' + sig_cat
            for sig in os.listdir(sig_cat_folder):
                if sig.endswith('.root'): sig_list.append(sig_cat_folder + '/' + sig)
        for bkg_cat in self.train_dd_background:
            bkg_cat_folder = self._inputFolder + '/' + bkg_cat
            for bkg in os.listdir(bkg_cat_folder):
                if bkg.endswith('.root'): bkg_dd_list.append(bkg_cat_folder + '/' + bkg)
        for bkg_cat in self.train_mc_background:
            bkg_cat_folder = self._inputFolder + '/' + bkg_cat
            for bkg in os.listdir(bkg_cat_folder):
                if bkg.endswith('.root'): bkg_mc_list.append(bkg_cat_folder + '/' + bkg)
        for bkg_cat in self.train_data_background:
            bkg_cat_folder = self._inputFolder + '/' + bkg_cat
            for bkg in os.listdir(bkg_cat_folder):
                if bkg.endswith('.root'): bkg_data_list.append(bkg_cat_folder + '/' + bkg)

        print('-------------------------------------------------')
        for sig in sig_list: print('XGB INFO: Adding signal sample: ', sig)
        #TODO put this to the config
        for filename in tqdm(sorted(sig_list), desc='XGB INFO: Loading training signals', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
            file = uproot.open(filename)
            for data in file[self.inputTree].iterate(self._sig_branches, library='pd', step_size=self._chunksize):
                data = self.preselect(data, 'signal')
                self.m_data_sig = pd.concat([self.m_data_sig, data], ignore_index=True)

        print('----------------------------------------------------------')
        if bkg_mc_list:
            for bkg in bkg_mc_list: print('XGB INFO: Adding mc background sample: ', bkg)
            for bkg in tqdm(sorted(bkg_mc_list), desc='XGB INFO: Loading training backgrounds', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
                #TODO put this to the config
                if "DY" in bkg:
                    branches = self._mc_branches #+ ["n_iso_photons"]
                else: 
                    branches = self._mc_branches
                file = uproot.open(bkg)
                for data in file[self.inputTree].iterate(branches, library='pd', step_size=self._chunksize):
                    if "DY" in bkg:
                        data = self.preselect(data, 'DYJetsToLL')
                    else:
                        data = self.preselect(data, 'mc_background')
                    self.m_data_bkg = pd.concat([self.m_data_bkg, data], ignore_index=True)
        
        print('----------------------------------------------------------')
        if bkg_dd_list:
            for bkg in bkg_dd_list: print('XGB INFO: Adding Data-driven background sample: ', bkg)
            #TODO put this to the config
            for filename in tqdm(sorted(bkg_dd_list), desc='XGB INFO: Loading training backgrounds', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
                file = uproot.open(filename)
                for data in file[self.inputTree].iterate(self._branches, library='pd', step_size=self._chunksize):
                    data = self.preselect(data, 'background')
                    self.m_data_bkg = pd.concat([self.m_data_bkg, data], ignore_index=True)

        print('----------------------------------------------------------')
        if bkg_data_list:
            for bkg in bkg_data_list: print('XGB INFO: Adding data side band background sample: ', bkg)
            #TODO put this to the config
            for filename in tqdm(sorted(bkg_data_list), desc='XGB INFO: Loading training backgrounds', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
                file = uproot.open(filename)
                for data in file[self.inputTree].iterate(self._data_branches, library='pd', step_size=self._chunksize):
                    data = self.preselect(data, 'data')
                    self.m_data_bkg = pd.concat([self.m_data_bkg, data], ignore_index=True)

    def plot_corr(self, data_type):
        if not os.path.isdir("plots/corr/"):
            os.makedirs("plots/corr/")
        if data_type == "sig":
            data = self.m_data_sig
        else:
            data = self.m_data_bkg
        columns = list(data.columns)
        for rem in ("is_center", "weight", "event", "gamma_mvaID_WP80", "gamma_mvaID_WPL", "n_iso_photons", "n_b_jets", "n_jets"):
            if rem in columns:
                columns.remove(rem)
                data = data.drop(rem, axis=1)
        print(columns)

        data = data.corr() * 100
        # data = data.dropna(axis=0, how='all').dropna(axis=1, how='all')

        # 将相关性矩阵转换为按相关性值排序的数据框
        sorted_corr = data.unstack().sort_values(ascending=False)

        # 选择相关性大于0.4的变量对
        strong_corr = sorted_corr[(sorted_corr.abs() > 40) & (sorted_corr < 100)]
        strong_corr = strong_corr[strong_corr.index.map(lambda x: x[0] < x[1])]
        max_var_length = max([len(var) for var in columns])
        # 输出结果
        for (var1, var2), corr in strong_corr.items():
            print(f"Correlation between {var1:<{max_var_length}} and {var2:<{max_var_length}} is {corr:.2f}")
        # print(data)

        plt.figure(figsize=(12, 9), dpi=300)
        plt.imshow(data)
        plt.colorbar()
        for i in range(0, len(columns)):
            line = list(data[columns[i]])
            for j in range(0, len(line)):
                plt.text(i, j, int(line[j]), verticalalignment='center', horizontalalignment='center', fontsize=8)
        plt.xticks(np.arange(0, len(columns)), columns, rotation=-90, fontsize=12)
        plt.yticks(np.arange(0, len(columns)), columns, fontsize=12)
        plt.title(self._region)
        plt.tight_layout()
        plt.savefig("plots/corr/corr_%s_%s.pdf" % (self._region, data_type))

    def reweightSignal(self):
        min_mass, max_mass = 100, 180
        if (self.m_data_bkg["H_mass"] > min_mass).all() and (self.m_data_bkg["H_mass"] < max_mass).all():
            n_bins = 100
            bin_edges = np.linspace(min_mass, max_mass, n_bins + 1)
            
            sig_hist, _ = np.histogram(self.m_data_sig['H_mass'], bins=bin_edges, weights=self.m_data_sig['weight'], density=True)
            bkg_hist, _ = np.histogram(self.m_data_bkg['H_mass'], bins=bin_edges, weights=self.m_data_bkg['weight'], density=True)

            # calculate ratio of signal to background
            ratio_hist = np.divide(bkg_hist, sig_hist, out=np.zeros_like(bkg_hist), where=sig_hist != 0)

            # refreche the weights
            bin_indices = np.digitize(self.m_data_sig['H_mass'], bins=bin_edges) - 1
            self.m_data_sig['weight'] = self.m_data_sig['weight'] * np.where((bin_indices >= 0) & (bin_indices < n_bins), ratio_hist[bin_indices], 1)


    def prepareData(self, fold=0):
        # training, validation, test split
        print('----------------------------------------------------------')
        print('XGB INFO: Splitting samples to training, validation, and test...')
        test_sig = self.m_data_sig[self.m_data_sig[self.randomIndex]%4 == fold]
        test_bkg = self.m_data_bkg[self.m_data_bkg[self.randomIndex]%4 == fold]
        val_sig = self.m_data_sig[(self.m_data_sig[self.randomIndex]-1)%4 == fold]
        val_bkg = self.m_data_bkg[(self.m_data_bkg[self.randomIndex]-1)%4 == fold]
        train_sig = self.m_data_sig[((self.m_data_sig[self.randomIndex]-2)%4 == fold) | ((self.m_data_sig[self.randomIndex]-3)%4 == fold)]
        train_bkg = self.m_data_bkg[((self.m_data_bkg[self.randomIndex]-2)%4 == fold) | ((self.m_data_bkg[self.randomIndex]-3)%4 == fold)]

        headers = ['Sample', 'Total', 'Training', 'Validation']
        sample_size_table = [
            ['Signal'    , len(train_sig)+len(val_sig), len(train_sig), len(val_sig)],
            ['Background', len(train_bkg)+len(val_bkg), len(train_bkg), len(val_bkg)],
        ]

        print(tabulate(sample_size_table, headers=headers, tablefmt='simple'))

        # setup the training data set
        print('XGB INFO: Setting the training arrays...')
        x_test_sig = test_sig[self.train_variables]
        x_test_bkg = test_bkg[self.train_variables]
        x_val_sig = val_sig[self.train_variables]
        x_val_bkg = val_bkg[self.train_variables]
        x_train_sig = train_sig[self.train_variables]
        x_train_bkg = train_bkg[self.train_variables]

        x_train = pd.concat([x_train_sig, x_train_bkg])
        x_val = pd.concat([x_val_sig, x_val_bkg])
        x_test = pd.concat([x_test_sig, x_test_bkg])

        # setup the weights
        print('XGB INFO: Setting the event weights...')

        if self.SF == -1: self.SF = 1.*train_sig.shape[0]/train_bkg.shape[0]
        train_sig_wt = train_sig[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_sig[self.weight].mean() * (1.+self.SF) * train_sig.shape[0] )
        train_bkg_wt = train_bkg[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_bkg[self.weight].mean() * (1.+self.SF) * train_bkg.shape[0] )
        val_sig_wt = val_sig[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_sig[self.weight].mean() * (1.+self.SF) * train_sig.shape[0] )
        val_bkg_wt = val_bkg[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_bkg[self.weight].mean() * (1.+self.SF) * train_bkg.shape[0] )
        test_sig_wt = test_sig[[self.weight]]
        test_bkg_wt = test_bkg[[self.weight]]

        self.m_train_wt[fold] = pd.concat([train_sig_wt, train_bkg_wt]).to_numpy()
        self.m_train_wt[fold][self.m_train_wt[fold] < 0] = 0
        self.m_val_wt[fold] = pd.concat([val_sig_wt, val_bkg_wt]).to_numpy()
        self.m_val_wt[fold][self.m_val_wt[fold] < 0] = 0
        self.m_test_wt[fold] = pd.concat([test_sig_wt, test_bkg_wt]).to_numpy()
        self.m_test_wt[fold][self.m_test_wt[fold] < 0] = 0

        # setup the truth labels
        print('XGB INFO: Signal labeled as one; background labeled as zero.')
        self.m_y_train[fold] = np.concatenate((np.ones(len(train_sig), dtype=np.uint8), np.zeros(len(train_bkg), dtype=np.uint8)))
        self.m_y_val[fold]   = np.concatenate((np.ones(len(val_sig)  , dtype=np.uint8), np.zeros(len(val_bkg)  , dtype=np.uint8)))
        self.m_y_test[fold]  = np.concatenate((np.ones(len(test_sig)  , dtype=np.uint8), np.zeros(len(test_bkg)  , dtype=np.uint8)))
        
        # construct DMatrix
        print('XGB INFO: Constucting D-Matrix...')
        self.m_dTrain[fold] = xgb.DMatrix(x_train, label=self.m_y_train[fold], weight=self.m_train_wt[fold])
        self.m_dVal[fold] = xgb.DMatrix(x_val, label=self.m_y_val[fold], weight=self.m_val_wt[fold])
        self.m_dTest[fold] = xgb.DMatrix(x_test)
        self.m_dTest_sig[fold] = xgb.DMatrix(x_test_sig)
        self.m_dTest_bkg[fold] = xgb.DMatrix(x_test_bkg)

    def trainModel(self, fold, param=None):

        if not param:
            print('XGB INFO: Setting the hyperparameters...')
            param = self.params[0 if len(self.params) == 1 else fold]
        print('param: ', param)

        # finally start training!!!
        print('XGB INFO: Start training!!!')
        evallist  = [(self.m_dTrain[fold], 'train'), (self.m_dVal[fold], 'eval')]
        evals_result = {}
        eval_result_history = []
        try:
            self.m_bst[fold] = xgb.train(param, self.m_dTrain[fold], self.numRound, evals=evallist, early_stopping_rounds=self.early_stopping_rounds, evals_result=evals_result)
        except KeyboardInterrupt:
            print('Finishing on SIGINT.')

    def testModel(self, fold):
        # get scores
        print('XGB INFO: Computing scores for different sample sets...')
        self.m_score_val[fold]   = self.m_bst[fold].predict(self.m_dVal[fold])
        self.m_score_test[fold]  = self.m_bst[fold].predict(self.m_dTest[fold])
        self.m_score_train[fold] = self.m_bst[fold].predict(self.m_dTrain[fold])
        self.m_score_test_sig[fold]  = self.m_bst[fold].predict(self.m_dTest_sig[fold])
        self.m_score_test_bkg[fold]  = self.m_bst[fold].predict(self.m_dTest_bkg[fold])

    def plotScore(self, fold=0, sample_set='val'):

        if sample_set == 'train':
            plt.hist(self.m_score_train[fold], bins='auto')
        elif sample_set == 'val':
            plt.hist(self.m_score_val[fold], bins='auto')
        elif sample_set == 'test':
            plt.hist(self.m_score_test[fold], bins='auto')
        plt.title('Score distribution %s' % sample_set)
        plt.show()
        if sample_set == 'test':
            plt.hist(self.m_score_test_sig[fold], bins='auto')
            plt.title('Score distribution test signal')
            plt.show()
            plt.hist(self.m_score_test_bkg[fold], bins='auto')
            plt.title('Score distribution test background')
            plt.show()

    def plotFeaturesImportance(self, fold=0, save=True, show=True, type='weight'):
        """Plot feature importance. Type can be 'weight', 'gain' or 'cover'"""

        xgb.plot_importance(booster=self.m_bst[fold], importance_type=type, show_values=show, values_format='{v:.2f}')
        plt.tight_layout()

        if save:
            # create output directory
            if not os.path.isdir('plots/feature_importance'):
                os.makedirs('plots/feature_importance')
            # save figure
            plt.savefig('plots/feature_importance/%d_BDT_%s_%d_%s.png' % (self._shield+1, self._region, fold, type))

        if show: plt.show()

        plt.clf()
        return

    def plotROC(self, fold=0, save=True, show=False):
    
        fpr_train, tpr_train, _ = roc_curve(self.m_y_train[fold], self.m_score_train[fold], sample_weight=self.m_train_wt[fold])
        tpr_train, fpr_train = np.array(list(zip(*sorted(zip(tpr_train, fpr_train)))))
        roc_auc_train = 1 - auc(tpr_train, fpr_train)
        fpr_val, tpr_val, _ = roc_curve(self.m_y_val[fold], self.m_score_val[fold], sample_weight=self.m_val_wt[fold])
        tpr_val, fpr_val = np.array(list(zip(*sorted(zip(tpr_val, fpr_val)))))
        roc_auc_val = 1 - auc(tpr_val, fpr_val)
        fpr_test, tpr_test, _ = roc_curve(self.m_y_test[fold], self.m_score_test[fold], sample_weight=self.m_test_wt[fold])
        tpr_test, fpr_test = np.array(list(zip(*sorted(zip(tpr_test, fpr_test)))))
        roc_auc_test = 1 - auc(tpr_test, fpr_test)

        fnr_train = 1.0 - fpr_train
        fnr_val = 1.0 - fpr_val
        fnr_test = 1.0 - fpr_test

        plt.grid(color='gray', linestyle='--', linewidth=1)
        plt.plot(tpr_train, fnr_train, label='Train set, area = %0.6f' % roc_auc_train, color='black', linestyle='dotted')
        plt.plot(tpr_val, fnr_val, label='Val set, area = %0.6f' % roc_auc_val, color='blue', linestyle='dashdot')
        plt.plot(tpr_test, fnr_test, label='Test set, area = %0.6f' % roc_auc_test, color='red', linestyle='dashed')
        plt.plot([0, 1], [1, 0], linestyle='--', color='black', label='Luck')
        plt.xlabel('Signal acceptance')
        plt.ylabel('Background rejection')
        plt.title('Receiver operating characteristic')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xticks(np.arange(0, 1, 0.1))
        plt.yticks(np.arange(0, 1, 0.1))
        plt.legend(loc='lower left', framealpha=1.0)
        plt.tight_layout()

        if save:
            # create output directory
            if not os.path.isdir('plots/roc_curve'):
                os.makedirs('plots/roc_curve')
            # save figure
            plt.savefig('plots/roc_curve/%d_BDT_%s_%d.pdf' % (self._shield+1, self._region, fold))

        if show: plt.show()

        plt.clf()

        return

    def getAUC(self, fold=0, sample_set='val'):

        if sample_set == 'train':
            fpr, tpr, _ = roc_curve(self.m_y_train[fold], self.m_score_train[fold], sample_weight=self.m_train_wt[fold])
        elif sample_set == 'val':
            fpr, tpr, _ = roc_curve(self.m_y_val[fold], self.m_score_val[fold], sample_weight=self.m_val_wt[fold])
        elif sample_set == 'test':
            fpr, tpr, _ = roc_curve(self.m_y_test[fold], self.m_score_test[fold], sample_weight=self.m_test_wt[fold])

        tpr, fpr = np.array(list(zip(*sorted(zip(tpr, fpr)))))
        roc_auc = 1 - auc(tpr, fpr)

        return self.params[0 if len(self.params) == 1 else fold], roc_auc

    def transformScore(self, fold=0, sample='sig'):

        print(f'XGB INFO: transforming scores based on {sample}')
        # transform the scores
        self.m_tsf[fold] = QuantileTransformer(n_quantiles=1000, output_distribution='uniform', subsample=1000000000, random_state=0)
        #plt.hist(score_test_sig, bins='auto')
        #plt.show()
        if sample == 'sig': self.m_tsf[fold].fit(self.m_score_test_sig[fold].reshape(-1, 1))
        elif sample == 'bkg': self.m_tsf[fold].fit(self.m_score_test_bkg[fold].reshape(-1, 1))
        #score_test_sig_t=tsf.transform(self.m_score_test_sig[fold].reshape(-1, 1)).reshape(-1)
        #plt.hist(score_test_sig_t, bins='auto')
        #plt.show()

    def save(self, fold=0):

        # create output directory
        if not os.path.isdir(self._outputFolder):
            os.makedirs(self._outputFolder)

        # save BDT model
        if fold in self.m_bst.keys(): self.m_bst[fold].save_model('%s/BDT_%s_%d.h5' % (self._outputFolder, self._region, fold))

        # save score transformer
        if fold in self.m_tsf.keys(): 
            with open('%s/BDT_tsf_%s_%d.pkl' %(self._outputFolder, self._region, fold), 'wb') as f:
                pickle.dump(self.m_tsf[fold], f, -1)

    def optunaHP(self, fold=0):
        import optuna
        from xgboost import XGBClassifier
        import optuna.visualization as vis

        def objective(trial):
            params = {
                'max_depth': trial.suggest_int('max_depth', 3, 50),
                'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.05, log=True),
                'subsample': trial.suggest_float('subsample', 0.4, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.4, 1.0),
                'gamma': trial.suggest_float('gamma', 1e-6, 10, log=True),
                'min_child_weight': trial.suggest_int('min_child_weight', 1, 50),
                'max_delta_step': trial.suggest_int('max_delta_step', 10, 100),
            }
            
            self.setParams(params, fold)
            self.trainModel(fold, self.params[fold])
            self.testModel(fold)
            eval_auc = self.getAUC(fold)[-1]
            train_auc = self.getAUC(fold, 'train')[-1]
            
            metrics = {
                'eval_auc': eval_auc,
                'sqrt_eval_auc_minus_train_auc': np.sqrt(train_auc * (2 * eval_auc - train_auc)),
                'eval_auc_minus_train_auc': eval_auc * 2 - train_auc,
                'eval_auc_over_train_auc': eval_auc ** 2 / ((eval_auc + train_auc) / 2)
            }
            return metrics.get(self.optuna_metric, eval_auc)
        
        exp_dir = f'models/optuna_{self._region}/'
        if not os.path.exists(exp_dir):
            os.makedirs(exp_dir)

        study_name = f"study_fold_{fold}"  # Unique study name
        storage_name = f"sqlite:///{exp_dir}optuna_study_fold_{fold}.sqlite3"  # Storage name for SQLite
        if not self.continue_optuna and os.path.exists(f"{exp_dir}optuna_study_fold_{fold}.sqlite3"):
            os.remove(f"{exp_dir}optuna_study_fold_{fold}.sqlite3")
        self.study = optuna.create_study(
                sampler=optuna.samplers.TPESampler(n_startup_trials=10, multivariate=True, group=True),
                study_name=study_name, 
                storage=storage_name, 
                load_if_exists=True, 
                direction='maximize'
            )
        self.study.optimize(objective, n_trials=self.n_calls)
        
        best_params_path = f"{exp_dir}BDT_region_{self._region}_fold{fold}.json"
     
        self.setParams(self.study.best_params, fold)
        with open(best_params_path, 'w') as f:
            json.dump(self.params[fold], f)

def main():

    args=getArgs()
    
    configPath = args.config
    xgb_model = XGBoostHandler(configPath, args.region)

    if args.inputFolder: xgb_model.setInputFolder(args.inputFolder)
    if args.outputFolder: xgb_model.setOutputFolder(args.outputFolder)
    if args.params: xgb_model.setParams(args.params)

    if args.hyperparams_path:
        for fold in range(4):
            path = "{}/BDT_region_{}_fold{}.json".format(args.hyperparams_path, args.region, fold)
            stream = open(path, 'r')
            xgb_model.params[fold] = json.load(stream)
            xgb_model.setParams(xgb_model.params[fold], fold)
        # xgb_model.setParams(hyperparameters)
        
    xgb_model.readData()

    if args.corr:
        xgb_model.plot_corr("sig")
        xgb_model.plot_corr("bkg")

    # xgb_model.reweightSignal()

    # looping over the "4 folds"
    for i in args.fold:

        print('===================================================')
        print(f' The #{i} fold of training')
        print('===================================================')

        #xgb_model.setParams({'eval_metric': ['auc', 'logloss']}, i)
        xgb_model.set_early_stopping_rounds(20)

        xgb_model.prepareData(i)

        if args.optuna:
            try:
                xgb_model.optunaHP(i)
            except KeyboardInterrupt:
                print('Finishing on SIGINT.')
            continue

        xgb_model.trainModel(i)

        xgb_model.testModel(i)
        print("param: %s, Val AUC: %f" % xgb_model.getAUC(i))

        #xgb_model.plotScore(i, 'test')
        if args.importance:
            xgb_model.plotFeaturesImportance(i, type='weight')
            xgb_model.plotFeaturesImportance(i, type='gain')
        if args.roc:
            xgb_model.plotROC(i)

        xgb_model.transformScore(i)

        if args.save: xgb_model.save(i)

    print('------------------------------------------------------------------------------')
    print('Finished training.')

    return

if __name__ == '__main__':
    main()
