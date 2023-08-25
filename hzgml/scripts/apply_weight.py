#!/usr/bin/env python
import os
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd
from root_pandas import *
import pickle
from tqdm import tqdm
import logging
import uproot
import awkward as ak
from pdb import set_trace
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError + 1
pd.options.mode.chained_assignment = None


def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-i', '--inputFolder', action='store', default='skimmed_ntuples', help='directory of training inputs')
    parser.add_argument('-o', '--outputFolder', action='store', default='skimmed_ntuples', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH', 'inclusive'], default='zero_jet', help='Region to process')

    return parser.parse_args()


class ApplyWeightHandler(object):
    "Class for applying XGBoost"

    def __init__(self, region=''):

        print('==================================')
        print('  ApplyWeigthHandler initialized  ')
        print('==================================')

        args=getArgs()

        self._region = region if region else 'inclusive'
        self._inputFolder = ''
        self._inputTree = region if region else 'inclusive'
        self._outputFolder = ''
        self._branches = []
        self._outbranches = []

        self.data_med = pd.DataFrame()
        self.data_fake = pd.DataFrame()
        self.mc_true = pd.DataFrame()
        self.mc_med = pd.DataFrame()

        self.pt_bin = {
            "ele":{'isBarrel':[15, 17, 20, 25, 30, 40, 5000],'isEndcap':[15, 17, 20, 25, 30, 40, 5000]},
            'mu':{'isBarrel':[15, 17, 20, 25, 30, 40, 5000],'isEndcap':[15, 17, 20, 25, 30, 40, 5000]}
            }
        self.fake_frac = {
            "ele":{'isBarrel':[0.586, 0.450, 0.365, 0.256, 0.158, 0.180],'isEndcap':[0.594, 0.417, 0.454, 0.369, 0.242, 0.113]},
            "mu":{'isBarrel':[0.527, 0.465, 0.346, 0.308, 0.165, 0.106],'isEndcap':[0.646, 0.537, 0.455, 0.259, 0.219, 0.145]}
            }#zpt
        self.fake_frac_err = {
            'ele':{'isBarrel':[0.019, 0.016, 0.015, 0.021, 0.019, 0.024],'isEndcap':[0.034, 0.029, 0.026, 0.032, 0.031, 0.036]},
            'mu':{'isBarrel':[0.016, 0.014, 0.013, 0.018, 0.016, 0.015],'isEndcap':[0.028, 0.027, 0.021, 0.026, 0.027, 0.029]}
            }#zpt
        self.fakeWeight = {
            "ele":{"isBarrel":{"w": [], "w_err":[], "w_down": []}, "isEndcap":{"w": [], "w_err":[], "w_down": []}},
            "mu":{"isBarrel":{"w": [], "w_err":[], "w_down": []}, "isEndcap":{"w": [], "w_err":[], "w_down": []}}
            }
        self.trueWeight = {
            "ele":{"isBarrel":{"w": [], "w_err":[], "w_down": []}, "isEndcap":{"w": [], "w_err":[], "w_down": []}},
            "mu":{"isBarrel":{"w": [], "w_err":[], "w_down": []}, "isEndcap":{"w": [], "w_err":[], "w_down": []}}
            }

        self.setInputFolder(args.inputFolder)
        self.setOutputFolder(args.outputFolder)

        print('===============================')
        print('     Finish initialization     ')
        print('===============================')

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def eta_sel(self, data, type):
        if type == 1:
            # return (data.loc[:, 'gamma_eta']>-1.4442) & (data.loc[:, 'gamma_eta']<1.4442)
            return (data.loc[:, 'gamma_eta']>-1.5) & (data.loc[:, 'gamma_eta']<1.5)
        if type == 0:
            # return ((data.loc[:, 'gamma_eta']>-2.5) & (data.loc[:, 'gamma_eta']<-1.566)) | ((data.loc[:, 'gamma_eta']>1.566) & (data.loc[:, 'gamma_eta']<2.5))
            return (data.loc[:, 'gamma_eta']<-1.5) | (data.loc[:, 'gamma_eta']>1.5)

    def readRoot(self, type, mc_tagger=0, sb_tagger=0, tagger=0):
        input_path = self._inputFolder + '/%s/2017.root' % type
        print("Reading file: %s" % input_path)
        data = uproot.open(input_path+":inclusive").arrays(library='pd')
        
        cut_com = (data.loc[:, 'Z_sublead_lepton_pt'] > 15)
        if "DY" in type:
            cut_com = cut_com & (data.loc[:, 'n_iso_photons'] == 0)
        if sb_tagger:
            cut = cut_com & ((self.eta_sel(data, 1) & (data.loc[:,'gamma_mvaID'] < 0.42) & (data.loc[:,'gamma_mvaID'] > -0.97) & (data['gamma_chiso']*data['gamma_pt'] < 1.141)) | (self.eta_sel(data, 0) & (data.loc[:,'gamma_mvaID'] < 0.14) & (data.loc[:,'gamma_mvaID'] > -0.97) & (data['gamma_chiso']*data['gamma_pt'] < 1.051)))
        else:
            cut = cut_com & ((self.eta_sel(data, 1) & (data.loc[:,'gamma_mvaID'] > 0.42)) | (self.eta_sel(data, 0) & (data.loc[:,'gamma_mvaID'] > 0.14)))

        if mc_tagger:
            cut = cut & (data.loc[:,'gamma_genPartFlav'] == 1)
        data_sel = data.loc[cut].copy().reset_index(drop=True)
        data_sel['tagger'] = (data_sel['gamma_pt']>0)*tagger

        if mc_tagger and sb_tagger:
            self.mc_true = pd.concat([self.mc_true, data_sel], ignore_index=True, sort=False)
        if not sb_tagger and not mc_tagger:
            self.data_med = pd.concat([self.data_med, data_sel], ignore_index=True, sort=False)
        if not mc_tagger and sb_tagger:
            self.data_fake = pd.concat([self.data_fake, data_sel], ignore_index=True, sort=False)
        if mc_tagger and not sb_tagger:
            self.mc_med = pd.concat([self.mc_med, data_sel], ignore_index=True, sort=False)

        del data_sel, data, cut_com

    def storeOutput(self, data, path):
        print("Storing ", path)
        data.to_root(path, key='inclusive', index=False)
        data_zero_jet = data[(data.n_jets == 0) & (data.n_leptons < 3)]
        data_one_jet = data[(data.n_jets == 1) & (data.n_leptons < 3)]
        data_two_jet = data[(data.n_jets >= 2) & (data.n_leptons < 3)]
        data_zero_to_one_jet = data[(data.n_jets < 2) & (data.n_leptons < 3)]
        data_VH_ttH =  data[data.n_leptons >= 3]
        data_zero_jet.to_root(path, key='zero_jet', mode='a', index=False)
        data_one_jet.to_root(path, key='one_jet', mode='a', index=False)
        data_zero_to_one_jet.to_root(path, key='zero_to_one_jet', mode='a', index=False)
        data_two_jet.to_root(path, key='two_jet', mode='a', index=False)
        data_VH_ttH.to_root(path, key='VH_ttH', mode='a', index=False)

        del data_zero_jet, data_one_jet, data_two_jet, data_zero_to_one_jet, data_VH_ttH

    def readFiles(self):
        self.readRoot('data',0,0,0)
        self.readRoot('data',0,1,0)
        for tagger, type in enumerate(['ZGToLLG','ZG2JToG2L2J', 'DYJetsToLL','TTGJets','TGJets']):
            self.readRoot(type,1,1,tagger+1)
            self.readRoot(type,1,0,tagger+1)

        self.storeOutput(self.data_med, "skimmed_ntuples/data_med/2017.root")

        print('==================================')
        print('Finish data reading and preprocess')
        print('==================================')

    def calculateWeight(self):
        for part in ['ele', 'mu']:
            if part == 'ele':
                part_cut = 11
            else: 
                part_cut = 13
            for isbarrel in ['isBarrel', 'isEndcap']:
                if isbarrel == 'isBarrel':
                    eta_cut = 1
                else: 
                    eta_cut = 0
                for i in range(len(self.pt_bin[part][isbarrel])-1):
                    print(part, " ", isbarrel, " ", self.pt_bin[part][isbarrel][i], " to ", self.pt_bin[part][isbarrel][i+1])
                    data_med_cut = self.eta_sel(self.data_med, eta_cut) & (self.data_med.loc[:,'Z_lead_lepton_id']==part_cut) & (self.data_med.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (self.data_med.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1]) #& ((self.data_med.loc[:,'H_mass']<120) | (self.data_med.loc[:,'H_mass']>130))
                    data_med_sum = self.data_med.loc[data_med_cut, 'weight'].sum() 
                    data_med_sum_err = (self.data_med.loc[data_med_cut, 'weight']*self.data_med.loc[data_med_cut, 'weight']).sum()

                    data_fake_cut = self.eta_sel(self.data_fake, eta_cut) & (self.data_fake.loc[:,'Z_lead_lepton_id']==part_cut) & (self.data_fake.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (self.data_fake.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])# & ((self.data_fake.loc[:,'H_mass']<120) | (self.data_fake.loc[:,'H_mass']>130))
                    data_fake_sum = self.data_fake.loc[data_fake_cut, 'weight'].sum() 
                    data_fake_sum_err = (self.data_fake.loc[data_fake_cut, 'weight']*self.data_fake.loc[data_fake_cut, 'weight']).sum()

                    mc_true_cut = self.eta_sel(self.mc_true, eta_cut) & (self.mc_true.loc[:,'Z_lead_lepton_id']==part_cut) & (self.mc_true.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (self.mc_true.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1]) #& ((self.mc_true.loc[:,'H_mass']<120) | (self.mc_true.loc[:,'H_mass']>130))
                    mc_true_sum = self.mc_true.loc[mc_true_cut,'weight'].sum()
                    mc_true_sum_err = (self.mc_true.loc[mc_true_cut,'weight']*self.mc_true.loc[mc_true_cut,'weight']).sum()

                    mc_med_cut = self.eta_sel(self.mc_med, eta_cut) & (self.mc_med.loc[:,'Z_lead_lepton_id']==part_cut) & (self.mc_med.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (self.mc_med.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1]) #& ((self.mc_med.loc[:,'H_mass']<120) | (self.mc_med.loc[:,'H_mass']>130)) 
                    mc_med_sum = self.mc_med.loc[mc_med_cut,'weight'].sum()
                    mc_med_sum_err = (self.mc_med.loc[mc_med_cut,'weight']*self.mc_med.loc[mc_med_cut,'weight']).sum()

                    fakeYield = data_fake_sum - mc_true_sum
                    fakeYield_err = data_fake_sum_err + mc_true_sum_err
                    trueWeight = (1. - self.fake_frac[part][isbarrel][i])
                    # self.weight[part][isbarrel]['w_up'].append(ratio * (self.fake_frac[part][isbarrel][i] + self.fake_frac_err[part][isbarrel][i]))
                    # self.weight[part][isbarrel]['w_down'].append(ratio * (self.fake_frac[part][isbarrel][i] - self.fake_frac_err[part][isbarrel][i]))
                    fw = data_med_sum / fakeYield * self.fake_frac[part][isbarrel][i]
                    fw_err = np.sqrt(
                        # (self.fake_frac[part][isbarrel][i]*data_med_sum_err/fakeYield)**2+
                        # (self.fake_frac[part][isbarrel][i]*data_med_sum*fakeYield_err/fakeYield**2)**2+
                        (self.fake_frac_err[part][isbarrel][i]*data_med_sum/fakeYield)**2)
                    tw = data_med_sum / mc_med_sum * trueWeight
                    tw_err = np.sqrt(
                        # (trueWeight*data_med_sum_err/mc_med_sum)**2+
                        # (trueWeight*data_med_sum*mc_med_sum_err/mc_med_sum**2)**2+
                        (self.fake_frac_err[part][isbarrel][i]*data_med_sum/mc_med_sum)**2)
                    print(fw, "+-", fw_err)
                    print(tw, "+-", tw_err)
                    self.fakeWeight[part][isbarrel]['w'].append(fw)
                    self.trueWeight[part][isbarrel]['w'].append(tw)
                    self.fakeWeight[part][isbarrel]['w_err'].append(fw_err)
                    self.trueWeight[part][isbarrel]['w_err'].append(tw_err)

        # print(self.fakeWeight)
        # print(self.trueWeight)

        print('==================================')
        print('  Finish fake weight calculation  ')
        print('==================================')

    def applyFakeWeight(self, type):
        output_path = self._outputFolder + '/%s/2017.root' % type
        if not os.path.isdir(self._outputFolder): os.makedirs(self._outputFolder)
        if os.path.isfile(output_path): os.remove(output_path)
        print("Target path: ", output_path)

        if type == "data_fake":
            in_data = self.data_fake.copy()
            minus = 1
        if type == "mc_true":
            in_data = self.mc_true.copy()
            minus = -1
        data_o = pd.DataFrame(data={"fake_weight":np.zeros_like(in_data['gamma_pt']),"fake_weight_err":np.zeros_like(in_data['gamma_pt'])})#,"fake_weight_up":np.zeros_like(in_data['gamma_pt']),"fake_weight_down":np.zeros_like(in_data['gamma_pt'])})

        for part in ['ele', 'mu']:
            if part == 'ele':
                part_cut = 11
            else: 
                part_cut = 13
            for isbarrel in ['isBarrel', 'isEndcap']:
                if isbarrel == 'isBarrel':
                    eta_cut = 1
                else: 
                    eta_cut = 0
                for i in range(len(self.pt_bin[part][isbarrel])-1):
                    data_o['fake_weight'] += (self.eta_sel(in_data, eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.fakeWeight[part][isbarrel]['w'][i] * minus
                    data_o['fake_weight_err'] += (self.eta_sel(in_data, eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.fakeWeight[part][isbarrel]['w_err'][i] * minus
                    # data_o['fake_weight_up'] += ((in_data.loc[:,'photon_is_barrel']==eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.weight[part][isbarrel]['w_up'][i] * minus
                    # data_o['fake_weight_down'] += ((in_data.loc[:,'photon_is_barrel']==eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.weight[part][isbarrel]['w_down'][i] * minus
        # print(data_o)
        data_o['weight'] = in_data['weight'] * data_o['fake_weight']
        data_o['weight_err'] = in_data['weight'] * data_o['fake_weight_err']
        in_data.pop('weight')
        in_data = pd.concat([in_data, data_o], sort=False, axis=1)
        self.storeOutput(in_data, output_path)

        del in_data, data_o

    def applyTrueWeight(self, type):
        output_path = self._outputFolder + '/%s/2017.root' % type
        if not os.path.isdir(self._outputFolder): os.makedirs(self._outputFolder)
        if os.path.isfile(output_path): os.remove(output_path)
        print("Target path: ", output_path)

        if type == "mc_med":
            in_data = self.mc_med.copy()
            minus = 1
        else:
            return 0
        data_o = pd.DataFrame(data={"true_weight":np.zeros_like(in_data['gamma_pt']), "true_weight_err":np.zeros_like(in_data['gamma_pt'])})#,"fake_weight_up":np.zeros_like(in_data['gamma_pt']),"fake_weight_down":np.zeros_like(in_data['gamma_pt'])})

        for part in ['ele', 'mu']:
            if part == 'ele':
                part_cut = 11
            else: 
                part_cut = 13
            for isbarrel in ['isBarrel', 'isEndcap']:
                if isbarrel == 'isBarrel':
                    eta_cut = 1
                else: 
                    eta_cut = 0
                for i in range(len(self.pt_bin[part][isbarrel])-1):
                    data_o['true_weight'] += (self.eta_sel(in_data, eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.trueWeight[part][isbarrel]['w'][i] * minus
                    data_o['true_weight_err'] += (self.eta_sel(in_data, eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.trueWeight[part][isbarrel]['w_err'][i] * minus
                    # data_o['fake_weight_up'] += ((in_data.loc[:,'photon_is_barrel']==eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.weight[part][isbarrel]['w_up'][i] * minus
                    # data_o['fake_weight_down'] += ((in_data.loc[:,'photon_is_barrel']==eta_cut) & (in_data.loc[:,'Z_lead_lepton_id']==part_cut) & (in_data.loc[:,'gamma_pt']>=self.pt_bin[part][isbarrel][i]) & (in_data.loc[:,'gamma_pt']<self.pt_bin[part][isbarrel][i+1])) * self.weight[part][isbarrel]['w_down'][i] * minus
        # print(data_o)
        data_o['weight'] = in_data['weight'] * data_o['true_weight']
        data_o['weight_err'] = in_data['weight'] * data_o['true_weight_err']
        in_data.pop('weight')
        in_data = pd.concat([in_data, data_o], sort=False, axis=1)
        in_data.to_root(output_path, key='inclusive', index=False)
        self.storeOutput(in_data, output_path)

        del in_data, data_o

def main():

    args=getArgs()
    
    fwa = ApplyWeightHandler(args.region)
    fwa.readFiles()
    fwa.calculateWeight()
    fwa.applyFakeWeight("data_fake")
    fwa.applyFakeWeight("mc_true")
    fwa.applyTrueWeight("mc_med")
    print('==================================')
    print('  Finish fake weight application  ')
    print('==================================')

    return 0

if __name__ == '__main__':
    main()
