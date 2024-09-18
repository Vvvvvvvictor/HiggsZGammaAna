import os
import sys
from ROOT import *
import json

class Analyzer_Config:
    def __init__(self, channel, year, region='',treename=''):
        self.channel            = channel
        self.year               = year
        self.version            = 'UL'
        self.region             = region
        self.sample_loc         = 'NONE'
        self.out_dir            = 'NONE'
        self.out_region_name    = 'NONE'
        self.root_output_name   = 'NONE'
        self.plot_output_path   = 'NONE'
        self.BDT_filename       = 'NONE'
        self.treename           = treename
        self.mvaCut             = {}
        self.sig_names          = []
        self.bkg_names          = []
        self.samp_names         = []
        self.sys_names          = []
        self.norm_SFs           = 1.0
        
        self.reweight_variables = ['Z_cos_theta', 'H_relpt', 'gamma_mvaID', 'gamma_relpt', 'gamma_ptRelErr']

        self.Config_Analyzer()

    def Config_Analyzer(self):
        if self.channel == 'inclusive' or self.channel == 'ggH' or self.channel == 'VBFH' or self.channel == 'WH' or self.channel == 'ZH' or self.channel == 'ttH':

            if self.region == "SR":
                self.out_region_name = self.version + '_SR'
            elif self.region == "CR":
                self.out_region_name = self.version + '_CR'
            else:
                self.out_region_name = self.version + '_full'

            if self.year == '2016':
                self.sample_loc       = '/eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples'
                self.out_dir          = 'plots_16UL'
                self.root_output_name = "ALP_plot_data16_{0}_{1}.root".format(self.out_region_name,self.treename)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
                #mvaCut = 0.8675
            elif self.year == '-2016':
                self.sample_loc       = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/16APV'
                self.out_dir          = 'plots_16APVUL'
                self.root_output_name = "ALP_plot_data16APV_{0}_{1}.root".format(self.out_region_name,self.treename)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
                #mvaCut = 0.8365
            elif self.year == '2017':
                self.sample_loc       = '/eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples'
                #self.sample_loc       = '/afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/Plotter/Data/2017'
                self.out_dir          = 'plots_17UL'
                self.root_output_name = "HZGamma_plot_data17_{0}_{1}.root".format(self.out_region_name,self.treename)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2017.pkl"
                #mvaCut = 0.8365
            elif self.year == '2018':
                self.sample_loc       = '/eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples'
                self.out_dir          = 'plots_18UL'
                self.root_output_name = "ALP_plot_data18_{0}_{1}.root".format(self.out_region_name,self.treename)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2018.pkl"
                #mvaCut = 0.9766
            elif 'run2' in self.year:
                self.sample_loc       = '/eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples_v1'
                self.out_dir          = 'plots_run2UL'
                self.root_output_name = "ALP_plot_run2_{0}_{1}.root".format(self.out_region_name,self.treename)
                self.BDT_filename     = "/afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/Plotter/Data/UL/BDT/model_HZGamma_BDT_zero_to_one_jet_fullrange.pkl"
                self.mvaCut           = [0.0, 0.4449999928474426, 0.6000000238418579, 0.7699999809265137, 0.8849999904632568, 1.]
            elif self.year == '2022':
                self.sample_loc       = '/eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples'
                #self.sample_loc       = '/afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/Plotter/Data/2017'
                self.out_dir          = 'plots_17UL'
                self.root_output_name = "HZGamma_plot_data17_{0}_{1}.root".format(self.out_region_name,self.treename)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2017.pkl"
            else:
                print('[[Analysis Configs]] do not included at 2016/2017/2018!')
                exit(0)

            if self.channel == 'inclusive':
                self.sig_names  = ['sig']
            else:
                self.sig_names  = [self.channel]

            self.bkg_names  = ['ZG2JToG2L2J', 'WZ',  'TTGJets', 'DYJetsToLL', 'WW', 'TGJets', 'ZGToLLG', 'TT', 'EWKZ2J', 'ttWJets', 'ttZJets', 'WGToLNuG', 'ZZ']
            self.samp_names = self.bkg_names + self.sig_names + ['data']

            self.plot_output_path = "{0}/plot_{1}_{2}".format(self.out_dir, self.out_region_name, self.treename)
            
            self.sys_names  = ['CMS_eff_g_up','CMS_eff_g_dn','CMS_pileup_up','CMS_pileup_dn','CMS_eff_lep_up','CMS_eff_lep_dn']

        else:
            print("channel is invalid: channel = %s" %self.channel)
            sys.exit()

    def Print_Config(self):
        print('Running analysis in channel: %s' %self.channel)
        print('getting ntuples from: %s' %self.sample_loc)
        print('using signals: ')
        print(self.sig_names)
        print('using backgrounds:')
        print(self.bkg_names)

    def Create_output_fiels(self):
        # Check the output files
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        if not os.path.exists(self.plot_output_path):
            os.makedirs(self.plot_output_path)

        out_file = TFile( self.out_dir + '/' + self.root_output_name , "RECREATE")

        return out_file

    def load_reweighting(self, json_name):
        with open(json_name,'r', encoding='UTF-8') as f:
            load_dict = json.load(f)

        return load_dict