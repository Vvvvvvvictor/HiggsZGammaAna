import os
import sys
from ROOT import *

class Analyzer_Config:
    def __init__(self, channel, year, region, mva):
        self.channel            = channel
        self.year               = year
        self.version            = 'UL'
        self.region             = region
        self.mva                = mva
        self.sample_loc         = 'NONE'
        self.out_dir            = 'NONE'
        self.out_region_name    = 'NONE'
        self.root_output_name   = 'NONE'
        self.plot_output_path   = 'NONE'
        self.BDT_filename       = 'NONE'
        self.mvaCut             = {}
        self.sig_names          = []
        self.bkg_names          = []
        self.samp_names         = []
        self.sys_names          = []

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
                self.sample_loc       = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/16'
                self.out_dir          = 'plots_16UL'
                self.root_output_name = "ALP_plot_data16_{0}.root".format(self.out_region_name)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
                #mvaCut = 0.8675
            elif self.year == '-2016':
                self.sample_loc       = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/16APV'
                self.out_dir          = 'plots_16APVUL'
                self.root_output_name = "ALP_plot_data16APV_{0}}.root".format(self.out_region_name)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
                #mvaCut = 0.8365
            elif self.year == '2017':
                self.sample_loc       = '/afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/Plot/Data/2017'
                self.out_dir          = 'plots_17UL'
                self.root_output_name = "HZGamma_plot_data17_{0}.root".format(self.out_region_name)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2017.pkl"
                #mvaCut = 0.8365
            elif self.year == '2018':
                self.sample_loc       = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/18'
                self.out_dir          = 'plots_18UL'
                self.root_output_name = "ALP_plot_data18_{0}.root".format(self.out_region_name)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2018.pkl"
                #mvaCut = 0.9766
            elif self.year == 'run2':
                self.sample_loc       = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/run2'
                self.out_dir          = 'plots_run2UL'
                self.root_output_name = "ALP_plot_run2_{0}.root".format(self.out_region_name)
                self.BDT_filename     = "/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param.pkl"
                self.mvaCut           = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}
            else:
                print('do not included at 2016/2017/2018!')
                exit(0)

            if self.channel == 'inclusive':
                self.sig_names  = ['sig']
            else:
                self.sig_names  = [self.channel]

            self.bkg_names  = ['DYJetsToLL', 'ZGToLLG', 'TT', 'ZG2JToG2L2J', 'TGJets', 'TTGJets', 'WW', 'WZ', 'LLAJJ']
            self.samp_names = self.bkg_names + self.sig_names + ['data']

            self.plot_output_path = "{0}/plot_{1}".format(self.out_dir, self.out_region_name)
            
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