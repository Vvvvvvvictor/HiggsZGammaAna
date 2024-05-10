import os
from Plot_Helper import MakeStack
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *
import numpy as np

def getMassSigma(ana_cfg):
    sigma_low = {}
    sigma_hig = {}
    for sample in ana_cfg.sig_names:
        files = TFile(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
        filesTree = files.Get("passedEvents")
        filesTree.Draw("{0}>>tree{0}".format("H_m"),"factor*pho1SFs*pho1SFs*({0})".format("H_m>-90&&passChaHadIso&&passNeuHadIso&&passdR_gl&&passHOverE&&H_m>115&&H_m<135"))
        Hist = gDirectory.Get("tree{0}".format("H_m"))
        sigma_bin = 0
        for i in range(Hist.GetNbinsX()/2):
            if Hist.Integral(50-i,51+i)/Hist.Integral()>0.683: 
                sigma_bin = i
                break
        sigma_low[sample] = Hist.GetBinCenter(50-sigma_bin)-125.0-Hist.GetBinWidth(50-sigma_bin)/2.0
        sigma_hig[sample] = Hist.GetBinCenter(51+sigma_bin)-125.0+Hist.GetBinWidth(51+sigma_bin)/2.0
    
    return sigma_low, sigma_hig

def reweighting(ana_cfg, json, values):
    
    x_bins = [-1., -0.4, 0.4, 1]
    y_bins = [0., 0.05, 0.1, 1]
    z_bins = [0.4, 0.8, 0.9, 1.]
    w_bins = [0., 0.15, 0.2, 0.4]
    v_bins = [0.01, 0.02, 0.03, 0.08]

    
    
    bins_list = [x_bins, y_bins, z_bins, w_bins, v_bins]

    inbins = 1

    for i in range(len(values)):
        if (values[i] >= max(bins_list[i])) or (values[i] <= min(bins_list[i])):
            inbins = 0

    if inbins:
        bin_indices = np.array([np.digitize(values[i], bins) for i, bins in enumerate(bins_list)])
        bin_number = np.ravel_multi_index(bin_indices - 1, [len(bins) - 1 for bins in bins_list])

        weight = json['bin_{}'.format(bin_number)][1] * ana_cfg.norm_SFs
        #print('nbin:', bin_number, json['bin_{}'.format(bin_number)][0], weight)
    else:
        weight = ana_cfg.norm_SFs

    #print(weight)

    return weight