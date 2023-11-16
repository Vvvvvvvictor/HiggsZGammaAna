import os
from Plot_Helper import MakeStack
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *

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