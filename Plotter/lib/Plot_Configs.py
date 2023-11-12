import os
import sys
from ROOT import *

class Plot_Config:

    def __init__(self, ana_cfg, year):
        self.ana_cfg = ana_cfg
        self.colors  = {}
        self.variables_map  = {}
        self.logY    = False
        self.year    = year
        if year == '2016':
            self.lumi    = '16.81'
        elif year == '-2016':
            self.lumi    = '19.52'
        elif year == '2017':
            self.lumi    = '41.48'
        elif year == '2018':
            self.lumi    = '59.83'
        elif year == 'run2Rereco':
            self.lumi    = '137.24'
        elif year == 'run2':
            self.lumi    = '138'
        else:
            print('do not in 2016/2017/2018!')
            exit(0)
        self.sig_scale  = 5
        self.sig_weight = 0.4133

        self.LoadColors()
        self.LoadVariables()
        self.InitializeHistos()



    def LoadColors(self):
        self.colors["Data"] = kBlack
        self.colors["sig"] = kRed
        self.colors["ggH"] = kRed
        self.colors["VBFH"] = kRed
        self.colors["WH"] = kRed
        self.colors["ZH"] = kRed
        self.colors["ttH"] = kRed
        

        self.colors["DYJetsToLL"]  =  kAzure + 7
        self.colors["ZGToLLG"]  =  kOrange
        self.colors["TT"]  =  kGreen
        self.colors["ZG2JToG2L2J"]  =  kViolet
        self.colors["TGJets"]  =  kSpring
        self.colors["TTGJets"]  =  kYellow
        self.colors["WW"]  =  kPink
        self.colors["WZ"]  =  kRed -7
        self.colors["ZZ"]  =  kYellow -7
        self.colors["LLAJJ"]  =  kGreen -7

    def LoadVariables(self):
        self.variables_map = {
            'Z_pt':                 [r"p_{T, Z}",               50, 0., 50.],
            'Z_eta':                [r"\eta_{Z}",               50, -2.5, 2.5],
            'Z_phi':                [r"\phi_{Z}",               50, -4., 4.],
            'Z_mass':               [r"m_{Z}",                  50, 50., 120.],
            'H_pt':                 [r"p_{T, \ell\ell\gamma}",  50, 0., 80.],
            'H_eta':                [r"\eta_{\ell\ell\gamma}",  50, -2.5, 2.5],
            'H_phi':                [r"\phi_{\ell\ell\gamma}",  50, -4., 4.],
            'H_mass':               [r"m_{\ell\ell\gamma}",     65, 105., 170.],
            'gamma_pt':             [r"p_{T,\gamma}",           50, 0., 50.],
            'gamma_eta':            [r"\eta_{\gamma}",          50, -2.5, 2.5],
            'gamma_phi':            [r"\phi_{\gamma}",          50, -4., 4.],
            'Z_lead_lepton_pt':     [r"p_{T,\ell 1}",           50, 0., 50.],
            'Z_lead_lepton_eta':    [r"\eta_{T,\ell 1}",        50, -2.5, 2.5],
            'Z_lead_lepton_phi':    [r"\phi_{T,\ell 1}",        50, -4., 4.],
            'Z_lead_lepton_mass':   [r"m_{\ell 1}",             50, 0., 1.],
            'Z_sublead_lepton_pt':  [r"p_{T,\ell 2}",           50, 0., 50.],
            'Z_sublead_lepton_eta': [r"\eta_{T,\ell 2}",        50, -2.5, 2.5],
            'Z_sublead_lepton_phi': [r"\phi_{T,\ell 2}",        50, -4., 4.],
            'Z_sublead_lepton_mass':[r"m_{\ell 2}",             50, 0., 1.]
        }

    def InitializeHistos(self):
        histos = {}

        for var_name in self.variables_map.keys():
            histos[var_name] = {}

            nbins = self.variables_map[var_name][1]
            x_min = self.variables_map[var_name][2]
            x_max = self.variables_map[var_name][3]
            for sample in self.ana_cfg.samp_names:
                histos[var_name][sample]    = TH1F( "{}_{}".format(var_name, sample), "{}_{}".format(var_name, sample), nbins, x_min, x_max)

        return histos


    def SetHistStyles(self, hist, sample):
        if sample == 'data':
            hist.SetMarkerStyle(20)
            hist.GetXaxis().SetTickLength(0.04)
            hist.GetXaxis().SetLabelSize(0)
            hist.GetXaxis().SetTitleOffset(0.95)
            hist.GetYaxis().SetLabelSize(0.04)
            hist.GetYaxis().SetTitleSize(0.05)
            hist.GetYaxis().SetTitleFont(42)
            hist.GetYaxis().SetTitleOffset(1.15)
        elif sample in self.ana_cfg.sig_names:
            hist.SetLineColor(self.colors[sample])
            hist.SetLineWidth(2)
            hist.SetFillStyle(0)
            #hist.SetFillColor(kGray)
        elif sample in self.ana_cfg.bkg_names:
            hist.SetFillColor(self.colors[sample])

    def SetHzaStyle():
        hzaStyle = TStyle("hzaPaperStyle","Hza Paper Style")

        hzaStyle.SetFrameFillColor(0)
        hzaStyle.SetStatColor(0)
        hzaStyle.SetOptStat(0)
        hzaStyle.SetTitleFillColor(0)
        hzaStyle.SetCanvasBorderMode(0)
        hzaStyle.SetPadBorderMode(0)
        hzaStyle.SetFrameBorderMode(0)
        hzaStyle.SetPadColor(kWhite)
        hzaStyle.SetCanvasColor(kWhite)

        hzaStyle.SetCanvasDefH(600) #Height of canvas
        hzaStyle.SetCanvasDefW(600) #Width of canvas
        hzaStyle.SetCanvasDefX(0)   #POsition on screen
        hzaStyle.SetCanvasDefY(0)

        hzaStyle.SetPadLeftMargin(0.13)
        hzaStyle.SetPadRightMargin(0.1)
        hzaStyle.SetPadTopMargin(0.085)
        hzaStyle.SetPadBottomMargin(0.12)

        # For hgg axis titles:
        hzaStyle.SetTitleColor(1, "XYZ")
        hzaStyle.SetTitleFont(42, "XYZ")
        hzaStyle.SetTitleSize(0.05, "XYZ")
        hzaStyle.SetTitleXOffset(0.95)#//0.9)
        hzaStyle.SetTitleYOffset(1.15)# // => 1.15 if exponents

        # For hgg axis labels:
        hzaStyle.SetLabelColor(1, "XYZ")
        hzaStyle.SetLabelFont(42, "XYZ")
        hzaStyle.SetLabelOffset(0.007, "XYZ")
        hzaStyle.SetLabelSize(0.04, "XYZ")

        # Legends
        hzaStyle.SetLegendBorderSize(0)
        hzaStyle.SetLegendFillColor(kWhite)
        hzaStyle.SetLegendFont(42)

        hzaStyle.SetFillColor(10)
        # Nothing for now
        hzaStyle.SetTextFont(42)
        hzaStyle.SetTextSize(0.03)
        hzaStyle.cd()
                