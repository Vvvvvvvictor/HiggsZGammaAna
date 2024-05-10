####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc
from Analyzer_Helper import reweighting
import Analyzer_Configs as AC
import Plot_Configs     as PC

import tdrstyle

import argparse
parser = argparse.ArgumentParser(description="Plotter for HZGamma analysis")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument("--region", dest="region", default="full", help="CR for sideband region, SR for signal region, full for full region")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
parser.add_argument('--reweighting', dest='reweighting', action='store_true', default=False, help='Apply reweighting?')
args = parser.parse_args()



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def main():

    # Setup the plot configuration
    cat = 'zero_to_one_jet'
    #cat = 'two_jet'
    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year, args.region, cat)
    
    if args.reweighting:
        analyzer_cfg.out_dir = analyzer_cfg.out_dir+"_reweighting"
        analyzer_cfg.plot_output_path = "{0}/plot_{1}_{2}".format(analyzer_cfg.out_dir, analyzer_cfg.out_region_name, analyzer_cfg.treename)
        reweighting_json = analyzer_cfg.load_reweighting('/afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/Plotter/Data/reweight_run2UL_{}.json'.format(cat))

    if args.ele:
        analyzer_cfg.out_dir = analyzer_cfg.out_dir+"_ele"
        analyzer_cfg.plot_output_path = "{0}/plot_{1}_{2}".format(analyzer_cfg.out_dir, analyzer_cfg.out_region_name, analyzer_cfg.treename)
    if args.mu:
        analyzer_cfg.out_dir = analyzer_cfg.out_dir+"_mu"
        analyzer_cfg.plot_output_path = "{0}/plot_{1}_{2}".format(analyzer_cfg.out_dir, analyzer_cfg.out_region_name, analyzer_cfg.treename)

    analyzer_cfg.Print_Config()
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    # creat output files
    out_file = analyzer_cfg.Create_output_fiels()
    # load Ntuples
    ntuples = LoadNtuples(analyzer_cfg)
    # load variables
    var_names = plot_cfg.variables_map.keys()
    # declare histograms
    histos = plot_cfg.InitializeHistos()

    # loop over samples and events
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print('\n\nOn sample: %s' %sample)
        print('total events: %d' %ntup.GetEntries())

        for iEvt in range( ntup.GetEntries() ):
            ntup.GetEvent(iEvt)

            if (iEvt % 100000 == 1):
                print("looking at event %d" %iEvt)
            #if (iEvt == 500): break

            
            if args.ele:
                if abs(ntup.Z_lead_lepton_id) == 13: 
                    continue
            if args.mu:
                if abs(ntup.Z_lead_lepton_id) == 11: 
                    continue
            
            weight = ntup.weight_central

            if  args.region == "SR" and (ntup.H_mass>130. or ntup.H_mass<120.): continue
            if  args.region == "CR" and (ntup.H_mass<130. and ntup.H_mass>120.): continue

            # remove DY duplicate
            if ("DY" in sample) and (ntup.n_iso_photons > 0):
                continue

            if ("ZGToLLG" in sample) and (ntup.n_iso_photons < 1):
                continue

            if args.reweighting:
                reweight = getattr(ntup, 'reweight')
                '''
                if sample in analyzer_cfg.bkg_names+analyzer_cfg.sig_names:
                    value_reweight = []
                    for val in analyzer_cfg.reweight_variables:
                        value_reweight.append(getattr(ntup,val))
                    #print(analyzer_cfg.norm_SFs)
                    reweight = reweighting(analyzer_cfg, reweighting_json, value_reweight)
                '''
            else:
                reweight = 1.0
                

            for var_name in var_names:

                if  args.blind and (sample == 'data') and ('H_mass' in var_name):
                    if getattr(ntup,var_name) > 130 or getattr(ntup,var_name) < 120 :
                        histos[var_name][sample].Fill( getattr(ntup,var_name), weight * reweight )
                else:
                    histos[var_name][sample].Fill( getattr(ntup,var_name), weight * reweight )

        # End of for iEvt in range( ntup.GetEntries() )
    # End of for sample in analyzer_cfg.samp_names
    
    # save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        hist_bkg = stacks['bkg'].GetStack().Last()
        hist_bkg.SetName(var_name+"_bkg")
        hist_bkg.Write()
        for sample in analyzer_cfg.samp_names:
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    # save stack plots and make ratio plots
    out_file.cd()
    scaled_sig = {}
    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        #scaled_sig = 0
        #ratio_plot = 0
        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        if args.ln:
            canv_log = CreateCanvas(var_name+'_log')
            DrawOnCanv(canv_log, var_name, plot_cfg, analyzer_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, logY=True)

            canv_log.Write()
            SaveCanvPic(canv_log, analyzer_cfg.plot_output_path, var_name+'_log')
        else:
            canv = CreateCanvas(var_name)
            DrawOnCanv(canv, var_name, plot_cfg, analyzer_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, logY=False)

            canv.Write()
            SaveCanvPic(canv, analyzer_cfg.plot_output_path, var_name)


    
    print('\n\n')
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()
    
    print('Done')


main()
