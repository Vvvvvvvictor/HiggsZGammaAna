####################################################
####################################################

import os
import sys
import numpy as np

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc
from Analyzer_Helper import getMassSigma
import Analyzer_Configs as AC
import Plot_Configs     as PC

import tdrstyle

import argparse
parser = argparse.ArgumentParser(description="Plotter for HZGamma analysis")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument("--region", dest="region", default="full", help="CR for sideband region, SR for signal region, full for full region")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
args = parser.parse_args()



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def main():

    # Setup the plot configuration
    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year, args.region, args.mva)
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
            #if (iEvt == 100): break

            
            if args.ele:
                if abs(ntup.Z_lead_lepton_id) == 13: 
                    continue
            if args.mu:
                if abs(ntup.Z_lead_lepton_id) == 11: 
                    continue
            
            weight = ntup.weight_central
            #weight = ntup.weight

            if  args.region == "SR" and (ntup.H_m>135. or ntup.H_m<115.): continue
            if  args.region == "CR" and (ntup.H_m<135. and ntup.H_m>115.): continue

            #print(getattr(ntup,'Z_pt'))

            if not (ntup.n_jets >=2 and ntup.n_leptons < 3):
                continue

            for var_name in var_names:
                histos[var_name][sample].Fill( getattr(ntup,var_name), weight )

        # End of for iEvt in range( ntup.GetEntries() )
    # End of for sample in analyzer_cfg.samp_names
    
    # save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
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
