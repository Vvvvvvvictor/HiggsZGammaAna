####################################################
####################################################

import os
import sys

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend
from Analyzer_Helper import MuPair, DimuCandidates, IsBtagLoose, IsBtagMed, IsBtagTight, CountYield
import Analyzer_Configs as AC
import Plot_Configs     as PC

gROOT.SetBatch(True)


def main():
    out_name = "ttH_exercise_out.root"    
    if not os.path.exists('plots'):
	os.makedirs('plots')
    out_file = TFile( "plots/" + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive')
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    var_names = ['mu1_lepMVA', 'mu2_lepMVA', 'dimu_mass', 'dimu_pt', 'mu3_lepMVA', 'ele_lepMVA']
    plot_cfg = PC.Plot_Config(analyzer_cfg)

    histos = {}
    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
        histos['mu1_lepMVA'][sample] = TH1F('mu1_lepMVA' + '_' + sample, 'mu1_lepMVA' + '_' + sample, 20,  -1,   1)
        histos['mu2_lepMVA'][sample] = TH1F('mu2_lepMVA' + '_' + sample, 'mu2_lepMVA' + '_' + sample, 20,  -1,   1)
        histos['dimu_mass'][sample]  = TH1F('dimu_mass'  + '_' + sample, 'dimu_mass'  + '_' + sample, 40, 110, 150)
        histos['dimu_pt'][sample]    = TH1F('dimu_pt'    + '_' + sample, 'dimu_pt'    + '_' + sample, 20,  20, 420)
	histos['mu3_lepMVA'][sample] = TH1F('mu3_lepMVA' + '_' + sample, 'mu3_lepMVA' + '_' + sample, 20,  -1,   1)
	histos['ele_lepMVA'][sample] = TH1F('ele_lepMVA' + '_' + sample, 'ele_lepMVA' + '_' + sample, 20,  -1,   1)

    for sample in analyzer_cfg.samp_names:
	ntup = ntuples[sample] # just a short name
	print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()
        for iEvt in range( ntup.GetEntries() ):
	    ntup.GetEvent(iEvt)
	    #if (iEvt % 50 != 1): continue
	    if (iEvt % 100000 == 1):
	        print "looking at event %d" %iEvt


	    ##################   baseline ttH selections   ##############
	    InttHChannel = True
	    if ntup.nMuons == 2 and ntup.nElectrons == 0:  IntHChannel = False
	    nLooseB = 0
	    nMedB   = 0
	    for jet_idx in range(ntup.nJets):
		if IsBtagLoose('2018', ntup.jet_DeepCSV[jet_idx]): nLooseB += 1
		if IsBtagMed  ('2018', ntup.jet_DeepCSV[jet_idx]): nMedB   += 1
	    #if nLooseB > 1 or nMedB > 0:  InttHChannel = False
	    if nLooseB < 2 and nMedB < 1:  InttHChannel = False
	    if not InttHChannel: continue
	    ##################### baseline ttH selections ###############


	    #####################  select muon pair  ###################
	    muPairs = DimuCandidates(ntup)
	    ## Veto Z-mass pair 
	    for muPair in muPairs:
	        if abs(muPair.mass - 91) < 5: InttHChannel = False
	    if not InttHChannel: continue	    

	    ## Select highest-pt muPair in [110, 150] window
	    dimu = None	
	    for muPair in muPairs:
		if ntup.muon_lepMVA[muPair.mu1_idx] < 0.4 or ntup.muon_lepMVA[muPair.mu2_idx] < 0.4: continue
		if muPair.mass > 110 and muPair.mass < 150:
		    if dimu == None:            dimu = muPair
		    elif muPair.pt > dimu.pt:   dimu = muPair
	    if dimu == None: continue
	    ####################  found muon pair  #######################

	    #################  get one extra lepton  #####################
	    InttHChannel = False
	    if ntup.nMuons > 2:
	      for mu_idx in range(ntup.nMuons):
		if mu_idx != dimu.mu1_idx and mu_idx != dimu.mu2_idx and ntup.muon_lepMVA[mu_idx] >= 0.4:
		  histos['mu3_lepMVA'][sample] .Fill( ntup.muon_lepMVA[mu_idx], ntup.xsec_norm * ntup.event_wgt)
		  InttHChannel = True
		  break ## for now, just look at one extra muon per events
	    if ntup.nElectrons > 0:
	      for ele_idx in range(ntup.nElectrons):
		if ntup.electron_lepMVA[ele_idx] >= 0.4:
		  histos['ele_lepMVA'][sample] .Fill( ntup.electron_lepMVA[0], ntup.xsec_norm * ntup.event_wgt)	
		  InttHChannel = True
		  break ## similarly, only look at one electron
	    if not InttHChannel: continue
	    #################  get one extra lepton  #####################

	    histos['mu1_lepMVA'][sample] .Fill( ntup.muon_lepMVA[dimu.mu1_idx], ntup.xsec_norm * ntup.event_wgt)
	    histos['mu2_lepMVA'][sample] .Fill( ntup.muon_lepMVA[dimu.mu2_idx], ntup.xsec_norm * ntup.event_wgt)
	    if sample != 'data' or abs(dimu.mass - 125) > 5:  # blind data in signal region
	        histos['dimu_mass'][sample].Fill( dimu.mass, ntup.xsec_norm * ntup.event_wgt)
	    histos['dimu_pt'][sample].Fill( dimu.pt        , ntup.xsec_norm * ntup.event_wgt)

    out_file.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
	    plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    lumi_label = MakeLumiLabel(plot_cfg.lumi)
    cms_label  = MakeCMSDASLabel()

    for var_name in var_names:
	stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
	scaled_sig = ScaleSignal(plot_cfg, stacks['sig'], var_name)
	ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
	legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

 	canv = CreateCanvas(var_name)
	DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label)
	canv.Write()

    print '\n\n'
    CountYield(analyzer_cfg, histos['dimu_mass'])
    out_file.Close()

main()











	
