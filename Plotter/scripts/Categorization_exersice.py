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

MU_MASS = 0.105658367

def main():
    out_name = "Cat_exersice_out.root"    
    if not os.path.exists('plots'):
        os.makedirs('plots')
    out_file = TFile( "plots/" + out_name , "RECREATE")

    analyzer_cfg = AC.Analyzer_Config('inclusive')
    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)


    var_names = ['mu1_pt', 'mu2_pt', 'dimu_mass', 'dimu_pt', 'dijet_mass', 'dimu_mass_2l', 'dimu_mass_nl']
    plot_cfg = PC.Plot_Config(analyzer_cfg)

    ### declare histograms
    histos = {}
    for var_name in var_names:
        histos[var_name] = {}
    for sample in analyzer_cfg.samp_names:
	histos['mu1_pt'][sample]     = TH1F('mu1_pt'     + '_' + sample, 'mu1_pt'     + '_' + sample, 20,  20, 420)
	histos['mu2_pt'][sample]     = TH1F('mu2_pt'     + '_' + sample, 'mu2_pt'     + '_' + sample, 20,  20, 320)
	histos['dimu_mass'][sample]  = TH1F('dimu_mass'  + '_' + sample, 'dimu_mass'  + '_' + sample, 40, 110, 150)
	histos['dimu_pt'][sample]    = TH1F('dimu_pt'    + '_' + sample, 'dimu_pt'    + '_' + sample, 20,  20, 420)
	histos['dijet_mass'][sample] = TH1F('dijet_mass' + '_' + sample, 'dijet_mass' + '_' + sample, 20, 100, 900)

	histos['dimu_mass_2l'][sample] = TH1F('dimu_mass_2l' + '_' + sample, 'dimu_mass_2l' + '_' + sample, 40, 110, 150)
	histos['dimu_mass_nl'][sample] = TH1F('dimu_mass_nl' + '_' + sample, 'dimu_mass_nl' + '_' + sample, 40, 110, 150)

    ### loop over samples and events
    for sample in analyzer_cfg.samp_names:
	ntup = ntuples[sample] # just a short name
	print '\n\nOn sample: %s' %sample
        print 'total events: %d' %ntup.GetEntries()
        for iEvt in range( ntup.GetEntries() ):
	    ntup.GetEvent(iEvt)
	    #if (iEvt % 50 != 1): continue
	    if (iEvt % 10000 == 1):
	        print "looking at event %d" %iEvt

	    ### find a oppositely charged muon pair
	    ### select highest-pt muon in each charge
	    muP_vec = TLorentzVector(0,0,0,0)
	    muN_vec = TLorentzVector(0,0,0,0)
	    for mu_idx in range(ntup.nMuons):
		if ntup.muon_charge[mu_idx] == 1  and ntup.muon_pt[mu_idx] > muP_vec.Pt(): 
		    muP_vec.SetPtEtaPhiM(ntup.muon_pt[mu_idx], ntup.muon_eta[mu_idx], ntup.muon_phi[mu_idx], MU_MASS)
		if ntup.muon_charge[mu_idx] == -1 and ntup.muon_pt[mu_idx] > muN_vec.Pt(): 
                    muN_vec.SetPtEtaPhiM(ntup.muon_pt[mu_idx], ntup.muon_eta[mu_idx], ntup.muon_phi[mu_idx], MU_MASS)

	    dimu_vec = muP_vec + muN_vec
	    ### found a oppositely charged muon pair


	    ### if at least 2 jets in the event, calculate the highest dijet mass
	    dijet_mass = 0
	    if ntup.nJets >= 2:
	      j1_vec = TLorentzVector(0,0,0,0)
	      j2_vec = TLorentzVector(0,0,0,0)
	      for j1_idx in range(ntup.nJets):
		for j2_idx in range(ntup.nJets):
		  if j1_idx < j2_idx:
		    j1_vec.SetPtEtaPhiM(ntup.jet_pt[0], ntup.jet_eta[0], ntup.jet_phi[0], 0) ## NOTE: This is not actually correct. Jets do have masses.
		    j2_vec.SetPtEtaPhiM(ntup.jet_pt[1], ntup.jet_eta[1], ntup.jet_phi[1], 0) ##       But I did not save them in the ntuples......  -- XWZ
		    dijet_vec = j1_vec + j2_vec
		    if dijet_mass < dijet_vec.M(): dijet_mass = dijet_vec.M()
	    ### got dijet mass

	    histos['mu1_pt'][sample] .Fill( ntup.muon_pt[0], ntup.xsec_norm * ntup.event_wgt )
	    histos['mu2_pt'][sample] .Fill( ntup.muon_pt[1], ntup.xsec_norm * ntup.event_wgt )
	    if sample != 'data' or abs(dimu_vec.M() - 125) > 5:  # blind data in signal region
	        histos['dimu_mass'][sample].Fill( dimu_vec.M()   , ntup.xsec_norm * ntup.event_wgt )
		if ntup.nMuons == 2 and ntup.nElectrons == 0:
		    histos['dimu_mass_2l'][sample].Fill( dimu_vec.M()   , ntup.xsec_norm * ntup.event_wgt )
		else:
		    histos['dimu_mass_nl'][sample].Fill( dimu_vec.M()   , ntup.xsec_norm * ntup.event_wgt )
	    histos['dimu_pt'][sample].Fill( dimu_vec.Pt()  , ntup.xsec_norm * ntup.event_wgt )
	    histos['dijet_mass'][sample].Fill( dijet_mass  , ntup.xsec_norm * ntup.event_wgt )
	## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    ### save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
	    plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    ### save stack plots and make ratio plots
    out_file.cd()
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

    print '\n\nFor inclusive selection'
    CountYield(analyzer_cfg, histos['dimu_mass'])
    print '\n\nFor 2 lepton selection'
    CountYield(analyzer_cfg, histos['dimu_mass_2l'])
    print '\n\nFor multi-lepton selection'
    CountYield(analyzer_cfg, histos['dimu_mass_nl'])
    out_file.Close()


main()











	
