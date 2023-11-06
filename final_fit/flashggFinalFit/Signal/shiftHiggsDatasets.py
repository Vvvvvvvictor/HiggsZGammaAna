import ROOT
from ROOT import *

import sys
import os

inDir = sys.argv[1]
outDir = inDir 

mass = [250] # radion mass 
for m in mass:
    print "Looking at Radion mass = ", m
    values = [-5,0,5]
    higgs_mass = 125

    ws_name = 'tagsDumper/cms_hgg_13TeV'
    dataset_name = 'ggF_X250_WWgg_qqlnugg_13TeV_HHWWggTag_0'
    #temp_ws = TFile(inDir+'testWS.root').Get(ws_name)
    temp_ws = TFile(inDir+'/Signal_X'+str(m)+'.root').Get(ws_name)
    # temp_ws.Print()
    for value in values:
        shift = value + higgs_mass
        # dataset = (temp_ws.data(dataset_name)).Clone('h4g_'+str(shift)+'_13TeV_4photons')
        # dataset = (temp_ws.data(dataset_name)).Clone('HHWWgg_'+str(shift)+'_13TeV_')
        dataset = (temp_ws.data(dataset_name)).Clone('ggF_'+str(m)+'_13TeV_HHWWggTag_0_'+str(shift)) # includes process and category 
        dataset.Print()
        # dataset.changeObservableName('dZ_bdtVtx','dZ')
        # dataset.changeObservableName('dZ_zeroVtx','dZ')
        dataset.changeObservableName('CMS_hgg_mass','CMS_hgg_mass_old')
        higgs_old = dataset.get()['CMS_hgg_mass_old']
        higgs_new = RooFormulaVar( 'CMS_hgg_mass', 'CMS_hgg_mass', "(@0+%.1f)"%value,RooArgList(higgs_old) );
        dataset.addColumn(higgs_new).setRange(105,145)
        dataset.Print()

        #output = TFile('/afs/cern.ch/work/t/twamorka/ThesisAnalysis/CMSSW_10_5_0/src/flashgg/testshift.root','RECREATE')
        # output = TFile(inDir+'VtxZero/w_signal_'+str(m)+'_'+str(shift)+'.root','RECREATE')
        # output = TFile('/eos/user/t/twamorka/newCatalog_fixVtx_3Oct2019/hadd_WS/reducedWS/w_signal_'+str(m)+'_'+str(shift)+'.root','RECREATE')
        # output = TFile('/afs/cern.ch/work/a/atishelm/21JuneFlashgg/CMSSW_10_5_0/src/flashgg/fggfinalfit_files/X_signal_'+str(m)+'_'+str(shift)+'_HHWWgg_qqlnu.root','RECREATE')
        output = TFile(outDir + 'X_signal_'+str(m)+'_'+str(shift)+'_HHWWgg_qqlnu.root','RECREATE')
        output.mkdir("tagsDumper")
        output.cd("tagsDumper")
        # ws_new = ROOT.RooWorkspace("cms_h4g_13TeV_4photons")
        ws_new = ROOT.RooWorkspace("cms_hgg_13TeV")
        getattr(ws_new,'import')(dataset,RooCmdArg())
        ws_new.Write()
        output.Close()
