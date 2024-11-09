#!/bin/bash

cd /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine

eval `scramv1 runtime -sh`

text2workspace.py output_Datacard/pruned_datacard_2017.txt -o ./output_datacard_rootfile/Datacard_2017_mu_inclusive.root -m 125 higgsMassRange=122,128 