#!/bin/bash

cd /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/Combine_results_tf

eval `scramv1 runtime -sh`

text2workspace.py Datacard.txt -o Datacard_mu_inclusive.root -m 125 higgsMassRange=122,128 