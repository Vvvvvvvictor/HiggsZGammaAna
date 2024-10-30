#!/bin/bash

cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine

eval `scramv1 runtime -sh`

text2workspace.py /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_mu/15_pruned_datacard_16APV_mu.txt -o ./output_datacard_rootfile_mu/15_Datacard_16APV_mu_mu_inclusive.root -m 125 higgsMassRange=122,128 