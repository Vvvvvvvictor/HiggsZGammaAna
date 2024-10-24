#!/bin/bash
ulimit -s unlimited
set -e
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src
export SCRAM_ARCH=el9_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard
export PYTHONPATH=$PYTHONPATH:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/tools

python3 /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/makeYields.py --cat cat0 --procs ggH_16_ele,ggH_16APV_ele,ggH_17_ele,ggH_18_ele,ggH_16_mu,ggH_16APV_mu,ggH_17_mu,ggH_18_mu --ext HToZa --mass 125 --inputWSDirMap run2=/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/data --sigModelWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_ele/signalFit/output,/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_mu/signalFit/output --sigModelExt packaged --bkgModelWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_param_UL/fit_results_run2 --bkgModelExt multipdf  --doSystematics
