
source /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/tools


python3 ./RunYields.py --inputWSDirMap run2=/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/data --procs ggH_16_ele,ggH_16APV_ele,ggH_17_ele,ggH_18_ele,ggH_16_mu,ggH_16APV_mu,ggH_17_mu,ggH_18_mu --sigModelWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_ele/signalFit/output,/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_mu/signalFit/output --bkgModelWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_param_UL/fit_results_run2