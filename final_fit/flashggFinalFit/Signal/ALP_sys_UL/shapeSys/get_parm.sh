#!/bin/bash
massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

years=( 16 16APV 17 18 )

nMass=${#massList[@]}
nyear=${#years[@]}

for ((iBin=0; iBin<$nMass; iBin++))
    do
    for ((jBin=0; jBin<$nyear; jBin++))
        do
        
        if [[ $1 -eq 0 ]]; then
            cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_run2_ele/M${massList[$iBin]}/${years[$jBin]}/CMS-HGG_sigfit_data_ggh_cat0.root ./fit_parm_ele/paramDump_ggh_cat0_M${massList[$iBin]}_${years[$jBin]}.root
            cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_run2_ele/M${massList[$iBin]}/${years[$jBin]}/plots_ALP/paramDump_ggh_cat0txt ./fit_parm_ele/paramDump_ggh_cat0_M${massList[$iBin]}_${years[$jBin]}.txt
        elif [[ $1 -eq 1 ]]; then
            cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_run2_mu/M${massList[$iBin]}/${years[$jBin]}/CMS-HGG_sigfit_data_ggh_cat0.root ./fit_parm_mu/paramDump_ggh_cat0_M${massList[$iBin]}_${years[$jBin]}.root
            cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_run2_mu/M${massList[$iBin]}/${years[$jBin]}/plots_ALP/paramDump_ggh_cat0txt ./fit_parm_mu/paramDump_ggh_cat0_M${massList[$iBin]}_${years[$jBin]}.txt
        fi
    done
done