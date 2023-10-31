#!/bin/bash

year=( 16 16APV 17 18 )
nYear=${#year[@]}
massList=( 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 )
#massList=( 1 )
massList_sig=( 10 10 15 15 15 15 20 20 20 20 25 25 25 25 30 30 )
#massList=( 4 )
nMass=${#massList[@]}

cd ./Signal
mkdir ./ALP_SigModel_param_UL/fit_results_runII

for ((iBin=0; iBin<$nMass; iBin++))
    do

    if [ $# -ne 0 ];then
        if [[ ${massList[ $iBin ]} -ne $1 ]];then
          continue
        fi
        if [[ ${years[ $jBin ]} != $2 ]];then
          continue
        fi
    fi

    mkdir ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}
    #rm -rf ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/*
    
    #cp ./ALP_SigModel_param_UL/fit_results_run2_mu/M${massList_sig[$iBin]}/16/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_mu.root
    #cp ./ALP_SigModel_param_UL/fit_results_run2_mu/M${massList_sig[$iBin]}/16APV/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_mu.root
    #cp ./ALP_SigModel_param_UL/fit_results_run2_mu/M${massList_sig[$iBin]}/17/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_mu.root
    #cp ./ALP_SigModel_param_UL/fit_results_run2_mu/M${massList_sig[$iBin]}/18/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_mu.root
    
    
    #cp ./ALP_SigModel_param_UL/fit_results_run2_ele/M${massList_sig[$iBin]}/16/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_ele.root
    #cp ./ALP_SigModel_param_UL/fit_results_run2_ele/M${massList_sig[$iBin]}/16APV/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_ele.root
    #cp ./ALP_SigModel_param_UL/fit_results_run2_ele/M${massList_sig[$iBin]}/17/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_ele.root
    #cp ./ALP_SigModel_param_UL/fit_results_run2_ele/M${massList_sig[$iBin]}/18/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_ele.root
    
    #cp ../Background/ALP_BkgModel_param_UL/fit_results_run2/${massList[$iBin]}/CMS-HGG_mva_13TeV_multipdf.root ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}
    
    
    ###python makeParametricModelDatacard_ALP_UL_v2.py -i ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_mu.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_mu.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_mu.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_mu.root -s ./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16APV.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_17.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_18.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16APV.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_17.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_18.dat -o ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt --NormSys ./ALP_sys_UL/NormSys/NormalizationSys_16_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_16APV_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_17_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_18_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_16_mu.txt,./ALP_sys_UL/NormSys/NormalizationSys_16APV_mu.txt,./ALP_sys_UL/NormSys/NormalizationSys_17_mu.txt,./ALP_sys_UL/NormSys/NormalizationSys_18_mu.txt --BDTSys ./ALP_sys_UL/BDTSys/BDTSys_16_ele.txt,./ALP_sys_UL/BDTSys/BDTSys_16APV_ele.txt,./ALP_sys_UL/BDTSys/BDTSys_17_ele.txt,./ALP_sys_UL/BDTSys/BDTSys_18_ele.txt,./ALP_sys_UL/BDTSys/BDTSys_16_mu.txt,./ALP_sys_UL/BDTSys/BDTSys_16APV_mu.txt,./ALP_sys_UL/BDTSys/BDTSys_17_mu.txt,./ALP_sys_UL/BDTSys/BDTSys_18_mu.txt -p ggH_16_ele,ggH_16APV_ele,ggH_17_ele,ggH_18_ele,ggH_16_mu,ggH_16APV_mu,ggH_17_mu,ggH_18_mu -c cat0 --isMultiPdf -m 125 --mA ${massList[$iBin]} --interp ./ALP_sys_UL/eff.txt
    #python makeParametricModelDatacard_ALP_UL_v2.py -i ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_ele.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_mu.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_mu.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_mu.root,./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_mu.root -s ./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16APV.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_17.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_18.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_16APV.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_17.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/M${massList_sig[$iBin]}/ScaleSmear_m${massList_sig[$iBin]}_18.dat -o ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt --NormSys ./ALP_sys_UL/NormSys/NormalizationSys_16_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_16APV_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_17_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_18_ele.txt,./ALP_sys_UL/NormSys/NormalizationSys_16_mu.txt,./ALP_sys_UL/NormSys/NormalizationSys_16APV_mu.txt,./ALP_sys_UL/NormSys/NormalizationSys_17_mu.txt,./ALP_sys_UL/NormSys/NormalizationSys_18_mu.txt -p ggH_16_ele,ggH_16APV_ele,ggH_17_ele,ggH_18_ele,ggH_16_mu,ggH_16APV_mu,ggH_17_mu,ggH_18_mu -c cat0 --isMultiPdf -m 125 --mA ${massList[$iBin]} --interp ./ALP_sys_UL/eff.txt

    cd ./ALP_SigModel_param_UL/fit_results_runII/M${massList[$iBin]}
    pwd
    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125 --rAbsAcc 0.00000001 --rMin -200 --rMax 200 --freezeParameters MH
    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125 --freezeParameters MH

    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits -m 125 --freezeParameters MH -n _observed
    
    #text2workspace.py datacard_ALPmass${massList[$iBin]}.txt -m 125 -o datacard_ALPmass${massList[$iBin]}.root

    ### expected Impacts ####
    #rm -rf ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/higgsCombine_paramFit*
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 10 --robustFit 1 --doInitialFit -t -1 --expectSignal 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 10 --robustFit 1 --doFits -t -1 --expectSignal 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Expected.json
    #plotImpacts.py -i M${massList[$iBin]}_Expected.json -o M${massList[$iBin]}_Expected_Impact
    
    ### Observed Impacts ####
    #rm -rf ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/higgsCombine_paramFit*

    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -0.1 --rMax 1 --robustFit 1 --doInitialFit --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH #-v 3
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -0.1 --rMax 1 --robustFit 1 --doFits --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Observed.json
    #plotImpacts.py -i M${massList[$iBin]}_Observed.json  -o M${massList[$iBin]}_Observed_Impact
    
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -0.1 --rMax 1 --robustFit 1 --doInitialFit --cminDefaultMinimizerStrategy 0 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --robustFit 1 --doFits --cminDefaultMinimizerStrategy 0 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Observed.json
    #plotImpacts.py -i M${massList[$iBin]}_Observed.json --blind -o M${massList[$iBin]}_Observed_Impact

    ## start GoF
    #cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/runbatch_1.* .
    #sed -i "s/_MASS/${massList[$iBin]}/g" runbatch_1.sh

    #condor_submit ./runbatch_1.jdl

    #combine -M GoodnessOfFit datacard_ALPmass${massList[$iBin]}.txt --algo=saturated -m 125 --setParameters MH=125 -n _saturated
    #combine -M GoodnessOfFit datacard_ALPmass${massList[$iBin]}.txt --algo=saturated -t 1000 -s 12345 -m 125 --setParameters MH=125 -n _saturated

    #combine -M GoodnessOfFit datacard_ALPmass${massList[$iBin]}.txt --algo=KS -m 125 --setParameters MH=125 -n _KS
    #combine -M GoodnessOfFit datacard_ALPmass${massList[$iBin]}.txt --algo=KS -t 100 -s 12345 -m 125 --setParameters MH=125 -n _KS

    #python /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/runGoF.py > GoF.log
    ## end GoF



    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 10 --robustFit 1 --doInitialFit -t -1 --expectSignal 1 --cminDefaultMinimizerStrategy 0 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 10 --robustFit 1 --doFits -t -1 --expectSignal 1 --cminDefaultMinimizerStrategy 0 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Expected.json
    #plotImpacts.py -i M${massList[$iBin]}_Expected.json -o M${massList[$iBin]}_Expected_Impact
    cd ../../../
    
    
done

cd ./ALP_SigModel_param_UL/fit_results_runII
#mkdir plot
cd plot
cp ../../fit_results_runII_expect/plot/*.py ./
sed -i "33s/fit_results_run2_/fit_results_/" com_plot.py
python com_plot.py -fb -y run2 -c runII --interp --observe
gs -r600 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=Figure_006.png -dBATCH -dNOPAUSE ALP_xs_UpperLimit.eps