#!/bin/bash
version='UL'
year=( 16 16APV 17 18 )
nYear=${#year[@]}
#massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
#massList=( 1 2 3 4 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
#massList=( 3 4 5 6 7 8 )
massList=( 6 )
# 3 4 5 6
#massList=( 3 4 6 7 8 9 10 11 12 13 14 15 )
nMass=${#massList[@]}

cd ./Signal
mkdir ./ALP_SigModel_param_${version}/fit_results_runII

for ((iBin=0; iBin<$nMass; iBin++))
    do
    mkdir ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}

    
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_mu/M${massList[$iBin]}/16/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_mu.root
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_mu/M${massList[$iBin]}/16APV/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_mu.root
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_mu/M${massList[$iBin]}/17/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_mu.root
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_mu/M${massList[$iBin]}/18/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_mu.root
    
    
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_ele/M${massList[$iBin]}/16/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_ele.root
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_ele/M${massList[$iBin]}/16APV/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_ele.root
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_ele/M${massList[$iBin]}/17/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_ele.root
    #cp ./ALP_SigModel_param_${version}/fit_results_run2_ele/M${massList[$iBin]}/18/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_ele.root
    
    #cp ../Background/ALP_BkgModel_param_${version}/fit_results_run2/${massList[$iBin]}/CMS-HGG_mva_13TeV_multipdf.root ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}
        
    ####python makeParametricModelDatacard_ALP_UL_v2.py -i ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_mu.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_mu.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_mu.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_mu.root -s ./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16APV.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_17.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_18.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16APV.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_17.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_18.dat -o ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt --NormSys ./ALP_sys_${version}/NormSys/NormalizationSys_16_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_16APV_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_17_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_18_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_16_mu.txt,./ALP_sys_${version}/NormSys/NormalizationSys_16APV_mu.txt,./ALP_sys_${version}/NormSys/NormalizationSys_17_mu.txt,./ALP_sys_${version}/NormSys/NormalizationSys_18_mu.txt --BDTSys ./ALP_sys_${version}/BDTSys/BDTSys_16_ele.txt,./ALP_sys_${version}/BDTSys/BDTSys_16APV_ele.txt,./ALP_sys_${version}/BDTSys/BDTSys_17_ele.txt,./ALP_sys_${version}/BDTSys/BDTSys_18_ele.txt,./ALP_sys_${version}/BDTSys/BDTSys_16_mu.txt,./ALP_sys_${version}/BDTSys/BDTSys_16APV_mu.txt,./ALP_sys_${version}/BDTSys/BDTSys_17_mu.txt,./ALP_sys_${version}/BDTSys/BDTSys_18_mu.txt -p ggH_16_ele,ggH_16APV_ele,ggH_17_ele,ggH_18_ele,ggH_16_mu,ggH_16APV_mu,ggH_17_mu,ggH_18_mu -c cat0 --isMultiPdf -m 125 --mA ${massList[$iBin]}
    #python makeParametricModelDatacard_ALP_UL_v2.py -i ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_ele.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16_mu.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV_mu.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17_mu.root,./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18_mu.root -s ./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16APV.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_17.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_18.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16APV.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_17.dat,./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_18.dat -o ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt --NormSys ./ALP_sys_${version}/NormSys/NormalizationSys_16_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_16APV_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_17_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_18_ele.txt,./ALP_sys_${version}/NormSys/NormalizationSys_16_mu.txt,./ALP_sys_${version}/NormSys/NormalizationSys_16APV_mu.txt,./ALP_sys_${version}/NormSys/NormalizationSys_17_mu.txt,./ALP_sys_${version}/NormSys/NormalizationSys_18_mu.txt -p ggH_16_ele,ggH_16APV_ele,ggH_17_ele,ggH_18_ele,ggH_16_mu,ggH_16APV_mu,ggH_17_mu,ggH_18_mu -c cat0 --isMultiPdf -m 125 --mA ${massList[$iBin]}

    cd ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}
    pwd
    
    # expected limit
    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125 --rAbsAcc 0.00000001 --rMin -200 --rMax 200 --freezeParameters MH
    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125 --freezeParameters MH
    # observed limit
    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits -m 125 --freezeParameters MH -n _observed

    rm higgsCombine.Test.POINT*
    rm roostats*
    rm M${massList[$iBin]}.*.out
    rm M${massList[$iBin]}.*.err
    rm M${massList[$iBin]}.*.log

    #combineTool.py -d datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --singlePoint 0.01:0.1:0.0002 --saveToys --saveHybridResult -T 1000 --clsAcc 0 --rAbsAcc 0 -m 125 --freezeParameters MH  --job-mode condor --task-name M${massList[$iBin]} --sub-opts='+JobFlavour="testmatch"'
    combineTool.py -d datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --singlePoint 0.005:0.1:0.001 --saveToys --saveHybridResult -T 100000 --clsAcc 0 --rAbsAcc 0 -m 125 --freezeParameters MH  --job-mode condor --task-name M${massList[$iBin]} --sub-opts='+JobFlavour="testmatch"'
    
    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --freezeParameters MH -m 125 --expectedFromGrid=0.5 --name _M${massList[$iBin]}_0.5
    #hadd -f merged.root higgsCombine.Test.POINT*

    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.025 --name _M${massList[$iBin]}_0.025
    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.16 --name _M${massList[$iBin]}_0.16
    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.5 --name _M${massList[$iBin]}_0.5
    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.84 --name _M${massList[$iBin]}_0.84
    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.975 --name _M${massList[$iBin]}_0.975

    #combine datacard_ALPmass${massList[$iBin]}.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --name _M${massList[$iBin]}_Observed


    #text2workspace.py datacard_ALPmass${massList[$iBin]}.txt -m 125 -o datacard_ALPmass${massList[$iBin]}.root

    ### expected Impacts ####
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 1.5 --robustFit 1 --doInitialFit -t -1 --expectSignal 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.5 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 1.5 --robustFit 1 --doInitialFit -t -1 --expectSignal 1 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 1.5 --robustFit 1 --doFits -t -1 --expectSignal 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Expected.json
    #plotImpacts.py -i M${massList[$iBin]}_Expected.json -o M${massList[$iBin]}_Expected_Impact
    
    ### Observed Impacts ####
    #rm -rf ./ALP_SigModel_param_${version}/fit_results_runII/M${massList[$iBin]}/higgsCombine_paramFit*
    
    # M6
    #combine -M MultiDimFit datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin 0 --rMax 0.5 --algo grid --points 50 --robustFit 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 1 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH
    #python plot1DScan.py higgsCombineTest.MultiDimFit.mH125.root -o single_scan
    #python /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/CombineHarvester/CombineTools/scripts/plot1DScan.py higgsCombineTest.MultiDimFit.mH125.root -o single_scan
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin 0 --rMax 1 --robustFit 1 --doInitialFit --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 1 --cminDefaultMinimizerTolerance 0.5 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH #-v 3
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin 0 --rMax 1 --robustFit 1 --doFits --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.3 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Observed.json
    #plotImpacts.py -i M${massList[$iBin]}_Observed.json -o M${massList[$iBin]}_Observed_Impact
    
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -0.1 --rMax 1 --robustFit 1 --doInitialFit --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH #-v 3
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -0.1 --rMax 1 --robustFit 1 --doFits --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125,r=0.01 --freezeParameters MH
    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Observed.json
    #plotImpacts.py -i M${massList[$iBin]}_Observed.json -o M${massList[$iBin]}_Observed_Impact

    #combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --robustFit 1 --doInitialFit --setParameters MH=125 --freezeParameters MH #-v 3
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

    #####combine -M FitDiagnostics datacard_ALPmass${massList[$iBin]}.root -t -1 --expectSignal 0 -m 125
    #combine --method FitDiagnostics datacard_ALPmass${massList[$iBin]}.root  --mass 125 --seed 2 --toys -1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --expectSignal 0 --freezeParameters MH
    #combine --method FitDiagnostics datacard_ALPmass${massList[$iBin]}.root  --mass 125 --seed 2 --toys -1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --expectSignal 1 --freezeParameters MH
    #python /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnosticsTest.root -g plots.root > diffNuisances.log
    
    cd ../../../
    
    
done

#cd ./ALP_SigModel_param_${version}/fit_results_runII
#mkdir plot
#cd plot
#cp ../../../ALP_plot_param_UL_ele/*.py ./
#sed -i "33s/fit_results_run2_/fit_results_/" com_plot.py
#python com_plot.py -fb -y run2 -c runII
