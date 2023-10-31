#!/bin/bash

#lable='16'
#Lumi='35.9'
#lable='17'
#Lumi='41.5'
#lable='18'
#Lumi='56.9'
lable='run2'
channel='ele'
Lumi_run2='137.64'

years=( 16 16APV 17 18 )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

#nBins=( 120 120 120 120 120 120 120 130 120 300 400 120 120 120 120 ) ## Sum-n Gaus
nBins=( 50 50 75 240 80 180 450 200 200 180 200 90 270 200 180 200 90 400 250 350 200 120 120 300 550 50 100 150 200 300 190 190 180 300 270 200 450 50 50 50 120 200 120 200 250 60 200 200 100 200 50 150 150 100 200 120 ) ## Sum-n Gaus
#nBins=550
#nBins=( 200 120 120 120 130 120 300 400 120 120 120 120 ) ## DCB+Gaus
nMass=${#massList[@]}

###### background fit ######

cd ./Background/

dir_out_bkg="./ALP_BkgModel_param_UL"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

path_in_bkg="./ALP_input/fit_${lable}_UL_${channel}"
path_out_bkg="$dir_out_bkg/fit_results_$lable_${channel}"
mkdir $path_out_bkg

for ((iBin=0; iBin<$nMass; iBin++))
    do
    mkdir "$path_out_bkg/${massList[$iBin]}"
    path_bkg="$path_out_bkg/${massList[$iBin]}"

    #./bin/fTest_ALP -i $path_in_bkg/ALP_data_bkg_Am${massList[$iBin]}_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest -c 1 --isFlashgg 0 --isData 0 -f data17, --mhLow 110 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135 --channel ${channel}

    #./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 1 --mhVal 125.0 --mhLow 110 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0  --mhLowBlind 115 --mhHighBlind 135 --channel ${channel}
    #echo "./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 1 --mhVal 125.0 --mhLow 110 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0  --mhLowBlind 115 --mhHighBlind 135"
done

#./bin/makeBkgPlots_ALP -b ./ALP_BkgModel_param/fit_results_run2/1/CMS-HGG_mva_13TeV_multipdf.root -d ./ALP_BkgModel_param/fit_results_run2/1/BkgPlots -o ./ALP_BkgModel_param/fit_results_run2/1/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 1 --mhVal 125.0 --mhLow 110 --mhHigh 180 --intLumi 134.3 -c 0 --isFlashgg 0  --mhLowBlind 115 --mhHighBlind 135 --verbose 2


cd ../Signal/
pwd

dir_out_sig="./ALP_SigModel_param_UL"

if [ ! -d $dir_out_sig ];then
  mkdir $dir_out_sig
  echo "creat $dir_out_sig"
  mkdir "$dir_out_sig/fit_results_$lable_${channel}"
else
  echo "$dir_out_sig already exist"
fi


path_in_sig="../Background/ALP_input/fit_${lable}_UL_${channel}"
mkdir "$dir_out_sig/fit_results_${lable}_${channel}"

for ((iBin=0; iBin<$nMass; iBin++))
    do

    mkdir "$dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}"

    for ((jBin=0; jBin<$nYear; jBin++))
      do
      
      #### single mass point and single year
      if [ $# -ne 0 ];then
        if [[ ${massList[ $iBin ]} -ne $1 ]];then
          continue
        fi
        if [[ ${years[ $jBin ]} != $2 ]];then
          continue
        fi
      fi


      mkdir "$dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}/${years[$jBin]}"

      path_out_bkg="$dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}/${years[$jBin]}"

      #./bin/signalFTest_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_${years[$jBin]}_workspace.root -d $path_out_bkg/ALPmassInde_data_sig.dat -o $path_out_bkg/HZAmassInde_ftest -p data -f cat0 -m 125 --mhLowBlind 110 --mhHighBlind 140  --verbose 1

      #sed -i "s/data/ggh/" $path_out_bkg/ALPmassInde_data_sig.dat


      #echo "./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins ${nBins[$[iBin*nYear] + $jBin]} --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 --useDCBplusGaus 1 "
      #./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins ${nBins[$[iBin*nYear] + $jBin]} --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 #--verbose 2
      #./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins $nBins --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 #--verbose 2
      #./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s ./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_${channel}/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_${years[$jBin]}.dat --nBins ${nBins[$[iBin*nYear] + $jBin]} --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 #--verbose 2
      #./bin/makeParametricSignalModelPlots_ALP -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root  -o $path_out_bkg/SignalModel -p ggh -f cat0 -L 100 -H 135 -m 125 #--binning ${nBins[$[iBin*nYear] + $jBin]} #--useDCBplusGaus 1

      ###python makeParametricModelDatacardFLASHgg_ALP.py -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -o $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt -p ggh -c cat0 --photonCatScales empty.dat ---p --intLumi ${Lumis[$jBin]} -m 125
      ###sed -i "45,54d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
      ###sed -i "21,41d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
      done

    cp ./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/16/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_${lable}_${channel}/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16.root
    cp ./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/16APV/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_${lable}_${channel}/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_16APV.root
    cp ./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/17/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_${lable}_${channel}/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_17.root
    cp ./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/18/CMS-HGG_sigfit_data_ggh_cat0.root ./ALP_SigModel_param_UL/fit_results_${lable}_${channel}/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0_18.root

    python makeParametricModelDatacard_ALP.py -i ./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/16/CMS-HGG_sigfit_data_ggh_cat0.root,./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/16APV/CMS-HGG_sigfit_data_ggh_cat0.root,./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/17/CMS-HGG_sigfit_data_ggh_cat0.root,./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/18/CMS-HGG_sigfit_data_ggh_cat0.root -s ./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_${channel}/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_16.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_${channel}/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_17.dat,./ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_${channel}/ScaleSmear_plot/ScaleSmear_m${massList[$iBin]}_18.dat -o ./ALP_SigModel_param_UL/fit_results_run2_${channel}/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt --NormSys ./ALP_sys_UL/NormSys/NormalizationSys_16_${channel}.txt,./ALP_sys_UL/NormSys/NormalizationSys_16APV_${channel}.txt,./ALP_sys_UL/NormSys/NormalizationSys_17_${channel}.txt,./ALP_sys_UL/NormSys/NormalizationSys_18_${channel}.txt --BDTSys ./ALP_sys_UL/BDTSys/BDTSys_16_${channel}.txt,./ALP_sys_UL/BDTSys/BDTSys_16APV_${channel}.txt,./ALP_sys_UL/BDTSys/BDTSys_17_${channel}.txt,./ALP_sys_UL/BDTSys/BDTSys_18_${channel}.txt --channel ${channel} -p ggH_16,ggH_16APV,ggH_17,ggH_18 -c cat0 --isMultiPdf -m 125 --mA ${massList[$iBin]}

    cp ../Background/ALP_BkgModel_param_UL/fit_results_${channel}/${massList[$iBin]}/CMS-HGG_mva_13TeV_multipdf.root $dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}

    cd $dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}
    combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001

    text2workspace.py datacard_ALPmass${massList[$iBin]}.txt -m 125 -o datacard_ALPmass${massList[$iBin]}.root

    combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 10 --robustFit 1 --doInitialFit -t -1 --expectSignal 1 --cminDefaultMinimizerStrategy 0 --setParameters MH=125 --freezeParameters MH
    combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 --rMin -1 --rMax 10 --robustFit 1 --doFits -t -1 --expectSignal 1 --cminDefaultMinimizerStrategy 0 --setParameters MH=125 --freezeParameters MH
    combineTool.py -M Impacts -d datacard_ALPmass${massList[$iBin]}.root -m 125 -o M${massList[$iBin]}_Expected.json
    plotImpacts.py -i M${massList[$iBin]}_Expected.json -o M${massList[$iBin]}_Expected_Impact

    cd ../../../

      

done

cd ALP_plot_param_UL_${channel}
python com_plot.py -fb -y ${lable} -c ${channel}
mv ALP_xs_UpperLimit.png ALP_xs_UpperLimit_${lable}.png
mv ALP_Br_UpperLimit.png ALP_Br_UpperLimit_${lable}.png
