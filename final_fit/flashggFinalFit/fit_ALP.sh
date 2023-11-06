#!/bin/bash

#year='16'
#Lumi='36.9'
#year='17'
#Lumi='41.5'
year='18'
Lumi='56.9'

cd ./Background/

dir_out_bkg="./ALP_BkgModel"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

path_in_bkg="./ALP_input/fit_$year"
path_out_bkg="$dir_out_bkg/fit_results_$year"
mkdir $path_out_bkg

#./bin/fTest_ALP -i $path_in_bkg/ALP_data_bkg_workspace.root --saveMultiPdf $path_out_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_out_bkg/HZAmassInde_fTest -c 1 --isFlashgg 0 --isData 0 -f data17, --mhLow 0.8 --mhHigh 40 --mhLowBlind 5 --mhHighBlind 30

#./bin/makeBkgPlots_ALP -b $path_out_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_out_bkg/BkgPlots -o $path_out_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 1 --mhVal 15 --mhLow 0 --mhHigh 40 --intLumi $Lumi -c 0 --isFlashgg 0 --mhLowBlind 5 --mhHighBlind 30

cd ../Signal/
pwd

dir_out_sig="./ALP_SigModel"

if [ ! -d $dir_out_sig ];then
  mkdir $dir_out_sig
  echo "creat $dir_out_sig"
  mkdir "$dir_out_sig/fit_results_$year"
else
  echo "$dir_out_sig already exist"
fi


massList=( 5 15 30 )
massList_up=( 6 16 31 )
massList_dn=( 4 14 29 )
nBins=( 400 320 150)
nMass=${#massList[@]}


mhLowBlind=( 3 10 25 )
mhHighBlind=( 7 20 35 )

path_in_sig="../Background/ALP_input/fit_$year"

for ((iBin=0; iBin<$nMass; iBin++))
    do

    mkdir "$dir_out_sig/fit_results_$year"
    mkdir "$dir_out_sig/fit_results_$year/M${massList[$iBin]}"
    path_out_bkg="$dir_out_sig/fit_results_$year/M${massList[$iBin]}"

    #./bin/signalFTest_ALP -i $path_in_sig/ALP_data_sig_workspace_M${massList[$iBin]}.root -d $path_out_bkg/ALPmassInde_data_sig.dat -o $path_out_bkg/HZAmassInde_ftest -p data -f cat0 -m ${massList[$iBin]} --mhLowBlind ${mhLowBlind[$iBin]} --mhHighBlind ${mhHighBlind[$iBin]}

    #sed -i "s/data/ggh/" $path_out_bkg/ALPmassInde_data_sig.dat

    ./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_m${massList_dn[$iBin]}_workspace.root,$path_in_sig/ALP_data_sig_m${massList[$iBin]}_workspace.root,$path_in_sig/ALP_data_sig_m${massList_up[$iBin]}_workspace.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins ${nBins[$iBin]} --changeIntLumi $Lumi --useSSF 1 -L ${mhLowBlind[$iBin]} --mhHigh ${mhHighBlind[$iBin]} --massList ${massList_dn[$iBin]},${massList[$iBin]},${massList_up[$iBin]} -v 2
    #./bin/makeParametricSignalModelPlots_ALP -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root  -o $path_out_bkg/SignalModel -p ggh -f cat0 -L ${mhLowBlind[$iBin]} -H ${mhHighBlind[$iBin]} -m ${massList[$iBin]}

    #python makeParametricModelDatacardFLASHgg_ALP.py -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -o $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt -p ggh -c cat0 --photonCatScales empty.dat --isMultiPdf --intLumi $Lumi -m ${massList[$iBin]}
    #sed -i "45,54d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
    #sed -i "21,41d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt

    #cp ../Background/ALP_BkgModel/fit_results_$year/CMS-HGG_mva_13TeV_multipdf.root $path_out_bkg

    cd $path_out_bkg
    #combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m ${massList[$iBin]} --rAbsAcc 0.00000001

    #text2workspace.py datacard_ALPmass${massList[$iBin]}.txt -m ${massList[$iBin]} -o datacard_ALPmass${massList[$iBin]}.root

    cd ../../../

done

#cd ALP_plot
#python com_plot.py -fb
