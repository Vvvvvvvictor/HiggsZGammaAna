#!/bin/bash

lable='run2'
Lumi_run2='137.64'

years=( 16 16APV 17 18 )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 )
#massList=( 16 )
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

path_in_bkg="./ALP_input/fit_${lable}_UL"
path_out_bkg="$dir_out_bkg/fit_results_$lable"
mkdir $path_out_bkg

for ((iBin=0; iBin<$nMass; iBin++))
    do
    mkdir "$path_out_bkg/${massList[$iBin]}"
    path_bkg="$path_out_bkg/${massList[$iBin]}"

    #./bin/fTest_ALP_turnOn -i $path_in_bkg/ALP_data_bkg_Am${massList[$iBin]}_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135
    #./bin/fTest_ALP -i $path_in_bkg/ALP_data_bkg_Am${massList[$iBin]}_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest -c 1 --isFlashgg 0 --isData 0 -f data17, --mhLow 110 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135

    #./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 2 --mhVal 125.0 --mhLow 95 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0  --mhLowBlind 115 --mhHighBlind 135
    ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --unblind --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0
    #echo "./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 1 --mhVal 125.0 --mhLow 110 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0  --mhLowBlind 115 --mhHighBlind 135"
done

#./bin/makeBkgPlots_ALP -b ./ALP_BkgModel_param/fit_results_run2/1/CMS-HGG_mva_13TeV_multipdf.root -d ./ALP_BkgModel_param/fit_results_run2/1/BkgPlots -o ./ALP_BkgModel_param/fit_results_run2/1/BkgPlots.root -S 13 --isMultiPdf --useBinnedData  --doBands --massStep 1 --mhVal 125.0 --mhLow 110 --mhHigh 180 --intLumi 134.3 -c 0 --isFlashgg 0  --mhLowBlind 115 --mhHighBlind 135 --verbose 2



#cd ALP_plot_param_UL_${channel}
#python com_plot.py -fb -y ${lable} -c ${channel}
#mv ALP_xs_UpperLimit.png ALP_xs_UpperLimit_${lable}.png
#mv ALP_Br_UpperLimit.png ALP_Br_UpperLimit_${lable}.png
