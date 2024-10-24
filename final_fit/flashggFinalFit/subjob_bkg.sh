#!/bin/bash

cmsenv

lable='run2'
version='UL'
#channel=( ele mu )
channel='ele'
#Lumi_run2='137.64'
Lumi_run2='138'

years=( 16 16APV 17 18 )
#years=( 16 16APV )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# massList=( 10 )

nMass=${#massList[@]}

###### background fit ######

cd ./Background/

dir_out_bkg="./ALP_BkgModel_param_${version}"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

path_in_bkg="../MVAcut/${lable}_${version}/output/data"

path_out_bkg="$dir_out_bkg/fit_results_${lable}"
mkdir -p $path_out_bkg

path_run="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit"
# path_bkg="$path_out_bkg/2"
# ./bin/fTest_ALP_turnOn -i $path_in_bkg/ALP_data_bkg_Am2_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP 2 -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
# exit 


for ((iBin=0; iBin<$nMass; iBin++))
  do
  mkdir -p "$path_out_bkg/${massList[$iBin]}"
  path_bkg="$path_out_bkg/${massList[$iBin]}"

  hep_sub $path_run/runjob_bkg.sh -g cms -mem 2000 -wt mid -o $path_bkg/job${massList[$iBin]}.log -e $path_bkg/job${massList[$iBin]}.err -argu $path_in_bkg/ALP_data_bkg_Am${massList[$iBin]}_workspace.root $path_bkg/CMS-HGG_mva_13TeV_multipdf.root $path_bkg/HZAmassInde_fTest ${massList[$iBin]} 1 0 0 data,
  done