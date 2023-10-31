#!/bin/bash

lable='run2'
version='UL'
#channel=( ele mu )
channel='ele'
Lumi_run2='138'

years=( 16 16APV 17 18 )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

nBins=( 180 200 50 220 220 190 200 200 200 200 150 500 200 150 200 200 300 200 200 150 160 180 50 200 400 200 150 320 180 250 250 200 200 300 250 200 200 150 250 200 350 200 200 250 150 200 130 150 220 230 110 200 200 200 200 200 ) ## Sum-n Gaus

nMass=${#massList[@]}


# prepare signal and background workspace
cd ./InputData/

#python makeWorkspace_data_cats.py -d ./two_jet/data.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 
#python makeWorkspace_sig_cats.py -d ./two_jet/sig.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 
cd ../

###### background fit ######

cd ./Background/

dir_out_bkg="./HZGamma_BkgModel_${version}"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

cat='cat3'
path_in_bkg="../InputData/two_jet_data/HZGamma_data_bkg_workspace_$cat.root"
path_out_bkg="$dir_out_bkg/fit_results_${lable}_$cat"
mkdir $path_out_bkg


#./bin/fTest_ALP_turnOn -i $path_in_bkg --saveMultiPdf $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -D $path_out_bkg/HZGmassInde_fTest -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 105 --mhHigh 170  --mhLowBlind 122 --mhHighBlind 128
#./bin/makeBkgPlots_ALP -b $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -d $path_out_bkg/BkgPlots -o $path_out_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --mhLow 105 --mhHigh 170 --mhLowBlind 122 --mhHighBlind 128 --intLumi 41.48 -c 0 --isFlashgg 0

cd ../Signal/
pwd

dir_out_sig="./HZGamma_SigModel_${version}"

if [ ! -d $dir_out_sig ];then
  mkdir $dir_out_sig
  echo "creat $dir_out_sig"
else
  echo "$dir_out_sig already exist"
fi

path_out_sig="$dir_out_sig/fit_results_${lable}_$cat"
mkdir $path_out_sig
path_in_sig="../InputData/two_jet_data"

#./bin/signalFTest_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_${years[$jBin]}_workspace_${channel}.root -d $path_out_bkg/ALPmassInde_data_sig.dat -o $path_out_bkg/HZAmassInde_ftest -p data -f cat0 -m 125 --mhLowBlind 110 --mhHighBlind 140 # --verbose 1

#sed -i "s/data/ggh/" $path_out_bkg/ALPmassInde_data_sig.dat

#cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_run2_ele/M1/16/ALPmassInde_data_sig.dat $path_out_sig
#./bin/SignalFit_ALP -i $path_in_sig/HZGamma_data_sig_Hm120_workspace_$cat.root,$path_in_sig/HZGamma_data_sig_Hm125_workspace_$cat.root,$path_in_sig/HZGamma_data_sig_Hm130_workspace_$cat.root -o $path_out_sig/CMS-HGG_sigfit_data_ggh_$cat.root -p $path_out_sig/plots_ALP --procs ggh -f cat0 -d $path_out_sig/ALPmassInde_data_sig.dat -s empty.dat --nBins 65 --changeIntLumi 41.48 --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 --verbose 2
#./bin/makeParametricSignalModelPlots_ALP -i $path_out_sig/CMS-HGG_sigfit_data_ggh_$cat.root  -o $path_out_sig/SignalModel -p ggh -f cat0 -L 105 -H 140 -m 125 --binning 70 

###python makeParametricModelDatacardFLASHgg_ALP.py -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -o $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt -p ggh -c cat0 --photonCatScales empty.dat ---p --intLumi ${Lumis[$jBin]} -m 125
###sed -i "45,54d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
###sed -i "21,41d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt


cats=( 'cat0' 'cat1' 'cat2' 'cat3' )
ncats=${#cats[@]}

mkdir "$dir_out_sig/Combine_results"

#for ((iCat=0; iCat<$ncats; iCat++))
#    do

#    cp $dir_out_sig/fit_results_${lable}_${cats[$iCat]}/CMS-HGG_sigfit_data_ggh_${cats[$iCat]}.root $dir_out_sig/Combine_results
#    cp ../Background/HZGamma_BkgModel_${version}/fit_results_${lable}_${cats[$iCat]}/CMS-HGG_mva_13TeV_multipdf_${cats[$iCat]}.root $dir_out_sig/Combine_results
#done

#cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_runII/M1/datacard_ALPmass1.txt $dir_out_sig/Combine_results

cd $dir_out_sig/Combine_results

#for ((iCat=0; iCat<$ncats; iCat++))
#    do

#    text2workspace.py datacard_${cats[$iCat]}.txt -m 125 -o datacard_${cats[$iCat]}.root

#    combine datacard_${cats[$iCat]}.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n ${cats[$iCat]}
#done

combine datacard_$cat.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat
#combine datacard_$cat.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat

#combineCards.py cat0=datacard_cat0.txt cat1=datacard_cat1.txt cat2=datacard_cat2.txt cat3=datacard_cat3.txt &> datacard_allCats.txt
#combine datacard_allCats.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n allCats

#combine -M Significance datacard_allCats.txt -m 200 --rMin -1 --rMax 2