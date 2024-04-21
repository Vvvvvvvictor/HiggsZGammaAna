#!/bin/bash

lable='run2'
#version='UL_run2_01jet'
version='UL_run2_2jet_v2'
# version='UL_xingchen_nll_SB'
#channel=( ele mu )
channel='ele'
Lumi_run2='138'

years=( 16 16APV 17 18 )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )


# prepare signal and background workspace
cd ./InputData/
echo "Preparing workspace!!!!"

#sync xingchen
#cd ./InputData/xingchen/
#python makeWorkspace_data_cats.py 

#python makeWorkspace_data_cats.py -d ./two_jet/data.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 
#python makeWorkspace_sig_cats.py -d ./two_jet/sig.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 
#python makeWorkspace_bkg_cats.py -d ./two_jet/bkgmc.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 

#python makeWorkspace_data_cats.py -d ./finalfit/zero_to_one_jet/data.root -j ./finalfit/significance/bin_binaries_zero_to_one_jet.txt -o ./two_jet_data_01jet
#python makeWorkspace_sig_cats.py -d ./finalfit/zero_to_one_jet/sig.root -j ./finalfit/significance/bin_binaries_zero_to_one_jet.txt -o ./two_jet_data_01jet 

# python makeWorkspace_data_cats.py -t test -d ./jet2_run2/two_jet/data.root -j ./jet2_run2/significances/bin_binaries_1D_two_jet.txt -o ./jet2_run2 
# python makeWorkspace_sig_cats_vbf.py -t test -d ./jet2_run2/two_jet/sig.root -j ./jet2_run2/significances/bin_binaries_1D_two_jet.txt -o ./jet2_run2 

#python makeWorkspace_data_cats.py -t zero_to_one_jet -d /eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples/data/zero_to_one_jet_run2.root -j ./significances/bin_binaries_01_jet_run2.txt -o ./jet01_run2 
#python makeWorkspace_sig_cats.py -t zero_to_one_jet -d /eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples/sig/zero_to_one_jet_run2.root -j ./significances/bin_binaries_01_jet_run2.txt -o ./jet01_run2 
#python makeWorkspace_sig_cats_refit.py -t zero_to_one_jet -d /eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples/sig/zero_to_one_jet_run2_refit.root -j ./significances/bin_binaries_01_jet_run2.txt -o ./jet01_run2_refit -c ele

cd ../

###### background fit ######

cd ./Background/
# make clean; make -j16;

echo "Doing background fitting!!!"

dir_out_bkg="./HZGamma_BkgModel_${version}"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

cat_name='vbf'
# cat='cat3'
#cats=( 'cat0_mu' 'cat1_mu' 'cat2_mu' 'cat3_mu' 'cat0_mu_refit' 'cat1_mu_refit' 'cat2_mu_refit' 'cat3_mu_refit' ) 
cats=( 'cat0' 'cat1' 'cat2' 'cat3' ) 
#cats=( 'cat0' ) 
ncats=${#cats[@]}

#path_in_bkg="../InputData/two_jet_data_dijet/HZGamma_data_bkg_workspace_$cat.root"
# sync
# path_in_bkg="../InputData/xingchen/output_file/HZGamma_data_bkg_workspace_$cat.root"
# 01 jet
#path_in_bkg="../InputData/jet01_run2/HZGamma_data_bkg_workspace_$cat.root"
# 2 jets

for ((iCat=0; iCat<$ncats; iCat++))
  do
  cat=${cats[$iCat]}
  path_in_bkg="../InputData/jet2_run2/HZGamma_data_bkg_workspace_$cat.root"
  path_out_bkg="$dir_out_bkg/fit_results_${lable}_$cat"
  mkdir $path_out_bkg
  echo "Path in bkg: " $path_in_bkg
  echo "Path out bkg: " $path_out_bkg

  ./bin/fTest_ALP_turnOn -i $path_in_bkg --saveMultiPdf $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -D $path_out_bkg/HZGmassInde_fTest -c $(($iCat+1)) -n 1 --channel ${cat_name}_${cat} --isFlashgg 0 --isData 0 -f data, --mhLow 105 --mhHigh 170  --mhLowBlind 122 --mhHighBlind 128 --sidebandOnly #--chi2fit #--verbose 2 #--runFtestCheckWithToys --sidebandOnly

  # ./bin/fTest_ALP_turnOn -i $path_in_bkg --saveMultiPdf $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -D $path_out_bkg/HZGmassInde_fTest -c 1 --channel ${cat_name}_${cat} --isFlashgg 0 --isData 0 -f data, --mhLow 105 --mhHigh 170  --mhLowBlind 122 --mhHighBlind 128 #--chi2fit #--verbose 2 #--runFtestCheckWithToys #--sidebandOnly
  ./bin/makeBkgPlots_ALP -b $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -d $path_out_bkg/BkgPlots -o $path_out_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --mhVal 125.0 --mhLow 105 --mhHigh 170 --mhLowBlind 122 --mhHighBlind 128 --intLumi 138 -c 0 --channel ${cat_name}_${cat} --isFlashgg 0 --doBands --massStep 2 #--unblind
done

cd ../Signal/
# make clean; make -j16;
pwd

dir_out_sig="./HZGamma_SigModel_${version}"

if [ ! -d $dir_out_sig ];then
  mkdir $dir_out_sig
  echo "creat $dir_out_sig"
else
  echo "$dir_out_sig already exist"
fi

echo "dir out sig: " $dir_out_sig

for ((iCat=0; iCat<$ncats; iCat++))
  do
  cat=${cats[$iCat]}
  path_out_sig="$dir_out_sig/fit_results_${lable}_$cat"
  mkdir $path_out_sig
  path_in_sig="../InputData/jet2_run2"

  #./bin/signalFTest_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_${years[$jBin]}_workspace_${channel}.root -d $path_out_bkg/ALPmassInde_data_sig.dat -o $path_out_bkg/HZAmassInde_ftest -p data -f cat0 -m 125 --mhLowBlind 110 --mhHighBlind 140 # --verbose 1

  #sed -i "s/data/ggh/" $path_out_bkg/ALPmassInde_data_sig.dat

  cp /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/Signal/HZGamma_SigModel_UL/fit_results_run2_cat0/ALPmassInde_data_sig.dat $path_out_sig
  sed -i "s/cat0/${cat_name}_${cat}/" $path_out_sig/ALPmassInde_data_sig.dat

  ./bin/SignalFit_ALP -i $path_in_sig/HZGamma_data_sig_Hm120_workspace_$cat.root,$path_in_sig/HZGamma_data_sig_Hm125_workspace_$cat.root,$path_in_sig/HZGamma_data_sig_Hm130_workspace_$cat.root -o $path_out_sig/CMS-HGG_sigfit_data_ggh_$cat.root -p $path_out_sig/plots_ALP --procs ggh -f ${cat_name}_${cat} -d $path_out_sig/ALPmassInde_data_sig.dat -s empty.dat --nBins 65 --changeIntLumi 138 --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 --verbose 2 #--checkYields 1 #--useDCBplusGaus 1 --verbose 2
  ./bin/makeParametricSignalModelPlots_ALP -i $path_out_sig/CMS-HGG_sigfit_data_ggh_$cat.root  -o $path_out_sig/SignalModel -p ggh -f ${cat_name}_${cat} -L 105 -H 140 -m 125 --binning 70 

  ###python makeParametricModelDatacardFLASHgg_ALP.py -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -o $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt -p ggh -c cat0 --photonCatScales empty.dat ---p --intLumi ${Lumis[$jBin]} -m 125
  ###sed -i "45,54d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
  ###sed -i "21,41d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
done

#cats=( 'cat0' 'cat1' 'cat2' 'cat3' )
# cats=( 'cat1' )
# ncats=${#cats[@]}

mkdir "$dir_out_sig/Combine_results"

for ((iCat=0; iCat<$ncats; iCat++))
   do

   cp $dir_out_sig/fit_results_${lable}_${cats[$iCat]}/CMS-HGG_sigfit_data_ggh_${cats[$iCat]}.root $dir_out_sig/Combine_results
   cp ../Background/HZGamma_BkgModel_${version}/fit_results_${lable}_${cats[$iCat]}/CMS-HGG_mva_13TeV_multipdf_${cats[$iCat]}.root $dir_out_sig/Combine_results
done

echo "Using combine to get the fitting result"
cd $dir_out_sig/Combine_results

for ((iCat=0; iCat<$ncats; iCat++))
    do

    text2workspace.py datacard_${cats[$iCat]}.txt -m 125 -o datacard_${cats[$iCat]}.root

    combine datacard_${cats[$iCat]}.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n ${cats[$iCat]}
    combine datacard_${cats[$iCat]}.txt -M Significance -t -1 --expectSignal=1 -m 125.0 -n ${cats[$iCat]} #--freezeParameters allConstrainedNuisances

    combine -M MultiDimFit datacard_${cats[$iCat]}.root -m 125 -t -1 --expectSignal=1 --rMin -5 --rMax 5 --algo grid --points 50 --robustFit 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 1 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125 --freezeParameters MH -n ${cats[$iCat]}
    plot1DScan.py higgsCombine${cats[$iCat]}.MultiDimFit.mH125.root -o single_scan_${cats[$iCat]} --main-label Expected

    # # start GoF
    # combine -M GoodnessOfFit datacard_${cats[$iCat]}.txt --algo=saturated -m 125 --setParameters MH=125 -n _${cats[$iCat]}_saturated
    # combine -M GoodnessOfFit datacard_${cats[$iCat]}.txt --algo=saturated -t 500 -s 12345 -m 125 --setParameters MH=125 -n _${cats[$iCat]}_saturated

    # combine -M GoodnessOfFit datacard_${cats[$iCat]}.txt --algo=KS -m 125 --setParameters MH=125 -n _${cats[$iCat]}_KS
    # combine -M GoodnessOfFit datacard_${cats[$iCat]}.txt --algo=KS -t 100 -s 12345 -m 125 --setParameters MH=125 -n _${cats[$iCat]}_KS

    # python /afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/runGoF.py > GoF.log
done

combine datacard_$cat.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat
combine datacard_$cat.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat

combineCards.py cat0=datacard_cat0.txt cat1=datacard_cat1.txt cat2=datacard_cat2.txt cat3=datacard_cat3.txt &> datacard_allCats.txt
combineCards.py cat0=datacard_cat2.txt cat1=datacard_cat3.txt cat2=datacard_cat0.txt cat3=datacard_cat1.txt &> datacard_allCats.txt
text2workspace.py datacard_allCats.txt -m 125 -o datacard_allCats.root
combine datacard_allCats.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n allCats #--freezeParameters allConstrainedNuisances
combine datacard_allCats.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n allCats #--freezeParameters allConstrainedNuisances #-v 2
# combine -M MultiDimFit datacard_allCats.root -m 125 -t -1 --expectSignal=1 --rMin -5 --rMax 5 --algo grid --points 50 --robustFit 1 --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 1 --cminDefaultMinimizerTolerance 0.1 --cminFallbackAlgo Minuit2,0:0.2 --cminFallbackAlgo Minuit2,0:0.4 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --setParameters MH=125 --freezeParameters allConstrainedNuisances -n allCats
# plot1DScan.py higgsCombineallCats.MultiDimFit.mH125.root -o single_scan_allCats --main-label Expected

# text2workspace.py datacard_allCats.txt -m 125 -o datacard_allCats.root
# combineCards.py cat0=datacard_cat1.txt cat1=datacard_cat2.txt cat2=datacard_cat3.txt cat3=datacard_cat0.txt  &> datacard_allCats_test.txt
# combine datacard_allCats_test.txt -M Significance -t -1 --expectSignal=1 -m 125.0 -n test
# combine  datacard_allCats.root -M FitDiagnostics -t -1 --expectSignal 1

# combineTool.py -d datacard_allCats.txt -M HybridNew --LHCmode LHC-significance --saveToys --fullBToys --expectSignal 1 --saveHybridResult -T 1000 -m 125 --freezeParameters MH --job-mode condor --task-name HZGamma --sub-opts='+JobFlavour="testmatch"'
# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-significance --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.5 --name _HZG_0.5
# combine -M HybridNew datacard_allCats.txt --LHCmode LHC-significance  --saveToys --fullBToys --expectSignal 1 --saveHybridResult -m 125 -T 10 -i 1 -s 123456
# combine -M Significance/ datacard_allCats.txt -m 200 --rMin -1 --rMax 2


### plot limits

#python ../../python/com_plot.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/HZGamma_SigModel_UL/Combine_results -o /afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/HZGamma_SigModel_UL/Combine_results