#!/bin/bash
basic_path=`pwd`
lable='run2'
version='UL'
#channel=( ele mu )
channel='ele'
Lumi_run2='138'
cat='cat7'

years=( 16 16APV 17 18 )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

nBins=( 180 200 50 220 220 190 200 200 200 200 150 500 200 150 200 200 300 200 200 150 160 180 50 200 400 200 150 320 180 250 250 200 200 300 250 200 200 150 250 200 350 200 200 250 150 200 130 150 220 230 110 200 200 200 200 200 ) ## Sum-n Gaus

nMass=${#massList[@]}


# prepare signal and background workspace
cd $basic_path/InputData/

# python makeWorkspace_data_cats.py -d ./two_jet/data.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 
# python makeWorkspace_sig_cats.py -d ./two_jet/sig.root -j ./significances/bin_binaries_two_jet.txt -o ./two_jet_data 

# python makeWorkspace_data_cats.py -d /eos/home-j/jiehan/root/outputs/two_jet/Data.root -j /eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt -o ./two_jet_data 
# python makeWorkspace_sig_cats.py -d /eos/home-j/jiehan/root/outputs/two_jet/sig.root -j /eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt -o ./two_jet_data 

# python makeWorkspace_data_cats.py -d /eos/home-j/jiehan/root/outputs/two_jet/Data.root -j /eos/home-j/jiehan/root/outputs/significances/bin_binaries_2D_two_jet.txt -o ./two_jet_data -mb True
# python makeWorkspace_sig_cats.py -d /eos/home-j/jiehan/root/outputs/two_jet/sig.root -j /eos/home-j/jiehan/root/outputs/significances/bin_binaries_2D_two_jet.txt -o ./two_jet_data -mb True
cd $basic_path/

###### background fit ######

cd $basic_path/Background/

dir_out_bkg="./HZGamma_BkgModel_${version}"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

cats=( 'cat0' 'cat1' 'cat2' 'cat3' 'cat4' 'cat5' 'cat6' 'cat7' )
ncats=${#cats[@]}

# make clean; make;

# for ((iCat=0; iCat<$ncats; iCat++))
# do
# cat=${cats[$iCat]}

# # TODO: for this part, you must run one cat each time and tune the parameters
# cat="cat2"
# path_in_bkg="../InputData/two_jet_data/HZGamma_data_bkg_workspace_$cat.root"
# path_out_bkg="$dir_out_bkg/fit_results_${lable}_$cat"
# mkdir $path_out_bkg

# ./bin/fTest_ALP_turnOn -i $path_in_bkg --saveMultiPdf $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -D $path_out_bkg/HZGmassInde_fTest -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 100 --mhHigh 180  --mhLowBlind 122 --mhHighBlind 128 #--verbose 2

# ./bin/makeBkgPlots_ALP -b $path_out_bkg/CMS-HGG_mva_13TeV_multipdf_$cat.root -d $path_out_bkg/BkgPlots -o $path_out_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --mhLow 100 --mhHigh 180 --mhLowBlind 122 --mhHighBlind 128 --intLumi 137 -c 0 --isFlashgg 0

# done

cd $basic_path/Signal/

dir_out_sig="./HZGamma_SigModel_${version}"

if [ ! -d $dir_out_sig ];then
  mkdir $dir_out_sig
  echo "creat $dir_out_sig"
else
  echo "$dir_out_sig already exist"
fi

# make clean; make;

# for ((iCat=0; iCat<$ncats; iCat++))
# do
# cat=${cats[$iCat]}

# path_out_sig="$dir_out_sig/fit_results_${lable}_$cat"
# mkdir $path_out_sig
# path_in_sig="../InputData/two_jet_data"

# # ./bin/signalFTest_ALP -i $path_in_sig/ALP_data_sig_Am_${years[$jBin]}_workspace_${channel}.root -d $path_out_bkg/allCatsInde_data_sig.dat -o $path_out_bkg/HZAmassInde_ftest -p data -f cat0 -m 125 --mhLowBlind 100 --mhHighBlind 180 # --verbose 1

# #sed -i "s/data/ggh/" $path_out_bkg/allCatsInde_data_sig.dat

# # cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_run2_ele/M1/16/allCatsInde_data_sig.dat $path_out_sig
# # cp /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/Signal/HZGamma_SigModel_UL/fit_results_run2_cat1/allCatsInde_data_sig.dat $path_out_sig
# ./bin/SignalFit_ALP -i $path_in_sig/HZGamma_data_sig_Hm120_workspace_$cat.root,$path_in_sig/HZGamma_data_sig_Hm125_workspace_$cat.root,$path_in_sig/HZGamma_data_sig_Hm130_workspace_$cat.root -o $path_out_sig/CMS-HGG_sigfit_data_ggh_$cat.root -p $path_out_sig/plots_ALP --procs ggh -f cat0 -d $path_out_sig/allCatsInde_data_sig.dat -s empty.dat --nBins 80 --changeIntLumi 137 --useSSF 1 -L 100 --mhHigh 180 --massList 120,125,130 #--useDCBplusGaus 1 --verbose 2
# ./bin/makeParametricSignalModelPlots_ALP -i $path_out_sig/CMS-HGG_sigfit_data_ggh_$cat.root  -o $path_out_sig/SignalModel -p ggh -f cat0 -L 100 -H 180 -m 125 --binning 80 

# ###python makeParametricModelDatacardFLASHgg_ALP.py -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -o $path_out_bkg/datacard_allCats.txt -p ggh -c cat0 --photonCatScales empty.dat ---p --intLumi ${Lumis[$jBin]} -m 125
# ###sed -i "45,54d" $path_out_bkg/datacard_allCats.txt
# ###sed -i "21,41d" $path_out_bkg/datacard_allCats.txt
# done

mkdir "$dir_out_sig/Combine_results"

# for ((iCat=0; iCat<$ncats; iCat++))
#    do

#    cp $dir_out_sig/fit_results_${lable}_${cats[$iCat]}/CMS-HGG_sigfit_data_ggh_${cats[$iCat]}.root $dir_out_sig/Combine_results
#    cp $basic_path/Background/HZGamma_BkgModel_${version}/fit_results_${lable}_${cats[$iCat]}/CMS-HGG_mva_13TeV_multipdf_${cats[$iCat]}.root $dir_out_sig/Combine_results
# done

# cp $basic_path/Background/HZGamma_BkgModel_${version}/fit_results_${lable}_$cat/CMS-HGG_mva_13TeV_multipdf_$cat.root $dir_out_sig/Combine_results

# cp /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_runII/M1/datacard_allCats1.txt $dir_out_sig/Combine_results

cd $dir_out_sig/Combine_results
# cd ../Combine_results

# for ((iCat=0; iCat<$ncats; iCat++))
#   do

#    text2workspace.py datacard_${cats[$iCat]}.txt -m 125.0 -o datacard_${cats[$iCat]}.root

#    combine datacard_${cats[$iCat]}.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n ${cats[$iCat]}
#    combine -M Significance datacard_${cats[$iCat]}.txt -t -1 --expectSignal=1 -m 125.0
#    combine -d datacard_${cats[$iCat]}.root -n .${cats[$iCat]}_bestfit -M MultiDimFit --algo grid --points 20 --noMCbonly 1  --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100  --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_MaxCalls=9999999  --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH  --X-rtd MINIMIZER_freezeDisassociatedParams  --X-rtd OPTIMIZE_BOUNDS=0  -v 3  --setParameterRanges r=-10,10 -m 125.0  --floatOtherPOIs 1 --X-rtd NO_INITIAL_SNAP  --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2  --saveWorkspace --alignEdges 1 -t -1 --expectSignal=1 #--toysFrequentist

#    plot1DScan.py higgsCombine.${cats[$iCat]}_bestfit.MultiDimFit.mH125.root -o single_scan_${cats[$iCat]} --y-max 2.5 --main-label Expected #--others
# done

# text2workspace.py datacard_$cat.txt -m 125.0 -o datacard_$cat.root
# combine -d datacard_$cat.root -n .cat0_bestfit -M MultiDimFit --algo grid --points 20 --noMCbonly 1  --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100  --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_MaxCalls=9999999  --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH  --X-rtd MINIMIZER_freezeDisassociatedParams  --X-rtd OPTIMIZE_BOUNDS=0  -v 3  --setParameterRanges r=-10,10 -m 125.0  --floatOtherPOIs 1 --X-rtd NO_INITIAL_SNAP  --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2  --saveWorkspace --alignEdges 1 -t -1 --expectSignal=1 #--toysFrequentist

# plot1DScan.py higgsCombine.${cat}_bestfit.MultiDimFit.mH125.root -o single_scan_cat0 --y-max 2.5 --main-label Expected #--others

# combine datacard_$cat.txt -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat
# combine -M Significance datacard_$cat.txt -t -1 --expectSignal=1 -m 125.0
#combine datacard_$cat.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat

combineCards.py cat0=datacard_cat0.txt cat1=datacard_cat1.txt cat2=datacard_cat2.txt cat3=datacard_cat3.txt cat4=datacard_cat4.txt cat5=datacard_cat5.txt cat6=datacard_cat6.txt cat7=datacard_cat7.txt &> datacard_allCats.txt
combine datacard_allCats.txt -M AsymptoticLimits --run=blind --rAbsAcc 0.00000001 -n allCats -m 125.0

text2workspace.py datacard_allCats.txt -m 125.0 -o datacard_allCats.root

# combineTool.py -d datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --singlePoint 0.01:0.2:0.002 --saveToys --saveHybridResult -T 100 --clsAcc 0 --rAbsAcc 0.000000001 -m 125 --freezeParameters MH # --job-mode condor --task-name M --sub-opts='+JobFlavour="testmatch"'
# combineTool.py -d datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --singlePoint 0.07:0.6:0.001 --saveToys --saveHybridResult -T 1000 --clsAcc 0 --rAbsAcc 0 -m 125 --freezeParameters MH  --job-mode condor --task-name M --sub-opts='+JobFlavour="testmatch"'
    
# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --freezeParameters MH -m 125 --expectedFromGrid=0.5 --name _M_0.5
# hadd -f merged.root higgsCombine.Test.POINT*

# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.025 --name _M_0.025
# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.16 --name _M_0.16
# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.5 --name _M_0.5
# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.84 --name _M_0.84
# combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --expectedFromGrid=0.975 --name _M_0.975

   #combine datacard_allCats.txt -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged.root --freezeParameters MH -m 125 --name _M_Observed

combine -d datacard_allCats.root -n .bestfit -M MultiDimFit --algo grid --points 31 --noMCbonly 1  --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100  --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_MaxCalls=9999999  --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH  --X-rtd MINIMIZER_freezeDisassociatedParams  --X-rtd OPTIMIZE_BOUNDS=0  -v 3  --setParameterRanges r=-0.5,2.5 -m 125  --floatOtherPOIs 1 --X-rtd NO_INITIAL_SNAP  --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2  --saveWorkspace --alignEdges 1 -t -1 --expectSignal=1 #--toysFrequentist

plot1DScan.py higgsCombine.bestfit.MultiDimFit.mH125.root -o single_scan --y-max 2.5 --main-label Expected #--others

# combine datacard_allCats.txt -M Significance -m 125.0 --rMin 0.1 --rMax 3
combine -M Significance datacard_allCats.txt -t -1 --expectSignal=1 -m 125.0
# # combine -M Significance datacard_cat0.txt -t -1 --expectSignal=2 -m 125.0

cd $basic_path/
# cd ../
