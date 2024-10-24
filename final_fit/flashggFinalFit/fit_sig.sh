#!/bin/bash

cmsenv

lable='run2'
version='UL'
#channel=( ele mu )
channel='ele'
#Lumi_run2='137.64'
Lumi_run2='138'

Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

ChannelList=( ele mu lep )
nChannel=${#ChannelList[@]}

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
nALPmass=${#ALPmassList[@]}

YearsList=( 16 16APV 17 18 )
nYear=${#YearsList[@]}

HmassList=( 120 125 130 )
nHmass=${#HmassList[@]}

#massList=( $1 )

#massList=( 1 )

# nbins 5GeV 3.7 1GeV 3,5

#nBins=( 120 120 120 120 120 120 120 130 120 300 400 120 120 120 120 ) ## Sum-n Gaus
#nBins=( 80 350 100 300 250 190 200 200 200 200 150 50 180 200 150 200 200 200 200 150 200 180 155 200 80 250 150 320 120 230 200 200 200 300 250 200 295 150 250 250 200 200 200 200 300 200 300 200 220 230 110 140 200 150 150 200 ) ## Sum-n Gaus
#nBins=( 200 200 150 300 250 190 200 200 200 200 150 50 180 200 150 200 200 200 200 150 200 180 155 200 80 250 150 320 120 230 200 200 200 300 250 200 295 150 250 250 200 200 200 200 300 200 300 200 220 230 110 140 200 150 150 200 ) ## Sum-n Gaus v1

#nBins=( 180 200 150 220 250 190 200 200 200 200 150 150 200 150 200 200 300 200 200 150 160 180 155 200 400 200 150 320 180 250 250 200 200 300 250 200 200 150 250 200 350 200 200 250 150 200 130 150 220 230 110 200 200 200 150 180 ) ## Sum-n Gaus
nBins=( 180 200 50 220 220 190 200 200 200 200 150 500 200 150 200 200 300 200 200 150 160 180 50 200 400 200 150 320 180 250 250 200 200 300 250 200 200 150 250 200 350 200 200 250 150 200 130 150 220 230 110 200 200 200 200 200 ) ## Sum-n Gaus

#nBins=150
#nBins=( 200 120 120 120 130 120 300 400 120 120 120 120 ) ## DCB+Gaus
nMass=${#massList[@]}

dir_sig="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal"
source /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/tools


# python3 ./fTest.py --mass_ALP 1 --year 16 --mass 125 --channel ele --ext 16_ele --inputWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig/ele
# exit

# loop structure: channel(3) -> ALP_mass(14) -> year(4) -> H mass(3)
# file structure: channel(3) -> 1_fTest_16_lep_Hm125.pdf
# file structure: channel(3) -> 1_nGauss_16_lep_Hm125.pdf

# for ((iChannel=0; iChannel<$nChannel; iChannel++))
#   do
#   for ((iALPmass=0; iALPmass<$nALPmass; iALPmass++))
#     do
#     for ((iYear=0; iYear<$nYear; iYear++))
#       do
#         for ((iHmass=0; iHmass<$nHmass; iHmass++))
#         do
#         python3 $dir_sig/fTest.py --mass_ALP ${ALPmassList[$iALPmass]} --year ${YearsList[$iYear]} --mass ${HmassList[$iHmass]} --channel ${ChannelList[$iChannel]} --skipWV --nProcsToFTest -1 --doPlots --inputWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig/${ChannelList[$iChannel]}
#         done
#       done
#     done
#   done

# python3 $dir_sig/scripts/signalFit.py --mass_ALP 1 --year 16 --channel ele --inputWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig/ele

python3 $dir_sig/RunPlotter.py --mass_ALP 20 --years 16 --channel mu

for ((iBin=0; iBin<$nMass; iBin++))
    do

    # mkdir -p "$dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}"

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

      #if [[ ${years[ $jBin ]} != $1 ]];then
      #    continue
      #fi

      # mkdir -p "$dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}/${years[$jBin]}"

      # path_out_bkg="$dir_out_sig/fit_results_${lable}_${channel}/M${massList[$iBin]}/${years[$jBin]}"

      # python3 ./fTest.py --mass_ALP 1 --year 16 --mass 125 --channel lep --procs GG2H --cat cat0 --ext 16_lep --skipWV --nProcsToFTest -1 --doPlots --inputWSDir /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig/lep

    #   ./bin/signalFTest_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_${years[$jBin]}_workspace_${channel}.root -d $path_out_bkg/ALPmassInde_data_sig.dat -o $path_out_bkg/HZAmassInde_ftest -p data -f cat0 -m 125 --mhLowBlind 110 --mhHighBlind 140 # --verbose 1
    #   sed -i "s/data/ggh/" $path_out_bkg/ALPmassInde_data_sig.dat


      #echo "./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins ${nBins[$[iBin*nYear] + $jBin]} --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 --useDCBplusGaus 1 "
      # ./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace_${channel}.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace_${channel}.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace_${channel}.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins ${nBins[$[iBin*nYear] + $jBin]} --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 #--verbose 2
      #./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace_${channel}.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace_${channel}.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace_${channel}.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s empty.dat --nBins $nBins --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 #--verbose 2
      #./bin/SignalFit_ALP -i $path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm120_${years[$jBin]}_workspace_${channel}.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm125_${years[$jBin]}_workspace_${channel}.root,$path_in_sig/ALP_data_sig_Am${massList[$iBin]}_Hm130_${years[$jBin]}_workspace_${channel}.root -o $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -p $path_out_bkg/plots_ALP --procs ggh -f cat0 -d $path_out_bkg/ALPmassInde_data_sig.dat -s ./ALP_sys_${version}/shapeSys/fit_run2_UL_ScaleSmear_${channel}/ScaleSmear_plot/M${massList[$iBin]}/ScaleSmear_m${massList[$iBin]}_${years[$jBin]}.dat --nBins ${nBins[$[iBin*nYear] + $jBin]} --changeIntLumi ${Lumis[$jBin]} --useSSF 1 -L 110 --mhHigh 140 --massList 120,125,130 #--useDCBplusGaus 1 #--verbose 2
      # ./bin/makeParametricSignalModelPlots_ALP -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root  -o $path_out_bkg/SignalModel -p ggh -f cat0 -L 100 -H 135 -m 125 --binning 70 # ${nBins[$[iBin*nYear] + $jBin]} #--useDCBplusGaus 1

      ###python makeParametricModelDatacardFLASHgg_ALP.py -i $path_out_bkg/CMS-HGG_sigfit_data_ggh_cat0.root -o $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt -p ggh -c cat0 --photonCatScales empty.dat ---p --intLumi ${Lumis[$jBin]} -m 125
      ###sed -i "45,54d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
      ###sed -i "21,41d" $path_out_bkg/datacard_ALPmass${massList[$iBin]}.txt
      done
done
