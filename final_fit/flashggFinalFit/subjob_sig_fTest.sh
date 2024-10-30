#!/bin/bash

cmsenv

lable='run2'
version='UL'
#channel=( ele mu )
# channel='ele'
#Lumi_run2='137.64'
Lumi_run2='138'

years=( 16 16APV 17 18 )
#years=( 16 16APV )
nYear=${#years[@]}
Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

ChannelList=( ele mu )
# ChannelList=( mu )
nChannel=${#ChannelList[@]}

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# ALPmassList=( 1 )
nALPmass=${#ALPmassList[@]}

YearsList=( 16 16APV 17 18 )
# YearsList=( 16APV )
nYear=${#YearsList[@]}

HmassList=( 125 )
nHmass=${#HmassList[@]}

#massList=( $1 )

# massList=( 3 7 )

# nbins 5GeV 3.7 1GeV 3,5

#nBins=( 120 120 120 120 120 120 120 130 120 300 400 120 120 120 120 ) ## Sum-n Gaus
#nBins=( 80 350 100 300 250 190 200 200 200 200 150 50 180 200 150 200 200 200 200 150 200 180 155 200 80 250 150 320 120 230 200 200 200 300 250 200 295 150 250 250 200 200 200 200 300 200 300 200 220 230 110 140 200 150 150 200 ) ## Sum-n Gaus
#nBins=( 200 200 150 300 250 190 200 200 200 200 150 50 180 200 150 200 200 200 200 150 200 180 155 200 80 250 150 320 120 230 200 200 200 300 250 200 295 150 250 250 200 200 200 200 300 200 300 200 220 230 110 140 200 150 150 200 ) ## Sum-n Gaus v1

#nBins=( 180 200 150 220 250 190 200 200 200 200 150 150 200 150 200 200 300 200 200 150 160 180 155 200 400 200 150 320 180 250 250 200 200 300 250 200 200 150 250 200 350 200 200 250 150 200 130 150 220 230 110 200 200 200 150 180 ) ## Sum-n Gaus
nBins=( 180 200 50 220 220 190 200 200 200 200 150 500 200 150 200 200 300 200 200 150 160 180 50 200 400 200 150 320 180 250 250 200 200 300 250 200 200 150 250 200 350 200 200 250 150 200 130 150 220 230 110 200 200 200 200 200 ) ## Sum-n Gaus

#nBins=150
#nBins=( 200 120 120 120 130 120 300 400 120 120 120 120 ) ## DCB+Gaus
nMass=${#massList[@]}

###### background fit ######

cd ./Signal/

dir_out_bkg="./ALP_BkgModel_param_${version}"

if [ ! -d $dir_out_bkg ];then
  # mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

# path_in_bkg="./ALP_input/fit_${lable}_${version}"
path_in_bkg="../MVAcut/${lable}_${version}/output/data"

path_out_bkg="$dir_out_bkg/fit_results_${lable}"
# mkdir -p $path_out_bkg

# path_in_sig="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig"

path_run="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit"

path_sig_log="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/log_files"

# loop structure: channel(3) -> ALP_mass(14) -> year(4) -> H mass(3)
# file structure: channel(3) -> 1_fTest_16_lep_Hm125.pdf

for ((iChannel=0; iChannel<$nChannel; iChannel++))
  do
  for ((iALPmass=0; iALPmass<$nALPmass; iALPmass++))
    do
    for ((iYear=0; iYear<$nYear; iYear++))
      do
        for ((iHmass=0; iHmass<$nHmass; iHmass++))
        do
          hep_sub $path_run/runjob_sig_fTest.sh -g cms -mem 2000 -wt mid -o $path_sig_log/${ChannelList[$iChannel]}/${ALPmassList[$iALPmass]}_fTest_job_${YearsList[$iYear]}_${ChannelList[$iChannel]}_Hm${HmassList[$iHmass]}.log -e $path_sig_log/${ChannelList[$iChannel]}/${ALPmassList[$iALPmass]}_fTest_job_${YearsList[$iYear]}_${ChannelList[$iChannel]}_Hm${HmassList[$iHmass]}.err -argu ${ALPmassList[$iALPmass]} ${YearsList[$iYear]} ${HmassList[$iHmass]} ${ChannelList[$iChannel]} "/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig/${ChannelList[$iChannel]}"
        done
      done
    done
  done
