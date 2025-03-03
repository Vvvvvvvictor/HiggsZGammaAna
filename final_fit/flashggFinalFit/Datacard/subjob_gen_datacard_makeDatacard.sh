#!/bin/bash

cmsenv

#!/bin/bash
ulimit -s unlimited
set -e
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src
export SCRAM_ARCH=el9_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard
export PYTHONPATH=$PYTHONPATH:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/tools

ChannelList=( ele mu )
# ChannelList=( ele )
nChannel=${#ChannelList[@]}

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# ALPmassList=( 1 )
nALPmass=${#ALPmassList[@]}

YearsList=( 16 16APV 17 18 )
# YearsList=( 16 )
nYear=${#YearsList[@]}

path_makeDatacard="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard"
path_inputWSDir="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output"

path_datacard_log="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/logfiles"

# python3 $path_makeYields/makeYields.py --inputWSDirMap 16=$path_inputWSDir --mass_ALP 1 --year 16 --channel ele


for ((iChannel=0; iChannel<$nChannel; iChannel++))
  do
  for ((iALPmass=0; iALPmass<$nALPmass; iALPmass++))
    do
    for ((iYear=0; iYear<$nYear; iYear++))
      do
        hep_sub $path_makeDatacard/runjob_gen_datacard_makeDatacard.sh -g cms -mem 2000 -wt mid -o $path_datacard_log/${ChannelList[$iChannel]}/${ALPmassList[$iALPmass]}_makeDatacard_job_${YearsList[$iYear]}_${ChannelList[$iChannel]}.log -e $path_datacard_log/${ChannelList[$iChannel]}/${ALPmassList[$iALPmass]}_makeDatacard_job_${YearsList[$iYear]}_${ChannelList[$iChannel]}.err -argu ${ALPmassList[$iALPmass]} ${YearsList[$iYear]} ${ChannelList[$iChannel]}
      done
    done
  done