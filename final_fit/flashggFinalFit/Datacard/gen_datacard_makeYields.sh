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
# YearsList=( 16APV )
nYear=${#YearsList[@]}

path_makeYields="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard"
path_inputWSDir="/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output"

python3 $path_makeYields/makeYields.py --inputWSDirMap 17=$path_inputWSDir --mass_ALP 30 --year 17 --channel mu


# for ((iChannel=0; iChannel<$nChannel; iChannel++))
#   do
#   for ((iALPmass=0; iALPmass<$nALPmass; iALPmass++))
#     do
#     for ((iYear=0; iYear<$nYear; iYear++))
#       do
#         python3 $path_makeYields/makeYields.py --inputWSDirMap ${YearsList[$iYear]}=$path_inputWSDir --mass_ALP ${ALPmassList[$iALPmass]} --year ${YearsList[$iYear]} --channel ${ChannelList[$iChannel]}
#       done
#     done
#   done


