#!/bin/bash

cmsenv

ChannelList=( ele mu )
# ChannelList=( ele )
nChannel=${#ChannelList[@]}

ALPmassList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# ALPmassList=( 1 )
nALPmass=${#ALPmassList[@]}

YearsList=( 16 16APV 17 18 )
# YearsList=( 16 )
nYear=${#YearsList[@]}

# Expected (--doObserved False)
python3 RunFits.py --inputJson inputs_HToZa.json --mode mu_inclusive --doObserved False --mass_ALP 2 --year 16 --channel ele

# Data (--doObserved True)



# for ((iChannel=0; iChannel<$nChannel; iChannel++))
#   do
#   for ((iALPmass=0; iALPmass<$nALPmass; iALPmass++))
#     do
#     for ((iYear=0; iYear<$nYear; iYear++))
#       do
#       python3 RunText2Workspace.py --mode mu_inclusive --mass_ALP ${ALPmassList[$iALPmass]} --year ${YearsList[$iYear]} --channel ${ChannelList[$iChannel]}
#       done
#     done
#   done

