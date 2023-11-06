#!/bin/bash
echo "Starting job"
echo "copying proxy file to /tmp area"
echo "start running"
cd /afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/
cmsenv
#cd ./flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_runII/M_MASS
#combine -M GoodnessOfFit datacard_ALPmass_MASS.txt --algo=saturated -m 125 --setParameters MH=125 -n _saturated
#combine -M GoodnessOfFit datacard_ALPmass_MASS.txt --algo=saturated -t 500 -s 12345 -m 125 --setParameters MH=125 -n _saturated

#combine -M GoodnessOfFit datacard_ALPmass_MASS.txt --algo=KS -m 125 --setParameters MH=125 -n _KS
#combine -M GoodnessOfFit datacard_ALPmass_MASS.txt --algo=KS -t 500 -s 12345 -m 125 --setParameters MH=125 -n _KS

cd ./flashggFinalFit
./runII_combine_UL.sh $1

echo "running done"
