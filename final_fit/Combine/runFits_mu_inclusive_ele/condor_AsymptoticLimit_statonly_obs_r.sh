#!/bin/sh
ulimit -s unlimited
set -e
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src
export SCRAM_ARCH=el9_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/runFits_mu_inclusive_ele

if [ $1 -eq 0 ]; then
  combine --redefineSignalPOI r --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M AsymptoticLimits -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/2_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 -n _AsymptoticLimit_statonly_obs_r
fi

