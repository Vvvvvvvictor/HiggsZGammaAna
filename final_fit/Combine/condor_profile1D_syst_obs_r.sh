#!/bin/sh
ulimit -s unlimited
set -e
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src
export SCRAM_ARCH=el9_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine

if [ $1 -eq 0 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 0 --lastPoint 0 -n _profile1D_syst_obs_r.POINTS.0.0
fi
if [ $1 -eq 1 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 1 --lastPoint 1 -n _profile1D_syst_obs_r.POINTS.1.1
fi
if [ $1 -eq 2 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 2 --lastPoint 2 -n _profile1D_syst_obs_r.POINTS.2.2
fi
if [ $1 -eq 3 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 3 --lastPoint 3 -n _profile1D_syst_obs_r.POINTS.3.3
fi
if [ $1 -eq 4 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 4 --lastPoint 4 -n _profile1D_syst_obs_r.POINTS.4.4
fi
if [ $1 -eq 5 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 5 --lastPoint 5 -n _profile1D_syst_obs_r.POINTS.5.5
fi
if [ $1 -eq 6 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 6 --lastPoint 6 -n _profile1D_syst_obs_r.POINTS.6.6
fi
if [ $1 -eq 7 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 7 --lastPoint 7 -n _profile1D_syst_obs_r.POINTS.7.7
fi
if [ $1 -eq 8 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 8 --lastPoint 8 -n _profile1D_syst_obs_r.POINTS.8.8
fi
if [ $1 -eq 9 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 9 --lastPoint 9 -n _profile1D_syst_obs_r.POINTS.9.9
fi
if [ $1 -eq 10 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 10 --lastPoint 10 -n _profile1D_syst_obs_r.POINTS.10.10
fi
if [ $1 -eq 11 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 11 --lastPoint 11 -n _profile1D_syst_obs_r.POINTS.11.11
fi
if [ $1 -eq 12 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 12 --lastPoint 12 -n _profile1D_syst_obs_r.POINTS.12.12
fi
if [ $1 -eq 13 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 13 --lastPoint 13 -n _profile1D_syst_obs_r.POINTS.13.13
fi
if [ $1 -eq 14 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 14 --lastPoint 14 -n _profile1D_syst_obs_r.POINTS.14.14
fi
if [ $1 -eq 15 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 15 --lastPoint 15 -n _profile1D_syst_obs_r.POINTS.15.15
fi
if [ $1 -eq 16 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 16 --lastPoint 16 -n _profile1D_syst_obs_r.POINTS.16.16
fi
if [ $1 -eq 17 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 17 --lastPoint 17 -n _profile1D_syst_obs_r.POINTS.17.17
fi
if [ $1 -eq 18 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 18 --lastPoint 18 -n _profile1D_syst_obs_r.POINTS.18.18
fi
if [ $1 -eq 19 ]; then
  combine --floatOtherPOIs 1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.38 -d /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_datacard_rootfile_ele/1_Datacard_16_ele_mu_inclusive.root --setParameterRanges r=0,2 --points 20 --firstPoint 19 --lastPoint 19 -n _profile1D_syst_obs_r.POINTS.19.19
fi

