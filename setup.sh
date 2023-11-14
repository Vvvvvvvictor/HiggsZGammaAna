#TODO: some problems here. Seems break the environment after "cmsenv"
cd /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA
conda deactivate
conda activate higgs-zg-ana

## python 3.8
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc9-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_100/ROOT/v6.24.00/x86_64-centos7-gcc9-opt/bin/thisroot.sh

## python 3.9
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc12-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_104/ROOT/6.28.04/x86_64-centos7-gcc12-opt/bin/thisroot.sh

cd /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/CMSSW_11_3_4/src/
cmsenv
cd /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/