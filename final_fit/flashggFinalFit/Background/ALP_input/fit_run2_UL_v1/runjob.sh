#!/bin/bash
/bin/hostname
gcc -v
pwd
export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
source /cvmfs/cms.cern.ch/cmsset_default.sh

#python makeWorkspace_data.py -m $1 --interp
#python makeWorkspace_sig.py -m $1 --ele
#python makeWorkspace_sig_fTest.py -m $1 --ele
#python makeWorkspace_sig.py -m $1 --mu
#python makeWorkspace_sig_fTest.py -m $1 --mu

python makeWorkspace_data.py -m 30
python makeWorkspace_sig.py -m 30 --ele
python makeWorkspace_sig.py -m 30 --mu