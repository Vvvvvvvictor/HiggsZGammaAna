#!/bin/bash
# import ROOT
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_102/ROOT/6.26.04/x86_64-centos9-gcc11-opt/bin/thisroot.sh
source /eos/user/${USER::1}/$USER/hzgmlenv/bin/activate

export PATH="`pwd`:${PATH}"
export PYTHONPATH="`pwd`:${PYTHONPATH}"
export THEANO_FLAGS="gcc.cxxflags='-march=core2'"

export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"

export PATH="`pwd`/hzgml:$PATH"
export PYTHONPATH="`pwd`/hzgml:$PYTHONPATH"

export PYTHONPATH="/eos/user/${USER::1}/$USER/hzgmlenv/lib/python3.9/site-packages/:$PYTHONPATH"
