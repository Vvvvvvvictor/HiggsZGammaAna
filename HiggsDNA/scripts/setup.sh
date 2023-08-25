#!/bin/bash
# include ROOT
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc9-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_100/ROOT/v6.24.00/x86_64-centos7-gcc9-opt/bin/thisroot.sh
source "`pwd`/higgsdnaenv/bin/activate"

export PATH="`pwd`:${PATH}"
export PYTHONPATH="`pwd`:${PYTHONPATH}"
export THEANO_FLAGS="gcc.cxxflags='-march=core2'"

export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"

export PATH="`pwd`/higgs-dna:$PATH"
export PYTHONPATH="`pwd`/higgs-dna:$PYTHONPATH"
