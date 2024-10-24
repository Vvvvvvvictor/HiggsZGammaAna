# Add tools dir to PYTHONPATH
eval `scramv1 runtime -sh`

export PYTHONPATH=$PYTHONPATH:${CMSSW_BASE}/src/flashggFinalFit/tools
export PYTHONPATH=$PYTHONPATH:${CMSSW_BASE}/src/flashggFinalFit/Trees2WS
export PYTHONPATH="/eos/user/${USER::1}/$USER/env4finalfit/lib/python3.9/site-packages/:$PYTHONPATH"
