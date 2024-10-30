# do these after config the CMSSW env. and run the code "cmsenv"
python3 -m venv /eos/user/${USER::1}/$USER/env4finalfit
source /eos/user/${USER::1}/$USER/env4finalfit/bin/activate
export PYTHONPATH="/eos/user/${USER::1}/$USER/env4finalfit/lib/python3.9/site-packages/:$PYTHONPATH"
pip3 install numpy pandas uproot --ignore-installed
pip3 install commonTools