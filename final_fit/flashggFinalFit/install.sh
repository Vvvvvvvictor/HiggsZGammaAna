# do these after config the CMSSW env. and run the code "cmsenv"
python3 -m venv /eos/user/${USER::1}/$USER/env4finalfit
source /eos/user/${USER::1}/$USER/env4finalfit/bin/activate
pip3 install numpy pandas uproot --ignore-installed
pip3 install commonTools