# include ROOT

# source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
# source /cvmfs/sft.cern.ch/lcg/releases/LCG_104/ROOT/6.28.04/x86_64-centos7-gcc11-opt/bin/thisroot.sh

source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_102/ROOT/6.26.04/x86_64-centos9-gcc11-opt/bin/thisroot.sh

python3 -m venv /eos/user/${USER::1}/$USER/hzgmlenv/
source /eos/user/${USER::1}/$USER/hzgmlenv/bin/activate

export PYTHONPATH="/eos/user/${USER::1}/$USER/hzgmlenv/lib/python3.9/site-packages/:$PYTHONPATH"

# pip install --upgrade pip
# pip install -r requirements.txt
pip install scikit-optimize==0.9
pip install xgboost==2.0.3
pip install optuna==3.6.1
pip install optuna_dashboard==0.15.1
pip install plotly==5.22.0 --ignore-installed 
pip install pyarrow==14.0.0 
pip install fastparquet==0.8.2

# pip install -U --ignore-installed scikit-learn tensorboardx servicex numba httpstan coffea

# # setting python path
# export PATH="`pwd`/scripts:${PATH}"
# export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"
# export PATH="`pwd`/hzgml:$PATH"
# export PYTHONPATH="`pwd`/hzgml:$PYTHONPATH"
