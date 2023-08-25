# include ROOT
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc9-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_100/ROOT/v6.24.00/x86_64-centos7-gcc9-opt/bin/thisroot.sh
python3 -m venv higgsdnaenv
source higgsdnaenv/bin/activate
pip install --upgrade pip
# pip install -r requirements.txt

# setting python path
export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"
export PATH="`pwd`/higgs-dna:$PATH"
export PYTHONPATH="`pwd`/higgs-dna:$PYTHONPATH"
