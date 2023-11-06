# include ROOT
source /eos/user/j/jiehan/hzgmlenv_test/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# setting python path
export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"
export PATH="`pwd`/hzgml:$PATH"
export PYTHONPATH="`pwd`/hzgml:$PYTHONPATH"
