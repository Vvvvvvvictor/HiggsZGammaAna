#!/bin/bash
/bin/hostname
gcc -v
pwd
export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/

# source /cvmfs/cms.cern.ch/cmsset_default.sh

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/publicfs/cms/user/laipeizhu/anaconda/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/publicfs/cms/user/laipeizhu/anaconda/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/publicfs/cms/user/laipeizhu/anaconda/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/publicfs/cms/user/laipeizhu/anaconda/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda init

conda activate python_Three

# python makeWorkspace_data.py -m $1

python makeWorkspace_sig.py -m $1

# python makeWorkspace_data.py -m 1

# python makeWorkspace_data.py -m $1 --interp
python makeWorkspace_sig.py -m $1 --ele
#python makeWorkspace_sig_fTest.py -m $1 --ele
python makeWorkspace_sig.py -m $1 --mu
#python makeWorkspace_sig_fTest.py -m $1 --mu

#python makeWorkspace_data.py -m 30
#python makeWorkspace_sig.py -m 30 --ele
#python makeWorkspace_sig.py -m 30 --mu