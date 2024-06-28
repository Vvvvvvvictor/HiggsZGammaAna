#!/bin/bash
echo "Starting job"
echo "copying proxy file to /tmp area"
echo "start running"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/eos/user/z/zewang/anaconda/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/eos/user/z/zewang/anaconda/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/eos/user/z/zewang/anaconda/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/eos/user/z/zewang/anaconda/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


cd /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/Plotter
conda activate hzgml

#python scripts/HZGamma_ZReFit.py --doReFit --splitJobs --nJobs 10 --iJob $1

bash scripts/run_skim_ntuples.sh $1

echo "running done"