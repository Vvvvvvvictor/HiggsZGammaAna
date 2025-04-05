#!/bin/bash
set -e

# Set up the environment (copied from install.sh)
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_102/ROOT/6.26.04/x86_64-centos9-gcc11-opt/bin/thisroot.sh

# Activate the virtual environment (replace <your_username> with your actual username)
source /eos/home-p/pelai/HZgamma/Machine_Learning/hzg-pku-hzgml/bin/activate

# Set PYTHONPATH
export PYTHONPATH="/eos/home-p/pelai/HZgamma/Machine_Learning/hzg-pku-hzgml/lib/python3.9/site-packages/:$PYTHONPATH"
export QT_QPA_PLATFORM=offscreen

# Get fold and region from command-line arguments passed by HTCondor
FOLD=$1
REGION=$2

echo "======== Starting script ========"
echo "FOLD: $FOLD"
echo "REGION: $REGION"
echo "Current working directory: $(pwd)"
echo "Python version: $(python --version)"
echo "Which python: $(which python)"
echo "VIRTUAL_ENV: $VIRTUAL_ENV"

# Run the appropriate Python command based on the region
if [ "$REGION" == "zero_to_one_jet" ]; then
    echo "Running: python scripts/train_bdt.py ..."
    python /afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/hzgml/scripts/train_bdt.py -r zero_to_one_jet --optuna --n-calls 40 --fold $FOLD --continue-optuna 0 --optuna_metric "eval_auc" --inputFolder "/eos/home-p/pelai/HZgamma/Root_Dataset/run2/NanoV9/Mix_Sig_WO_Systematic"
elif [ "$REGION" == "two_jet" ]; then
    echo "Running: python scripts/train_bdt.py ..."
    python /afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/hzgml/scripts/train_bdt.py -r two_jet --optuna --n-calls 20 --fold $FOLD --continue-optuna 0 --optuna_metric "eval_auc" --inputFolder "/eos/home-p/pelai/HZgamma/Root_Dataset/run2/NanoV9/Mix_Sig_WO_Systematic"
fi