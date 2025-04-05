#!/bin/bash
set -e

# Set up the environment
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_102/ROOT/6.26.04/x86_64-centos9-gcc11-opt/bin/thisroot.sh

# Activate the virtual environment
source /eos/home-p/pelai/HZgamma/Machine_Learning/hzg-pku-hzgml/bin/activate

# Set PYTHONPATH and QT environment variables
export PYTHONPATH="/eos/home-p/pelai/HZgamma/Machine_Learning/hzg-pku-hzgml/lib/python3.9/site-packages/:$PYTHONPATH"
export QT_QPA_PLATFORM=offscreen

# Define the base path
BASE_PATH="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/hzgml"

# Argument to determine which configuration to run
CONFIG=$1

# Run the appropriate Python command based on the argument
if [ "$CONFIG" == "zero_to_one_jet" ]; then
    python "${BASE_PATH}/scripts/train_bdt.py" -r zero_to_one_jet --save \
        --hyperparams_path "${BASE_PATH}/models/optuna_zero_to_one_jet" \
        --inputFolder "/eos/home-p/pelai/HZgamma/Root_Dataset/run2/NanoV9/Mix_Sig_WO_Systematic"
elif [ "$CONFIG" == "two_jet" ]; then
    python "${BASE_PATH}/scripts/train_bdt.py" -r two_jet --save \
        --hyperparams_path "${BASE_PATH}/models/optuna_two_jet" \
        --inputFolder "/eos/home-p/pelai/HZgamma/Root_Dataset/run2/NanoV9/Mix_Sig_WO_Systematic"
else
    echo "Invalid configuration: $CONFIG"
    exit 1
fi