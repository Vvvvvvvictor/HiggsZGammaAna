#!/bin/bash

# initdir=$1
region=$1
fold=$2

# cd $initdir
# source scripts/setup.sh
for region in two_jet;do # zero_jet one_jet two_jet VH_ttH VBF; do
# region=two_jet
# for fold in {0..3};do
echo python scripts/train_bdt.py -r $region --optuna
python scripts/train_bdt.py -r $region --optuna --n-calls 100
# echo python scripts/train_bdt.py -r $region --skopt --skopt-plot
# python scripts/train_bdt.py -r $region  --skopt  --skopt-plot --n-calls 12
# python scripts/train_bdt.py -r $region  --skopt  --skopt-plot --hyperparams_path "models/skopt" --n-calls 36
# python scripts/train_bdt.py -r $region  --skopt  --skopt-plot --hyperparams_path "models/skopt" --n-calls 64
# done
done
        
