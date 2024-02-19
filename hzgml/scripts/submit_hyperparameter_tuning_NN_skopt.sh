#!/bin/bash

initdir=$1
region=$2
fold=$3

cd $initdir
source scripts/setup.sh
echo python scripts/train_NN.py -r $region -f $fold --skopt
python scripts/train_NN.py -r $region -f $fold --skopt
        