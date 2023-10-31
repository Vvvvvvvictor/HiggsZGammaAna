#!/bin/bash

seconds=$((600))

for i in {1..50000}
do

    ./runII_combine_UL.sh
    echo "Delete condor output file, sleep 10 mins"
    sleep ${seconds}

done