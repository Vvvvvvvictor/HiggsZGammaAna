#!/bin/bash

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
# massList=( 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 )

nMass=${#massList[@]}

for ((iBin=0; iBin<$nMass; iBin++))
    do
    hep_sub runjob.sh -g cms -mem 4000 -wt mid -o ./log_file/job${massList[$iBin]}.log -e ./log_file/job${massList[$iBin]}.err -argu ${massList[$iBin]}
    done