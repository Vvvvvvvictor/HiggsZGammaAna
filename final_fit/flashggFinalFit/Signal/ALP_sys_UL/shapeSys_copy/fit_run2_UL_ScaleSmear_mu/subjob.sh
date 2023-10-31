#!/bin/bash

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )

nMass=${#massList[@]}

for ((iBin=0; iBin<$nMass; iBin++))
    do
    hep_sub runjob.sh -g cms -o job${massList[$iBin]}.out -e job${massList[$iBin]}.err -argu ${massList[$iBin]}
done