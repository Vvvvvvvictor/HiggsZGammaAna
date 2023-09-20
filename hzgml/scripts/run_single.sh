#!/bin/bash
#########################################################################
#                                                                       #
#    The example wrapper for the training and categorization task.      #
#                                                                       #
#    Please contact ~jay.chan@cern.ch~ in case there is any issue.      #
#                                                                       #
#########################################################################

echo $#
if [ $# == 1 ]; then
    T=$1; A=-1; S=-1;
elif [ $# == 2 ]; then
    A=$2; S=-1; T=$1
else 
    exit 0
fi
echo "==================================================="
echo "region: $T, num: $2"
echo "==================================================="
echo "Shielded parameter is: $S . Added variables is: $A ."
############################
#  Training the BDT models
############################
# python scripts/train_bdt.py -r $T --save -s $S -a $A
# python scripts/train_bdt.py -r VBF --save -s $S -a $A

###########################################
#  Applying the BDT models to all samples
###########################################
# python scripts/apply_bdt.py -r $T -s $S -a $A

###########################################################
#  Optimizing the BDT boundaries for zero-jet and two-jet
###########################################################
python scripts/categorization_1D.py -r $T -b 3 -s $S -a $A

##############################################
#  Optimizing the BDT boundaries for two-jet
##############################################
# python scripts/categorization_2D.py -r two_jet -b 4 -v 4 -s $S

# root -l -q 'scripts/draw_bdtDis.cpp("'$T'", 1, 1000000)'
