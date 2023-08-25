#!/bin/bash
#########################################################################
#                                                                       #
#    The example wrapper for the training and categorization task.      #
#                                                                       #
#    Please contact ~jay.chan@cern.ch~ in case there is any issue.      #
#                                                                       #
#########################################################################

echo $#
if [ $# == 1 ]
then
    A=$1; S=-1;
else 
    S=-1; A=-1;
fi
echo "==================================================="
echo "num: $2"
echo "==================================================="
echo "Shielded parameter is: $S . Added variables is: $A ."
############################
#  Training the BDT models
############################
python scripts/train_NN.py -r zero_jet --save -s $S -a $A
python scripts/train_NN.py -r one_jet --save -s $S -a $A
python scripts/train_NN.py -r two_jet --save -s $S -a $A
python scripts/train_NN.py -r VH_ttH --save -s $S -a $A

###########################################
#  Applying the BDT models to all samples
###########################################
python scripts/apply_NN.py -r zero_jet -s $S -a $A
python scripts/apply_NN.py -r one_jet -s $S -a $A
python scripts/apply_NN.py -r two_jet -s $S -a $A
python scripts/apply_NN.py -r VH_ttH -s $S -a $A

###########################################################
#  Optimizing the BDT boundaries for zero-jet and two-jet
###########################################################
python scripts/categorization_1D.py -r zero_jet -b 7 -s $S -a $A
python scripts/categorization_1D.py -r one_jet -b 7 -s $S -a $A
python scripts/categorization_1D.py -r two_jet -b 7 -s $S -a $A
python scripts/categorization_1D.py -r VH_ttH -b 7 -s $S -a $A

##############################################
#  Optimizing the BDT boundaries for two-jet
##############################################
# python scripts/categorization_2D.py -r two_jet -b 4 -v 4 -s $S

# root -l -q 'scripts/draw_bdtDis.cpp("zero_jet", 1, 1000000)'
# root -l -q 'scripts/draw_bdtDis.cpp("one_jet", 1, 1000000)'
# root -l -q 'scripts/draw_bdtDis.cpp("two_jet", 1, 1000000)'
# root -l -q 'scripts/draw_bdtDis.cpp("VH_ttH", 1, 1000000)'
