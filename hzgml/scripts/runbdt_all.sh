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

# python scripts/train_bdt.py -r two_jet --optuna --n-calls 20 --continue-optuna 0
# for fold in {0..3};do
# python scripts/train_bdt.py -r two_jet --optuna --n-calls 80 --fold $fold --continue-optuna 1 
# done
python scripts/train_bdt.py -r two_jet --save  --hyperparams_path "models/optuna_new_better"
python scripts/apply_bdt.py -r two_jet
python scripts/categorization_1D.py -r two_jet -b 4 --minN 10 --floatB
# python scripts/Generate_fake_photon_template.py -r two_jet
# python ../plot_python/make_data_driven_two_jet.py 
# python ../plot_python/compare_sig_bkg_bdt_sosb.py
python ../SSTest/Generate_fake_photon_template_2J.py
# rm ../plot_python/pic/2J/*
# python ../plot_python/compare_data_driven_hmass_two_jet.py


# # python scripts/train_bdt.py -r two_jet --skopt-plot --params '{"silent": 1, "eval_metric": ["logloss", "auc"], "grow_policy": "lossguide", "nthread": 4, "objective": "binary:logistic", "tree_method": "hist", "booster": "gbtree", "alpha": 0.5144783323380544, "colsample_bytree": 0.9588358073169332, "gamma": 3.3537213020169725, "max_delta_step": 19.5023193765768, "min_child_weight": 77.0, "subsample": 0.9644739775053346, "eta": 0.01681114970710191, "max_bin": 330.0, "max_depth": 6.0}'

# # python scripts/reweight.py > log
# # source scripts/submit_hyperparameter_tuning_bdt_skopt.sh

# ############################
# #  Training the BDT models
# ############################
# # python scripts/train_bdt.py -r zero_jet --save -s $S -a $A
# # python scripts/train_bdt.py -r one_jet --save -s $S -a $A
# python scripts/train_bdt.py -r two_jet --optuna --n-calls 10 --continue-optuna 0
# for fold in {0..3};do
# python scripts/train_bdt.py -r two_jet --optuna --n-calls 40 --fold $fold --continue-optuna 1 
# done
# # python scripts/train_bdt.py -r two_jet --optuna --n-calls 50 --fold 3 --continue-optuna 0
# python scripts/train_bdt.py -r two_jet --save --hyperparams_path "models/optuna"
# # python scripts/train_bdt.py -r two_jet --save -s $S -a $A --hyperparams_path "models/skopt_for_norm"
# # python scripts/train_bdt.py -r VBF --save -s $S -a $A --hyperparams_path "models/skopt"
# # python scripts/train_bdt.py -r VH_ttH --save -s $S -a $A

# ###########################################
# #  Applying the BDT models to all samples
# ###########################################
# # python scripts/apply_bdt.py -r zero_jet -s $S -a $A
# # python scripts/apply_bdt.py -r one_jet -s $S -a $A
# # python scripts/apply_bdt.py -r zero_to_one_jet
# python scripts/apply_bdt.py -r two_jet
# # python scripts/apply_bdt.py -r VBF -s $S -a $A
# # python scripts/apply_bdt.py -r VH_ttH -s $S -a $A

# # python scripts/Generate_fake_photon_template.py -r two_jet
# # # python scripts/Generate_fake_photon_template_zebing.py -r zero_to_one_jet

# ###########################################################
# #  Optimizing the BDT boundaries for zero-jet and two-jet
# ###########################################################     
# # python scripts/categorization_1D.py -r zero_jet -b 4 -s $S -a $A --minN 300
# # python scripts/categorization_1D.py -r one_jet -b 4 -s $S -a $A --minN 300
# python scripts/categorization_1D.py -r two_jet -b 4 --minN 10 --floatB
# # python scripts/categorization_1D.py -r two_jet -b 4 --minN 200
# # python scripts/categorization_1D.py -r VH_ttH -b 2 -s $S -a $A --minN 20

# ##############################################
# #  Optimizing the BDT boundaries for two-jet
# ##############################################
# # python scripts/categorization_2D.py -r two_jet -b 4 -v 3 --minN 10
# # python scripts/categorization_2D_vbf_2j.py -r two_jet -b 4 -v 4 --minN 10

# python ../SSTest/Generate_template.py
# rm ../plot_python/pic/2J/*
# python ../plot_python/compare_data_driven_hmass_two_jet.py
# # python ../plot_python/compare_data_driven_hmass_zero_to_one_jet.py
# # python ../plot_python/compare_bkg_hmass.py
# # python ../plot_python/compare_sig_bkg_bdt_sosb.py 

# # zero_jet one_jet two_jet VH_ttH

# # for i in zero_jet one_jet two_jet VH_ttH;
# # doo
# # root -l -q 'scripts/draw_bdtDis.cpp("'$i'", 0, 1)';
# # root -l -q 'scripts/draw_HMassDis.cpp("'$i'", 1)';
# # root -l -q 'scripts/draw_HMassDis_sum.cpp("'$i'", 1)';
# # done
