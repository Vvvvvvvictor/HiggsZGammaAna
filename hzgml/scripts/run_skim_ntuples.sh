#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

################
# Signal samples
################

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/ggH/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/VBF/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/WminusH/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/WplusH/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZH/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/ttH/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/ggH_M125_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/ggH/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/VBFH_M125_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/VBF/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/WminusH_M125_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/WminusH/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/WplusH_M125_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/WplusH/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/ZH_M125_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZH/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/ttH_M125_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/ttH/2017.root

##############
# Data samples
##############

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/data/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/data/Data_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/data/2017.root

###################
# Prompt MC samples
###################

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZGToLLG/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZG2JToG2L2J/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/TGJets/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/TTGJets/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/ZGToLLG_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZGToLLG/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/ZG2JToG2L2J_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZG2JToG2L2J/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/TTGJets_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/TTGJets/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/TGJets_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/TGJets/2017.root

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/DYJetsToLL/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/dy/DYJetsToLL_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/DYJetsToLL/2017.root

# Use fake photon background estimation with data-driven

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_med/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_fake/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_true/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_med/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/apply_weight.py

# ######################
# Non prompt MC samples
# ######################

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/LLAJJ/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/TT/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/WW/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/WZ/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZZ/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/LLAJJ_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/LLAJJ/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/TT_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/TT/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/WW_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/WW/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/WZ_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/WZ/2017.root
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/ZZ_2017/merged_nominal.parquet -o /eos/home-j/jiehan/root/2017/skimmed_ntuples/ZZ/2017.root

echo "==============FINISHED==========="
