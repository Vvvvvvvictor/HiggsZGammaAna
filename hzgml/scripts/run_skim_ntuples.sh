#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

################
# Signal samples
################

mkdir -p /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ggH/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/VBF/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WminusH/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WplusH/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZH/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ttH/
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/ggH_M125_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ggH/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/VBFH_M125_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/VBF/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/WminusH_M125_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WminusH/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/WplusH_M125_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WplusH/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/ZH_M125_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZH/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/signal/ttH_M125_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ttH/2017.root

##############
# Data samples
##############

mkdir -p /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/data/
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/data/Data_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/data/2017.root

###################
# Prompt MC samples
###################

mkdir -p /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZGToLLG/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZG2JToG2L2J/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TGJets/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TTGJets/
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/ZGToLLG_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZGToLLG/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/ZG2JToG2L2J_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZG2JToG2L2J/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/TTGJets_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TTGJets/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/TGJets_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TGJets/2017.root

mkdir -p /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/DYJetsToLL/
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/dy/DYJetsToLL_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/DYJetsToLL/2017.root

Use fake photon background estimation with data-driven

mkdir -p /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/data_med/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/data_fake/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/mc_true/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/mc_med/
python scripts/apply_weight.py

# ######################
# Non prompt MC samples
# ######################

mkdir -p /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/LLAJJ/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TT/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WW/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WZ/ /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZZ/
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/LLAJJ_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/LLAJJ/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/TT_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TT/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/WW_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WW/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/WZ_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/WZ/2017.root
python scripts/skim_ntuples.py -i /eos/home-j/jiehan/parquet/2017/mva_based/background/ZZ_2017/merged_nominal.parquet -o /afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZZ/2017.root

echo "==============FINISHED==========="
