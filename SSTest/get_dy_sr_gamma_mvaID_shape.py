import utils as ut
import pandas as pd
import os
import json

# Load the data
path = "/eos/user/j/jiehan/data_for_norm_float_v1/DYJets_01J"
selection = ["regions == 0", "H_mass > 100", "H_mass < 180", "H_mass + Z_mass > 185", "gamma_eta < 10", "gamma_mvaID > -1"]
bins = 100
tree = "zero_to_one_jet"

data_B, data_E = pd.DataFrame(), pd.DataFrame()
if os.path.isdir(path):
    for file in os.listdir(path):
        if file.endswith(".root"):
            data = ut.ReadFile(f"{path}/{file}", tree=tree, selections=selection)
            data_B = pd.concat([data_B, data.query("(gamma_eta<1.5) & (gamma_eta>-1.5)")])
            data_E = pd.concat([data_E, data.query("(gamma_eta>1.5) | (gamma_eta<-1.5)")])

data_B.to_pickle("/eos/user/j/jiehan/data_for_norm_float_v1/gamma_mvaID_01J_B.pkl")
data_E.to_pickle("/eos/user/j/jiehan/data_for_norm_float_v1/gamma_mvaID_01J_E.pkl")
    