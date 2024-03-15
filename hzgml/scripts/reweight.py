import uproot
import pandas as pd
import numpy as np
import os

def combine(path, my_path, process):
    tree_name = "two_jet"
    file_path = path +  process
    if os.path.isdir(file_path):
        if process == "Data": process = "data"
        output_path = my_path + process + "/"
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        output_file = uproot.recreate(output_path+"two_jet_run2.root")
        data = pd.DataFrame()

        for root_file in os.listdir(file_path):
            if "run2" in root_file:
                continue
            print("processing {}...".format(file_path+"/"+root_file))
            input_file = uproot.open(file_path+"/"+root_file)

            data_o = input_file[tree_name].arrays(library='pd').copy()
            print(data_o[data_o["n_jets"]>=2]["weight"].sum())
            data = pd.concat([data, data_o], ignore_index=True, sort=False)

            input_file.close()
        output_file[tree_name] = pd.concat([data], ignore_index=True, sort=False)
        output_file.close()

def reweight(path, my_path, process):
    file_path = path +  process
    if os.path.isdir(file_path):
        output_path = my_path + process + "/"
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        for root_file in os.listdir(file_path):
            if "two_jet" not in root_file:
                continue
            output_file = uproot.recreate(output_path+root_file)
            print("processing {}...".format(file_path+"/"+root_file))
            input_file = uproot.open(file_path+"/"+root_file)

            for tree_name in input_file.keys():
                tree_name = tree_name.split(";")[0]
                data = input_file[tree_name].arrays(library='pd').copy()
                print(data[data["n_jets"]>=2]["weight"].sum())
                if "reweight" in data.keys():
                    print(data["weight"].sum())
                    data["weight"] = data["weight"] * data["reweight"]
                print(data["weight"].sum())
                output_file[tree_name] = pd.concat([data], ignore_index=True, sort=False)
                
            input_file.close()
            output_file.close()


bkg_list = ["data", "ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]

my_path = "/eos/user/j/jiehan/root/skimmed_ntuples_two_jet/"

path = "/eos/user/z/zewang/HZGamma_data/run2UL/skimmed_ntuples/"
process_list = os.listdir(path)
for process in process_list:
    reweight(path, my_path, process)
    # if process in bkg_list:
    #     reweight(path, my_path, process)

# path = "/eos/user/j/jiehan/root/skimmed_ntuples/"
# process_list = os.listdir(path)
# for process in process_list:
#     combine(path, my_path, process)
#     # if process not in bkg_list:
#     #     combine(path, my_path, process)