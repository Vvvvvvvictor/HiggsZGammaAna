import uproot
import os
from pdb import set_trace

def process_file(input_file, output_file):
    with uproot.open(input_file) as file:
        df = file['two_jet'].arrays(library="pd")
        df = df.reset_index(drop=True)
        df["event"] = df.index
        with uproot.recreate(output_file) as f:
            f["two_jet"] = df


input_path = "/eos/user/j/jiehan/root/skimmed_ntuples_two_jet/"
output_path = "/eos/user/j/jiehan/root/skimmed_ntuples_two_jet_reindex/"

for dir in os.listdir(input_path):
    if os.path.isdir(input_path + dir) == False:
        continue
    if not os.path.exists(output_path + dir):
        os.makedirs(output_path + dir)
    for file in os.listdir(input_path + dir):
        input_file = input_path + dir + "/" + file
        output_file = output_path + dir + "/" + file
        process_file(input_file, output_file)
        print("Reindexed", input_file, "to", output_file)