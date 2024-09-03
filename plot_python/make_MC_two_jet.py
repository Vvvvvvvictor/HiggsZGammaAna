import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

# Load data
branches = ["H_mass", "weight", "bdt_score_t", "regions"]

boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt", "r").readlines()[2].split()[1:]]
print(boundaries)

################################################################
# Compare true photon samples plus data-driven photon samples with data in signal region
################################################################


bkg_list = ["DYJetsToLL", "ZGToLLG", "ZG2JToG2L2J", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
def apply(bkg_list, bounds):
    bkg_samples = []
    for bkg in bkg_list:
        bkg_sample = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
        bkg_sample = bkg_sample.query("bdt_score_t > @bounds[0] and bdt_score_t < @bounds[1]")
        bkg_samples.append(bkg_sample)
        
    data_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/data.root")["test"].arrays(branches, library="pd")
    data_sr = data_sr[(data_sr["bdt_score_t"] > bounds[0]) & (data_sr["bdt_score_t"] < bounds[1])]
    
    num_data_side_band = np.sum(data_sr[(data_sr["H_mass"] < 120) | (data_sr["H_mass"] > 130)]["weight"])
    num_bkg_side_band = np.sum([np.sum(bkg_sample[(bkg_sample["H_mass"] < 120) | (bkg_sample["H_mass"] > 130)]["weight"]) for bkg_sample in bkg_samples])
    sf = num_data_side_band / num_bkg_side_band

    # apply sf to bkg_samples
    for i, bkg_sample in enumerate(bkg_samples):
        bkg_sample['weight_cor'] = bkg_sample['weight'] * sf
        bkg_sample['weight_sf'] = sf
        bkg_sample['tag'] = np.ones(len(bkg_sample)) * i
        
    print(f"bin {bounds[0]}-{bounds[1]}: sf: {sf}")

    output = pd.concat(bkg_samples)

    return output

# Apply
with uproot.recreate("/eos/home-j/jiehan/root/outputs/two_jet/MC_bkg_v3.root") as output_file:
    output_dataframe = pd.DataFrame()
    for i in range(0, len(boundaries)-1):
        partial_output = apply(bkg_list, [boundaries[i], boundaries[i+1]])
        output_dataframe = pd.concat([output_dataframe, partial_output])
    output_file['two_jet'] = output_dataframe
