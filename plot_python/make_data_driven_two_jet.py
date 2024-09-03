import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

# Load data
branches = ["H_mass", "weight", "bdt_score_t", "regions"]
data_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/Data_fake.root")["two_jet"].arrays(branches, library="pd")
zg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZGToLLG_fake.root")["two_jet"].arrays(branches, library="pd")
ewkzg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZG2JToG2L2J_fake.root")["two_jet"].arrays(branches, library="pd")
dy_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/DYJetsToLL.root")["test"].arrays(branches, library="pd")

boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt", "r").readlines()[2].split()[1:]]
boundaries = [0.0] + boundaries
print(boundaries)

################################################################
# Compare true photon samples plus data-driven photon samples with data in signal region
################################################################


bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]#, "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
fake_bkg_list = ["EWKZ2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets", "DYJetsToLL"]
def apply(data_cr, zg_cr, ewkzg_cr, bkg_list, fake_bkg_list, bounds):
    data_cr, zg_cr = data_cr[(data_cr["bdt_score_t"] > bounds[0]) & (data_cr["bdt_score_t"] < bounds[1])], zg_cr[(zg_cr["bdt_score_t"] > bounds[0]) & (zg_cr["bdt_score_t"] < bounds[1])]
    ewkzg_cr = ewkzg_cr[(ewkzg_cr["bdt_score_t"] > bounds[0]) & (ewkzg_cr["bdt_score_t"] < bounds[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=80, weights=data_cr["weight"], range=[100, 180])
    ewkzg_cr_hist, bin_edges = np.histogram(ewkzg_cr["H_mass"], bins=80, weights=ewkzg_cr["weight"], range=[100, 180])
    zg_cr_hist, bin_edges = np.histogram(zg_cr["H_mass"], bins=80, weights=zg_cr["weight"], range=[100, 180])

    cr_hist = data_cr_hist - ewkzg_cr_hist - zg_cr_hist

    fake_bkg_sum = 0
    for fake_bkg in fake_bkg_list:
        fake_bkg_sample = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{fake_bkg}.root")["test"].arrays(branches, library="pd")
        fake_bkg_sample = fake_bkg_sample[(fake_bkg_sample["bdt_score_t"] > bounds[0]) & (fake_bkg_sample["bdt_score_t"] < bounds[1])]
        fake_bkg_sum += np.sum(fake_bkg_sample["weight"])

    data_driven_sf = fake_bkg_sum / np.sum(cr_hist)
    
    bkg_samples = []
    for bkg in bkg_list:
        bkg_sample = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
        bkg_sample = bkg_sample.query("bdt_score_t > @bounds[0] and bdt_score_t < @bounds[1]")
        bkg_samples.append(bkg_sample)
        
    data_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/data.root")["test"].arrays(branches, library="pd")
    data_sr = data_sr[(data_sr["bdt_score_t"] > bounds[0]) & (data_sr["bdt_score_t"] < bounds[1])]
    
    num_data_side_band = np.sum(data_sr[(data_sr["H_mass"] < 120) | (data_sr["H_mass"] > 130)]["weight"])
    num_cr_side_band = np.sum(cr_hist[np.where((bin_edges[1:] < 120) | (bin_edges[:-1] > 130))]) * data_driven_sf
    num_bkg_side_band = np.sum([np.sum(bkg_sample[(bkg_sample["H_mass"] < 120) | (bkg_sample["H_mass"] > 130)]["weight"]) for bkg_sample in bkg_samples])
    sf = num_data_side_band / (num_cr_side_band + num_bkg_side_band)

    data_cr['weight_cor'] = data_cr['weight'] * sf * data_driven_sf
    data_cr['weight_sf'] = sf * data_driven_sf
    data_cr['tag'] = np.zeros(len(data_cr))
    
    zg_cr['weight_cor'] = zg_cr['weight'] * sf * data_driven_sf * -1
    zg_cr['weight_sf'] = sf * data_driven_sf * -1
    zg_cr['tag'] = np.ones(len(zg_cr))
    
    ewkzg_cr['weight_cor'] = ewkzg_cr['weight'] * sf * data_driven_sf * -1
    ewkzg_cr['weight_sf'] = sf * data_driven_sf * -1
    ewkzg_cr['tag'] = np.ones(len(ewkzg_cr)) * 2
    
    # apply sf to bkg_samples
    for i, bkg_sample in enumerate(bkg_samples):
        bkg_sample['weight_cor'] = bkg_sample['weight'] * sf
        bkg_sample['weight_sf'] = sf
        bkg_sample['tag'] = np.ones(len(bkg_sample)) * (3 + i)

    output = pd.concat([data_cr, zg_cr, ewkzg_cr]+bkg_samples)

    return output

# Apply
with uproot.recreate("/eos/home-j/jiehan/root/outputs/two_jet/data_driven_bkg_v3.root") as output_file:
    output_dataframe = pd.DataFrame()
    for i in range(0, len(boundaries)-1):
        partial_output = apply(data_cr, zg_cr, ewkzg_cr, bkg_list, fake_bkg_list, [boundaries[i], boundaries[i+1]])
        output_dataframe = pd.concat([output_dataframe, partial_output])
    output_file['test'] = output_dataframe
