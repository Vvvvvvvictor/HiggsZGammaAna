import pandas as pd
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace


VAR = "gamma_pt"
XLABLE = r"photon pt"
BINS = 40
RMIN = 15
RMAX = 55
WEIGHT = "weight_central"



data_list = [
    "/eos/user/j/jiehan/parquet/nanov9/data/Data_2016/merged_nominal.parquet"
]
bkg_mc_list = [
    "/eos/user/j/jiehan/parquet/nanov9/bkgmc/DYJetsToLL_2016/merged_nominal.parquet",
    "/eos/user/j/jiehan/parquet/nanov9/bkgmc/ZGToLLG_2016/merged_nominal.parquet"
]

bkg_mc_hist_list, data_hist_list = [], []

def convert_parquet_to_hist(file_list, hist_list):
    for data in file_list:
        samples = pd.read_parquet(data)
        if "weight_central_no_lumi" in samples.keys():
            print(data, ":", samples["weight_central"][1]/samples["weight_central_no_lumi"][1], sep="\n")
        # if "DY" in data:
            # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
            # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
            # samples = samples[samples["n_iso_photons"]==0]
            # samples[WEIGHT]*=3
        # if "ZG" in data:
            # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
            # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
            # samples = samples[samples["n_iso_photons"]>0]
            # samples[WEIGHT]*=3
        hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
        hist_list.append(hist)

convert_parquet_to_hist(data_list, data_hist_list)
convert_parquet_to_hist(bkg_mc_list, bkg_mc_hist_list)

hep.style.use("CMS")
fig, ax = plt.subplots(figsize=(10, 10), dpi=400)

bins = np.linspace(RMIN, RMAX, BINS+1)
hep.histplot(
    data_hist_list, 
    bins=bins, 
    stack=True, 
    histtype='errorbar', 
    color='r',
    ax=ax,
    label=["Data"]
)
hep.histplot(
    bkg_mc_hist_list, 
    bins=bins, 
    histtype='fill', 
    color=['b','g'],
    ax=ax,
    label=['DYJets','SM ZG']
)

ax.set_xlabel("{} (GeV)".format(XLABLE), fontsize=25)
ax.set_ylabel("Events / {:.0f} GeV".format((RMAX-RMIN)/BINS), fontsize=25)
ax.set_xlim(RMIN, RMAX)
ax.legend()

plt.savefig("pic/a.png")
