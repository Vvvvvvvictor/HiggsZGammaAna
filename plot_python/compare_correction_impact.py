import pandas as pd
import numpy as np
import os

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

def get_nominal_data(year, sample):
    data = pd.read_parquet(os.path.join(path, f"{sample}_{year}/merged_nominal.parquet"))
    return data

def get_correction_data(year, sample, correction_name, correction_type):
    data = pd.read_parquet(os.path.join(path, f"{sample}_{year}/merged_{correction_name}_{correction_type}.parquet"))
    return data

path = "/eos/home-j/jiehan/parquet/nanov9/data_for_norm_v2/"
photon_correction_name = ["fnuf", "material", "scale", "smear"]
jet_correction_name = ["JER", "JES"]
muon_correction_name = ["Muon_pt"]
met_correction_name = ["MET_Unclustered", "MET_JES"]
years = ["2016preVFP", "2016postVFP", "2017", "2018"]

# section for signal comparison
signal_samples = ["ggH_M125", "VBF_M125"]

for year in years:
    for sample in signal_samples:
        nominal_data = get_nominal_data(year, sample)
        
        for correction_name in photon_correction_name:
            up_data = get_correction_data(year, sample, correction_name, "up")
            down_data = get_correction_data(year, sample, correction_name, "down")
            
            # y-axis ratio of two subplots is 4:1
            fig, ax = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": [4, 1]}, sharex=True)
            
            nominal_hist = np.histogram(nominal_data["gamma_pt"], bins=100, range=(15, 95), weights=nominal_data["weight_central"])     
            nominal_error_hist = np.histogram(nominal_data["gamma_pt"], bins=100, range=(15, 95), weights=nominal_data["weight_central"]**2)
            up_hist = np.histogram(up_data["gamma_pt"], bins=100, range=(15, 95), weights=up_data["weight_central"])
            down_hist = np.histogram(down_data["gamma_pt"], bins=100, range=(15, 95), weights=down_data["weight_central"])
            
            ax[0].errorbar(nominal_hist[1][:-1], nominal_hist[0], yerr=np.sqrt(nominal_error_hist[0]), fmt="o", label="nominal", color="black")
            ax[0].fill_between(up_hist[1][:-1], up_hist[0], nominal_hist[0], alpha=0.5, color="red", label="up")
            ax[0].fill_between(down_hist[1][:-1], down_hist[0], nominal_hist[0], alpha=0.5, color="blue", label="down")
            
            ax[0].set_title(f"{sample} {year} photon pt {correction_name}")
            ax[0].set_xlabel("")
            ax[0].set_ylabel("Events")
            ax[0].legend()
            ax[0].set_ylim(0, 1.05*np.max(nominal_hist[0]))
            ax[0].set_xlim(15, 95)
            
            # ratio plot
            ratio_up = up_hist[0]/nominal_hist[0]
            ratio_down = down_hist[0]/nominal_hist[0]
            
            # ax[1].errorbar(nominal_hist[1][:-1], np.zeros_like(nominal_hist[0]), yerr=np.sqrt(nominal_error_hist[0])/nominal_hist[0], fmt="o", color="black")
            ax[1].fill_between(up_hist[1][:-1], ratio_up, np.ones_like(nominal_hist[0]), alpha=0.5, color="red", label="up")
            ax[1].fill_between(down_hist[1][:-1], ratio_down, np.ones_like(nominal_hist[0]), alpha=0.5, color="blue", label="down")
            
            ax[1].set_xlim(15, 95)
            
            ax[1].set_xlabel("photon pt")
            ax[1].set_ylabel("ratio")
            
            # save figure
            fig.tight_layout()
            fig.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{sample}_{year}_gamma_pt_{correction_name}.png")
            
            plt.clf()
            
        for correction_name in jet_correction_name:
            up_data = get_correction_data(year, sample, correction_name, "up")
            down_data = get_correction_data(year, sample, correction_name, "down")
            
            fig, ax = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": [4, 1]}, sharex=True)
            
            nominal_hist = np.histogram(nominal_data["jet_1_pt"], bins=100, range=(30, 100), weights=nominal_data["weight_central"])     
            nominal_error_hist = np.histogram(nominal_data["jet_1_pt"], bins=100, range=(30, 100), weights=nominal_data["weight_central"]**2)
            up_hist = np.histogram(up_data["jet_1_pt"], bins=100, range=(30, 100), weights=up_data["weight_central"])
            down_hist = np.histogram(down_data["jet_1_pt"], bins=100, range=(30, 100), weights=down_data["weight_central"])
            
            ax[0].errorbar(nominal_hist[1][:-1], nominal_hist[0], yerr=np.sqrt(nominal_error_hist[0]), fmt="o", label="nominal", color="black")
            ax[0].fill_between(up_hist[1][:-1], up_hist[0], nominal_hist[0], alpha=0.5, color="red", label="up")
            ax[0].fill_between(down_hist[1][:-1], down_hist[0], nominal_hist[0], alpha=0.5, color="blue", label="down")
            
            ax[0].set_title(f"{sample} {year} leading jet pt {correction_name}")
            ax[0].set_xlabel("")
            ax[0].set_ylabel("Events")
            ax[0].legend()
            ax[0].set_ylim(0, 1.05*np.max(nominal_hist[0]))
            ax[0].set_xlim(30, 100)
            
            # ratio plot
            ratio_up = up_hist[0]/nominal_hist[0]
            ratio_down = down_hist[0]/nominal_hist[0]
            
            # ax[1].errorbar(nominal_hist[1][:-1], np.zeros_like(nominal_hist[0]), yerr=np.sqrt(nominal_error_hist[0]), fmt="o", label="nominal", color="black")
            ax[1].fill_between(up_hist[1][:-1], ratio_up, np.ones_like(nominal_hist[0]), alpha=0.5, color="red", label="up")
            ax[1].fill_between(down_hist[1][:-1], ratio_down, np.ones_like(nominal_hist[0]), alpha=0.5, color="blue", label="down")
            
            ax[1].set_xlim(30, 100)
            
            ax[1].set_xlabel("leading jet pt")
            ax[1].set_ylabel("ratio")
            
            # save figure
            fig.tight_layout()
            fig.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{sample}_{year}_jet_1_pt_{correction_name}.png")
            
            plt.clf()
            
            fig, ax = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": [4, 1]}, sharex=True)
            
            nominal_hist = np.histogram(nominal_data["jet_2_pt"], bins=100, range=(30, 100), weights=nominal_data["weight_central"])     
            nominal_error_hist = np.histogram(nominal_data["jet_2_pt"], bins=100, range=(30, 100), weights=nominal_data["weight_central"]**2)
            up_hist = np.histogram(up_data["jet_2_pt"], bins=100, range=(30, 100), weights=up_data["weight_central"])
            down_hist = np.histogram(down_data["jet_2_pt"], bins=100, range=(30, 100), weights=down_data["weight_central"])
            
            ax[0].errorbar(nominal_hist[1][:-1], nominal_hist[0], yerr=np.sqrt(nominal_error_hist[0]), fmt="o", label="nominal", color="black")
            ax[0].fill_between(up_hist[1][:-1], up_hist[0], nominal_hist[0], alpha=0.5, color="red", label="up")
            ax[0].fill_between(down_hist[1][:-1], down_hist[0], nominal_hist[0], alpha=0.5, color="blue", label="down")
            
            ax[0].set_title(f"{sample} {year} subleading jet pt {correction_name}")
            ax[0].set_xlabel("")
            ax[0].set_ylabel("Events")
            ax[0].legend()
            ax[0].set_ylim(0, 1.05*np.max(nominal_hist[0]))
            ax[0].set_xlim(30, 100)
            
            # ratio plot
            ratio_up = up_hist[0]/nominal_hist[0]
            ratio_down = down_hist[0]/nominal_hist[0]
            
            # ax[1].errorbar(nominal_hist[1][:-1], np.zeros_like(nominal_hist[0]), yerr=np.sqrt(nominal_error_hist[0]), fmt="o", label="nominal", color="black")
            ax[1].fill_between(up_hist[1][:-1], ratio_up, np.ones_like(nominal_hist[0]), alpha=0.5, color="red", label="up")
            ax[1].fill_between(down_hist[1][:-1], ratio_down, np.ones_like(nominal_hist[0]), alpha=0.5, color="blue", label="down")
            
            ax[1].set_xlim(30, 100)
            
            ax[1].set_xlabel("subleading jet pt")
            ax[1].set_ylabel("ratio")
            
            # save figure
            fig.tight_layout()
            fig.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{sample}_{year}_jet_2_pt_{correction_name}.png")
            