# sw_lumiXyear
# slly_m
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import differential_evolution

plt.style.use(hep.style.CMS)

signal_names = ["ggH_M125", "VBF_M125"]
background_names = ["DYJetsToLL", "ZGToLLG"]

filepath = "/eos/user/j/jiehan/root/outputs/two_jet/"

# Load data from ROOT files into pandas dataframes
sig_data, bkg_data = [], []
for name in signal_names:
    sig_data.append(uproot.open(filepath + name + ".root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "w_lumiXyear", "lly_m"], library="pd"))
for name in background_names:
    bkg_data.append(uproot.open(filepath + name + ".root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "w_lumiXyear", "lly_m"], library="pd"))

# Create 2D histograms for signal and background
bins = 100
sig_hist, bkg_hist, full_bkg_hist = np.zeros((bins, bins)), np.zeros((bins, bins)), np.zeros((bins, bins))
for sig in sig_data:
    sig = sig.query("lly_m > 120 & lly_m < 130")
    sig_hist += np.histogram2d(sig["bdt_score_t"], sig["vbf_score_t"], bins=bins, weights=sig["w_lumiXyear"])[0]
for bkg in bkg_data:
    bkg = bkg.query("lly_m > 100")
    full_bkg_hist += np.histogram2d(bkg["bdt_score_t"], bkg["vbf_score_t"], bins=bins, weights=bkg["w_lumiXyear"])[0]
    bkg = bkg.query("lly_m > 120 & lly_m < 130")
    bkg_hist += np.histogram2d(bkg["bdt_score_t"], bkg["vbf_score_t"], bins=bins, weights=bkg["w_lumiXyear"])[0]

# Function to calculate the sum of the squares of s/b
def calculate_significance(params):
    n, m = int(params[0]), int(params[1])
    bin_edges_x = [0, n, bins]
    bin_edges_y = [0, m, bins]
    
    significance = 0
    for i in range(len(bin_edges_x)-1):
        for j in range(len(bin_edges_y)-1):
            s = np.sum(sig_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            b = np.sum(bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            full_b = np.sum(full_bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            if full_b < 10 or b < 2:
                significance += -10
            else:
                significance += 2 * ( ( s + b ) * np.log( 1 + s / b ) - s )
    
    return -significance  # We minimize the negative to maximize the positive

# Use differential evolution to find the best boundaries
bounds = [(1, bins), (1, bins)]
result = differential_evolution(calculate_significance, bounds, strategy='best1bin', maxiter=1000, popsize=15, tol=0.001)

best_n, best_m = int(result.x[0]), int(result.x[1])
best_significance = np.sqrt(-result.fun)

print(f"The best boundaries are n={best_n}, m={best_m} with a sum of squares of significance: {best_significance}")

# print significance in each category
bin_edges_x = [0, best_n, bins]
bin_edges_y = [0, best_m, bins]
for i in range(len(bin_edges_x)-1):
    for j in range(len(bin_edges_y)-1):
        s = np.sum(sig_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
        b = np.sum(bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
        print(f"Category {j * 2**(j) + i }: s={s}, b={b}, significance={np.sqrt(2 * ( ( s + b ) * np.log( 1 + s / b ) - s ))}")
        
with open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/best_boundaries.txt", "w") as f:
    f.write(f"{best_n/100:.2f},{best_m/100:.2f}")