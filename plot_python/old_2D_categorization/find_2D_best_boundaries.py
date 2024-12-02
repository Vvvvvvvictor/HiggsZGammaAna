# sw_lumiXyear
# slly_m
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import differential_evolution
from pdb import set_trace

plt.style.use(hep.style.CMS)

signal_names = ['ggH_M125','VBF_M125']
background_names = ["ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J"]
type = "fullSimrw" # "fullSim" or "fullSimrw"
if type == "fullSimrw":
    data_names = ["Data"] 

filepath = "/eos/user/j/jiehan/root/outputs/two_jet/"

# Load data from ROOT files into pandas dataframes
sig_data, bkg_data = [], []
for name in signal_names:
    sig_data.append(uproot.open(filepath + name + ".root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight", "H_mass"], library="pd"))
for name in background_names:
    bkg_data.append(uproot.open(filepath + name + ".root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight", "H_mass"], library="pd"))
if type == "fullSimrw":
    data = uproot.open(filepath + "Data.root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight", "H_mass"], library="pd")

# Create 2D histograms for signal and background
bins = 100
sig_hist, bkg_hist, full_bkg_hist = np.zeros((bins, bins)), np.zeros((bins, bins)), np.zeros((bins, bins))
for sig in sig_data:
    sig = sig.query("H_mass > 120 & H_mass < 130")
    sig_hist += np.histogram2d(sig["bdt_score_t"], sig["vbf_score_t"], bins=bins, weights=sig["weight"])[0]
for bkg in bkg_data:
    bkg = bkg.query("H_mass > 100")
    full_bkg_hist += np.histogram2d(bkg["bdt_score_t"], bkg["vbf_score_t"], bins=bins, weights=bkg["weight"])[0]
    bkg = bkg.query("H_mass > 120 & H_mass < 130")
    bkg_hist += np.histogram2d(bkg["bdt_score_t"], bkg["vbf_score_t"], bins=bins, weights=bkg["weight"])[0]
if type == "fullSimrw":
    data = data.query("H_mass < 120 | H_mass > 130")
    data_hist = np.histogram2d(data["bdt_score_t"], data["vbf_score_t"], bins=bins, weights=data["weight"])[0]

# Function to calculate the sum of the squares of s/b
def calculate_significance(params):
    bin_edges_x = [0, bins]
    bin_edges_y = [0, bins]
    for p in params[:-1]:
        bin_edges_x.insert(-1, int(p))
    bin_edges_y.insert(-1, int(params[-1]))
    if not (bin_edges_x == sorted(bin_edges_x) and bin_edges_y == sorted(bin_edges_y)):
        return 0
    
    significance = 0
    for i in range(len(bin_edges_x)-1):
        for j in range(len(bin_edges_y)-1):
            s = np.sum(sig_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            b = np.sum(bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            full_b = np.sum(full_bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            if (full_b-b) < 10 or b < 2 or full_b < b:
                significance += 0
            else:
                if type == "fullSimrw":
                    d = np.sum(data_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
                    b_cor = b * (d / (full_b - b))
                    if d < 10 or b_cor < 2:
                        significance += 0
                    else:
                        significance += 2 * ( ( s + b ) * np.log( 1 + s / b ) - s )
                else:
                    significance += 2 * ( ( s + b ) * np.log( 1 + s / b ) - s )
               
    return -significance  # We minimize the negative to maximize the positive

# Use differential evolution to find the best boundaries
nx, ny = 1, 1
bounds = [(1, bins)] * (nx + ny)
result = differential_evolution(calculate_significance, bounds, strategy='best1bin', maxiter=2000, popsize=15, tol=0.0001, mutation=(0.5, 1), recombination=0.7, seed=42, polish=True, disp=True, updating='deferred')

best_p= list(map(int, list(result.x)))
best_significance = np.sqrt(-result.fun)

print(f"The best boundaries are {best_p} with a significance of {best_significance}")

# print significance in each category
bin_edges_x = [0, bins]
bin_edges_y = [0, bins]
for p in best_p[:-1]:
    bin_edges_x.insert(-1, p)
bin_edges_y.insert(-1, best_p[-1])
significances = []
for i in range(len(bin_edges_x)-1):
    for j in range(len(bin_edges_y)-1):
        s = np.sum(sig_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
        b = np.sum(bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
        if type == "fullSimrw":
            raw_b = b
            full_b = np.sum(full_bkg_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            d = np.sum(data_hist[int(bin_edges_x[i]):int(bin_edges_x[i+1]), int(bin_edges_y[j]):int(bin_edges_y[j+1])])
            b = b * (d / (full_b - b))
            significance = np.sqrt(2 * ( ( s + b ) * np.log( 1 + s / b ) - s ))
            print(f"Category {(len(bin_edges_x) -1) * j + i}: s={s}, b={b}, raw_b={raw_b}, d={d}, bsb={full_b-b}, significance={np.sqrt(2 * ( ( s + b ) * np.log( 1 + s / b ) - s ))}")
        else:
            significance = np.sqrt(2 * ( ( s + b ) * np.log( 1 + s / b ) - s ))
            print(f"Category {(len(bin_edges_x) -1) * j + i}: s={s}, b={b}, significance={np.sqrt(2 * ( ( s + b ) * np.log( 1 + s / b ) - s ))}")
        if b < 2:
            significance = 0
        if type == "fullSimrw":
            if d < 10 or (full_b-b) < 10 or full_b < b:
                significance = 0
        significances.append(significance)
            
with open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/best_boundaries.txt", "w") as f:
    str = ",".join([f"{p/100:2}" for p in best_p])
    f.write(str)
    f.write("\n")
    str = ",".join([f"{s:.4f}" for s in significances])
    f.write(str)
    