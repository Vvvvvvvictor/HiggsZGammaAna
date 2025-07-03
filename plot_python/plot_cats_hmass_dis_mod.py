import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import os
import argparse # Import argparse
plt.style.use(hep.style.CMS)
from pdb import set_trace

# Define output directory
output_dir = "/eos/home-j/jiehan/root/outputs/"

# Add color dictionary
color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBSZ+$\gamma$": "#92dadd", 
              r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", 
              "multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2"}

# Add background label mapping
background_labels = {
    "ZGToLLG": r"Z$+\gamma$", 
    "DYJetsToLL": "Z+Fake Photon", 
    "EWKZ2J": r"VBSZ+$\gamma$"
}

def read_config(file_name):
    with open(file_name, "r") as f:
        return list(map(float, f.readline().strip().split(" ")))
  
def get_hist(file_name, tree_name, var_name, bin_boundaries, ranges):
    data = uproot.open(file_name)[tree_name].arrays([var_name, "weight", "bdt_score_t", "n_jets"], library="pd").query(f"bdt_score_t>{bin_boundaries[0]} & bdt_score_t<{bin_boundaries[1]}")
    weight = data["weight"]
    data = data[var_name]
    return np.histogram(data, weights=weight, range=ranges[1:], bins=ranges[0])

def get_err_hist(file_name, tree_name, var_name, bin_boundaries, ranges):
    data = uproot.open(file_name)[tree_name].arrays([var_name, "weight", "bdt_score_t", "n_jets"], library="pd").query(f"bdt_score_t>{bin_boundaries[0]} & bdt_score_t<{bin_boundaries[1]}")
    weight = data["weight"]
    data = data[var_name]
    return np.histogram(data, weights=weight**2, range=ranges[1:], bins=ranges[0])

# --- Add command line argument parsing ---
parser = argparse.ArgumentParser(description="Plot H_mass distributions for different categories.")
parser.add_argument("--channel", type=str, default="zero_to_one_jet", help="The analysis channel (e.g., two_jet, one_jet).")
args = parser.parse_args()
channel = args.channel
# --- End of new code ---

boundaries = read_config(f"{output_dir}test/significances/bin_boundaries_1D_{channel}.txt")
backgrounds = ["ZGToLLG", "DYJetsToLL", "EWKZ2J"]
signal = ["ggH_M125", "VBF_M125"] #, "ZH_M125", "WplusH_M125", "WminusH_M125", "ttH_M125"]
data = ["Data"]

output_plot_dir = f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{channel}/"
if not os.path.exists(output_plot_dir):
    os.makedirs(output_plot_dir)

folder = "test"  # Can be changed to "test" to switch folders
if boundaries[0] > 0.01:
    boundaries.insert(0, 0)

for j in range(-5, 6):
    for k in range(1, len(boundaries)-1):
        boundaries = read_config(f"{output_dir}{folder}/significances/bin_boundaries_1D_{channel}.txt")
        if boundaries[0] > 0.01:
            boundaries.insert(0, 0)
        boundaries[k] = boundaries[k] + j * 0.01
        for i in range(len(boundaries)-1):
            if ((i != k) and (i != k-1)): continue
            if boundaries[i] == boundaries[i+1]: continue
            bkg_hist = []
            bkg_err_hist = []
            for bkg in backgrounds:
                hist, bins = get_hist(f"{output_dir}{folder}/{channel}/{bkg}.root", channel, "H_mass", boundaries[i:], [85, 95, 180])
                bkg_hist.append(hist)
                # Calculate error for each background
                err_hist, _ = get_err_hist(f"{output_dir}{folder}/{channel}/{bkg}.root", channel, "H_mass", boundaries[i:], [85, 95, 180])
                bkg_err_hist.append(np.sqrt(err_hist))
                
            sig_hist = np.zeros(len(bins)-1)
            VBF_hist, bins = get_hist(f"{output_dir}{folder}/{channel}/VBF_M125.root", channel, "H_mass", boundaries[i:], [85, 95, 180])
            # sig_err_hist = np.zeros(len(bins)-1)
            for sig in signal:
                hist, bins = get_hist(f"{output_dir}{folder}/{channel}/{sig}.root", channel, "H_mass", boundaries[i:], [85, 95, 180])
                sig_hist += hist
                # sig_err_hist += get_err_hist(f"{output_dir}{folder}/{channel}/{sig}.root", channel, "H_mass", boundaries[i:], [85, 95, 180])[0]
                
            data_hist = np.zeros(len(bins)-1)
            for dat in data:
                hist, bins = get_hist(f"{output_dir}{folder}/{channel}/{dat}.root", channel, "H_mass", boundaries[i:], [85, 95, 180])
                data_hist += hist

            sf = np.sum(data_hist) / np.sum(bkg_hist)
            bkg_hist = [bkg * sf for bkg in bkg_hist]
            bkg_err_hist = [err * sf for err in bkg_err_hist]
            
            # Calculate total background and its error
            total_bkg = np.sum(bkg_hist, axis=0)
            total_bkg_err = np.sqrt(np.sum([err**2 for err in bkg_err_hist], axis=0))
            
            # Create subplot layout
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1})
            
            # Main plot
            labels = [background_labels[bkg] for bkg in backgrounds]
            colors = [color_dict[background_labels[bkg]] for bkg in backgrounds]
            hep.histplot(bkg_hist, bins, stack=True, label=labels, histtype="fill", color=colors, ax=ax1)   

            mask = (bins[:-1] + bins[1:]) / 2
            data_hist_masked = np.where((mask > 120) & (mask < 130), 0, data_hist)
            ax1.errorbar(mask, data_hist_masked, yerr=np.sqrt(data_hist_masked), fmt="o", label="data", color="black")
            
            # Add error band for total background
            bin_centers = (bins[:-1] + bins[1:]) / 2
            ax1.fill_between(bins, np.concatenate([[total_bkg[0] - total_bkg_err[0]], total_bkg - total_bkg_err]), 
                            np.concatenate([[total_bkg[0] + total_bkg_err[0]], total_bkg + total_bkg_err]), 
                            step='pre', color='gray', alpha=0.3, hatch='///', label='MC uncertainty', zorder=10)

            ax1.plot((bins[:-1] + bins[1:]) / 2, VBF_hist, label="VBF", color="yellow", linestyle="--", linewidth=2)
            ax1.plot((bins[:-1] + bins[1:]) / 2, sig_hist, label="signal", color="red", linewidth=3)
            ax1.text(1, 1, f"SF: {sf:.2f}", transform=ax1.transAxes, fontsize=16, horizontalalignment='right', verticalalignment='top')
            ax1.set_xlim(95, 180)
            print(f"Plotting category {len(boundaries)-1-i} with boundaries {boundaries[i:]} for k={k}, j={j}")
            print(f"Max data: {data_hist.max()}, Max background: {np.sum(bkg_hist, axis=0).max()}")
            ax1.set_ylim(0, max(np.sum(bkg_hist, axis=0).max(), data_hist.max()) * 1.1)
            ax1.set_ylabel("Events")
            ax1.set_title(f"Events in category {len(boundaries)-1-i} of {channel.replace('_', ' ')} channel")
            ax1.legend(ncol=2, fontsize=16, loc='best')
            ax1.set_xticklabels([])  # Remove x-axis labels from main plot
            
            # Ratio plot
            ratio = np.divide(data_hist_masked, total_bkg, out=np.zeros_like(data_hist_masked), where=total_bkg!=0)
            ratio_err = np.divide(np.sqrt(data_hist_masked), total_bkg, out=np.zeros_like(data_hist_masked), where=total_bkg!=0)
            
            # Background uncertainty in ratio plot
            bkg_ratio_err = np.divide(total_bkg_err, total_bkg, out=np.zeros_like(total_bkg_err), where=total_bkg!=0)
            
            ax2.fill_between(bins, np.concatenate([[1 - bkg_ratio_err[0]], 1 - bkg_ratio_err]), 
                            np.concatenate([[1 + bkg_ratio_err[0]], 1 + bkg_ratio_err]), 
                            step='pre', color='gray', alpha=0.3, hatch='///')
            ax2.axhline(y=1, color='black', linestyle='-', alpha=0.8)
            ax2.errorbar(mask, ratio, yerr=ratio_err, fmt="o", color="black")
            
            ax2.set_xlim(95, 180)
            ax2.set_ylim(0.5, 1.5)
            ax2.set_xlabel(r"$m_{ll\gamma}$ [GeV/c$^2$]")
            ax2.set_ylabel("Data/MC")
            ax2.grid(True, alpha=0.3)
            
            # plt.savefig(f"{output_plot_dir}hmass_dis_all.png")
            plt.savefig(f"{output_plot_dir}hmass_dis_{len(boundaries)-1-i}_k{k}_j{j}.png")
            plt.clf()
