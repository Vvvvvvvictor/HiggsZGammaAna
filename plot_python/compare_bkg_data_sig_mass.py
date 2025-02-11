import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

config_dict = {
    "H_mass": {"range": (100, 180), "bins": 80, "title": r"$m_{ll\gamma}(GeV/c^{2})$"}
}

TREE = "two_jet"
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run2/"

sig = {"sig": ["ggH_M125", "VBF_M125", "WminusH_M125", "WplusH_M125", "ZH_M125", "ttH_M125"], "ggH": ["ggH_M125"], "VBF": ["VBF_M125"]}
data = {"data": ["Data"]}
bkg = {r"Z$+\gamma$": ["ZGToLLG"], "Z+Fake Photon": ["DYJetsToLL", "EWKZ2J"], r"VBSZ+$\gamma$": ["ZG2JToG2L2J"], r"t$\bar{t}$": ["TT"], r"t$\gamma$/t$\bar{t}\gamma$": ["TTGJets", "TGJets"], "multiboson": ["WW", "WZ", "ZZ"], r"t$\bar{t}$+X": ["ttZJets", "ttWJets"]}
color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBSZ+$\gamma$": "#92dadd", r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", "multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2", "data": "black", "sig": "red", "ggH": "magenta", "VBF": "green"}

def convert_root_to_hist(file_dict, selection=None):
    mass_hist = np.zeros(80)
    error = np.zeros(BINS)
    hists = []
    for file in file_dict.values():
        type_hist = np.zeros(BINS)
        for f in file:
            for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
                print("Reading", PATH+f+"/"+year+f".root:{TREE}", "...")
                if "H_mass" not in VAR:
                    try:
                        samples = uproot.open(PATH+f+"/"+year+f".root:{TREE}").arrays([VAR, WEIGHT, "H_mass"], library="pd")
                    except:
                        continue
                else:
                    try:
                        samples = uproot.open(PATH+f+"/"+year+f".root:{TREE}").arrays([VAR, WEIGHT], library="pd")
                    except:
                        continue
                if selection != None:
                    samples = samples.query(selection)
                hist, _ = np.histogram(samples.query("H_mass<120 | H_mass>130")["H_mass"], bins=80, range=[100, 180], weights=samples.query("H_mass<120 | H_mass>130")[WEIGHT])
                mass_hist = mass_hist + hist
                hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
                hist[0] += np.sum(samples[WEIGHT][samples[VAR] < RMIN])
                hist[-1] += np.sum(samples[WEIGHT][samples[VAR] > RMAX])
                type_hist = type_hist + hist
                hist, _ = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT]**2)
                hist[0] += np.sum(samples[WEIGHT][samples[VAR] < RMIN])**2
                hist[-1] += np.sum(samples[WEIGHT][samples[VAR] > RMAX])**2
                error = error + hist
        hists.append(type_hist)
    return hists, np.sqrt(error), bins, mass_hist

for i in config_dict:
    # if i != 4 :
    #     continue
    VAR = i
    XLABLE = config_dict[VAR]["title"]
    RMIN = config_dict[VAR]["range"][0]
    RMAX = config_dict[VAR]["range"][1]
    BINS = config_dict[VAR]["bins"]
    
    print("\n\n", VAR, RMIN, RMAX, "\n\n")

    hist2, hist2_err, bins, mass2_hist = convert_root_to_hist(bkg, selection="((H_mass<120) | (H_mass>130)) & (H_mass>100) & (H_mass<180)" if "H_mass" not in VAR else None)
    hist1, hist1_err, _, mass1_hist = convert_root_to_hist(data, selection="((H_mass<120) | (H_mass>130)) & (H_mass>100) & (H_mass<180)")
    hist3, hist3_err, _, mass3_hist = convert_root_to_hist(sig, selection="(H_mass>100) & (H_mass<180)")
    sf = sum(mass1_hist)/sum(mass2_hist)
    hist2 = [i*sf for i in hist2]

    points = (bins[:-1] + bins[1:]) / 2

    # split figure into 2x1
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=200, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
    ax1.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False, left=True, right=True, labelright=False)
    ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)

    if len(hist1) == 1:
        ax1.errorbar(points, hist1[0], yerr=hist1_err, xerr=(bins[1:] - bins[:-1]) / 2, fmt="o", label=f"Data[N={np.sum(hist1):.1f}]", color="black", markersize=7, linewidth=3)
    else:
        raise ValueError("Data should be only 1")
    colors = [color_dict[i] for i in bkg]
    labels = [i+f"[N={np.sum(hist2[j]):.1f}]" for j, i in enumerate(bkg)]
    hep.histplot(hist2, bins, color=colors, label=labels, stack=True, histtype="fill", ax=ax1)
    for j, i in enumerate(hist3):
        hep.histplot(i*50, bins, color=color_dict[list(sig.keys())[j]], label=list(sig.keys())[j]+f"x50[N={np.sum(hist3[j])*50:.1f}]", histtype="step", ax=ax1, linewidth=3)

    bin_width = (RMAX - RMIN) / BINS
    if bin_width < 0.01:
        y_label = f"Events/{bin_width:.4f}"
    elif bin_width < 0.1:
        y_label = f"Events/{bin_width:.3f}"
    elif bin_width < 1:
        y_label = f"Events/{bin_width:.2f}"
    else:
        y_label = f"Events/{bin_width:.1f}"
    ax1.set_ylabel(y_label, fontsize=24)
    ax1.legend(fontsize=14, ncol=2, handletextpad=0.4, columnspacing=0.5)
    ax1.grid()
    ax1.set_xlim(RMIN, RMAX)
    ax1.set_ylim(0, 1.7*max(np.sum(hist1, axis=0).max(), np.sum(hist2, axis=0).max()))

    ax1.annotate(rf"L=137.2 fb$^{{-1}}$, sf={sf:.2f}", xy=(1, 1.01), xycoords='axes fraction', fontsize=16, ha="right")

    data_sum, bkg_sum = np.where(np.sum(hist1, axis=0) == 0, 1e-8, np.sum(hist1, axis=0)), np.where(np.sum(hist2, axis=0) == 0, 1e-8, np.sum(hist2, axis=0))
    ratio = hist1[0]/bkg_sum
    ratio_err = hist1_err/data_sum
    err = hist2_err/bkg_sum

    mask = (ratio > 0) & (ratio < 2)
    ax2.errorbar(points[mask], ratio[mask], yerr=ratio_err[mask], xerr=((bins[1:] - bins[:-1]) / 2)[mask], fmt="o", color="black", markersize=7, linewidth=3)
    ax2.bar(points, 2 * err, bottom=1 - err, color="#ffedcd", width=(RMAX-RMIN)/BINS, label="Bkg Uncertainty")
    ax2.set_xlim(RMIN, RMAX)
    ax2.set_ylim(0, 2)
    ax2.set_ylabel("Data/Bkg", fontsize=24)
    ax2.axhline(1, color="black", linestyle="--")
    ax2.set_yticks(np.arange(0, 2.1, 0.5))
    ax2.set_yticklabels([f"{i}" if i in [0.5, 1, 1.5] else "" for i in np.arange(0, 2.1, 0.5)])

    ax2.set_xlabel(XLABLE, fontsize=24)

    plt.tight_layout()
    if os.path.exists("pic/dataVbkg") == False:
        os.makedirs("pic/dataVbkg")
    plt.savefig(f"pic/dataVbkg/{VAR}.pdf")
    plt.clf()