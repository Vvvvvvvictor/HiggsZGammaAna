import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

config_dict = {
    # "mass_jj": {"range": (0, 300), "title": r"$m_{jj}(GeV/c^{2})$"},
    # "Z_relpt": {"range": (0, 1), "title": r"${p_T^{ll}\cdot c}/{m_{ll\gamma}}$"},
    # "gamma_relpt": {"range": (0.1, 1), "title": r"${p_T^{\gamma}\cdot c}/{m_{ll\gamma}}$"},
    # "pt_balance_0j": {"range": (0, 1), "title": r"$Zeppenfeld_{\gamma}(0j)$"},
    # "H_mass": {"range": (100, 180), "title": r"$m_{ll\gamma}(GeV/c^{2})$"},
    "delta_eta_jj": {"range": (0, 8), "bins": 40, "title": r"$\Delta\eta_{jj}$"},
    "delta_phi_jj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{jj}$"},
    "delta_phi_zgjj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{ll\gamma,jj}$"},
    "gamma_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta_{\gamma}$"},
    "gamma_mvaID": {"range": (0.14, 1), "bins": 43, "title": "photon MVA"},
    # "H_ptt": {"range": (0, 160), "title": r"$p_{T_{t}}^{ll\gamma}$"},
    "jet_1_pt": {"range": (30, 330), "bins": 50, "title": r"$p_{T_{j1}}(GeV/c)$"},
    "jet_2_pt": {"range": (30, 150), "bins": 50, "title": r"$p_{T_{j2}}(GeV/c)$"},
    "jet1G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(\gamma,j1)$"},
    "jet2G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(\gamma,j2)$"},
    # "jet_1_btagDeepFlavB": {"range": (0, 0.1), "title": "j1 btag"},
    # "jet_2_btagDeepFlavB": {"range": (0, 0.15), "title": "j2 btag"},
    "l1g_deltaR": {"range": (0.3, 4.3), "bins": 40, "title": r"max($\Delta R(l,\gamma)$)"},
    "l2g_deltaR": {"range": (0.3, 3.3), "bins": 40, "title": r"min($\Delta R(l,\gamma)$)"},
    "lep_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\theta$"},
    "lep_phi": {"range": (-3.2, 3.2), "bins": 40, "title": r"$\phi$"},
    "photon_zeppenfeld": {"range": (0, 5), "bins": 50, "title": "Zeppenfeld $\gamma$"},
    "pt_balance": {"range": (0, 1), "bins": 50, "title": "system balance"},
    "Z_cos_theta": {"range": (-1, 1), "bins": 50, "title": r"$\cos\Theta$"},
    "Z_lead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta_{l1}$"},
    "Z_sublead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta_{l2}$"},
    "H_relpt": {"range": (0, 3), "bins": 50, "title": r"${p_{T_{ll\gamma}}\cdot c}/{m_{ll\gamma}}$"},
    "gamma_ptRelErr": {"range": (0.01, 0.11), "bins": 50, "title": r"$\sigma_{p_T^{\gamma}}/p_T^{\gamma}$"},
    # "HZ_deltaRap": {"range": (-0.7, 0.7), "title": r"$\Delta y(ll,ll\gamma)$"}
}

TREE = "two_jet"
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run2/"

data = {"data": ["Data"]}
bkg = {r"Z$+\gamma$": ["ZGToLLG"], "Z+Fake Photon": ["DYJetsToLL", "EWKZ2J"], r"VBSZ+$\gamma$": ["ZG2JToG2L2J"], r"t$\bar{t}$": ["TT"], r"t$\gamma$/t$\bar{t}\gamma$": ["TTGJets", "TGJets"], "multiboson": ["WW", "WZ", "ZZ"], r"t$\bar{t}$+X": ["ttZJets", "ttWJets"]}
color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBSZ+$\gamma$": "#92dadd", r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", "multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2"}

def convert_root_to_hist(file_dict, selection=""):
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
                samples = samples.query(selection)
                hist, _ = np.histogram(samples["H_mass"], bins=80, range=[100, 180], weights=samples[WEIGHT])
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

    hist2, hist2_err, bins, mass2_hist = convert_root_to_hist(bkg, selection="(H_mass<120) | (H_mass>130)" if "H_mass" not in VAR else "")

    hist1, hist1_err, _, mass1_hist = convert_root_to_hist(data, selection="(H_mass<120) | (H_mass>130)")
    sf = sum(mass1_hist)/sum(mass2_hist)
    hist2 = [i*sf for i in hist2]

    points = (bins[:-1] + bins[1:]) / 2

    # split figure into 2x1
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=200, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
    ax1.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False, left=True, right=True, labelright=False)
    ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)

    # plt.figure(figsize=(8, 8), dpi=200)
    # plt.plot(bins, hist1/sum(hist1), "o-", label="Signal")
    if len(hist1) == 1:
        # ax1.errorbar(points, hist1[0], yerr=hist1_err, fmt="o", label=f"Data[N={np.sum(hist1):.1f}]", color="black", markersize=7, linewidth=3)
        ax1.errorbar(points, hist1[0], yerr=hist1_err, xerr=(bins[1:] - bins[:-1]) / 2, fmt="o", label=f"Data[N={np.sum(hist1):.1f}]", color="black", markersize=7, linewidth=3)
    else:
        raise ValueError("Data should be only 1")
    colors = [color_dict[i] for i in bkg]
    labels = [i+f"[N={np.sum(hist2[j]):.1f}]" for j, i in enumerate(bkg)]
    hep.histplot(hist2, bins, color=colors, label=labels, stack=True, histtype="fill", ax=ax1)
    # plt.errorbar(bins, hist2/sum(hist2), yerr=hist2_err/sum(hist2), fmt="s-", label="Bkg", color="blue", markersize=7, linewidth=3)
    
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
    ax1.set_ylim(0, 1.6*max(np.sum(hist1, axis=0).max(), np.sum(hist2, axis=0).max()))

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