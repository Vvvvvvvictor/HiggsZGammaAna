import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

config_dict = {
    "mass_jj": {"range": (0, 300), "title": r"$m_{jj}(GeV/c^{2})$"},
    # "Z_relpt": {"range": (0, 1), "title": r"${p_T^{ll}\cdot c}/{m_{ll\gamma}}$"},
    # "gamma_relpt": {"range": (0.1, 1), "title": r"${p_T^{\gamma}\cdot c}/{m_{ll\gamma}}$"},
    # "pt_balance_0j": {"range": (0, 1), "title": r"$Zeppenfeld_{\gamma}(0j)$"},
    # "H_mass": {"range": (100, 180), "title": r"$m_{ll\gamma}(GeV/c^{2})$"},
    # "H_ptt": {"range": (0, 160), "title": r"$p_{T_{t}}^{ll\gamma}$"},
    # "jet_1_btagDeepFlavB": {"range": (0, 0.1), "title": "j1 btag"},
    # "jet_2_btagDeepFlavB": {"range": (0, 0.15), "title": "j2 btag"},
    # "jet_1_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta_{j1}$"},
    # "jet_2_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta_{j2}$"},
    # "HZ_deltaRap": {"range": (-0.7, 0.7), "title": r"$\Delta y(ll,ll\gamma)$"},
    # "gamma_pt": {"range": (15, 100), "bins": 50, "title": r"$p_{T_{\gamma}}(GeV/c)$"},
    # "Z_pt": {"range": (0, 200), "bins": 50, "title": r"$p_{T_{ll}}(GeV/c)$"},
    # "Z_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta_{ll}$"},
    # "H_pt": {"range": (0, 250), "bins": 50, "title": r"$p_{T_{ll\gamma}}(GeV/c)$"},
    # "H_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta_{ll\gamma}$"},
    # "jet_pair_pt": {"range": (0, 300), "bins": 50, "title": r"$p_{T_{jj}}(GeV/c)$"},
    # "system_pt": {"range": (0, 250), "bins": 50, "title": r"$p_{T_{ll\gamma jj}}(GeV/c)$"},
    # "delta_eta_zgjj": {"range": (0, 8), "bins": 40, "title": r"$\Delta\eta_{zg,jj}$"},

    "llphoton_hmiss_photon_dphi": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{ll\gamma}$"},

    "delta_eta_jj": {"range": (0, 8), "bins": 40, "title": r"$\Delta\eta_{jj}$"},
    "delta_phi_jj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{jj}$"},
    "delta_phi_zgjj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{ll\gamma,jj}$"},
    "gamma_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta_{\gamma}$"},
    "gamma_mvaID": {"range": (0.14, 1), "bins": 43, "title": "photon MVA"},
    "jet_1_pt": {"range": (30, 330), "bins": 50, "title": r"$p_{T_{j1}}(GeV/c)$"},
    "jet_2_pt": {"range": (30, 150), "bins": 50, "title": r"$p_{T_{j2}}(GeV/c)$"},
    "jet1G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(\gamma,j1)$"},
    "jet2G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(\gamma,j2)$"},
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
    "gamma_ptRelErr": {"range": (0.01, 0.11), "bins": 50, "title": r"$\sigma_{p_T^{\gamma}}/p_T^{\gamma}$"}
}

TREE = "two_jet"
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run2/"
years = ["2018"]

data = {"data": ["Data"]}
bkg = {r"Z$+\gamma$": ["ZGToLLG"], "Z+Fake Photon": ["DYJetsToLL", "EWKZ2J"], r"VBSZ+$\gamma$": ["ZG2JToG2L2J"], r"t$\bar{t}$": ["TT"], r"t$\gamma$/t$\bar{t}\gamma$": ["TTGJets", "TGJets"], "multiboson": ["WW", "WZ", "ZZ", "WWG", "WZG", "ZZG",], r"t$\bar{t}$+X": ["ttZJets", "ttWJets"]}
color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBSZ+$\gamma$": "#92dadd", r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", "multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2"}

def get_sf(bkg_dict, data_dict, sf_var, range, selection=""):
    var_set = [sf_var, WEIGHT, "H_mass"]
    for var in selection.replace("|", "&").split("&"):
        var = var.strip().strip("()").split("<")[0].split(">")[0].split("=")[0].strip()
        if var not in var_set:
            var_set.append(var)
    bins = np.linspace(range[0], range[1], range[2]+1)
    bins[0] = -999
    bins[-1] = np.inf
    print(bins)
    sf, bkg = np.zeros(len(bins)-1), np.zeros(len(bins)-1)
    data_list = []
    for key, datasets in data_dict.items():
        data_list += datasets
    for dataset in data_list:
        for year in years:
            print(f"Reading {PATH}{dataset}/{year}.root:{TREE} ...")
            samples = uproot.open(f"{PATH}{dataset}/{year}.root:{TREE}").arrays(var_set, library="pd")
            samples = samples.query(selection)
            hist, _ = np.histogram(samples[sf_var], bins=bins, weights=samples[WEIGHT])
            sf += hist
    bkg_list = []
    for key, datasets in bkg_dict.items():
        bkg_list += datasets
    for dataset in bkg_list:
        for year in years:
            print(f"Reading {PATH}{dataset}/{year}.root:{TREE} ...")
            samples = uproot.open(f"{PATH}{dataset}/{year}.root:{TREE}").arrays(var_set, library="pd")
            samples = samples.query(selection)
            hist, _ = np.histogram(samples[sf_var], bins=bins, weights=samples[WEIGHT])
            bkg += hist
    sf = sf / bkg
    return bins, sf


def convert_root_to_hist(file_dict, sf_var=None, sf_bins=None, sf=None, selection=""):
    mass_hist = np.zeros(80)
    error = np.zeros(BINS)
    hists = []
    for file in file_dict.values():
        type_hist = np.zeros(BINS)
        for f in file:
            for year in years:
                print("Reading", PATH+f+"/"+year+f".root:{TREE}", "...")
                var_set = [VAR, WEIGHT]
                if sf_var is not None and sf_var not in var_set:
                    var_set.append(sf_var)
                if "H_mass" not in VAR:
                    var_set.append("H_mass")
                for var in selection.replace("|", "&").split("&"):
                    var = var.strip().strip("()").split("<")[0].split(">")[0].split("=")[0].strip()
                    if var not in var_set:
                        var_set.append(var)
                try:
                    samples = uproot.open(PATH+f+"/"+year+f".root:{TREE}").arrays(var_set, library="pd")
                except:
                    continue
                samples = samples.query(selection)
                weight = samples[WEIGHT]
                if sf_var is not None and sf_var in samples and sf_bins is not None and sf is not None:
                    sf_weight = np.ones(len(samples))
                    for i, (low, high) in enumerate(zip(sf_bins[:-1], sf_bins[1:])):
                        sf_weight[(samples[sf_var] > low) & (samples[sf_var] < high)] = sf[i]
                    weight = weight * sf_weight
                samples[WEIGHT] = weight
                hist, _ = np.histogram(samples.query("H_mass<120 | H_mass>130")["H_mass"], bins=80, range=[100, 180], weights=samples.query("H_mass<120 | H_mass>130")[WEIGHT])
                mass_hist = mass_hist + hist
                hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
                hist[0] += np.sum(samples[WEIGHT][samples[VAR] < RMIN])
                hist[-1] += np.sum(samples[WEIGHT][samples[VAR] > RMAX])
                type_hist = type_hist + hist
                hist, _ = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT]**2)
                hist[0] += np.sum(samples[WEIGHT][samples[VAR] < RMIN]**2)
                hist[-1] += np.sum(samples[WEIGHT][samples[VAR] > RMAX]**2)
                error = error + hist
        hists.append(type_hist)
    print(np.sqrt(error))
    print(sum(hists))
    return hists, np.sqrt(error), bins, mass_hist

mass_selection = "((H_mass<120) | (H_mass>130))"
selection = "jet_2_pt < 50"
# sf_var = "delta_phi_zgjj"
# sf_bins, sf = get_sf(bkg, data, sf_var, (2, 3.14, 20), selection=f"{mass_selection} & {selection}")

for i in config_dict:
    # if i != 4 :
    #     continue
    VAR = i
    XLABLE = config_dict[VAR]["title"]
    RMIN = config_dict[VAR]["range"][0]
    RMAX = config_dict[VAR]["range"][1]
    BINS = config_dict[VAR]["bins"]
    
    print("\n\n", VAR, RMIN, RMAX, "\n\n")

    hist1, hist1_err, _, mass1_hist = convert_root_to_hist(data, selection=f"{mass_selection} & {selection}")
    hist2, hist2_err, bins, mass2_hist = convert_root_to_hist(bkg, selection=f"{mass_selection} & {selection}" if "H_mass" not in VAR else selection)

    mass_sf = sum(mass1_hist)/sum(mass2_hist)
    hist2 = [i*mass_sf for i in hist2]

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

    ax1.annotate(rf"L=41.48 fb$^{{-1}}$, sf={mass_sf:.2f}", xy=(1, 1.01), xycoords='axes fraction', fontsize=16, ha="right")

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