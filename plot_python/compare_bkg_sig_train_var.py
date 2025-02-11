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

sig = {"sig": ["ggH_M125", "VBF_M125"], "ggH": ["ggH_M125"], "VBF": ["VBF_M125"]}
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

for var in config_dict:
    # if i != 4 :
    #     continue
    VAR = var
    XLABLE = config_dict[VAR]["title"]
    RMIN = config_dict[VAR]["range"][0]
    RMAX = config_dict[VAR]["range"][1]
    BINS = config_dict[VAR]["bins"]
    
    print("\n\n", VAR, RMIN, RMAX, "\n\n")

    hist2, hist2_err, bins, mass2_hist = convert_root_to_hist(bkg, selection="(H_mass>100) & (H_mass<180)")
    hist1, hist1_err, _, mass1_hist = convert_root_to_hist(sig, selection="(H_mass>100) & (H_mass<180)")

    sig_yields = [np.sum(i) for i in hist1]
    bkg_yields = [np.sum(i) for i in hist2]
    hist1 = [i/sig_yields[0] for i in hist1]
    hist1_err = [i/sig_yields[0] for i in hist1_err]
    hist2 = [i/np.sum(bkg_yields) for i in hist2]

    points = (bins[:-1] + bins[1:]) / 2

    fig, ax1 = plt.subplots(1, 1, figsize=(8, 8), dpi=200)
    ax1.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)

    colors = [color_dict[i] for i in bkg]
    labels = [i+f"[N={bkg_yields[j]:.1f}]" for j, i in enumerate(bkg.keys())]
    hep.histplot(hist2, bins, color=colors, label=labels, stack=True, histtype="fill", ax=ax1)
    for j, i in enumerate(hist1[0:]):
        hep.histplot(i, bins, color=color_dict[list(sig.keys())[j]], label=list(sig.keys())[j]+f"[N={sig_yields[j]:.1f}]", histtype="step", ax=ax1, linewidth=3)
    
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
    ax1.set_ylim(0, 1.6*max(hist1[0].max(), np.sum(hist2, axis=0).max()))
    ax1.annotate(rf"L=137.2 fb$^{{-1}}$", xy=(1, 1.01), xycoords='axes fraction', fontsize=16, ha="right")
    ax1.set_xlabel(XLABLE, fontsize=24)

    plt.tight_layout()
    if os.path.exists("pic/sigVbkg") == False:
        os.makedirs("pic/sigVbkg")
    plt.savefig(f"pic/sigVbkg/{VAR}.pdf")
    plt.clf()