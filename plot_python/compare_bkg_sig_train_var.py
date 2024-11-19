import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

config_dict = {
    # # "mass_jj": {"range": (0, 300), "title": r"$m_{jj}(GeV/c^{2})$"},
    # # "Z_relpt": {"range": (0, 1), "title": r"${p_T^{ll}\cdot c}/{m_{ll\gamma}}$"},
    # # "gamma_relpt": {"range": (0.1, 1), "title": r"${p_T^{\gamma}\cdot c}/{m_{ll\gamma}}$"},
    # # "pt_balance_0j": {"range": (0, 1), "title": r"$Zeppenfeld_{\gamma}(0j)$"},
    # # "H_mass": {"range": (100, 180), "title": r"$m_{ll\gamma}(GeV/c^{2})$"},
    # "delta_eta_jj": {"range": (0, 9), "title": r"$\Delta\eta_{jj}$"},
    # "delta_phi_jj": {"range": (0, 3.2), "title": r"$\Delta\phi_{jj}$"},
    # "delta_phi_zgjj": {"range": (0, 3.2), "title": r"$\Delta\phi_{ll\gamma,jj}$"},
    # "gamma_eta": {"range": (-2.5, 2.5), "title": r"$\eta_{\gamma}$"},
    # "gamma_mvaID": {"range": (0.14, 1), "title": "photon MVA"},
    # # "H_ptt": {"range": (0, 160), "title": r"$p_{T_{t}}^{ll\gamma}$"},
    # "jet_1_pt": {"range": (30, 330), "title": r"$p_{T_{j1}}(GeV/c)$"},
    # "jet_2_pt": {"range": (30, 150), "title": r"$p_{T_{j2}}(GeV/c)$"},
    # "jet1G_deltaR": {"range": (0.4, 6.4), "title": r"$\Delta R(\gamma,j1)$"},
    # "jet2G_deltaR": {"range": (0.4, 6.4), "title": r"$\Delta R(\gamma,j2)$"},
    # # "jet_1_btagDeepFlavB": {"range": (0, 0.1), "title": "j1 btag"},
    # # "jet_2_btagDeepFlavB": {"range": (0, 0.15), "title": "j2 btag"},
    # "l1g_deltaR": {"range": (0.4, 4.8), "title": r"max($\Delta R(l,\gamma)$)"},
    # "l2g_deltaR": {"range": (0.4, 3.4), "title": r"min($\Delta R(l,\gamma)$)"},
    # "lep_cos_theta": {"range": (-1, 1), "title": r"$\cos\theta$"},
    # "lep_phi": {"range": (-3.2, 3.2), "title": r"$\phi$"},
    # "photon_zeppenfeld": {"range": (0, 5), "title": "Zeppenfeld $\gamma$"},
    # "pt_balance": {"range": (0, 1), "title": "system balance"},
    # "Z_cos_theta": {"range": (-1, 1), "title": r"$\cos\Theta$"},
    # "Z_lead_lepton_eta": {"range": (-2.5, 2.5), "title": r"$\eta_{l1}$"},
    # "Z_sublead_lepton_eta": {"range": (-2.5, 2.5), "title": r"$\eta_{l2}$"},
    # "H_relpt": {"range": (0, 3), "title": r"${p_{T_{ll\gamma}}\cdot c}/{m_{ll\gamma}}$"},
    "gamma_ptRelErr": {"range": (0, 0.1), "title": r"$\sigma_{p_T^{\gamma}}/p_T^{\gamma}$"},
    # "HZ_deltaRap": {"range": (-0.7, 0.7), "title": r"$\Delta y(ll,ll\gamma)$"}
}

TREE = "two_jet"

BINS = 25
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run2/"

sig = ["ggH_mix", "VBF_mix"]
sig_st = ["ggH_M125", "VBF_M125", "WminusH_M125", "WplusH_M125", "ZH_M125", "ttH_M125"]
bkg= ["ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]

def convert_root_to_hist(file, target=None):
    type_hist = np.zeros(BINS)
    mass_hist = np.zeros(100)
    error = np.zeros(BINS)
    if target is not None:
        for f in file:
            for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
                print("Reading", PATH+f+"/"+year+f".root:{TREE}", "...")
                if "H_mass" not in VAR:
                    try:
                        samples = uproot.open(PATH+f+"/"+year+f".root:{TREE}").arrays([VAR, "weight", "H_mass"], library="pd")
                    except:
                        continue
                else:
                    try:
                        samples = uproot.open(PATH+f+"/"+year+f".root:{TREE}").arrays([VAR, "weight"], library="pd")
                    except:
                        continue
                hist, _ = np.histogram(samples["H_mass"], bins=100, range=[100, 180], weights=samples[WEIGHT])
                mass_hist = mass_hist + hist
                hist, _ = np.histogram(samples["H_mass"], bins=100, range=[100, 180], weights=samples[WEIGHT]**2)
                error = error + hist
        ratio_hist = np.divide(target, mass_hist, out=np.zeros_like(target), where=mass_hist != 0)
        error = np.sqrt(np.divide(error, mass_hist**2, out=np.zeros_like(error), where=mass_hist != 0))
        print(ratio_hist)

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
            hist, _ = np.histogram(samples["H_mass"], bins=100, range=[100, 180], weights=samples[WEIGHT])
            mass_hist = mass_hist + hist
            if target is not None:
                bin_indices = np.digitize(samples['H_mass'], bins=np.linspace(100, 180, 101)) - 1
                samples['weight'] = samples[WEIGHT] * np.where((bin_indices >= 0) & (bin_indices < 100), ratio_hist[bin_indices], 1)
            hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
            type_hist = type_hist + hist
            hist, _ = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT]**2)
            error = error + hist
    return type_hist, np.sqrt(error), bins, mass_hist

for i in config_dict:
    # if i != 4 :
    #     continue
    VAR = i
    XLABLE = config_dict[VAR]["title"]
    RMIN = config_dict[VAR]["range"][0]
    RMAX = config_dict[VAR]["range"][1]
    
    print("\n\n", VAR, RMIN, RMAX, "\n\n")

    hist2, hist2_err, bins, mass_hist = convert_root_to_hist(bkg)
    # hist1, bins, _ = convert_root_to_hist(sig, mass_hist)
    # hist1, bins, _ = convert_root_to_hist(sig)
    hist1_1, hist1_1_error, _, _ = convert_root_to_hist(sig_st)

    bins = (bins[:-1] + bins[1:]) / 2

    plt.figure(figsize=(8, 8), dpi=200)
    # plt.plot(bins, hist1/sum(hist1), "o-", label="Signal")
    plt.errorbar(bins, hist1_1/sum(hist1_1), yerr=hist1_1_error/sum(hist1_1), fmt="o-", label="Sig", color="red", markersize=7, linewidth=3)
    plt.errorbar(bins, hist2/sum(hist2), yerr=hist2_err/sum(hist2), fmt="s-", label="Bkg", color="blue", markersize=7, linewidth=3)    
    
    plt.xlabel(XLABLE, fontsize=28)
    plt.ylabel("Events", fontsize=28)
    plt.legend(fontsize=24)
    plt.grid()
    plt.tight_layout()
    plt.title(f"Signal Background comparison", fontsize=32)
    if os.path.exists("pic/train_var") == False:
        os.makedirs("pic/train_var")
    plt.savefig(f"pic/train_var/{VAR}.png")
    plt.clf()