import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mplhep as hep
import re
import os
from pdb import set_trace
plt.style.use(hep.style.CMS)

config_dict = {
    "mass_jj": {"range": (0, 300, 60), "title": r"$m_{jj}(GeV/c^{2})$", "match_var": "jj_m"},
    "delta_eta_jj": {"range": (0, 9, 45), "title": r"$\Delta\eta_{jj}$", "match_var": "jj_deta"},
    "delta_phi_jj": {"range": (0, 3.2, 40), "title": r"$\Delta\phi_{jj}$", "match_var": "jj_dphi"},
    "delta_phi_zgjj": {"range": (0, 3.2, 40), "title": r"$\Delta\phi_{ll\gamma,jj}$", "match_var": "llyjj_dphi"},
    "gamma_eta": {"range": (-2.5, 2.5, 50), "title": r"$\eta_{\gamma}$", "match_var": "y_eta"},
    "gamma_mvaID": {"range": (0.14, 1, 50), "title": "photon MVA", "match_var": "y_mva"},
    "H_ptt": {"range": (0, 160, 40), "title": r"$p_{T_{t}}^{ll\gamma}$", "match_var": "lly_ptt"},
    "jet_1_pt": {"range": (30, 330, 50), "title": r"$p_{T_{j1}}(GeV/c)$", "match_var": "j1_pt"},
    "jet_2_pt": {"range": (30, 150, 40), "title": r"$p_{T_{j2}}(GeV/c)$", "match_var": "j2_pt"},
    "jet1G_deltaR": {"range": (0.4, 6.4, 50), "title": r"$\Delta R(\gamma,j1)$", "match_var": "yj1_dr"},
    "jet2G_deltaR": {"range": (0.4, 6.4, 50), "title": r"$\Delta R(\gamma,j2)$", "match_var": "yj2_dr"},
    "l1g_deltaR": {"range": (0.4, 4.9, 50), "title": r"max($\Delta R(l,\gamma)$)", "match_var": "yl_drmax"},
    "l2g_deltaR": {"range": (0.4, 3.4, 50), "title": r"min($\Delta R(l,\gamma)$)", "match_var": "yl_drmin"},
    "lep_cos_theta": {"range": (-1, 1, 50), "title": r"$\cos\theta$", "match_var": "costheta"},
    "lep_phi": {"range": (-3.2, 3.2, 40), "title": r"$\phi$", "match_var": "phi"},
    "photon_zeppenfeld": {"range": (0, 5, 50), "title": "Zeppenfeld $\gamma$", "match_var": "llyjj_zep"},
    "pt_balance": {"range": (0, 1, 50), "title": "system balance", "match_var": "llyjj_ptbal"},
    "Z_cos_theta": {"range": (-1, 1, 50), "title": r"$\cos\Theta$", "match_var": "cosTheta"},
    "Z_lead_lepton_eta": {"range": (-2.5, 2.5, 50), "title": r"$\eta_{l1}$", "match_var": "l1_eta"},
    "Z_sublead_lepton_eta": {"range": (-2.5, 2.5, 50), "title": r"$\eta_{l2}$", "match_var": "l2_eta"},
    "H_relpt": {"range": (0, 3, 50), "title": r"${p_{T_{ll\gamma}}\cdot c}/{m_{ll\gamma}}$", "match_var": "lly_ptmass"},
    "gamma_ptRelErr": {"range": (0, 0.1, 50), "title": r"$\sigma_{p_T^{\gamma}}/p_T^{\gamma}$", "match_var": "y_res"},
    "n_jets": {"range": (2, 6, 5), "title": "n_jets", "match_var": "njet"}
}

for i in config_dict.keys():
    print(f"{config_dict[i]['match_var']} to {i}")
exit()

var_map = {
    "H_mass": "lly_m",
    "weight": "w_lumiXyear",
    "Z_mass": "leplep_m",
    "gamma_relpt": "y_ptmass"
}

def extract_unit(label):
    match = re.search(r"\((.*?)\)", label)
    return rf"${match.group(1)}$" if match else ""

def plot_2_hist(data1, data2, var, range, xlabel, save_path):
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(8, 8), dpi=200, gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    
    hist1, bins = np.histogram(data1[var], bins=range[2], range=(range[0], range[1]), weights=data1["weight"])
    hist1_err, _ = np.histogram(data1[var], bins=range[2], range=(range[0], range[1]), weights=data1["weight"]**2)
    hist2, _ = np.histogram(data2[config_dict[var]["match_var"]], bins=range[2], range=(range[0], range[1]), weights=data2[var_map["weight"]])
    hist2_err, _ = np.histogram(data2[config_dict[var]["match_var"]], bins=range[2], range=(range[0], range[1]), weights=data2[var_map["weight"]]**2)
    
    hep.histplot(hist1, bins, histtype="errorbar", yerr=np.sqrt(hist1_err), label="Mine", ax=ax1, color="red")
    hep.histplot(hist2, bins, histtype="errorbar", yerr=np.sqrt(hist2_err), label="Rui's", ax=ax1, color="blue")
    ax1.set_ylabel(f"Events/{(range[1]-range[0])/range[2]:.1}"+extract_unit(xlabel))
    ax1.legend()
    
    ratio = hist1/hist2
    ratio_err = np.sqrt(hist1_err / np.where(hist2 != 0, hist2**2, 1) + hist2_err * np.where(hist2 != 0, hist1**2 / hist2**4, 0))
    hep.histplot(ratio, bins, histtype="errorbar", yerr=ratio_err, ax=ax2, color="black")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Mine/Rui's")   
    ax2.axhline(1, color='gray', linestyle='--')
    
    plt.tight_layout()
    plt.savefig(save_path)

dataset_name = ["ggH_M125", "VBF_M125", "DYJetsToLL", "ZGToLLG", "ZG2JToG2L2J"]
filepaths = ["/eos/user/j/jiehan/root/skimmed_ntuples_run2/", "/eos/home-j/jiehan/root/skimmed_ntuples_rui/"]
for name in dataset_name:
    data_mine, data_rui = pd.DataFrame(), pd.DataFrame()
    input_branches_mine = list(config_dict.keys()) + list(var_map.keys())
    input_branches_rui = [config_dict[var]["match_var"] for var in config_dict.keys()] + [var_map[var] for var in var_map.keys()]
    for file in os.listdir(f"{filepaths[0]}{name}"):
        if ".root" in file:
            data = uproot.open(f"{filepaths[0]}{name}/{file}")["two_jet"].arrays(input_branches_mine, library="pd")
            data_mine = pd.concat([data_mine, data])
    for file in os.listdir(f"{filepaths[1]}{name}"):
        if ".root" in file:
            data = uproot.open(f"{filepaths[1]}{name}/{file}")["tree"].arrays(input_branches_rui, library="pd") 
            data_rui = pd.concat([data_rui, data])
    data_rui = data_rui.query("lly_m > 100 & lly_m < 180 & leplep_m > 80 & leplep_m < 100 & leplep_m + lly_m > 185 & y_ptmass > 0.1363636 & njet >= 2")
    for var in config_dict.keys():
        plot_2_hist(data_mine, data_rui, var, config_dict[var]["range"], config_dict[var]["title"], f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/var_com/{name}_{var}.png")
