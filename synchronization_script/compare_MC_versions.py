import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

config_dict = {
    # "Z_relpt": {"range": (0, 1), "title": r"${p_T^{ll}\cdot c}/{m_{ll\gamma}}$"},
    # "gamma_relpt": {"range": (0.1, 1), "title": r"${p_T^{\gamma}\cdot c}/{m_{ll\gamma}}$"},
    # "pt_balance_0j": {"range": (0, 1), "title": r"$Zeppenfeld_{\gamma}(0j)$"},
    # "H_mass": {"range": (100, 180), "title": r"$m_{ll\gamma}(GeV/c^{2})$"},
    # "H_ptt": {"range": (0, 160), "title": r"$p_{T_{t}}^{ll\gamma}$"},
    # "jet_1_btagDeepFlavB": {"range": (0, 0.1), "title": "j1 btag"},
    # "jet_2_btagDeepFlavB": {"range": (0, 0.15), "title": "j2 btag"},
    # "HZ_deltaRap": {"range": (-0.7, 0.7), "title": r"$\Delta y(ll,ll\gamma)$"},

    # "jet_1_mass": {"range": (0, 40), "bins": 50, "title": r"$m_{j}(GeV/c^{2})$"},
    # "jet_1_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta_{j1}$"},

    "mass_jj": {"range": (0, 1000), "bins": 50, "title": r"$m_{jj}(GeV/c^{2})$"},
    "delta_eta_jj": {"range": (0, 8), "bins": 40, "title": r"$\Delta\eta_{j}$"},
    "delta_phi_jj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{j}$"},
    "jet_2_pt": {"range": (30, 150), "bins": 50, "title": r"$p_{T_{j2}}(GeV/c)$"},
    "jet2G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(\gamma,j2)$"},

    "n_jets": {"range": (2, 6), "bins": 4, "title": r"$N_{jets}$"},

    "llphoton_hmiss_photon_dphi": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi(ll\gamma-rev, photon)$"},
    "delta_phi_zgjj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi_{ll\gamma,j}$"},
    "gamma_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta_{\gamma}$"},
    "gamma_mvaID": {"range": (0.14, 1), "bins": 43, "title": "photon MVA"},
    "jet_1_pt": {"range": (30, 330), "bins": 50, "title": r"$p_{T_{j1}}(GeV/c)$"},
    "jet1G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(\gamma,j1)$"},
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
WEIGHT = "weight_corr"
OLD_PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_rui/"
NEW_PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_rui_new/"

sig = {"sig": ["ggH_M125", "VBF_M125"], "ggH": ["ggH_M125"], "VBF": ["VBF_M125"]}
bkg = {r"Z$+\gamma$": ["ZGToLLG"], 
        "Z+Fake Photon": ["DYJetsToLL"], #, "EWKZ2J"], 
        r"VBSZ+$\gamma$": ["EWKZ2J"], #"ZG2JToG2L2J"], 
        # r"t$\bar{t}$": ["TT"], r"t$\gamma$/t$\bar{t}\gamma$": ["TTGJets", "TGJets"], "multiboson": ["WW", "WZ", "ZZ", "WWG", "WZG", "ZZG"], r"t$\bar{t}$+X": ["ttZJets", "ttWJets"]
        }
color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBSZ+$\gamma$": "#92dadd", r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", "multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2", "data": "black", "sig": "red", "ggH": "magenta", "VBF": "green"}

def convert_root_to_hist(file_dict, path, selection=None):
    mass_hist = np.zeros(80)
    error = np.zeros(BINS)
    hists = []
    for file in file_dict.values():
        type_hist = np.zeros(BINS)
        for f in file:
            for year in ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]:
                print("Reading", path+f+"/"+year+f".root:{TREE}", "...")
                if "H_mass" not in VAR:
                    try:
                        samples = uproot.open(path+f+"/"+year+f".root:{TREE}").arrays([VAR, WEIGHT, "H_mass"], library="pd")
                    except:
                        print("Fail to read", path+f+"/"+year+f".root:{TREE}")
                        continue
                else:
                    try:
                        samples = uproot.open(path+f+"/"+year+f".root:{TREE}").arrays([VAR, WEIGHT], library="pd")
                    except:
                        print("Fail to read", path+f+"/"+year+f".root:{TREE}")
                        continue
                if selection != None:
                    samples = samples.query(selection)
                hist, _ = np.histogram(samples["H_mass"], bins=80, range=[100, 180], weights=samples[WEIGHT])
                mass_hist = mass_hist + hist
                hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
                hist[0] += np.sum(samples[WEIGHT][(samples[VAR] < RMIN) & (samples[VAR] > -900)])
                hist[-1] += np.sum(samples[WEIGHT][(samples[VAR] > RMAX) & (samples[VAR] < 900)])
                type_hist = type_hist + hist
                hist, _ = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT]**2)
                hist[0] += np.sum(samples[WEIGHT][(samples[VAR] < RMIN) & (samples[VAR] > -900)])**2
                hist[-1] += np.sum(samples[WEIGHT][(samples[VAR] > RMAX) & (samples[VAR] < 900)])**2
                error = error + hist
                # set_trace()
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
    
    selection = "(H_mass>120) & (H_mass<130)"
    
    old_bkg_hists, _, bins, _ = convert_root_to_hist(bkg, OLD_PATH, selection=selection)
    new_bkg_hists, _, _, _ = convert_root_to_hist(bkg, NEW_PATH, selection=selection)
    old_sig_hists, _, _, _ = convert_root_to_hist(sig, OLD_PATH, selection=selection)
    new_sig_hists, _, _, _ = convert_root_to_hist(sig, NEW_PATH, selection=selection)

    old_bkg_yields = [np.sum(i) for i in old_bkg_hists]
    new_bkg_yields = [np.sum(i) for i in new_bkg_hists]
    old_sig_yields = [np.sum(i) for i in old_sig_hists]
    new_sig_yields = [np.sum(i) for i in new_sig_hists]

    # Normalize to unit area to compare shapes
    old_bkg_hists_norm = [h / np.sum(old_bkg_yields) if np.sum(old_bkg_yields) > 0 else h for h in old_bkg_hists]
    new_bkg_hists_norm = [h / np.sum(new_bkg_yields) if np.sum(new_bkg_yields) > 0 else h for h in new_bkg_hists]
    old_sig_hists_norm = [h / np.sum(old_sig_yields) if np.sum(old_sig_yields) > 0 else h for h in old_sig_hists]
    new_sig_hists_norm = [h / np.sum(new_sig_yields) if np.sum(new_sig_yields) > 0 else h for h in new_sig_hists]

    points = (bins[:-1] + bins[1:]) / 2

    # Create figure with two subplots, ratio 7:3
    # Background plot
    fig_bkg = plt.figure(figsize=(8, 8), dpi=200)
    gs_bkg = fig_bkg.add_gridspec(2, 1, height_ratios=[7, 3], hspace=0)
    ax1_bkg = fig_bkg.add_subplot(gs_bkg[0])
    ax2_bkg = fig_bkg.add_subplot(gs_bkg[1])

    # Upper plot for background
    ax1_bkg.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False, left=True, right=True, labelright=False)
    ax1_bkg.tick_params(axis='x', labelbottom=False)  # Remove x tick labels

    bkg_colors = [color_dict[i] for i in bkg]
    bkg_labels_new = [i+f" (new) [N={new_bkg_yields[j]:.1f}]" for j, i in enumerate(bkg.keys())]
    
    total_old_bkg_norm = np.sum(old_bkg_hists_norm, axis=0)
    old_bkg_total_yield = np.sum(old_bkg_yields)

    hep.histplot(new_bkg_hists_norm, bins, color=bkg_colors, label=bkg_labels_new, stack=True, histtype="fill", ax=ax1_bkg)
    
    # Use bar for old histogram
    ax1_bkg.bar(bins[:-1], total_old_bkg_norm, width=np.diff(bins), align='edge', facecolor='gray', alpha=0.5, label=f"Old Bkg (total) [N={old_bkg_total_yield:.1f}]")
    
    bin_width = (RMAX - RMIN) / BINS
    if bin_width < 0.01:
        y_label = f"Normalized Events/{bin_width:.4f}"
    elif bin_width < 0.1:
        y_label = f"Normalized Events/{bin_width:.3f}"
    elif bin_width < 1:
        y_label = f"Normalized Events/{bin_width:.2f}"
    else:
        y_label = f"Normalized Events/{bin_width:.1f}"
    ax1_bkg.set_ylabel(y_label, fontsize=24)
    ax1_bkg.legend(fontsize=12, ncol=2, handletextpad=0.4, columnspacing=0.5)
    ax1_bkg.grid()
    ax1_bkg.set_xlim(RMIN, RMAX)
    
    total_new_bkg = np.sum(new_bkg_hists_norm, axis=0)
    ax1_bkg.set_ylim(0, 1.3 * max(np.max(total_new_bkg), np.max(total_old_bkg_norm)))

    ax1_bkg.annotate(rf"L=137.61 fb$^{{-1}}$(13TeV)+62.32 fb$^{{-1}}$(13.6TeV)", xy=(1, 1.01), xycoords='axes fraction', fontsize=16, ha="right")
    ax1_bkg.text(-0.1, 1.01, "CMS preliminary", transform=ax1_bkg.transAxes, fontsize=16, fontweight='bold', 
             verticalalignment='bottom', horizontalalignment='left')

    # Lower plot - background ratio plot
    total_old_bkg = np.sum(old_bkg_hists_norm, axis=0)
    bkg_ratio = np.divide(total_new_bkg, total_old_bkg, out=np.ones_like(total_new_bkg), where=total_old_bkg!=0)
    
    hep.histplot(bkg_ratio, bins, color='blue', label='Bkg', histtype="step", ax=ax2_bkg, linewidth=2)
    
    ax2_bkg.set_ylabel("New/Old", fontsize=24)
    ax2_bkg.set_xlabel(XLABLE, fontsize=24)
    ax2_bkg.tick_params(axis='both', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)
    ax2_bkg.grid()
    ax2_bkg.set_xlim(RMIN, RMAX)
    ax2_bkg.set_ylim(0.5, 1.5)
    ax2_bkg.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    # ax2_bkg.legend()

    plt.tight_layout()
    if not os.path.exists(f"pic/MC_comparison/{TREE}"):
        os.makedirs(f"pic/MC_comparison/{TREE}")
    plt.savefig(f"pic/MC_comparison/{TREE}/{VAR}_bkg.png")
    plt.clf()

    # Signal plot
    fig_sig = plt.figure(figsize=(8, 8), dpi=200)
    gs_sig = fig_sig.add_gridspec(2, 1, height_ratios=[7, 3], hspace=0)
    ax1_sig = fig_sig.add_subplot(gs_sig[0])
    ax2_sig = fig_sig.add_subplot(gs_sig[1])

    # Upper plot for signal
    ax1_sig.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False, left=True, right=True, labelright=False)
    ax1_sig.tick_params(axis='x', labelbottom=False)  # Remove x tick labels

    sig_colors = [color_dict[i] for i in sig]
    sig_labels_new = [i+f" (new) [N={new_sig_yields[j]:.1f}]" for j, i in enumerate(sig.keys())]
    
    total_old_sig_norm = np.sum(old_sig_hists_norm, axis=0)
    old_sig_total_yield = np.sum(old_sig_yields)

    hep.histplot(new_sig_hists_norm, bins, color=sig_colors, label=sig_labels_new, stack=True, histtype="fill", ax=ax1_sig)
    
    # Use bar for old histogram
    ax1_sig.bar(bins[:-1], total_old_sig_norm, width=np.diff(bins), align='edge', facecolor='gray', alpha=0.8, label=f"Old Sig (total) [N={old_sig_total_yield:.1f}]")
    
    ax1_sig.set_ylabel(y_label, fontsize=24)
    ax1_sig.legend(fontsize=12, ncol=2, handletextpad=0.4, columnspacing=0.5)
    ax1_sig.grid()
    ax1_sig.set_xlim(RMIN, RMAX)
    
    total_new_sig = np.sum(new_sig_hists_norm, axis=0)
    ax1_sig.set_ylim(0, 1.3 * max(np.max(total_new_sig), np.max(total_old_sig_norm)))

    ax1_sig.annotate(rf"L=137.61 fb$^{{-1}}$(13TeV)+62.32 fb$^{{-1}}$(13.6TeV)", xy=(1, 1.01), xycoords='axes fraction', fontsize=16, ha="right")
    ax1_sig.text(-0.1, 1.01, "CMS preliminary", transform=ax1_sig.transAxes, fontsize=16, fontweight='bold', 
             verticalalignment='bottom', horizontalalignment='left')

    # Lower plot - signal ratio plot
    total_old_sig = np.sum(old_sig_hists_norm, axis=0)
    sig_ratio = np.divide(total_new_sig, total_old_sig, out=np.ones_like(total_new_sig), where=total_old_sig!=0)
    
    hep.histplot(sig_ratio, bins, color='red', label='Sig', histtype="step", ax=ax2_sig, linewidth=2)
    
    ax2_sig.set_ylabel("New/Old", fontsize=24)
    ax2_sig.set_xlabel(XLABLE, fontsize=24)
    ax2_sig.tick_params(axis='both', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)
    ax2_sig.grid()
    ax2_sig.set_xlim(RMIN, RMAX)
    ax2_sig.set_ylim(0.5, 1.5)
    ax2_sig.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    # ax2_sig.legend()

    plt.tight_layout()
    if not os.path.exists(f"pic/MC_comparison/{TREE}"):
        os.makedirs(f"pic/MC_comparison/{TREE}")
    plt.savefig(f"pic/MC_comparison/{TREE}/{VAR}_sig.png")
    plt.clf()
