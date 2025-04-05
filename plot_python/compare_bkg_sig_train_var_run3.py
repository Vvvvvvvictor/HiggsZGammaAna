import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

Tree_dict = {"ggF": "zero_to_one_jet", "VBF": "two_jet"}
name_dict_ggF = {"H_relpt": "1_Higgs_Relative_pT", "gamma_mvaID": "2_Photon_ID_MVA", "gamma_ptRelErr": "3_Photon_Energy_Resolution", "l1g_deltaR": "4_dR_Lep1_gamma", "l2g_deltaR": "5_dR_Lep2_gamma", "Z_cos_theta": "6_cosTheta", "lep_cos_theta": "7_costheta", "lep_phi": "8_phi", "Z_lead_lepton_eta": "9_eta_Lep1", "Z_sublead_lepton_eta": "10_eta_Lep2", "gamma_eta": "11_eta_gamma"}
name_dict_VBF = {"H_relpt": "1_Higgs_Relative_pT", "gamma_mvaID": "2_Photon_ID_MVA", "gamma_ptRelErr": "3_Photon_Energy_Resolution", "l1g_deltaR": "4_dR_Lep1_gamma", "l2g_deltaR": "5_dR_Lep2_gamma", "Z_cos_theta": "6_cosTheta", "lep_cos_theta": "7_costheta", "lep_phi": "8_phi", "photon_zeppenfeld": "9_gamma_Zeppenfeld", "jet1G_deltaR": "10_dR_Jet1_gamma", "jet2G_deltaR": "11_dR_Jet2_gamma", "Z_lead_lepton_eta": "12_eta_Lep1", "Z_sublead_lepton_eta": "13_eta_Lep2", "gamma_eta": "14_eta_gamma", "delta_eta_jj": "15_dEta_Jet1_Jet2", "delta_phi_jj": "16_dPhi_Jet1_Jet2", "jet_1_pt": "17_Jet1_pt", "jet_2_pt": "18_Jet2_pt", "pt_balance": "19_System_pT_Balance", "delta_phi_zgjj": "20_dPhi_jj_llgamma"}

for primary_process_name, tree_name in Tree_dict.items(): 
        # "mass_jj": {"range": (0, 300), "title": r"$m_{jj}(GeV/c^{2})$"},
        # "Z_relpt": {"range": (0, 1), "title": r"${p_T^{ll}\cdot c}/{m_{ll\gamma}}$"},
        # "gamma_relpt": {"range": (0.1, 1), "title": r"${p_T^{\gamma}\cdot c}/{m_{ll\gamma}}$"},
        # "pt_balance_0j": {"range": (0, 1), "title": r"$Zeppenfeld_{\gamma}(0j)$"},
        # "H_mass": {"range": (100, 180), "title": r"$m_{ll\gamma}(GeV/c^{2})$"},
        # "HZ_deltaRap": {"range": (-0.7, 0.7), "title": r"$\Delta y(ll,ll\gamma)$"}
        # "H_ptt": {"range": (0, 160), "title": r"$p_{T_{t}}^{ll\gamma}$"},
        # "jet_1_btagDeepFlavB": {"range": (0, 0.1), "title": "j1 btag"},
        # "jet_2_btagDeepFlavB": {"range": (0, 0.15), "title": "j2 btag"},

    if primary_process_name == "ggF":    

        config_dict = {
        # ggF Training Variables
        #1
        "H_relpt": {"range": (0, 2.5), "bins": 50, "title": r"${p_{T_{ll\gamma}}}\;/\;{m_{ll\gamma}}$"},
        #2
        "gamma_mvaID": {"range": (0.14, 1), "bins": 43, "title": r"$\gamma$ ID MVA"},
        #3
        "gamma_ptRelErr": {"range": (0.01, 0.1), "bins": 50, "title": r"$\sigma_{\gamma}\;/\;E_{\gamma}$"},
        #4
        "l1g_deltaR": {"range": (0.3, 4.3), "bins": 40, "title": r"$\Delta R(Lep1,\;\gamma)$"},
        #5
        "l2g_deltaR": {"range": (0.3, 3.3), "bins": 40, "title": r"$\Delta R(Lep2,\;\gamma)$"},
        #6
        "Z_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\Theta$"},
        #7
        "lep_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\theta$"},
        #8
        "lep_phi": {"range": (-3.2, 3.2), "bins": 40, "title": r"$\phi$"},
        #9
        "Z_lead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(Lep1)$"},
        #10
        "Z_sublead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(Lep2)$"},
        #11
        "gamma_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(\gamma)$"}
        }

    elif primary_process_name == "VBF":     
            
        config_dict = {
        # VBF Training Variables
        #1
        "H_relpt": {"range": (0, 2.5), "bins": 50, "title": r"${p_{T_{ll\gamma}}}\;/\;{m_{ll\gamma}}$"},
        #2
        "gamma_mvaID": {"range": (0.14, 1), "bins": 43, "title": r"$\gamma$ ID MVA"},
        #3
        "gamma_ptRelErr": {"range": (0.01, 0.1), "bins": 50, "title": r"$\sigma_{\gamma}\;/\;E_{\gamma}$"},
        #4
        "l1g_deltaR": {"range": (0.3, 4.3), "bins": 40, "title": r"$\Delta R(Lep1,\;\gamma)$"},
        #5
        "l2g_deltaR": {"range": (0.3, 3.3), "bins": 40, "title": r"$\Delta R(Lep2,\;\gamma)$"},
        #6
        "Z_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\Theta$"},
        #7
        "lep_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\theta$"},
        #8
        "lep_phi": {"range": (-3.2, 3.2), "bins": 40, "title": r"$\phi$"},
        #9
        "photon_zeppenfeld": {"range": (0, 5), "bins": 50, "title": "$\gamma$ Zeppenfeld"},
        #10
        "jet1G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(Jet1,\;\gamma)$"},
        #11
        "jet2G_deltaR": {"range": (0.4, 6.4), "bins": 40, "title": r"$\Delta R(Jet2,\;\gamma)$"},
        #12
        "Z_lead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(Lep1)$"},
        #13
        "Z_sublead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(Lep2)$"},
        #14
        "gamma_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(\gamma)$"},
        #15
        "delta_eta_jj": {"range": (0, 8), "bins": 40, "title": r"$\Delta\eta(Jet1,\;Jet2)$"},
        #16
        "delta_phi_jj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi(Jet1,\;Jet2)$"},
        #17
        "jet_1_pt": {"range": (30, 280), "bins": 50, "title": r"$Jet1\;p_{T}\;[GeV]$"},
        #18
        "jet_2_pt": {"range": (30, 130), "bins": 50, "title": r"$Jet2\;p_{T}\;[GeV]$"},
        #19
        "pt_balance": {"range": (0, 1), "bins": 50, "title": r"System $p_{T}$ Balance"},
        #20
        "delta_phi_zgjj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi(jj,\;ll\gamma)$"}
        }


    TREE = tree_name
    WEIGHT = "weight"
    PATH = "/eos/home-p/pelai/HZgamma/Root_Dataset/run3/NanoV12/"

    sig = {"sig": ["ggH_M125", "VBF_M125"], "ggH": ["ggH_M125"], "VBF": ["VBF_M125"]}
    bkg = {"Multiboson": ["WW", "WZ", "ZZ"], r"t$\bar{t}$": ["TT", "TTtoLNu2Q"], r"VBS Z+$\gamma$": ["ZG2JToG2L2J"], r"Z$+\gamma$": ["DYGto2LG"], "Z+Fake Photon": ["DYJetsToLL"]}
    color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBS Z+$\gamma$": "#92dadd", r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", "Multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2","ggH": "#FF3B4F", "VBF": "#4EFFEF"}

    def convert_root_to_hist(file_dict, selection=None):
        mass_hist = np.zeros(80)
        error = np.zeros(BINS)
        hists = []
        for file in file_dict.values():
            type_hist = np.zeros(BINS)
            for f in file:
                for year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]:
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
                    hist[0] += np.sum(samples[WEIGHT][samples[VAR] < RMIN]**2)
                    hist[-1] += np.sum(samples[WEIGHT][samples[VAR] > RMAX]**2)
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

        hist1, hist1_err, _, mass1_hist = convert_root_to_hist(sig)
        hist2, hist2_err, bins, mass2_hist = convert_root_to_hist(bkg)

        sig_yields = [np.sum(i) for i in hist1]
        bkg_yields = [np.sum(i) for i in hist2]
        # Normalize to Sum of Two Signal MC 
        # hist1 = [i/sig_yields[0] for i in hist1]
        # Normalize to Each Signal 
        hist1 = [i / sig_yields[idx] if sig_yields[idx] != 0 else 0 for idx, i in enumerate(hist1)]
        hist2 = [i/np.sum(bkg_yields) for i in hist2]

        points = (bins[:-1] + bins[1:]) / 2

        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8), dpi=200)
        ax1.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)

        colors = [color_dict[i] for i in bkg]
        labels = [i+f" [N={bkg_yields[j]:.0f}]" for j, i in enumerate(bkg.keys())]
        hep.histplot(hist2, bins, color=colors, label=labels, stack=True, histtype="fill", ax=ax1)
        # ax1.errorbar(points, hist1[0], yerr=hist1_err, xerr=(bins[1:] - bins[:-1]) / 2, fmt="o", label=f"Sig [N={sig_yields[0]:.1f}]", color="black", markersize=7, linewidth=3)
        hep.histplot(hist2, bins, color="black", stack=True, histtype="step", lw=0.5, ax=ax1, zorder=1)
        for j, i in enumerate(hist1[1:]):
            hep.histplot(i, bins, color=color_dict[list(sig.keys())[j+1]], label=list(sig.keys())[j+1]+f" [N={sig_yields[j+1]:.0f}]", histtype="step", ax=ax1, linewidth=2.8, linestyle="-")
        
        bin_width = (RMAX - RMIN) / BINS
        if bin_width < 0.01:
            if VAR == "jet_1_pt" or VAR == "jet_2_pt":
                y_label = f"Events / {bin_width:.3f} GeV"
            else: 
                y_label = f"Events / {bin_width:.3f}"
        elif bin_width < 0.1:
            if VAR == "jet_1_pt" or VAR == "jet_2_pt":
                y_label = f"Events / {bin_width:.2f} GeV"
            else: 
                y_label = f"Events / {bin_width:.2f}"
        elif bin_width < 1:
            if VAR == "jet_1_pt" or VAR == "jet_2_pt":
                y_label = f"Events / {bin_width:.1f} GeV"
            else: 
                y_label = f"Events / {bin_width:.1f}"
        else:
            if VAR == "jet_1_pt" or VAR == "jet_2_pt":
                y_label = f"Events / {bin_width:.0f} GeV"
            else: 
                y_label = f"Events / {bin_width:.0f}"
        
        ax1.set_ylabel(y_label, fontsize=24)

        # ggF
        if primary_process_name == "ggF":
            # ax1.legend(fontsize=16, ncol=2, handletextpad=0.4, columnspacing=0.5, bbox_to_anchor=(1.02, 1.01),loc="upper right", frameon=False, handlelength=0.7, handleheight=1.2, labelspacing=0.1)
            ax1.legend(fontsize=17, ncol=2, loc="best", frameon=False, handletextpad=0.4, columnspacing=0.5, handlelength=0.7, handleheight=1.2, labelspacing=0.1)

        # VBF
        elif primary_process_name == "VBF":
            # ax1.legend(fontsize=16, ncol=2, handletextpad=0.4, columnspacing=0.5, bbox_to_anchor=(1.02, 1.01),loc="upper right", frameon=False, handlelength=0.7, handleheight=1.2, labelspacing=0.1)
            ax1.legend(fontsize=17, ncol=2, loc="best", frameon=False, handletextpad=0.4, columnspacing=0.5, handlelength=0.7, handleheight=1.2, labelspacing=0.1)

        # ax1.grid()
        ax1.set_xlim(RMIN, RMAX)
        ax1.set_xlabel(XLABLE, fontsize=24)
        # hist1 = Signal 
        # hist2 = Bkg
        max_bin_value = max(np.max(hist1, axis=1).max(), np.sum(hist2, axis=0).max())
        ax1.set_ylim(0, 1.5 * max_bin_value)

        ax1.annotate(rf"L=62.3 fb$^{{-1}}$", xy=(1, 1.01), xycoords='axes fraction', fontsize=16, ha="right")
        ax1.annotate(r"$\mathbf{CMS}\ \text{Preliminary}$",xy=(0, 1.01),xycoords='axes fraction',fontsize=22,ha="left")
        
        plt.tight_layout()
        if os.path.exists(f"pic/run3/sigVbkg_{primary_process_name}") == False:
            os.makedirs(f"pic/run3/sigVbkg_{primary_process_name}")
        
        if primary_process_name == "ggF":
            plt.savefig(f"pic/run3/sigVbkg_{primary_process_name}/{name_dict_ggF[VAR]}.pdf")
        elif primary_process_name == "VBF":
            plt.savefig(f"pic/run3/sigVbkg_{primary_process_name}/{name_dict_VBF[VAR]}.pdf")
        
        plt.close(fig)  # 关闭当前 figure
        plt.clf()