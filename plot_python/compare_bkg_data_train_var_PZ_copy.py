import pandas as pd
import os
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

Tree_dict = {"ggF": "zero_to_one_jet", "VBF": "two_jet"}
# Tree_dict = {"ggF": "zero_to_one_jet"}
# Tree_dict = {"ggF": "one_jet", "VBF": "two_jet"}
name_dict_ggF = {"H_relpt": "1_Higgs_Relative_pT", "gamma_mvaID": "2_Photon_ID_MVA", "gamma_ptRelErr": "3_Photon_Energy_Resolution", "l1g_deltaR": "4_dR_Lep1_gamma", "l2g_deltaR": "5_dR_Lep2_gamma", "Z_cos_theta": "6_cosTheta", "lep_cos_theta": "7_costheta", "lep_phi": "8_phi", "Z_lead_lepton_eta": "9_eta_Lep1", "Z_sublead_lepton_eta": "10_eta_Lep2", "gamma_eta": "11_eta_gamma", "jet_1_eta": "jet_1_eta", "jet_2_eta": "jet_2_eta", "n_jets": "n_jets", "jet_1_pt":"jet_1_pt", "jet_2_pt":"jet_2_pt"}
name_dict_VBF = {"H_relpt": "1_Higgs_Relative_pT", "gamma_mvaID": "2_Photon_ID_MVA", "gamma_ptRelErr": "3_Photon_Energy_Resolution", "l1g_deltaR": "4_dR_Lep1_gamma", "l2g_deltaR": "5_dR_Lep2_gamma", "Z_cos_theta": "6_cosTheta", "lep_cos_theta": "7_costheta", "lep_phi": "8_phi", "photon_zeppenfeld": "9_gamma_Zeppenfeld", "jet1G_deltaR": "10_dR_Jet1_gamma", "jet2G_deltaR": "11_dR_Jet2_gamma", "Z_lead_lepton_eta": "12_eta_Lep1", "Z_sublead_lepton_eta": "13_eta_Lep2", "gamma_eta": "14_eta_gamma", "delta_eta_jj": "15_dEta_Jet1_Jet2", "delta_phi_jj": "16_dPhi_Jet1_Jet2", "jet_1_pt": "17_Jet1_pt", "jet_2_pt": "18_Jet2_pt", "pt_balance": "19_System_pT_Balance", "delta_phi_zgjj": "20_dPhi_jj_llgamma", "jet_1_eta": "jet_1_eta", "jet_2_eta": "jet_2_eta", "n_jets": "n_jets"}

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
        # #1
        # "H_relpt": {"range": (0, 2.5), "bins": 50, "title": r"${p_{T_{ll\gamma}}}\;/\;{m_{ll\gamma}}$"},
        # #2
        # "gamma_mvaID": {"range": (0.14, 1), "bins": 43, "title": r"$\gamma$ ID MVA"},
        # #3
        # "gamma_ptRelErr": {"range": (0.01, 0.1), "bins": 50, "title": r"$\sigma_{\gamma}\;/\;E_{\gamma}$"},
        # #4
        # "l1g_deltaR": {"range": (0.3, 4.3), "bins": 40, "title": r"$\Delta R(Lep1,\;\gamma)$"},
        # #5
        # "l2g_deltaR": {"range": (0.3, 3.3), "bins": 40, "title": r"$\Delta R(Lep2,\;\gamma)$"},
        # #6
        # "Z_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\Theta$"},
        # #7
        # "lep_cos_theta": {"range": (-1, 1), "bins": 40, "title": r"$\cos\theta$"},
        # #8
        # "lep_phi": {"range": (-3.2, 3.2), "bins": 40, "title": r"$\phi$"},
        # #9
        # "Z_lead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(Lep1)$"},
        # #10
        # "Z_sublead_lepton_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(Lep2)$"},
        # #11
        # "gamma_eta": {"range": (-2.5, 2.5), "bins": 50, "title": r"$\eta(\gamma)$"},
        #jet_1_pt
        "jet_1_pt": {"range": (30, 280), "bins": 50, "title": r"$Jet1\;p_{T}\;[GeV]$"},
        #jet eta1
        "jet_1_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta(Jet1)$"},
        #njet
        "n_jets": {"range": (0, 6), "bins": 6, "title": "Number of Jet(s)"}
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
        "delta_phi_zgjj": {"range": (0, 3.2), "bins": 40, "title": r"$\Delta\phi(jj,\;ll\gamma)$"},
        #jet eta1
        "jet_1_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta(Jet1)$"},
        #jet eta2
        "jet_2_eta": {"range": (-5, 5), "bins": 50, "title": r"$\eta(Jet2)$"},
        #njet
        "n_jets": {"range": (0, 6), "bins": 6, "title": "Number of Jet(s)"}
        }

    # inclusive, zero_jet, one_jet, zero_to_one_jet, two_jet, VH_ttH, VH, ZH, ttH_had, ttH_lep
    # TREE = "zero_to_one_jet"
    TREE = tree_name
    WEIGHT = "weight"
    PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_rui_new/"
    # PATH = "/eos/home-p/pelai/HZgamma/Root_Dataset/run3/NanoV12/"
    # PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run2/"

    data = {"data": ["Data"]}
    # bkg = {r"Z$+\gamma$": ["ZGToLLG"], "Z+Fake Photon": ["DYJetsToLL", "EWKZ2J"], r"VBSZ+$\gamma$": ["ZG2JToG2L2J"], r"t$\bar{t}$": ["TT"], r"t$\gamma$/t$\bar{t}\gamma$": ["TTGJets", "TGJets"], "multiboson": ["WW", "WZ", "ZZ"], r"t$\bar{t}$+X": ["ttZJets", "ttWJets"]}
    # bkg = {"Multiboson": ["WW", "WZ", "ZZ"], r"t$\bar{t}$+X": ["ttZJets", "ttWJets"], r"VBS Z+$\gamma$": ["ZG2JToG2L2J"], r"t$\gamma$/t$\bar{t}\gamma$": ["TTGJets", "TGJets"], r"t$\bar{t}$": ["TT"], r"Z$+\gamma$": ["DYGto2LG"], "Z+Fake Photon": ["DYJetsToLL"]}
    # bkg = {"Multiboson": ["WW", "WZ", "ZZ"], r"t$\bar{t}$": ["TT", "TTtoLNu2Q"], r"VBS Z+$\gamma$": ["ZG2JToG2L2J"], r"Z$+\gamma$": ["DYGto2LG"], "Z+Fake Photon": ["DYJetsToLL"]}
    # bkg = {r"Z$+\gamma$": ["DYGto2LG"], "Z+Fake Photon": ["DYJetsToLL"]}
    # Rui's Sample
    bkg = {r"Z$+\gamma$": ["ZGToLLG"], "Z+Fake Photon": ["DYJetsToLL","EWKZ2J"]}
    color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBS Z+$\gamma$": "#92dadd", r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", "Multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2"}
    
    years = {"run2":["2016preVFP", "2016postVFP", "2017", "2018"],"run3":["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]}
    lumis = {"run2": 138, "run3": 62}  # fb^-1

    # 新增：動態建立 selection，若變數含 "jet" 則過濾 jet_1_pt != -999
    # 改進：依變數精準套用 jet 條件
    def build_selection(var: str) -> str:
        """
        根據欲畫的變量動態產生 selection 字串：
        - 非 H_mass 變量使用 H 質量 sideband
        - jet_1_* 相關變量要求 jet_1_pt!=-999 或 jet_1_eta>-5
        - jet_2_* 相關變量要求 jet_2_pt!=-999 或 jet_2_eta>-5
        - 需要兩個 jet 的組合變量同時要求兩個 jet 存在
        """
        parts = []
        # 重要：用雙層括號把 OR 群組起來，避免與後面的 & 優先序錯誤
        if var != "H_mass":
            parts.append("((H_mass<120) | (H_mass>130))")
        else:
            # 個別 jet_1 變數
            if "jet_1_pt" in var:
                parts.append("(jet_1_pt!=-999)")
            elif "jet_1_eta" in var:
                parts.append("(jet_1_eta!=-999)")

        return " & ".join(parts)

    def convert_root_to_hist(file_dict, selection="", years_to_run=None):
        mass_hist = np.zeros(80)
        error = np.zeros(BINS)
        hists = []
        bins = np.linspace(RMIN, RMAX, BINS + 1)
        # 若未提供，預設用 run3
        if years_to_run is None:
            years_to_run = years["run3"]
        for file in file_dict.values():
            type_hist = np.zeros(BINS)
            for f in file:
                for year in years_to_run:
                    print("Reading", PATH + f  + "/" + year + ".root:" + TREE, "...")
                    if "H_mass" not in VAR:
                        try:
                            samples = uproot.open(PATH + f + "/" + year + ".root:" + TREE).arrays([VAR, WEIGHT, "H_mass"], library="pd")
                        except:
                            continue
                    else:
                        try:
                            samples = uproot.open(PATH + f + "/" + year + ".root:" + TREE).arrays([VAR, WEIGHT], library="pd")
                        except:
                            continue
                    if selection:
                        samples = samples.query(selection)
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

        # 針對 run2 與 run3 各跑一遍
        for run_key in ["run2", "run3"]:
            selection_str = build_selection(VAR)
            # 簡單列印一次 selection 以便確認是否生效
            # print(f"[{run_key}] VAR={VAR} selection => {selection_str}")

            hist2, hist2_err, bin_edges, mass2_hist = convert_root_to_hist(
                bkg, selection=selection_str, years_to_run=years[run_key]
            )
            hist1, hist1_err, _, mass1_hist = convert_root_to_hist(
                data, selection=selection_str, years_to_run=years[run_key]
            )

            #--------------------------------------------
            s1 = np.sum(mass1_hist)
            s2 = np.sum(mass2_hist)
            if s2 > 0:
                sf = s1 / s2
            else:
                sf = 1.0
                print("[WARN] mass2_hist sum is zero; set sf=1.0 and skip rescaling.")
            #--------------------------------------------
            hist2 = [h*sf for h in hist2]

            points = (bin_edges[:-1] + bin_edges[1:]) / 2

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=200, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
            ax1.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False, left=True, right=True, labelright=False)
            ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, left=True, right=True, labelright=False)

            # plt.figure(figsize=(8, 8), dpi=200)
            # plt.plot(bins, hist1/sum(hist1), "o-", label="Signal")
            if len(hist1) == 1:
                # ax1.errorbar(points, hist1[0], yerr=hist1_err, fmt="o", label=f"Data[N={np.sum(hist1):.1f}]", color="black", markersize=7, linewidth=3)
                # ax1.errorbar(points, hist1[0], yerr=hist1_err, xerr=(bins[1:] - bins[:-1]) / 2, fmt="o", label=f"Data[N={np.sum(hist1):.1f}]", color="black", markersize=7, linewidth=3)
                ax1.errorbar(points, hist1[0], yerr=hist1_err, xerr=(bin_edges[1:] - bin_edges[:-1]) / 2, fmt="o", label=f"Data [N={np.sum(hist1):.0f}]", color="black", markersize=7, linewidth=3)
            else:
                raise ValueError("Data should be only 1")
            colors = [color_dict[i] for i in bkg]
            labels = [i+f" [N={np.sum(hist2[j]):.0f}]" for j, i in enumerate(bkg)]
            # print("hist2 shape:", np.array(hist2).shape)
            # print("bin_edges shape:", bin_edges.shape)
            hep.histplot(hist2, bin_edges, color=colors, label=labels, stack=True, histtype="fill", ax=ax1, zorder=1)
            
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

            # 若 hist2 / hist2_err 目前是 list，需要先轉成 np.array
            hist2_array = np.array(hist2)          # (N_bkgs, N_bins)
            hist2_err_array = np.array(hist2_err)  # (N_bkgs, N_bins)
            bkg_total = np.sum(hist2_array, axis=0)                   
            bkg_total_err = hist2_err_array  # 等于 hist2_err, 即已经是总误差

            # ----------------------------------------------------------
            # 准备 fill_between 不确定度带: 需要 N+1 长度
            upper = np.concatenate([
                bkg_total + bkg_total_err, 
                [bkg_total[-1] + bkg_total_err[-1]]
            ])  # shape: (N_bins+1,)
            lower = np.concatenate([
                bkg_total - bkg_total_err, 
                [bkg_total[-1] - bkg_total_err[-1]]
            ])  # shape: (N_bins+1,)

            # Statistic Uncertainty: Hatch Style
            # ax1.fill_between(bin_edges, lower, upper, step="post", facecolor="none", edgecolor="black", linewidth=0, hatch="///", alpha=1.0, label="Stat. Uncert.")
            ax1.fill_between(bin_edges, lower, upper, step="post", facecolor="none", edgecolor="black", linewidth=1.0, hatch="///", alpha=1.0, label="Stat. Uncert.")
            
            # Statistic Uncertainty: Gray Box Style 
            # ax1.fill_between(bin_edges, lower, upper, step="post", facecolor="gray", edgecolor="none", linewidth=0, alpha=0.4, label="Stat. Uncert.")

            # 画出总背景的阶梯线
            hep.histplot(hist2, bin_edges, color="black", stack=True, histtype="step", lw=0.8, ax=ax1, zorder=2)

            bin_width = (RMAX - RMIN) / BINS
            if bin_width < 0.01:
                if i == "jet_1_pt" or i == "jet_2_pt":
                    y_label = f"Events / {bin_width:.3f} GeV"
                else: 
                    y_label = f"Events / {bin_width:.3f}"
            elif bin_width < 0.1:
                if i == "jet_1_pt" or i == "jet_2_pt":
                    y_label = f"Events / {bin_width:.2f} GeV"
                else: 
                    y_label = f"Events / {bin_width:.2f}"
            elif bin_width < 1:
                if i == "jet_1_pt" or i == "jet_2_pt":
                    y_label = f"Events / {bin_width:.1f} GeV"
                else: 
                    y_label = f"Events / {bin_width:.1f}"
            else:
                if i == "jet_1_pt" or i == "jet_2_pt":
                    y_label = f"Events / {bin_width:.0f} GeV"
                else: 
                    y_label = f"Events / {bin_width:.0f}"
            
            ax1.set_ylabel(y_label, fontsize=24)

            # ggF
            if primary_process_name == "ggF":
                # ax1.legend(fontsize=16, ncol=2, handletextpad=0.4, columnspacing=0.5, bbox_to_anchor=(1.02, 1.01),loc="upper right", frameon=False, handlelength=0.7, handleheight=1.2, labelspacing=0.1)
                ax1.legend(fontsize=16, ncol=2, loc="best", frameon=False, handletextpad=0.4, columnspacing=0.5, handlelength=0.7, handleheight=1.2, labelspacing=0.1)

            # VBF
            elif primary_process_name == "VBF":
                # ax1.legend(fontsize=16, ncol=2, handletextpad=0.4, columnspacing=0.5, bbox_to_anchor=(1.02, 1.01),loc="upper right", frameon=False, handlelength=0.7, handleheight=1.2, labelspacing=0.1)
                ax1.legend(fontsize=17, ncol=2, loc="best", frameon=False, handletextpad=0.4, columnspacing=0.5, handlelength=0.7, handleheight=1.2, labelspacing=0.1)

            # ax1.grid()
            ax1.set_xlim(RMIN, RMAX)
            #---------------------------------
            # ax1.set_ylim(0, 1.7*max(np.sum(hist1, axis=0).max(), np.sum(hist2, axis=0).max()))
            ymax = max(np.sum(hist1, axis=0).max(),np.sum(hist2, axis=0).max())

            if ymax > 0:
                ax1.set_ylim(0, 1.7 * ymax)
            else:
                # 给一个固定的小范围，或者干脆不设 ylim
                ax1.set_ylim(0, 1.0)
                print("[WARN] Both hist1 and hist2 are empty; using default y-limits.")
            #---------------------------------
            # ax1.annotate(rf"L=137.2 fb$^{{-1}}$, SF={sf:.2f}", xy=(1, 1.01), xycoords='axes fraction', fontsize=18, ha="right")
            ax1.annotate(rf"L={lumis[run_key]} fb$^{{-1}}$, SF={sf:.2f}", xy=(1, 1.01), xycoords='axes fraction', fontsize=18, ha="right")
            ax1.annotate(r"$\mathbf{CMS}\ \text{Preliminary}$",xy=(0, 1.01),xycoords='axes fraction',fontsize=22,ha="left")

            data_sum, bkg_sum = np.where(np.sum(hist1, axis=0) == 0, 1e-8, np.sum(hist1, axis=0)), np.where(np.sum(hist2, axis=0) == 0, 1e-8, np.sum(hist2, axis=0))
            ratio = hist1[0]/bkg_sum
            ratio_err = hist1_err/data_sum
            err = hist2_err/bkg_sum

            mask = (ratio > 0) & (ratio < 2)
            # ax2.errorbar(points[mask], ratio[mask], yerr=ratio_err[mask], xerr=((bins[1:] - bins[:-1]) / 2)[mask], fmt="o", color="black", markersize=7, linewidth=3)
            ax2.errorbar(points[mask], ratio[mask], yerr=ratio_err[mask], xerr=((bin_edges[1:] - bin_edges[:-1]) / 2)[mask], fmt="o", color="black", markersize=7, linewidth=3)
            # ax2.bar(points, 2 * err, bottom=1 - err, color="#ffedcd", width=(RMAX-RMIN)/BINS, label="Bkg Uncertainty")
            ax2.bar(points, 2 * err, bottom=1 - err, color="#CCCCCC", width=(RMAX-RMIN)/BINS, label="Bkg Uncertainty")
            # ax2.bar(points, 2 * err, bottom=1 - err, color="none", hatch="///", width=(RMAX-RMIN)/BINS, label="Bkg Uncertainty", edgecolor="black")
            ax2.set_xlim(RMIN, RMAX)
            ax2.set_ylim(0, 2)
            ax2.set_ylabel(r"$\frac{\text{Data}}{\text{MC}}$", fontsize=28, rotation=90)
            ax2.yaxis.set_label_coords(-0.08, 0.72)
            ax2.axhline(1, color="black", linestyle="--")
            ax2.set_yticks(np.arange(0, 2.1, 0.5))
            ax2.set_yticklabels([f"{int(i)}" if i == 1 else f"{i}" if i in [0.5, 1.5] else "" for i in np.arange(0, 2.1, 0.5)])

            ax2.set_xlabel(XLABLE, fontsize=24)

            plt.tight_layout()

            # 依 run_key 輸出到 run2 或 run3 目錄
            out_dir = f"/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/plot_python/pic/{run_key}/dataVbkg_{primary_process_name}"
            if os.path.exists(out_dir) == False:
                os.makedirs(out_dir)
            if primary_process_name == "ggF":
                plt.savefig(f"{out_dir}/{name_dict_ggF[VAR]}.pdf")
            elif primary_process_name == "VBF":
                plt.savefig(f"{out_dir}/{name_dict_VBF[VAR]}.pdf")

            plt.clf()