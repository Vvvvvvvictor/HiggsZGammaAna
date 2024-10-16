import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace

var_names = [
    # "mass_jj",
    # "Z_relpt",
    # "gamma_relpt",
    # "pt_balance_0j",
    # "H_mass", 
    # "delta_eta_jj",
    # "delta_phi_jj",
    # "delta_phi_zgjj",
    # "gamma_eta",
    # "gamma_mvaID",
    # "H_ptt",
    # "jet_1_pt",
    # "jet_2_pt",
    # "jet1G_deltaR",
    # "jet2G_deltaR",
    "jet_1_eta"#,
    # "jet_2_eta",
    # "jet_1_btagDeepFlavB",
    # "jet_2_btagDeepFlavB",
    # "l1g_deltaR",
    # "l2g_deltaR",
    # "lep_cos_theta",
    # "lep_phi",
    # "photon_zeppenfeld",
    # "pt_balance",
    # "Z_cos_theta",
    # "Z_lead_lepton_eta",
    # "Z_sublead_lepton_eta",
    # "H_relpt",
    # "HZ_deltaRap"
    ]
x_ranges = [
    # (0, 300),
    # (0, 1),
    # (15/110, 0.6),
    # (0, 1),
    # (100, 180),
    # (0, 9),
    # (0, 3.2),
    # (0, 3.2),
    # (-2.5, 2.5),
    # (0.14, 1),
    # (0, 160),
    # (30, 330),
    # (30, 150),
    # (0.4, 6.4),
    # (0.4, 6.4),
    (-4.7, 4.7)#,
    # (-4.7, 4.7),
    # (0, 0.1),
    # (0, 0.15),
    # (0.4, 4.8),
    # (0.4, 3.4),
    # (-1, 1),
    # (-3.2, 3.2),
    # (0, 5),
    # (0, 1),
    # (-1, 1),
    # (-2.5, 2.5),
    # (-2.5, 2.5),
    # (0, 3),
    # (-0.7, 0.7)
    ]
nums_bin = [
    # 50,
    # 50,
    # 50,
    # 50,
    # 80,
    # 45,
    # 80,
    # 80,
    # 50,
    # 86,
    # 80,
    # 100,
    # 60,
    # 60,
    60#,
    # 94,
    # 94,
    # 50,
    # 60,
    # 80,
    # 60,
    # 50,
    # 64,
    # 50,
    # 50,
    # 50,
    # 50,
    # 50,
    # 50,
    # 70
]
x_titles = [
    # r"$m_{jj}(GeV/c^{2})$",
    # r"$\frac{p_T^{ll}}{p_{T_{Z\gamma}}}$",
    # r"$\frac{p_T^{\gamma}}{p_{T_{Z\gamma}}}$",
    # r"$Zeppenfeld_{\gamma}(0j)$",
    # r"$m_{ll\gamma}(GeV/c^{2})$", 
    # r"$\Delta\eta_{jj}$",
    # r"$\Delta\phi_{jj}$",
    # r"$\Delta\phi_{z\gamma,jj}$",
    # r"$\eta_{\gamma}$",
    # "photon MVA",
    # r"$p_{T_{t}}^{Z\gamma}$",
    # r"$p_{T_{j1}}(GeV/c)$",
    # r"$p_{T_{j2}}(GeV/c)$",
    # r"$\Delta R(\gamma,j1)$",
    # r"$\Delta R(\gamma,j2)$",
    r"$\eta_{j1}$"#,
    # r"$\eta_{j2}$",
    # "j1 btag",
    # "j2 btag",
    # r"max($\Delta R(l,\gamma)$)",
    # r"min($\Delta R(l,\gamma)$)",
    # r"$\cos\theta$",
    # r"$\phi$",
    # "Zeppenfeld $\gamma$",
    # "system balance",
    # r"$\cos\Theta$",
    # r"$\eta_{l1}$",
    # r"$\eta_{l2}$",
    # r"$\frac{p_{T_{ll\gamma}}\cdot c}{m_{ll\gamma}}$",
    # r"$\Delta y(ll,ll\gamma)$"
]

# VAR = "n_jets"
# XLABLE = r"$n_{j}$"
# BINS = 6
# RMIN = 0
# RMAX = 6

CHANNEL = "inclusive"
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run3/"

bkg_mc_file_list = [
    # ["ZGToLLG"],
    ["DYGto2LG_10to50", "DYGto2LG_50to100"],
    ["DYJetsToLL"],
    ["ZG2JToG2L2J"],
    # ["EWKZ2J"],
    # ["TT"],
    # ["TTGJets"],
    # ["ttWJets", "ttZJets"],
    # ["WW", "WZ", "ZZ"]
]
data_file_list = [
    ["Data"]
]

data_names = ["Data"]
bkg_mc_names = ["SM ZG", "DYJets", "EWK ZG"]#, "EWK Z+Jets", "TT", "TTG+Jets", "TTVJets", "Diboson"]
# "SM ZG", "DYJets", "EWK ZG", "EWK Z+Jets", "TT", "TTG+Jets", "TTVJets", "Diboson"


for i in range(len(var_names)):
    VAR = var_names[i]
    XLABLE = x_titles[i]
    RMIN, RMAX = x_ranges[i]
    BINS = nums_bin[i]
    
    bkg_mc_hist_list, data_hist_list, bkg_mc_hist_list_minus, data_hist_list_minus = [], [], [], []

    def convert_parquet_to_hist(file_list, hist_list):
        for data in file_list:
            samples = pd.read_parquet(data)
            # if "weight_central_no_lumi" in samples.keys():
            #     print(data, ":", samples["weight_central"][0]/samples["weight_central_no_lumi"][0], sep="\n")
            # if "DYJ" in data:
            #     # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
            #     # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
            #     samples = samples[samples["n_iso_photons"]==0]
            #     # samples[WEIGHT]*=3
            # if "DYG" in data:
            # #     # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
            # #     # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
            #     samples = samples[samples["n_iso_photons"]>0]
            # #     # samples[WEIGHT]*=3
            if "H_mass" in VAR and "data" in data:
                samples = samples[(samples["H_mass"]<122) | (samples["H_mass"]>128)]
            hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
            hist_list.append(hist)

    def convert_root_to_hist(file_list, hist_list):
        years = ["2022preEE", "2022postEE"] # ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
        error_hist = np.zeros(BINS)
        for type in file_list:
            type_hist = np.zeros(BINS)
            for data in type:
                for year in years:
                    try:
                        print("Reading", PATH+data+f"/{year}.root:{CHANNEL}", "...")
                    except:
                        set_trace()
                    branches = np.array([VAR, "weight", "jet_1_eta", "jet_1_pt", "H_mass"])
                    branches = np.unique(branches)
                    samples = uproot.open(PATH+data+f"/{year}.root:{CHANNEL}").arrays(branches, library="pd")
                    # veto = ((abs(samples["jet_1_eta"])>2.5) & (abs(samples["jet_1_eta"])<4.7)) | ((abs(samples["jet_2_eta"])>2.5) & (abs(samples["jet_2_eta"])<4.7)) | ((abs(samples["jet_3_eta"])>2.5) & (abs(samples["jet_3_eta"])<4.7)) | ((abs(samples["jet_4_eta"])>2.5) & (abs(samples["jet_4_eta"])<4.7))
                    # samples = samples[~veto]
                    if ("H_mass" in VAR and "Data" in data) or "H_mass" not in VAR:
                        samples = samples[(samples["H_mass"]<120) | (samples["H_mass"]>130)]
                    samples = samples.query("(jet_1_eta<2.5 & jet_1_eta>-2.5) | (jet_1_pt > 80)")
                    print(samples["weight"].sum())
                    hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
                    type_hist = type_hist + hist
                    hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT]**2)
                    error_hist = error_hist + hist
            hist_list.append(type_hist)
        return error_hist
    
    # convert_parquet_to_hist(data_list, data_hist_list)
    # convert_parquet_to_hist(bkg_mc_list, bkg_mc_hist_list)

    # PATH = "/eos/user/j/jiehan/root/skimmed_ntuples/"
    data_error = convert_root_to_hist(data_file_list, data_hist_list)
    bkg_error = convert_root_to_hist(bkg_mc_file_list, bkg_mc_hist_list)
    # PATH = "/eos/user/j/jiehan/root/skimmed_ntuples_run3/"
    # data_error_minus = convert_root_to_hist(data_file_list, data_hist_list_minus)
    # bkg_error_minus = convert_root_to_hist(bkg_mc_file_list, bkg_mc_hist_list_minus)

    hep.style.use("CMS")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), dpi=400, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})

    bins = np.linspace(RMIN, RMAX, BINS+1)
    # data_hist_list = [j-i for i, j in zip(data_hist_list, data_hist_list_minus)]
    # data_error = np.sqrt(data_error + data_error_minus)
    data_yield = np.sum(np.array(data_hist_list), axis=1)
    hep.histplot(
        data_hist_list, 
        yerr=np.sqrt(data_error),
        bins=bins, 
        histtype='errorbar', 
        color='black',
        ax=ax1,
        label=[f"{i}({j:.2f})" for i, j in zip(data_names, data_yield)]
    )
    # bkg_mc_hist_list = [j-i for i, j in zip(bkg_mc_hist_list, bkg_mc_hist_list_minus)]
    # bkg_error = np.sqrt(bkg_error + bkg_error_minus)
    bkg_yield = np.sum(np.array(bkg_mc_hist_list), axis=1)
    hep.histplot(
        bkg_mc_hist_list, 
        yerr=np.sqrt(bkg_error),
        bins=bins, 
        histtype='fill', 
        color=['b', 'g'],
        ax=ax1,
        label=[f"{i}({j:.2f})" for i, j in zip(bkg_mc_names, bkg_yield)],
        stack=True
    )

    # Adjust axis for the main plot
    ax1.set_xlabel("")
    ax1.set_ylabel("Events / {:.2f} {}".format((RMAX-RMIN)/BINS, XLABLE.split("(")[-1].split(")")[0] if "GeV" in XLABLE else ""), fontsize=25)
    ax1.set_xlim(RMIN, RMAX)
    ax1.set_ylim(0, 1.5*max(max(np.sum(np.array(data_hist_list), axis=0)+np.sqrt(data_error)), max(np.sum(np.array(bkg_mc_hist_list), axis=0)+np.sqrt(bkg_error))))
    ax1.set_xticklabels([])
    ax1.legend(loc="upper right", ncol=2)
    ax1.title.set_text("Data vs MC(Run2)")

    # Create ratio plot
    ratio = np.sum(np.array(data_hist_list), axis=0) / np.sum(np.array(bkg_mc_hist_list), axis=0)
    ratio_error = np.sqrt(data_error) / np.sum(np.array(bkg_mc_hist_list), axis=0)

    ax2.errorbar(bins[:-1] + 0.5 * (bins[1] - bins[0]), ratio, yerr=ratio_error, fmt='o', color='black')
    ax2.axhline(1.0, color='gray', linestyle='--')

    # Adjust axis for the ratio plot
    ax2.set_xlabel("{}".format(XLABLE), fontsize=25)
    ax2.set_ylabel("Data / MC", fontsize=20)
    ax2.set_xlim(RMIN, RMAX)
    ax2.set_ylim(0, 2)

    plt.tight_layout()

    plt.savefig("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/data_vs_mc_{}.png".format(VAR))