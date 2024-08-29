import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

var_names = [
    "mass_jj",
    "Z_relpt",
    "gamma_relpt",
    "pt_balance_0j",
    "H_mass", 
    "delta_eta_jj",
    "delta_phi_jj",
    "delta_phi_zgjj",
    "gamma_eta",
    "gamma_mvaID",
    "H_ptt",
    "jet_1_pt",
    "jet_2_pt",
    "jet1G_deltaR",
    "jet2G_deltaR",
    "jet_1_btagDeepFlavB",
    "jet_2_btagDeepFlavB",
    "l1g_deltaR",
    "l2g_deltaR",
    "lep_cos_theta",
    "lep_phi",
    "photon_zeppenfeld",
    "pt_balance",
    "Z_cos_theta",
    "Z_lead_lepton_eta",
    "Z_sublead_lepton_eta",
    "H_relpt",
    "HZ_deltaRap"
    ]
x_ranges = [
    (0, 300),
    (0, 1),
    (0.1, 1),
    (0, 1),
    (100, 180),
    (0, 9),
    (0, 3.2),
    (0, 3.2),
    (-2.5, 2.5),
    (0.14, 1),
    (0, 160),
    (30, 330),
    (30, 150),
    (0.4, 6.4),
    (0.4, 6.4),
    (0, 0.1),
    (0, 0.15),
    (0.4, 4.8),
    (0.4, 3.4),
    (-1, 1),
    (-3.2, 3.2),
    (0, 5),
    (0, 1),
    (-1, 1),
    (-2.5, 2.5),
    (-2.5, 2.5),
    (0, 3),
    (-0.7, 0.7)
    ]
x_titles = [
    r"$m_{jj}(GeV/c^{2})$",
    r"$\frac{p_T^{ll}}{p_{T_{Z\gamma}}}$",
    r"$\frac{p_T^{\gamma}}{p_{T_{Z\gamma}}}$",
    r"$Zeppenfeld_{\gamma}(0j)$",
    r"$m_{ll\gamma}(GeV/c^{2})$", 
    r"$\Delta\eta_{jj}$",
    r"$\Delta\phi_{jj}$",
    r"$\Delta\phi_{z\gamma,jj}$",
    r"$\eta_{\gamma}$",
    "photon MVA",
    r"$p_{T_{t}}^{Z\gamma}$",
    r"$p_{T_{j1}}(GeV/c)$",
    r"$p_{T_{j2}}(GeV/c)$",
    r"$\Delta R(\gamma,j1)$",
    r"$\Delta R(\gamma,j2)$",
    "j1 btag",
    "j2 btag",
    r"max($\Delta R(l,\gamma)$)",
    r"min($\Delta R(l,\gamma)$)",
    r"$\cos\theta$",
    r"$\phi$",
    "Zeppenfeld $\gamma$",
    "system balance",
    r"$\cos\Theta$",
    r"$\eta_{l1}$",
    r"$\eta_{l2}$",
    r"$\frac{p_{T_{ll\gamma}}\cdot c}{m_{ll\gamma}}$",
    r"$\Delta y(ll,ll\gamma)$"
]

BINS = 40
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples/"

sig = ["ggH_mix", "VBF_mix"]
sig_st = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
bkg= ["ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]

def convert_root_to_hist(file, target=None):
    type_hist = np.zeros(BINS)
    mass_hist = np.zeros(100)
    if target is not None:
        for f in file:
            for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
                print("Reading", PATH+f+"/"+year+".root:two_jet", "...")
                if "H_mass" not in VAR:
                    try:
                        samples = uproot.open(PATH+f+"/"+year+".root:two_jet").arrays([VAR, "weight", "H_mass"], library="pd")
                    except:
                        continue
                else:
                    try:
                        samples = uproot.open(PATH+f+"/"+year+".root:two_jet").arrays([VAR, "weight"], library="pd")
                    except:
                        continue
                hist, _ = np.histogram(samples["H_mass"], bins=100, range=[100, 180], weights=samples[WEIGHT])
                mass_hist = mass_hist + hist
        ratio_hist = np.divide(target, mass_hist, out=np.zeros_like(target), where=mass_hist != 0)
        print(ratio_hist)

    for f in file:
        for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
            print("Reading", PATH+f+"/"+year+".root:two_jet", "...")
            if "H_mass" not in VAR:
                try:
                    samples = uproot.open(PATH+f+"/"+year+".root:two_jet").arrays([VAR, WEIGHT, "H_mass"], library="pd")
                except:
                    continue
            else:
                try:
                    samples = uproot.open(PATH+f+"/"+year+".root:two_jet").arrays([VAR, WEIGHT], library="pd")
                except:
                    continue
            hist, _ = np.histogram(samples["H_mass"], bins=100, range=[100, 180], weights=samples[WEIGHT])
            mass_hist = mass_hist + hist
            if target is not None:
                bin_indices = np.digitize(samples['H_mass'], bins=np.linspace(100, 180, 101)) - 1
                samples['weight'] = samples[WEIGHT] * np.where((bin_indices >= 0) & (bin_indices < 100), ratio_hist[bin_indices], 1)
            hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
            type_hist = type_hist + hist
    return type_hist, bins, mass_hist

for i in range(len(var_names)):
    # if i != 4 :
    #     continue
    VAR = var_names[i]
    XLABLE = x_titles[i]
    RMIN = x_ranges[i][0]
    RMAX = x_ranges[i][1]

    print("\n\n", VAR, RMIN, RMAX, "\n\n")

    hist2, bins, mass_hist = convert_root_to_hist(bkg)
    hist1, bins, _ = convert_root_to_hist(sig, mass_hist)
    # hist1, bins, _ = convert_root_to_hist(sig)
    hist1_1, bins, _ = convert_root_to_hist(sig_st)

    bins = (bins[:-1] + bins[1:]) / 2

    plt.plot(bins, hist1/sum(hist1), "o-", label="Signal")
    plt.plot(bins, hist1_1/sum(hist1_1), "*-", label="Signal std")
    plt.plot(bins, hist2/sum(hist2), "x-", label="Background")

    plt.xlabel(XLABLE)
    plt.ylabel("Events")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.title(f"Signal Background comparison")
    plt.savefig(f"pic/compare_{VAR}.png")
    plt.clf()