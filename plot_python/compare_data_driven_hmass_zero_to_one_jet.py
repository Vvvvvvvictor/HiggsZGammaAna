import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

def apply_sf_to_weights(data, sf, bin_edges):
    """
    对 DataFrame 中的权重应用缩放因子。
    
    参数:
    - data: 包含 'gamma_pt' 和 'weight' 列的 DataFrame
    - sf: 缩放因子的数组
    - bin_edges: bin 的边界数组，用于划分 gamma pt
    
    返回:
    - 应用缩放因子后的 DataFrame
    """
    # 确定每个 Z_pt 所属的 bin
    bins = np.digitize(data['Z_pt'], bin_edges) - 1
    # 保证 bins 在有效范围内
    bins = np.clip(bins, 0, len(sf) - 1)
    # 应用缩放因子到 weight 列
    data['weight'] *= sf[bins]
    return data

# Load data
branches = ["H_mass", "Z_pt", "weight", "bdt_score", "gamma_relpt"]

data_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
zg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/ZGToLLG_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
ewkzg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/ZG2JToG2L2J_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
    
# boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_zero_to_one_jet.txt", "r").readlines()[2].split()[1:]]
boundaries = [0.28379330039024353, 0.4557725191116333, 0.5796570777893066, 0.7069960236549377, 1.]
print(boundaries)

# ewkzg_cr = ewkzg_cr[ewkzg_cr["H_mass"]<100]

################################################################
# Compare data-driven photon samples with fake photon samples in control region
################################################################

bkg_list = ["DYJets_01J"]
def draw(data_cr, zg_cr, ewkzg_cr, boundaries):
    data_cr, zg_cr = data_cr[(data_cr["bdt_score"]>boundaries[0]) & (data_cr["bdt_score"]<boundaries[1])], zg_cr[(zg_cr["bdt_score"]>boundaries[0]) & (zg_cr["bdt_score"]<boundaries[1])]
    ewkzg_cr = ewkzg_cr[(ewkzg_cr["bdt_score"]>boundaries[0]) & (ewkzg_cr["bdt_score"]<boundaries[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=80, weights=data_cr["weight"], range=[100, 180])
    zg_cr_hist, bin_edges = np.histogram(zg_cr["H_mass"], bins=80, weights=zg_cr["weight"], range=[100, 180])
    ewkzg_cr_hist, bin_edges = np.histogram(ewkzg_cr["H_mass"], bins=80, weights=ewkzg_cr["weight"], range=[100, 180])
    cr_hist = (data_cr_hist - zg_cr_hist - ewkzg_cr_hist)
    bkg_hists = []
    bkg_sum = 0
    for i in range(len(bkg_list)):
        bkg = bkg_list[i]
        bkg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}.root")["zero_to_one_jet"].arrays(branches, library="pd")
        bkg_cr = bkg_cr[(bkg_cr["bdt_score"]>boundaries[0]) & (bkg_cr["bdt_score"]<boundaries[1])]
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
        bkg_sum += np.sum(bkg_cr_hist)
        print(bkg_sum)
        print("sum of cr ", bkg, " in ", boundaries, ": ", len(bkg_cr["weight"]), "error:", np.sqrt(1/np.sum(1/bkg_cr["weight"]))*100)
    bkg_hists = [bkg_hist/bkg_sum for bkg_hist in bkg_hists]
    print("sum of cr data in ", boundaries, ": ", data_cr["weight"].sum(), "error:", np.sqrt(1/np.sum(1/data_cr["weight"]*data_cr["weight"].sum()/bkg_sum))*100)
    
    plt.figure()
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill")
    plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, cr_hist/np.sum(cr_hist), fmt='o', label="Data driven", yerr=np.sqrt(cr_hist)/np.sum(cr_hist), color='k')
    plt.legend(ncol=2, fontsize=14, loc="best")
    plt.title(f"Data-driven bkg. vs bkg: {boundaries[0]} < BDT < {boundaries[1]}")
    plt.xlim(100, 180)
    plt.grid()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_driven_hmass_{boundaries[0]}_{boundaries[1]}.png")
    
# Plot
# for i in range(0, len(boundaries)-1):
#     draw(data_cr, zg_cr, ewkzg_cr, [boundaries[i], boundaries[i+1]])

################################################################
# Compare data and MC samples in control region
################################################################

bkg_list = ["ZGToLLG", "ZG2JToG2L2J", "DYJets_01J"]

def draw(data_cr, bkg_list, bounds):
    data_cr = data_cr[(data_cr["bdt_score"]>bounds[0]) & (data_cr["bdt_score"]<bounds[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=80, weights=data_cr["weight"], range=[100, 180])
    
    bkg_hists = []
    for i in range(len(bkg_list)):
        bkg = bkg_list[i]
        bkg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
        bkg_cr = bkg_cr[(bkg_cr["bdt_score"]>bounds[0]) & (bkg_cr["bdt_score"]<bounds[1])]
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
    
    plt.figure()
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill")
    plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_cr_hist, fmt='o', label="Data", yerr=np.sqrt(data_cr_hist), color='k')
    plt.legend(fontsize=16, loc="best")
    plt.grid()
    plt.title(f"Data vs MC(CR): {bounds[0]} < BDT < {bounds[1]}")
    plt.xlabel("Higgs Mass")
    plt.ylabel("Events")
    plt.xlim(100, 180)
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_mc_hmass_{bounds[0]}_{bounds[1]}.png")
    
# Plot
# for i in range(0, len(boundaries)-1):
#     draw(data_cr, bkg_list, [boundaries[i], boundaries[i+1]])

################################################################
# Compare true photon samples plus data-driven photon samples with data in signal region
################################################################


bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]#, "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
def draw(data_cr, zg_cr, ewkzg_cr, bkg_list, bounds):
    data_cr, zg_cr = data_cr[(data_cr["bdt_score"]>bounds[0]) & (data_cr["bdt_score"]<bounds[1])], zg_cr[(zg_cr["bdt_score"]>bounds[0]) & (zg_cr["bdt_score"]<bounds[1])]
    ewkzg_cr = ewkzg_cr[(ewkzg_cr["bdt_score"]>bounds[0]) & (ewkzg_cr["bdt_score"]<bounds[1])]
    
    data_gpt_hist, bin_edges = np.histogram(data_cr["Z_pt"], bins=100, weights=data_cr["weight"], range=[0, 100])
    data_gpt_hist[-1] += np.sum(data_cr[data_cr["Z_pt"] > 100]["weight"])
    zg_gpt_hist, bin_edges = np.histogram(zg_cr["Z_pt"], bins=100, weights=zg_cr["weight"], range=[0, 100])
    zg_gpt_hist[-1] += np.sum(zg_cr[zg_cr["Z_pt"] > 100]["weight"])
    ewkzg_gpt_hist, bin_edges = np.histogram(ewkzg_cr["Z_pt"], bins=100, weights=ewkzg_cr["weight"], range=[0, 100])
    ewkzg_gpt_hist[-1] += np.sum(ewkzg_cr[ewkzg_cr["Z_pt"] > 100]["weight"])
    cr_hist = (data_gpt_hist - zg_gpt_hist - ewkzg_gpt_hist)
    
    bkg_samples = []
    for bkg in bkg_list:
        bkg_sample = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}.root")["zero_to_one_jet"].arrays(branches, library="pd")
        bkg_sample = bkg_sample[(bkg_sample["bdt_score"]>bounds[0]) & (bkg_sample["bdt_score"]<bounds[1])]
        bkg_samples.append(bkg_sample)
        
    data_sr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data.root")["zero_to_one_jet"].arrays(branches, library="pd")
    data_sr = data_sr[(data_sr["bdt_score"]>bounds[0]) & (data_sr["bdt_score"]<bounds[1])]
    
    # bkg_gpt_hists = np.zeros(100)
    # for bkg_sample in bkg_samples:
    #     bkg_gpt_hists += np.histogram(bkg_sample["Z_pt"], bins=100, weights=bkg_sample["weight"], range=[0, 100])[0]
    #     bkg_gpt_hists[-1] += np.sum(bkg_sample[bkg_sample["Z_pt"] > 100]["weight"])
    
    # data_gpt_hist_sr, bin_edges = np.histogram(data_sr["Z_pt"], bins=100, weights=data_sr["weight"], range=[0, 100])
    # data_gpt_hist_sr[-1] += np.sum(data_sr[data_sr["Z_pt"] > 100]["weight"])
    
    # # gpt_sf = np.nor(data_gpt_hist_sr-bkg_gpt_hists)/cr_hist
    # gpt_sf = (data_gpt_hist_sr-bkg_gpt_hists)/np.sum(data_gpt_hist_sr-bkg_gpt_hists)/cr_hist*np.sum(cr_hist)
    # bin_edges[-1] = 200
    # print(np.sum(data_cr["weight"]))
    # data_cr_mod = apply_sf_to_weights(data_cr, gpt_sf, bin_edges)
    # print(np.sum(data_cr_mod["weight"]))
    # zg_cr_mod = apply_sf_to_weights(zg_cr, gpt_sf, bin_edges)
    # ewkzg_cr_mod = apply_sf_to_weights(ewkzg_cr, gpt_sf, bin_edges)
    
    data_cr_mod, zg_cr_mod, ewkzg_cr_mod = data_cr, zg_cr, ewkzg_cr
    
    data_cr_hist, bin_edges = np.histogram(data_cr_mod["H_mass"], bins=80, weights=data_cr_mod["weight"], range=[100, 180])
    zg_cr_hist, bin_edges = np.histogram(zg_cr_mod["H_mass"], bins=80, weights=zg_cr_mod["weight"], range=[100, 180])
    ewkzg_cr_hist, bin_edges = np.histogram(ewkzg_cr_mod["H_mass"], bins=80, weights=ewkzg_cr_mod["weight"], range=[100, 180])
    
    cr_err_hist = (np.histogram(data_cr_mod["H_mass"], bins=80, weights=data_cr_mod["weight"]**2, range=[100, 180])[0] + np.histogram(zg_cr_mod["H_mass"], bins=80, weights=zg_cr_mod["weight"]**2, range=[100, 180])[0] + np.histogram(ewkzg_cr_mod["H_mass"], bins=80, weights=ewkzg_cr_mod["weight"]**2, range=[100, 180])[0]) ** 0.5

    # cr_hist = (data_cr_hist - zg_cr_hist)
    cr_hist = (data_cr_hist - zg_cr_hist - ewkzg_cr_hist)
  
    plt.figure()
    bkg_hists = []
    bkg_sum = 0
    bkg_err_hist = np.zeros(80)
    for bkg_cr in bkg_samples:
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
        bkg_err_hist += np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"]**2, range=[100, 180])[0]
        bkg_sum = bkg_sum + np.sum(bkg_cr_hist[np.where((bin_edges[1:]<120) | (bin_edges[:-1]>130))])
    bkg_err_hist = bkg_err_hist ** 0.5
    
    data_hist, bin_edges = np.histogram(data_sr[(data_sr["H_mass"]<120) | (data_sr["H_mass"]>130)]["H_mass"], bins=80, weights=data_sr[(data_sr["H_mass"]<120) | (data_sr["H_mass"]>130)]["weight"], range=[100, 180])
    cr_side_band = cr_hist[np.where((bin_edges[1:]<120) | (bin_edges[:-1]>130))]
    sf = (np.sum(data_hist)-bkg_sum)/np.sum(cr_side_band)
    # print(np.sum(data_hist)-bkg_sum, np.sum(cr_side_band), sf)
  
    # print(cr_hist)

    bkg_hist_sum = bkg_hists + [cr_hist*sf]
    ratio = (data_hist / np.sum(bkg_hist_sum, axis=0))
    ratio_test = (data_hist - np.sum(bkg_hists, axis=0)) / cr_hist / sf
    ratio_test_err = ((np.sqrt(data_hist) / cr_hist / sf)**2) ** 0.5# + (bkg_err_hist / cr_hist / sf)**2 + (cr_err_hist * (data_hist - np.sum(bkg_hists, axis=0)) / cr_hist**2 / sf**2)**2) ** 0.5

    # fit this ratio distribution with linear function
    from scipy.optimize import curve_fit
    def linear(x, a, b, c, d, e, f):
        return a*x**5 + b*x**4 + c*x**3 + d*x**2 + e*x + f
    
    x = bin_edges[:-1]+np.diff(bin_edges)/2
    y = ratio_test[np.where((x<120) | (x>130))]
    # yerr = 1/(data_hist - np.sum(bkg_hists, axis=0))[np.where((x<120) | (x>130))]
    yerr = ratio_test_err[np.where((x<120) | (x>130))]
    x = x[np.where((x<120) | (x>130))]
    popt, pcov = curve_fit(linear, x-100, y, sigma=yerr, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    print(popt, perr)

    # correct the bkg_hist_sum and ratio with the fitted function
    # bkg_hist_sum = [bkg_hist_sum[i] * linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt) for i in range(len(bkg_hist_sum))]
    bkg_hist_sum = [cr_hist*sf*linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt)] + bkg_hists

    # Add a ratio plot in the lower panel
    fig, [ax1, ax2, ax3] = plt.subplots(3, 1, gridspec_kw={"height_ratios": [4, 1, 1], "hspace": 0}, sharex='col')

    ax1.set_title(f"Data vs Bkg.: {bounds[0]:.2f} < BDT < {bounds[1]:.2f}") 
    hep.histplot(bkg_hist_sum, bins=bin_edges, label=["Fake photon"]+bkg_list, stack=True, histtype="fill", ax=ax1) 
    ax1.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_hist, fmt='o', label="Data", yerr=np.sqrt(data_hist), color='k')
    ax1.legend(ncol=1, fontsize=16)
    ax1.set_xlim(100, 180)
    ax1.grid()
    ax1.label_outer()

    ax2.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, yerr=np.sqrt(data_hist)/np.sum(bkg_hist_sum), fmt='o')
    # ax2.plot(bin_edges[:-1]+np.diff(bin_edges)/2, linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt)+ratio-ratio_test, 'r--')
    ax2.text(103, 1.45, "pre-fit")
    ax2.axhline(y=1, color='k', linestyle='--')
    ax2.set_ylim(0.2, 1.8)
    ax2.set_xlim(100, 180)
    ax2.set_ylabel("Dat. / Bkg.", fontsize=20)
    ax2.grid()
    ax2.label_outer()

    ratio = (data_hist / np.sum(bkg_hist_sum, axis=0))
    ax3.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, yerr=np.sqrt(data_hist)/np.sum(bkg_hist_sum), fmt='o')
    ax3.text(103, 1.45, "post-fit")
    ax3.axhline(y=1, color='k', linestyle='--')
    ax3.set_ylim(0.2, 1.8)
    ax3.set_xlim(100, 180)
    ax3.set_xlabel("Higgs Mass")
    ax3.set_ylabel("Dat. / Bkg.", fontsize=20)
    ax3.grid()
    ax3.label_outer()

    plt.tight_layout()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_driven_bkg_{bounds[0]}_{bounds[1]}.png")
    plt.clf()
    
    # plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio_test, yerr=np.sqrt(data_hist-np.sum(bkg_hists, axis=0))/cr_hist/sf, fmt='o', label="Data", color='k')
    plt.errorbar(x, y, yerr=yerr, fmt='o', label="Data", color='k')
    plt.plot(bin_edges[:-1]+np.diff(bin_edges)/2, linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt), 'r--', label="Fit")
    plt.axhline(y=1, color='k', linestyle='--')
    plt.legend()
    plt.ylim(-3, 4)
    plt.xlim(100, 180)
    plt.xlabel("Higgs Mass")
    plt.ylabel("Dat. / Fake Bkg.")
    plt.grid()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_driven_bkg_fit_{bounds[0]}_{bounds[1]}.png")

# Plot
for i in range(0, len(boundaries)-1):
    draw(data_cr, zg_cr, ewkzg_cr, bkg_list, [boundaries[i], boundaries[i+1]])

###############################################################################
# Compare data and MC in signal reigon
###############################################################################

bkg_list = ["ZGToLLG", "ZG2JToG2L2J", "DYJets_01J"]
data = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data.root")["zero_to_one_jet"].arrays(branches, library="pd")
bkgs = [uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}.root")["zero_to_one_jet"].arrays(branches, library="pd") for bkg in bkg_list]

def draw(data, bkgs, bounds):
    data = data[(data["bdt_score"]>bounds[0]) & (data["bdt_score"]<bounds[1])]
    data = data[(data["H_mass"]<122) | (data["H_mass"]>128)]
    data_hist, bin_edges = np.histogram(data["H_mass"], bins=80, weights=data["weight"], range=[100, 180])
    
    bkg_hists = []
    for bkg in bkgs:
        bkg = bkg[(bkg["bdt_score"]>bounds[0]) & (bkg["bdt_score"]<bounds[1])]
        bkg_hist, bin_edges = np.histogram(bkg["H_mass"], bins=80, weights=bkg["weight"], range=[100, 180])
        bkg_hists.append(bkg_hist)
    
    fig, [ax1, ax2] = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill", ax=ax1)
    ax1.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_hist, fmt='o', label="Data", yerr=np.sqrt(data_hist), color='k')
    ax1.legend(fontsize=16, loc="best")
    ax1.grid()
    ax1.set_title(f"Data vs MC(SR): {bounds[0]} < BDT < {bounds[1]}")
    ax1.set_xlabel("Higgs Mass")
    ax1.set_ylabel("Events")
    ax1.set_xlim(100, 180)
    
    ratio = (data_hist / np.sum(bkg_hists, axis=0))
    ax2.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, yerr=np.sqrt(ratio/np.sum(bkg_hists, axis=0)), fmt='o')
    ax2, ax1.axhline(y=1, color='k', linestyle='--')
    ax2.set_ylim(0, 2)
    ax2.set_xlim(100, 180)
    ax2.set_xlabel("Higgs Mass")
    ax2.set_ylabel("Data / MC")
    ax2.grid()
    
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_mc_sr_{bounds[0]}_{bounds[1]}.png")
    
# Plot
# for i in range(0, len(boundaries)-1):
#     draw(data, bkgs, [boundaries[i], boundaries[i+1]])

###############################################################################
# Compare data and MC in Control region
###############################################################################

bkg_list = ["ZGToLLG", "ZG2JToG2L2J", "DYJets_01J"]
data = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
bkgs = [uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}_fake.root")["zero_to_one_jet"].arrays(branches, library="pd") for bkg in bkg_list]

def draw(data, bkgs, bounds):
    data = data[(data["bdt_score"]>bounds[0]) & (data["bdt_score"]<bounds[1])]
    data_hist, bin_edges = np.histogram(data["H_mass"], bins=80, weights=data["weight"], range=[100, 180])
    
    bkg_hists = []
    for bkg in bkgs:
        bkg = bkg[(bkg["bdt_score"]>bounds[0]) & (bkg["bdt_score"]<bounds[1])]
        bkg_hist, bin_edges = np.histogram(bkg["H_mass"], bins=80, weights=bkg["weight"], range=[100, 180])
        bkg_hists.append(bkg_hist)
    
    fig, [ax1, ax2] = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill", ax=ax1)
    ax1.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_hist, fmt='o', label="Data", yerr=np.sqrt(data_hist), color='k')
    ax1.legend(fontsize=16, loc="best")
    ax1.grid()
    ax1.set_title(f"Data vs MC(CR): {bounds[0]} < BDT < {bounds[1]}")
    ax1.set_xlabel("Higgs Mass")
    ax1.set_ylabel("Events")
    ax1.set_xlim(100, 180)
    
    ratio = (data_hist / np.sum(bkg_hists, axis=0))
    ax2.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, yerr=np.sqrt(ratio/np.sum(bkg_hists, axis=0)), fmt='o')
    ax2, ax1.axhline(y=1, color='k', linestyle='--')
    ax2.set_ylim(0, 2)
    ax2.set_xlim(100, 180)
    ax2.set_xlabel("Higgs Mass")
    ax2.set_ylabel("Data / MC")
    ax2.grid()
    
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_mc_cr_{bounds[0]}_{bounds[1]}.png")
    
# Plot
# for i in range(0, len(boundaries)-1):
#     draw(data, bkgs, [boundaries[i], boundaries[i+1]])

###############################################################################
# Compare data-driven photon samples plus MC samples(SR) with data(SR) for different variables distributions
###############################################################################
variables = ["H_mass", "gamma_pt", "Z_pt", "gamma_mvaID"]
xlabels = [r"Higgs Mass [GeV/$c^2$]", "Photon pT [GeV/c]", "Z pT [GeV/c]", "Photon MVA ID"]
x_ranges = [(100, 180), (15, 75), (0, 100), (0, 1)]
xbin = 40
bkg_list = ["DYJets_01J"]
true_bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]
branches = ["weight", "bdt_score", "gamma_relpt"] + variables

data_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
zg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/ZGToLLG_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")
ewkzg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/ZG2JToG2L2J_fake.root")["zero_to_one_jet"].arrays(branches, library="pd")

for i in range(len(variables)):
    def draw(data_cr, zg_cr, ewkzg_cr, bkg_list, bounds, variable, xlabel, x_range):
        data_cr = data_cr[(data_cr["bdt_score"]>bounds[0]) & (data_cr["bdt_score"]<bounds[1]) & (data_cr["gamma_relpt"] > 0.3)]
        data_cr_hist, bin_edges = np.histogram(data_cr[variable], bins=xbin, weights=data_cr["weight"], range=x_range)
        
        zg_cr = zg_cr[(zg_cr["bdt_score"]>bounds[0]) & (zg_cr["bdt_score"]<bounds[1]) & (zg_cr["gamma_relpt"] > 0.3)]
        zg_sr_hist, bin_edges = np.histogram(zg_cr[variable], bins=xbin, weights=zg_cr["weight"], range=x_range)
        
        ewkzg_cr = ewkzg_cr[(ewkzg_cr["bdt_score"]>bounds[0]) & (ewkzg_cr["bdt_score"]<bounds[1]) & (ewkzg_cr["gamma_relpt"] > 0.3)]
        ewkzg_sr_hist, bin_edges = np.histogram(ewkzg_cr[variable], bins=xbin, weights=ewkzg_cr["weight"], range=x_range)
        
        bkg_hists = []
        for i in range(len(bkg_list)):
            bkg = bkg_list[i]
            bkg_sr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}.root")["zero_to_one_jet"].arrays(branches, library="pd")
            bkg_sr = bkg_sr[(bkg_sr["bdt_score"]>bounds[0]) & (bkg_sr["bdt_score"]<bounds[1]) & (bkg_sr["gamma_relpt"] > 0.3)]
            bkg_sr_hist, bin_edges = np.histogram(bkg_sr[variable], bins=xbin, weights=bkg_sr["weight"], range=x_range)
            bkg_hists.append(bkg_sr_hist)
            
        bkg_sum = np.sum(bkg_hists)
        bkg_hists = [bkg_hist/bkg_sum for bkg_hist in bkg_hists]
        data_sum = np.sum(data_cr_hist)
        data_cr_hist = data_cr_hist/data_sum
        
        plt.figure()
        hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill")
        plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_cr_hist, fmt='o', label="Data", yerr=np.sqrt(data_cr_hist/data_sum), color='k')
        plt.legend(fontsize=16, loc="best")
        plt.grid()
        plt.title(f"Data vs MC(SR): {bounds[0]} < BDT < {bounds[1]}")
        plt.xlabel(xlabel)
        plt.ylabel("Events")
        plt.xlim(x_range[0], x_range[1])
        plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_mc_{variable}_{bounds[0]}_{bounds[1]}.png")
    
    # Plot
    # for j in range(0, len(boundaries)-1):
    #     draw(data_cr, zg_cr, ewkzg_cr, bkg_list, [boundaries[j], boundaries[j+1]], variables[i], xlabels[i], x_ranges[i])
        
for i in range(len(variables)):
    def draw(bkg_list, true_bkg_list, bounds, variable, xlabel, x_range):
        data_sr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data.root")["zero_to_one_jet"].arrays(branches, library="pd")
        data_sr = data_sr[(data_sr["bdt_score"]>bounds[0]) & (data_sr["bdt_score"]<bounds[1])]
        data_sr_hist, bin_edges = np.histogram(data_sr[variable], bins=xbin, weights=data_sr["weight"], range=x_range)
        
        bkg_hists = []
        for i in range(len(bkg_list)):
            bkg = bkg_list[i]
            bkg_sr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}.root")["zero_to_one_jet"].arrays(branches, library="pd")
            bkg_sr = bkg_sr[(bkg_sr["bdt_score"]>bounds[0]) & (bkg_sr["bdt_score"]<bounds[1])]
            bkg_sr_hist, bin_edges = np.histogram(bkg_sr[variable], bins=xbin, weights=bkg_sr["weight"], range=x_range)
            bkg_hists.append(bkg_sr_hist)
        for i in range(len(true_bkg_list)):
            bkg = true_bkg_list[i]
            bkg_sr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/{bkg}.root")["zero_to_one_jet"].arrays(branches, library="pd")
            bkg_sr = bkg_sr[(bkg_sr["bdt_score"]>bounds[0]) & (bkg_sr["bdt_score"]<bounds[1])]
            bkg_sr_hist, bin_edges = np.histogram(bkg_sr[variable], bins=xbin, weights=bkg_sr["weight"], range=x_range)
            bkg_hists.append(bkg_sr_hist)
            
        fig, [ax1, ax2] = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})
        hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list+true_bkg_list, stack=True, histtype="fill", ax=ax1)
        ax1.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_sr_hist, fmt='o', label="Data", yerr=np.sqrt(data_sr_hist), color='k')
        ax1.legend(fontsize=16, loc="best")
        ax1.grid()
        ax1.set_title(f"Data vs MC(SR): {bounds[0]} < BDT < {bounds[1]}")
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel("Events")
        ax1.set_xlim(x_range[0], x_range[1])
        
        ratio = (data_sr_hist / np.sum(bkg_hists, axis=0))
        ax2.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, yerr=np.sqrt(ratio/np.sum(bkg_hists, axis=0)), fmt='o')
        ax2, ax1.axhline(y=1, color='k', linestyle='--')
        ax2.set_ylim(0, 2)
        ax2.set_xlim(x_range[0], x_range[1])
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel("Data / MC")
        ax2.grid()
        
        plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/01J/compare_data_mc_sr_{variable}_{bounds[0]}_{bounds[1]}.png")
        
    # Plot
    # for j in range(0, len(boundaries)-1):
    #     draw(bkg_list, true_bkg_list, [boundaries[j], boundaries[j+1]], variables[i], xlabels[i], x_ranges[i])
        
        