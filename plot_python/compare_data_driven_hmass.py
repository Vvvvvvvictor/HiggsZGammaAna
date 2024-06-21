import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

# Load data
branches = ["H_mass", "weight", "bdt_score_t"]
data_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/Data_fake.root")["two_jet"].arrays(branches, library="pd")
zg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZGToLLG_fake.root")["two_jet"].arrays(branches, library="pd")

dy_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/DYJetsToLL.root")["test"].arrays(branches, library="pd")

boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt", "r").readlines()[2].split()[1:]]

bkg_list = ["DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
def draw(data_cr, zg_cr, dy_sr, boundaries):
    data_cr, zg_cr = data_cr[(data_cr["bdt_score_t"]>boundaries[0]) & (data_cr["bdt_score_t"]<boundaries[1])], zg_cr[(zg_cr["bdt_score_t"]>boundaries[0]) & (zg_cr["bdt_score_t"]<boundaries[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=40, weights=data_cr["weight"], range=[100, 180])
    zg_cr_hist, bin_edges = np.histogram(zg_cr["H_mass"], bins=40, weights=zg_cr["weight"], range=[100, 180])
    cr_hist = (data_cr_hist - zg_cr_hist)
    bkg_hists = []
    bkg_sum = 0
    for i in range(len(bkg_list)):
        bkg = bkg_list[i]
        bkg_cr = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
        bkg_cr = bkg_cr[(bkg_cr["bdt_score_t"]>boundaries[0]) & (bkg_cr["bdt_score_t"]<boundaries[1])]
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=40, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
        bkg_sum += np.sum(bkg_cr_hist)
    bkg_hists = [bkg_hist/bkg_sum for bkg_hist in bkg_hists]
    
    plt.figure()
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill")
    plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, cr_hist/np.sum(cr_hist), fmt='o', label="Data driven", yerr=np.sqrt(cr_hist)/np.sum(cr_hist), color='k')
    plt.legend(ncol=2, fontsize=14, loc="best")
    plt.title(f"Data-driven bkg. vs bkg: {boundaries[0]} < BDT < {boundaries[1]}")
    plt.xlim(100, 180)
    plt.grid()
    plt.savefig(f"pic/compare_data_driven_hmass_{boundaries[0]}_{boundaries[1]}.png")
    
# Plot
for i in range(0, len(boundaries)-1):
    draw(data_cr, zg_cr, dy_sr, [boundaries[i], boundaries[i+1]])

# bkg_list = ["ZGToLLG"]#, "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
# def draw(data_cr, zg_cr, bkg_list, bounds):
#     data_cr, zg_cr = data_cr[(data_cr["bdt_score_t"]>bounds[0]) & (data_cr["bdt_score_t"]<bounds[1])], zg_cr[(zg_cr["bdt_score_t"]>bounds[0]) & (zg_cr["bdt_score_t"]<bounds[1])]
#     data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=40, weights=data_cr["weight"], range=[100, 180])
#     zg_cr_hist, bin_edges = np.histogram(zg_cr["H_mass"], bins=40, weights=zg_cr["weight"], range=[100, 180])
#     cr_hist = (data_cr_hist - zg_cr_hist)
    
#     plt.figure()
#     bkg_hists = []
#     bkg_sum = 0
#     for bkg in bkg_list:
#         bkg_cr = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
#         bkg_cr = bkg_cr[(bkg_cr["bdt_score_t"]>bounds[0]) & (bkg_cr["bdt_score_t"]<bounds[1])]
#         bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=40, weights=bkg_cr["weight"], range=[100, 180])
#         bkg_hists.append(bkg_cr_hist)
#         bkg_sum = bkg_sum + np.sum(bkg_cr_hist[np.where((bin_edges[1:]<120) | (bin_edges[:-1]>130))])
        
#     data_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/data.root")["test"].arrays(branches, library="pd")
#     data_sr = data_sr[(data_sr["bdt_score_t"]>bounds[0]) & (data_sr["bdt_score_t"]<bounds[1])]
#     data_hist, bin_edges = np.histogram(data_sr[(data_sr["H_mass"]<120) | (data_sr["H_mass"]>130)]["H_mass"], bins=40, weights=data_sr[(data_sr["H_mass"]<120) | (data_sr["H_mass"]>130)]["weight"], range=[100, 180])
    
#     sf = np.sum(data_hist)-bkg_sum
    
#     print(cr_hist)
#     # Add a ratio plot in the lower panel
#     plt.figure()
#     fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})
    
#     plt.subplot(2, 1, 1)
#     plt.title(f"Data vs Bkg.: {bounds[0]} < BDT < {bounds[1]}")
#     cr_side_band = cr_hist[np.where((bin_edges[1:]<120) | (bin_edges[:-1]>130))]
#     hep.histplot([cr_hist/np.sum(cr_hist)*sf]+bkg_hists, bins=bin_edges, label=["DY+jets"]+bkg_list, stack=True, histtype="fill") 
#     plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_hist, fmt='o', label="Data", yerr=np.sqrt(data_hist), color='k')
#     plt.legend(ncol=2, fontsize=16)
#     plt.xlim(100, 180)
#     plt.grid()

#     plt.subplot(5, 1, 5)
#     ratio = (data_hist / np.sum(bkg_hists+[cr_hist/np.sum(cr_hist)*sf], axis=0))
#     plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, xerr=np.diff(bin_edges)/2, fmt='o')
#     plt.axhline(y=1, color='k', linestyle='--')
#     plt.ylim(0, 2)
#     plt.xlim(100, 180)
#     plt.xlabel("Higgs Mass")
#     plt.ylabel("Data / Bkg.")
#     plt.grid()

#     plt.tight_layout()
#     plt.savefig(f"pic/compare_data_driven_bkg_{bounds[0]}_{bounds[1]}.png")

# # Plot
# for i in range(0, len(boundaries)-1):
#     draw(data_cr, zg_cr, bkg_list, [boundaries[i], boundaries[i+1]])