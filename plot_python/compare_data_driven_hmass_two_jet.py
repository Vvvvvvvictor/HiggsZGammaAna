import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

from pdb import set_trace

# Load data
branches = ["H_mass", "weight", "bdt_score_t", "regions"]
data_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/Data_fake.root")["two_jet"].arrays(branches, library="pd")
zg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZGToLLG_fake.root")["two_jet"].arrays(branches, library="pd")
ewkzg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZG2JToG2L2J_fake.root")["two_jet"].arrays(branches, library="pd")
dy_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/DYJetsToLL.root")["test"].arrays(branches, library="pd")

boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt", "r").readlines()[2].split()[1:]]
print(boundaries)

# ewkzg_cr = ewkzg_cr[ewkzg_cr["H_mass"]<100]

##############################################################################
# Compare data-driven photon samples with fake photon samples(SR)
##############################################################################

bkg_list = ["EWKZ2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets", "DYJetsToLL"]
def draw(data_cr, zg_cr, ewkzg_cr, dy_sr, bounds):
    data_cr, zg_cr = data_cr[(data_cr["bdt_score_t"]>bounds[0]) & (data_cr["bdt_score_t"]<bounds[1])], zg_cr[(zg_cr["bdt_score_t"]>bounds[0]) & (zg_cr["bdt_score_t"]<bounds[1])]
    ewkzg_cr = ewkzg_cr[(ewkzg_cr["bdt_score_t"]>bounds[0]) & (ewkzg_cr["bdt_score_t"]<bounds[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=80, weights=data_cr["weight"], range=[100, 180])
    ewkzg_cr_hist, bin_edges = np.histogram(ewkzg_cr["H_mass"], bins=80, weights=ewkzg_cr["weight"], range=[100, 180])
    zg_cr_hist, bin_edges = np.histogram(zg_cr["H_mass"], bins=80, weights=zg_cr["weight"], range=[100, 180])
    print("data error: ", np.sqrt(1/np.sum(1/data_cr["weight"]))*100)
    
    cr_hist = data_cr_hist - ewkzg_cr_hist - zg_cr_hist
    
    bkg_hists = []
    bkg_sum_err = 0
    for i in range(len(bkg_list)):
        bkg = bkg_list[i]
        bkg_cr = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
        bkg_cr = bkg_cr[(bkg_cr["bdt_score_t"]>bounds[0]) & (bkg_cr["bdt_score_t"]<bounds[1])]
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
        if np.sum(1/bkg_cr["weight"]) <= 0:
            continue
        bkg_sum_err += np.sqrt(1/np.sum(1/bkg_cr["weight"]))*np.sum(bkg_cr_hist)
    bkg_sum = np.sum(bkg_hists)
    print(f"fake photon MC error: {bkg_sum_err/bkg_sum*100}")
    bkg_hists = [bkg_hist/bkg_sum for bkg_hist in bkg_hists]
    
    plt.figure()
    cmap = plt.get_cmap('tab20')  # 'tab20' has 20 different colors
    colors = [cmap(i) for i in range(len(bkg_hists))]
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill", color=colors)
    plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, cr_hist/np.sum(cr_hist), fmt='o', label="Data driven", yerr=np.sqrt(cr_hist)/np.sum(cr_hist), color='k')
    plt.legend(ncol=2, fontsize=14, loc="best")
    plt.title(f"Data-driven bkg. vs bkg: {bounds[0]} < BDT < {bounds[1]}")
    plt.xlim(100, 180)
    plt.grid()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/2J/compare_data_driven_hmass_{bounds[0]}_{bounds[1]}.png")
    plt.clf()
    
# Plot
for i in range(0, len(boundaries)-1):
    draw(data_cr, zg_cr, ewkzg_cr, dy_sr, [boundaries[i], boundaries[i+1]])
    
##############################################################################
# Compare data-driven photon samples plus true photon samples(SR) with data(SR)
##############################################################################

bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]#, "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
fake_bkg_list = ["EWKZ2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets", "DYJetsToLL"]

def draw(data_cr, zg_cr, ewkzg_cr, bkg_list, fake_bkg_list, bounds):
    data_cr, zg_cr = data_cr[(data_cr["bdt_score_t"] > bounds[0]) & (data_cr["bdt_score_t"] < bounds[1])], zg_cr[(zg_cr["bdt_score_t"] > bounds[0]) & (zg_cr["bdt_score_t"] < bounds[1])]
    ewkzg_cr = ewkzg_cr[(ewkzg_cr["bdt_score_t"] > bounds[0]) & (ewkzg_cr["bdt_score_t"] < bounds[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=80, weights=data_cr["weight"], range=[100, 180])
    ewkzg_cr_hist, bin_edges = np.histogram(ewkzg_cr["H_mass"], bins=80, weights=ewkzg_cr["weight"], range=[100, 180])
    zg_cr_hist, bin_edges = np.histogram(zg_cr["H_mass"], bins=80, weights=zg_cr["weight"], range=[100, 180])

    cr_hist = data_cr_hist - ewkzg_cr_hist - zg_cr_hist

    fake_bkg_sum = 0
    for fake_bkg in fake_bkg_list:
        fake_bkg_sample = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{fake_bkg}.root")["test"].arrays(branches, library="pd")
        fake_bkg_sample = fake_bkg_sample[(fake_bkg_sample["bdt_score_t"] > bounds[0]) & (fake_bkg_sample["bdt_score_t"] < bounds[1])]
        fake_bkg_sum += np.sum(fake_bkg_sample["weight"])

    data_driven_sf = fake_bkg_sum / np.sum(cr_hist)

    plt.figure()
    bkg_hists = []
    bkg_sum = 0
    for bkg in bkg_list:
        bkg_cr = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
        bkg_cr = bkg_cr[(bkg_cr["bdt_score_t"] > bounds[0]) & (bkg_cr["bdt_score_t"] < bounds[1])]
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
        bkg_sum = bkg_sum + np.sum(bkg_cr_hist[np.where((bin_edges[1:] < 120) | (bin_edges[:-1] > 130))])

    data_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/data.root")["test"].arrays(branches, library="pd")
    data_sr = data_sr[(data_sr["bdt_score_t"] > bounds[0]) & (data_sr["bdt_score_t"] < bounds[1])]
    data_hist, bin_edges = np.histogram(data_sr[(data_sr["H_mass"] < 120) | (data_sr["H_mass"] > 130)]["H_mass"], bins=80, weights=data_sr[(data_sr["H_mass"] < 120) | (data_sr["H_mass"] > 130)]["weight"], range=[100, 180])
    cr_side_band = cr_hist[np.where((bin_edges[1:] < 120) | (bin_edges[:-1] > 130))]
    sf = np.sum(data_hist) / (bkg_sum + np.sum(cr_side_band) * data_driven_sf)

    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})

    ax1.set_title(f"Data vs Bkg.: {bounds[0]} < BDT < {bounds[1]}")

    final_bkg_hists = [bkg_hist * sf for bkg_hist in bkg_hists] + [cr_hist * sf * data_driven_sf]

    print(f"\nData driven scale factor: {data_driven_sf}, Data/MC scale factor: {sf}")
    # print(np.sum(final_bkg_hists, axis=0), end="\n\n")
    
    ratio = (data_hist / np.sum(final_bkg_hists, axis=0))

    if bounds[0] == 0.15:
        # fit this ratio distribution with linear function
        from scipy.optimize import curve_fit
        def linear(x, a, b):
            return a*x + b
        
        x = bin_edges[:-1]+np.diff(bin_edges)/2
        mask = np.where((data_hist>0) | (np.sum(final_bkg_hists, axis=0)>0) & (ratio>0) & (ratio<2))
        y = ratio[mask]
        yerr = 1/np.sqrt(data_hist[mask])
        x = x[mask]
        popt, pcov = curve_fit(linear, x-100, y, sigma=yerr, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        print(popt, perr)

        # correct the final_bkg_hists and ratio with the fitted function
        final_bkg_hists = [final_bkg_hists[i] * linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt) for i in range(len(final_bkg_hists))]
        ratio = (data_hist / np.sum(final_bkg_hists, axis=0))

    cmap = plt.get_cmap('tab20')  # 'tab20' has 20 different colors
    colors = [cmap(i) for i in range(len(bkg_list)+1)]
    hep.histplot(final_bkg_hists, bins=bin_edges, label=bkg_list+["Fake photon"], color=colors, stack=True, histtype="fill", ax=ax1)
    ax1.errorbar(bin_edges[:-1] + np.diff(bin_edges) / 2, data_hist, fmt='o', label="Data", yerr=np.sqrt(data_hist), color='k')
    ax1.legend(ncol=2, fontsize=16)
    ax1.set_xlim(100, 180)
    ax1.annotate(f"Scale fractor: {sf:.2f}", xy=(0.7, 0.8), xycoords='axes fraction', fontsize=20)
    ax1.grid()

    ax2.errorbar(bin_edges[:-1] + np.diff(bin_edges) / 2, ratio, yerr=np.sqrt(data_hist) / np.sum(final_bkg_hists, axis=0), fmt='o')
    # ax2.plot(bin_edges[:-1] + np.diff(bin_edges) / 2, linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt), label=f"Fit: {popt[0]:.2f}x + {popt[1]:.2f}", color='r')
    ax2.axhline(y=1, color='k', linestyle='--')
    ax2.set_ylim(0, 2)
    ax2.set_xlim(100, 180)
    ax2.set_xlabel("Higgs Mass")
    ax2.set_ylabel("Data / Bkg.")
    ax2.grid()

    plt.tight_layout()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/2J/compare_data_driven_bkg_{bounds[0]}_{bounds[1]}.png")
    plt.clf()
    
# Plot
for i in range(0, len(boundaries)-1):
    draw(data_cr, zg_cr, ewkzg_cr, bkg_list, fake_bkg_list, [boundaries[i], boundaries[i+1]])
    
##############################################################################
# Compare data-driven photon samples with MC samples in control region
##############################################################################

bkg_list = ["ZGToLLG", "ZG2JToG2L2J", "DYJets_2J"]
def draw(data_cr, bkg_list, bounds):
    data_cr= data_cr[(data_cr["bdt_score_t"]>bounds[0]) & (data_cr["bdt_score_t"]<bounds[1])]
    data_cr_hist, bin_edges = np.histogram(data_cr["H_mass"], bins=80, weights=data_cr["weight"], range=[100, 180])
    
    bkg_hist_list = []
    for bkg in bkg_list:
        bkg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/{bkg}_fake.root")["two_jet"].arrays(branches, library="pd")
        bkg_cr = bkg_cr[(bkg_cr["bdt_score_t"]>bounds[0]) & (bkg_cr["bdt_score_t"]<bounds[1])]
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hist_list.append(bkg_cr_hist)
    
    plt.figure()
    fig, [ax1, ax2] = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})
    # print(bkg_hist_list)
    hep.histplot(bkg_hist_list, bins=bin_edges, label=[bkg_list[i]+f"({sum(bkg_hist_list[i]):.2f})" for i in range(len(bkg_list))], stack=True, histtype="fill", ax=ax1)
    # print(data_cr_hist)
    ax1.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_cr_hist, fmt='o', label=f"Data({sum(data_cr_hist):.2f})", yerr=np.sqrt(data_cr_hist), color='k')
    ax1.legend(fontsize=16, loc="best")
    ax1.set_title(f"Data vs bkg: {bounds[0]} < BDT < {bounds[1]}")
    ax1.set_xlim(100, 180)
    ax1.set_ylabel("Events")
    ax1.grid()
    
    ax2.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_cr_hist/np.sum(bkg_hist_list, axis=0), fmt='o', yerr=np.sqrt(data_cr_hist)/np.sum(bkg_hist_list, axis=0))
    ax2.axhline(y=1, color='k', linestyle='--')
    ax2.set_ylim(0, 2)
    ax2.set_xlim(100, 180)
    ax2.set_xlabel("Higgs Mass")
    ax2.set_ylabel("Data / Bkg.")
    
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/2J/compare_data_mc_cr_hmass_{bounds[0]}_{bounds[1]}.png")
    plt.clf()
    
# Plot
# for i in range(0, len(boundaries)-1):
#     draw(data_cr, bkg_list, [boundaries[i], boundaries[i+1]])
draw(data_cr, bkg_list, [0, 1])

###############################################################################
# Compare data-driven photon samples plus MC samples(SR) with data(SR) for different variables distributions
###############################################################################
variables = ["H_mass", "gamma_pt", "Z_pt", "gamma_mvaID"]
xlabels = [r"Higgs Mass [GeV/$c^2$]", "Photon pT [GeV/c]", "Z pT [GeV/c]", "Photon MVA ID"]
x_ranges = [(100, 180), (15, 75), (0, 100), (0, 1)]
xbin = 40
bkg_list = ["EWKZ2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets", "DYJetsToLL"]
branches = ["weight", "bdt_score_t"] + variables

for i in range(len(variables)):
    data_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/Data_fake.root")["two_jet"].arrays(branches, library="pd")
    zg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZGToLLG_fake.root")["two_jet"].arrays(branches, library="pd")
    ewkzg_cr = uproot.open("/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/ZG2JToG2L2J_fake.root")["two_jet"].arrays(branches, library="pd")
    
    bkg_srs = []
    for bkg in bkg_list:
        bkg_sr = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
        bkg_srs.append(bkg_sr)
        
    for j in range(len(boundaries)-1):
        bounds = [boundaries[j], boundaries[j+1]]
        print("\n", bounds)
        
        data_cr_temp = data_cr[(data_cr["bdt_score_t"]>bounds[0]) & (data_cr["bdt_score_t"]<bounds[1])]
        zg_cr_temp = zg_cr[(zg_cr["bdt_score_t"]>bounds[0]) & (zg_cr["bdt_score_t"]<bounds[1])]
        ewkzg_cr_temp = ewkzg_cr[(ewkzg_cr["bdt_score_t"]>bounds[0]) & (ewkzg_cr["bdt_score_t"]<bounds[1])]
    
        data_cr_temp_hist, bin_edges = np.histogram(data_cr_temp[variables[i]], bins=xbin, weights=data_cr_temp["weight"], range=x_ranges[i])
        ewkzg_cr_temp_hist = np.histogram(ewkzg_cr_temp[variables[i]], bins=xbin, weights=ewkzg_cr_temp["weight"], range=x_ranges[i])[0]
        zg_cr_temp_hist = np.histogram(zg_cr_temp[variables[i]], bins=xbin, weights=zg_cr_temp["weight"], range=x_ranges[i])[0]
        
        data_driven_hist = data_cr_temp_hist - ewkzg_cr_temp_hist - zg_cr_temp_hist
        
        mc_sr_hists, mc_sr_nums = [], []
        for bkg_sr in bkg_srs:
            bkg_sr = bkg_sr[(bkg_sr["bdt_score_t"]>bounds[0]) & (bkg_sr["bdt_score_t"]<bounds[1])]
            mc_sr_hist = np.histogram(bkg_sr[variables[i]], bins=xbin, weights=bkg_sr["weight"], range=x_ranges[i])[0]
            mc_sr_hists.append(mc_sr_hist)
            mc_sr_nums.append(np.sum(mc_sr_hist))
            
        plt.figure()
        fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})
    
        plt.subplot(2, 1, 1)
        cmap = plt.get_cmap('tab20')  # 'tab20' has 20 different colors
        colors = [cmap(i) for i in range(len(bkg_list))]
        norm_mc_sr_hists = [mc_sr_hist/np.sum(mc_sr_hists) for mc_sr_hist in mc_sr_hists]
        hep.histplot(norm_mc_sr_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill", color=colors)
        plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, data_driven_hist/np.sum(data_driven_hist), fmt='o', label="Data driven", yerr=np.sqrt(data_driven_hist)/np.sum(data_driven_hist), color='k')
        plt.legend(ncol=2, fontsize=16, loc="best")
        plt.title(f"Data-driven bkg. vs bkg: {variables[i]}")
        plt.xlim(x_ranges[i][0], x_ranges[i][1])
        plt.grid()
        
        plt.subplot(5, 1, 5)
        ratio = (data_driven_hist / np.sum(data_driven_hist) / np.sum(norm_mc_sr_hists, axis=0))
        plt.errorbar(bin_edges[:-1]+np.diff(bin_edges)/2, ratio, yerr=np.sqrt(data_driven_hist)/np.sum(mc_sr_hists, axis=0), fmt='o')
        plt.axhline(y=1, color='k', linestyle='--')
        plt.ylim(0, 2)
        plt.xlim(x_ranges[i][0], x_ranges[i][1])
        plt.xlabel(xlabels[i])
        plt.ylabel("Data / Bkg.")
        plt.grid()
        
        plt.tight_layout()
        plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/2J/compare_data_driven_bkg_{variables[i]}_{bounds[0]}_{bounds[1]}.png")
        plt.clf()


###############################################################################
# Compare data and MC in signal region
###############################################################################

bkg_list = ["ZGToLLG", "ZG2JToG2L2J", "EWKZ2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets", "DYJetsToLL"]

data_sr = uproot.open("/eos/home-j/jiehan/root/outputs/two_jet/data.root")["test"].arrays(branches, library="pd")
data_sr = data_sr[(data_sr["bdt_score_t"]>boundaries[0]) & (data_sr["bdt_score_t"]<boundaries[-1])]

bkg_srs = []
for bkg in bkg_list:
    bkg_sr = uproot.open(f"/eos/home-j/jiehan/root/outputs/two_jet/{bkg}.root")["test"].arrays(branches, library="pd")
    bkg_srs.append(bkg_sr)
    
def draw(data_sr, bkg_srs, bkg_list, bounds):
    data_sr = data_sr[(data_sr["bdt_score_t"] > bounds[0]) & (data_sr["bdt_score_t"] < bounds[1])]
    data_sr = data_sr[(data_sr["H_mass"] < 122) | (data_sr["H_mass"] > 128)]
    data_hist, bin_edges = np.histogram(data_sr["H_mass"], bins=80, weights=data_sr["weight"], range=[100, 180])
    data_sr_num = data_sr.query("H_mass > 128 | H_mass < 122")["weight"].sum()

    bkg_hists = []
    bkg_sr_num = 0
    for bkg_sr in bkg_srs:
        bkg_sr = bkg_sr[(bkg_sr["bdt_score_t"] > bounds[0]) & (bkg_sr["bdt_score_t"] < bounds[1])]
        bkg_sr_num += bkg_sr.query("H_mass > 128 | H_mass < 122")["weight"].sum()
        bkg_hist, bin_edges = np.histogram(bkg_sr["H_mass"], bins=80, weights=bkg_sr["weight"], range=[100, 180])
        bkg_hists.append(bkg_hist)

    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [4, 1]})

    cmap = plt.get_cmap('tab20')  # 'tab20' has 20 different colors
    colors = [cmap(i) for i in range(len(bkg_list))]
    
    sf = data_sr_num / bkg_sr_num
    print(f"\n\n{sf}\n\n")
    bkg_hists = [bkg_hist * sf for bkg_hist in bkg_hists]
    hep.histplot(bkg_hists, bins=bin_edges, label=bkg_list, stack=True, histtype="fill", color=colors, ax=ax1)
    ax1.errorbar(bin_edges[:-1] + np.diff(bin_edges) / 2, data_hist, fmt='o', label="Data", yerr=np.sqrt(data_hist), color='k')
    ax1.legend(ncol=2, fontsize=16, loc="best")
    ax1.set_title(f"Data vs bkg(SR): {bounds[0]} < BDT < {bounds[1]}")
    ax1.set_xlim(100, 180)
    ax1.grid()
    ax1.annotate(f"Scale fractor: {sf:.2f}", xy=(0.7, 0.5), xycoords='axes fraction', fontsize=20)

    ratio = data_hist / np.where(np.sum(bkg_hists, axis=0) <= 0, np.nan, np.sum(bkg_hists, axis=0))
    ax2.errorbar(bin_edges[:-1] + np.diff(bin_edges) / 2, ratio, yerr=np.sqrt(data_hist) / np.sum(bkg_hists, axis=0), fmt='o')
    ax2.axhline(y=1, color='k', linestyle='--')
    ax2.set_ylim(0, 2)
    ax2.set_xlim(100, 180)
    ax2.set_xlabel("Higgs Mass")
    ax2.set_ylabel("Data / Bkg.")
    ax2.grid()

    plt.tight_layout()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/2J/compare_data_mc_sr_hmass_{bounds[0]}_{bounds[1]}.png")
    
# Plot
for i in range(0, len(boundaries)-1):
    draw(data_sr, bkg_srs, bkg_list, [boundaries[i], boundaries[i+1]])