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

data_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/Data_fake.root")["zero_to_one_jet"].arrays(library="pd")
zg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/ZGToLLG_fake.root")["zero_to_one_jet"].arrays(library="pd")
ewkzg_cr = uproot.open(f"/eos/home-j/jiehan/data_for_norm_float_v1/output/zero_to_one_jet/ZG2JToG2L2J_fake.root")["zero_to_one_jet"].arrays(library="pd")
    
# boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_zero_to_one_jet.txt", "r").readlines()[2].split()[1:]]
boundaries = [0.28379330039024353, 0.4557725191116333, 0.5796570777893066, 0.7069960236549377, 1.]
print(boundaries)

################################################################
# Compare true photon samples plus data-driven photon samples with data in signal region
################################################################


bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]#, "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
def apply(data_cr, zg_cr, ewkzg_cr, bkg_list, bounds):
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

    # cr_hist = (data_cr_hist - zg_cr_hist)
    cr_hist = (data_cr_hist - zg_cr_hist - ewkzg_cr_hist)
  
    plt.figure()
    bkg_hists = []
    bkg_sum = 0
    for bkg_cr in bkg_samples:
        bkg_cr_hist, bin_edges = np.histogram(bkg_cr["H_mass"], bins=80, weights=bkg_cr["weight"], range=[100, 180])
        bkg_hists.append(bkg_cr_hist)
        bkg_sum = bkg_sum + np.sum(bkg_cr_hist[np.where((bin_edges[1:]<120) | (bin_edges[:-1]>130))])
    
    data_hist, bin_edges = np.histogram(data_sr[(data_sr["H_mass"]<120) | (data_sr["H_mass"]>130)]["H_mass"], bins=80, weights=data_sr[(data_sr["H_mass"]<120) | (data_sr["H_mass"]>130)]["weight"], range=[100, 180])
    cr_side_band = cr_hist[np.where((bin_edges[1:]<120) | (bin_edges[:-1]>130))]
    sf = (np.sum(data_hist)-bkg_sum)/np.sum(cr_side_band)

    bkg_hist_sum = bkg_hists + [cr_hist*sf]
    ratio = (data_hist / np.sum(bkg_hist_sum, axis=0))

    # fit this ratio distribution with linear function
    from scipy.optimize import curve_fit
    def linear(x, a, b, c):
        return a*x**2 + b*x + c
    
    x = bin_edges[:-1]+np.diff(bin_edges)/2
    y = ratio[np.where((x<120) | (x>130))]
    yerr = 1/data_hist[np.where((x<120) | (x>130))]
    x = x[np.where((x<120) | (x>130))]
    popt, pcov = curve_fit(linear, x-100, y, sigma=yerr, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    print(popt, perr)

    # correct the bkg_hist_sum and ratio with the fitted function
    bkg_hist_sum = [bkg_hist_sum[i] * linear(bin_edges[:-1]+np.diff(bin_edges)/2-100, *popt) for i in range(len(bkg_hist_sum))]

    data_cr['weight'] = data_cr['weight'] * sf * linear(data_cr["H_mass"], *popt)
    zg_cr['weight'] = zg_cr['weight'] * sf * linear(zg_cr["H_mass"], *popt)
    ewkzg_cr['weight'] = ewkzg_cr['weight'] * sf * linear(ewkzg_cr["H_mass"], *popt)

    output = pd.concat([data_cr, zg_cr, ewkzg_cr])

    return output

# Apply
with uproot.recreate("/eos/home-j/jiehan/root/outputs/zero_to_one_jet/data_driven_bkg.root") as output_file:
    output_dataframe = pd.DataFrame()
    for i in range(0, len(boundaries)-1):
        partial_output = apply(data_cr, zg_cr, ewkzg_cr, bkg_list, [boundaries[i], boundaries[i+1]])
        output_dataframe = pd.concat([output_dataframe, partial_output])
    output_file['test'] = output_dataframe
