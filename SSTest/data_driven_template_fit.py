import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import mplhep as hep
plt.style.use(hep.style.CMS)

# 打开 ROOT 文件
file = ROOT.TFile("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root")
channel = "two_jet"

# 将 histograms 转换为 numpy 数组
def hist_to_array(hist):
    nbins = hist.GetNbinsX()
    values = np.array([hist.GetBinContent(i) for i in range(1, nbins + 1)])
    errors = np.array([hist.GetBinError(i) for i in range(1, nbins + 1)])
    return values, errors

for i in range(4):
    # 读取 histograms
    data_zg_hist = file.Get(f"data_zg_enriched_{channel}_cat{i}")
    fake_bkg_zg_hist = file.Get(f"fake_bkg_zg_enriched_{channel}_cat{i}")
    true_bkg_zg_hist = file.Get(f"true_bkg_zg_enriched_{channel}_cat{i}")
    
    data_dy_hist = file.Get(f"data_dy_enriched_{channel}_cat{i}")
    fake_bkg_dy_hist = file.Get(f"fake_bkg_dy_enriched_{channel}_cat{i}")
    true_bkg_dy_hist = file.Get(f"true_bkg_dy_enriched_{channel}_cat{i}")

    data_zg, data_zg_err = hist_to_array(data_zg_hist)
    fake_bkg_zg, fake_bkg_zg_err = hist_to_array(fake_bkg_zg_hist)
    true_bkg_zg, true_bkg_zg_err = hist_to_array(true_bkg_zg_hist)
    
    data_dy, data_dy_err = hist_to_array(data_dy_hist)
    fake_bkg_dy, fake_bkg_dy_err = hist_to_array(fake_bkg_dy_hist)
    true_bkg_dy, true_bkg_dy_err = hist_to_array(true_bkg_dy_hist)

    # 获取 bins 的范围
    bins = np.arange(data_zg_hist.GetNbinsX())

    # 定义拟合的区间
    fit_mask = ((bins < 22) | (bins >= 28))
    
    # 定义拟合函数，两个参数是bkg1和bkg2的归一化值
    def fit_function(x, norm1, norm2):
        return np.concatenate((norm1 * fake_bkg_zg[fit_mask] + norm2 * true_bkg_zg[fit_mask], norm1 * fake_bkg_dy[fit_mask] + norm2 * true_bkg_dy[fit_mask]))

    x = np.concatenate((bins[fit_mask], bins[fit_mask]))
    data = np.concatenate((data_zg[fit_mask], data_dy[fit_mask]))

    # 初始猜测的归一化值，可以根据实际情况调整
    initial_guess = [10, 10]

    # 使用curve_fit来拟合
    params, params_covariance = curve_fit(fit_function, x, data, p0=initial_guess, bounds=([0.0, 0.0], [np.inf, np.inf]))

    # 获取拟合的归一化值
    norm1, norm2 = params

    # 计算拟合的histogram
    fit_hist = fit_function(x, norm1, norm2)
    print(1-fit_hist/data)

    # 画图
    plt.figure()
    fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    x = bins + 100
    ax1.errorbar(x, data_zg+data_dy, label='Data', fmt='o', yerr=np.sqrt(data_zg_err**2 + data_dy_err**2))
    ax1.errorbar(x, (fake_bkg_zg+fake_bkg_dy) * norm1, label=f'fake BKG({norm1:.2f})', fmt='o', yerr=np.sqrt(fake_bkg_zg_err**2 + fake_bkg_dy_err**2)*norm1)
    ax1.errorbar(x, (true_bkg_zg+true_bkg_dy) * norm2, label=f'true BKG({norm2:.2f})', fmt='o', yerr=np.sqrt(true_bkg_zg_err**2 + true_bkg_dy_err**2)*norm2)
    ax1.errorbar(x, (fake_bkg_zg+fake_bkg_dy) * norm1 + (true_bkg_zg+true_bkg_dy) * norm2, label='Fitted', fmt='o', yerr=np.sqrt((fake_bkg_zg_err**2 + fake_bkg_dy_err**2)*norm1**2 + (true_bkg_zg_err**2 + true_bkg_dy_err**2)*norm2**2))
    ax1.legend()
    ax1.set_xlabel("H_mass")
    ax1.set_ylabel("Events")
    ax1.set_title(f"Category {i}")
    ax1.set_xlim(100, 180)
    
    ratio = (data_zg + data_dy) / ((fake_bkg_zg + fake_bkg_dy) * norm1 + (true_bkg_zg + true_bkg_dy) * norm2)
    ratio_err = np.sqrt((data_zg_err**2 + data_dy_err**2) / ((fake_bkg_zg + fake_bkg_dy) * norm1 + (true_bkg_zg + true_bkg_dy) * norm2) + ((data_zg + data_dy) * (fake_bkg_zg_err**2 + fake_bkg_dy_err**2) * norm1**2 + (data_zg + data_dy) * (true_bkg_zg_err**2 + true_bkg_dy_err**2) * norm2**2) / ((fake_bkg_zg + fake_bkg_dy) * norm1 + (true_bkg_zg + true_bkg_dy) * norm2)**2)
    ax2.errorbar(x, ratio, yerr=ratio_err, fmt='o')
    ax2.axhline(1, color='r', linestyle='--')
    ax2.set_xlabel("H_mass")
    ax2.set_ylabel("Data / Fitted")
    ax2.set_ylim(0., 2.)
    ax2.set_xlim(100, 180)
    
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/figure/fit_{i}.png")
    plt.clf()

    print(f"Fitted normalization for fake: {norm1}")
    print(f"Fitted normalization for true: {norm2}")

    # 画图：只包括zg部分
    plt.figure()
    fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    x = bins + 100
    ax1.errorbar(x, data_zg, label='Data', fmt='o', yerr=data_zg_err)
    ax1.errorbar(x, fake_bkg_zg * norm1, label=f'fake BKG({norm1:.2f})', fmt='o', yerr=fake_bkg_zg_err*norm1)
    ax1.errorbar(x, true_bkg_zg * norm2, label=f'true BKG({norm2:.2f})', fmt='o', yerr=true_bkg_zg_err*norm2)

    ax1.errorbar(x, fake_bkg_zg * norm1 + true_bkg_zg * norm2, label='Fitted', fmt='o', yerr=np.sqrt(fake_bkg_zg_err**2*norm1**2 + true_bkg_zg_err**2*norm2**2))
    ax1.legend()
    ax1.set_xlabel("H_mass")
    ax1.set_ylabel("Events")
    ax1.set_title(f"Category {i}")
    ax1.set_xlim(100, 180)

    ratio = data_zg / (fake_bkg_zg * norm1 + true_bkg_zg * norm2)
    ratio_err = np.sqrt(data_zg_err**2 / (fake_bkg_zg * norm1 + true_bkg_zg * norm2) + (data_zg * fake_bkg_zg_err**2 * norm1**2 + data_zg * true_bkg_zg_err**2 * norm2**2) / (fake_bkg_zg * norm1 + true_bkg_zg * norm2)**2)

    ax2.errorbar(x, ratio, yerr=ratio_err, fmt='o')
    ax2.axhline(1, color='r', linestyle='--')
    ax2.set_xlabel("H_mass")
    ax2.set_ylabel("Data / Fitted")
    ax2.set_ylim(0., 2.)
    ax2.set_xlim(100, 180)

    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/figure/fit_zg_{i}.png")
    plt.clf()

    print(f"Fitted normalization for fake: {norm1}")
    print(f"Fitted normalization for true: {norm2}")

    # 画图：只包括dy部分
    plt.figure()
    fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    x = bins + 100
    ax1.errorbar(x, data_dy, label='Data', fmt='o', yerr=data_dy_err)
    ax1.errorbar(x, fake_bkg_dy * norm1, label=f'fake BKG({norm1:.2f})', fmt='o', yerr=fake_bkg_dy_err*norm1)
    ax1.errorbar(x, true_bkg_dy * norm2, label=f'true BKG({norm2:.2f})', fmt='o', yerr=true_bkg_dy_err*norm2)
    ax1.errorbar(x, fake_bkg_dy * norm1 + true_bkg_dy * norm2, label='Fitted', fmt='o', yerr=np.sqrt(fake_bkg_dy_err**2*norm1**2 + true_bkg_dy_err**2*norm2**2))

    ax1.legend()
    ax1.set_xlabel("H_mass")
    ax1.set_ylabel("Events")
    ax1.set_title(f"Category {i}")
    ax1.set_xlim(100, 180)
    
    ratio = data_dy / (fake_bkg_dy * norm1 + true_bkg_dy * norm2)
    ax2.errorbar(x, ratio, fmt='o')
    ax2.axhline(1, color='r', linestyle='--')
    ax2.set_xlabel("H_mass")
    ax2.set_ylabel("Data / Fitted")
    ax2.set_ylim(0., 2.)
    ax2.set_xlim(100, 180)

    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/figure/fit_dy_{i}.png")
    plt.clf()

    print(f"Fitted normalization for fake: {norm1}")
    print(f"Fitted normalization for true: {norm2}")

# 关闭文件
file.Close()
