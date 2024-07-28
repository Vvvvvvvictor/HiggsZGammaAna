import numpy as np
from scipy.optimize import curve_fit
import utils as ut
import uproot
import ROOT
from pdb import set_trace

log_path = "/eos/home-j/jiehan/root/outputs/"
data_driven_path = "/eos/home-j/jiehan/data_for_norm_float_v1/output/"
true_bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]
data_driven_bkg_list = ["Data", "ZGToLLG", "ZG2JToG2L2J"]
sig_list = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
data_list = ["Data"]
channels = ["zero_to_one_jet"]
com_sel = []
target_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"

xmin, xmax, nbin = 100, 180, 80
f = ROOT.TFile(target_path, "RECREATE")

for chan in channels:
    result = ut.GetBDTBinary(log_path+"significances/bin_binaries_1D_{}.txt".format(chan))
    num_bin = int(result[0])
    binaries = result[1:num_bin+2]

    for cat in range(num_bin):
        selections = com_sel + ["(bdt_score>{}) & (bdt_score<{})".format(binaries[cat], binaries[cat+1])]
        
        name = "data_full_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree=chan, selections = selections)
            data_hist, yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)

        name = "data_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        data_hist.Sumw2()
        data_hist_py = np.zeros(nbin)
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree=chan, selections=selections + ["(H_mass<122) | (H_mass>128)"])
            data_hist, data_yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)
        data_hist_py = np.zeros(nbin)
        for i in range(1, nbin + 1):
            data_hist_py[i-1] = data_hist.GetBinContent(i)
        
        name = "data_zg_enriched_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree=chan, selections = selections + ["(H_mass<122) | (H_mass>128)", "gamma_relpt>0.3"])
            data_hist, data_yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)
        
        name = "data_dy_enriched_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree=chan, selections = selections + ["(H_mass<122) | (H_mass>128)", "gamma_relpt<=0.3"])
            data_hist, data_yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)
        
        name = "sig_{}_cat{}".format(chan, cat)
        sig_hist = ROOT.TH1D(name,name, 2*nbin, xmin, xmax)
        sig_hist.Sumw2()
        # for i, sig in enumerate(sig_list):
        #     arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, sig), selections = selections)
        #     sig_hist, yield_h = ut.AddHist(arrays, "H_mass", sig_hist)
        arrays = ut.ReadFile(log_path+"{}/signal.root".format(chan), tree="zero_to_one_jet", selections = selections)
        sig_hist, yield_h = ut.AddHist(arrays, "H_mass", sig_hist)
        sig_hist.Write(name)

        name = "bkg_{}_cat{}".format(chan, cat)
        fake_bkg_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        true_bkg_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        fake_bkg_hist.Sumw2()
        true_bkg_hist.Sumw2()
        
        true_bkg_sb = 0
        for i, bkg in enumerate(true_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), tree=chan, selections = selections)
            true_bkg_sb += np.sum(arrays.query("H_mass<122 | H_mass>128")["weight"])
            true_bkg_hist, true_yield_h = ut.AddHist(arrays, "H_mass", true_bkg_hist)
            
        fake_bkg_sb = 0
        for i, bkg in enumerate(data_driven_bkg_list):
            arrays = ut.ReadFile(data_driven_path+"{}/{}_fake.root".format(chan, bkg), tree=chan, selections = selections)
            if bkg != "Data":
                arrays["weight"] = arrays["weight"] * (-1)
            fake_bkg_sb += np.sum(arrays.query("H_mass<122 | H_mass>128")["weight"])
            fake_bkg_hist, fake_yield_h = ut.AddHist(arrays, "H_mass", fake_bkg_hist)
        fake_bkg_hist.Scale((data_yield_h - true_bkg_sb)/fake_bkg_sb)
        
        bkg_hist = true_bkg_hist.Clone()
        bkg_hist.Add(fake_bkg_hist)

        bkg_hist_py = np.zeros(nbin)
        for i in range(1, nbin + 1):
            bkg_hist_py[i-1] = bkg_hist.GetBinContent(i)

        ratio = data_hist_py / bkg_hist_py

        from scipy.optimize import curve_fit
        def linear(x, a, b, c):
            return a*x**2 + b*x + c
        
        bin_edges = np.linspace(xmin, xmax, nbin+1)
        x = bin_edges[:-1]+np.diff(bin_edges)/2
        mask = np.where((x<122) | (x>128))
        y = ratio[mask]
        yerr = 1/data_hist_py[mask]
        x = x[mask]
        popt, pcov = curve_fit(linear, x-100, y, sigma=yerr, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        print(popt, perr)

        x = bin_edges[:-1]+np.diff(bin_edges)/2
        bkg_hist_py_mod = bkg_hist_py * linear(x-100, *popt)
        for i in range(1, nbin + 1):
            bkg_hist.SetBinContent(i, bkg_hist_py_mod[i-1])
            
        bkg_hist.Write(name)
        
        true_bkg_hist.Scale(1./true_bkg_hist.Integral())
        true_bkg_hist.Write("true_"+name)
        fake_bkg_hist.Scale(1./fake_bkg_hist.Integral())
        fake_bkg_hist.Write("fake_"+name)
        
        name = "fake_bkg_zg_enriched_{}_cat{}".format(chan, cat)
        fake_bkg_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        fake_bkg_hist.Sumw2()
        for i, bkg in enumerate(data_driven_bkg_list):
            arrays = ut.ReadFile(data_driven_path+"{}/{}_fake.root".format(chan, bkg), tree=chan, selections = selections + ["gamma_relpt>0.3"])
            if bkg != "Data":
                arrays["weight"] = arrays["weight"] * (-1)
            fake_bkg_hist, fake_yield_h = ut.AddHist(arrays, "H_mass", fake_bkg_hist)
        fake_bkg_hist.Scale(1./fake_bkg_hist.Integral())
        fake_bkg_hist.Write(name)
        
        name = "fake_bkg_dy_enriched_{}_cat{}".format(chan, cat)
        fake_bkg_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        fake_bkg_hist.Sumw2()
        for i, bkg in enumerate(data_driven_bkg_list):
            arrays = ut.ReadFile(data_driven_path+"{}/{}_fake.root".format(chan, bkg), tree=chan, selections = selections + ["gamma_relpt<=0.3"])
            if bkg != "Data":
                arrays["weight"] = arrays["weight"] * (-1)
            fake_bkg_hist, fake_yield_h = ut.AddHist(arrays, "H_mass", fake_bkg_hist)
        fake_bkg_hist.Scale(1./fake_bkg_hist.Integral())
        fake_bkg_hist.Write(name)
        
        name = "true_bkg_zg_enriched_{}_cat{}".format(chan, cat)
        true_bkg_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        true_bkg_hist.Sumw2()
        for i, bkg in enumerate(true_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), tree=chan, selections = selections + ["gamma_relpt>0.3"])
            true_bkg_hist, true_yield_h = ut.AddHist(arrays, "H_mass", true_bkg_hist)
        true_bkg_hist.Scale(1./true_bkg_hist.Integral())
        true_bkg_hist.Write(name)
        
        name = "true_bkg_dy_enriched_{}_cat{}".format(chan, cat)
        true_bkg_hist = ROOT.TH1D(name,name, nbin, xmin, xmax)
        true_bkg_hist.Sumw2()
        for i, bkg in enumerate(true_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), tree=chan, selections = selections + ["gamma_relpt<=0.3"])
            true_bkg_hist, true_yield_h = ut.AddHist(arrays, "H_mass", true_bkg_hist)
        true_bkg_hist.Scale(1./true_bkg_hist.Integral())
        true_bkg_hist.Write(name)

print("We get background and signal template: {}".format(target_path))

