import utils as ut
import uproot
import numpy as np
import ROOT
from pdb import set_trace

log_path = "/eos/home-j/jiehan/root/outputs/"
data_driven_path = "/eos/home-j/jiehan/data_for_norm_float_v1/output/"
true_bkg_list = ["ZGToLLG", "ZG2JToG2L2J"]
data_driven_bkg_list = ["Data", "ZGToLLG", "ZG2JToG2L2J"]
fake_photon_bkg_list = ["DYJetsToLL", "WGToLNuG", "TTGJets", "TT", "TGJets", "WW", "WZ", "ZZ", "ttWJets", "ttZJets"]
sig_list = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
data_list = ["data"]
channels = ["two_jet"]
com_sel = []
target_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"

f = ROOT.TFile(target_path, "RECREATE")

for chan in channels:
    result = ut.GetBDTBinary(log_path+"significances/bin_binaries_1D_{}.txt".format(chan))
    num_bin = int(result[0])
    binaries = result[1:num_bin+2]

    for cat in range(num_bin):
        selections = com_sel + ["(bdt_score_t>{}) & (bdt_score_t<{})".format(binaries[cat], binaries[cat+1])]
        
        name = "sig_{}_cat{}".format(chan, cat)
        sig_hist = ROOT.TH1D(name,name, 320, 100, 180)
        sig_hist.Sumw2()
        for i, sig in enumerate(sig_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, sig), selections = selections)
            sig_hist, yield_h = ut.AddHist(arrays, "H_mass", sig_hist)
        sig_hist.Write(name)

        name = "data_full_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, 320, 100, 180)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), selections = selections)
            data_hist, yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)

        name = "data_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, 320, 100, 180)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), selections = selections + ["(H_mass<122) | (H_mass>128)"])
            data_hist, data_yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)
        
        name = "data_zg_enriched_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, 320, 100, 180)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), selections = selections + ["(H_mass<122) | (H_mass>128)", "gamma_relpt>0.4"])
            data_hist, data_yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)
        
        name = "data_dy_enriched_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, 320, 100, 180)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), selections = selections + ["(H_mass<122) | (H_mass>128)", "gamma_relpt<=0.4"])
            data_hist, data_yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)
        
        name = "bkg_{}_cat{}".format(chan, cat)
        fake_bkg_hist = ROOT.TH1D(name,name, 320, 100, 180)
        true_bkg_hist = ROOT.TH1D(name,name, 320, 100, 180)
        fake_bkg_hist.Sumw2()
        true_bkg_hist.Sumw2()
        
        true_bkg_sb = 0
        for i, bkg in enumerate(true_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg),  selections = selections)
            true_bkg_sb += np.sum(arrays.query("H_mass<122 | H_mass>128")["weight"])
            true_bkg_hist, true_yield_h = ut.AddHist(arrays, "H_mass", true_bkg_hist)
            
        fake_bkg = 0
        for i, bkg in enumerate(fake_photon_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), selections = selections)
            fake_bkg += np.sum(arrays["weight"])\
            
        data_driven_bkg_sb = 0
        data_driven_bkg = 0
        for i, bkg in enumerate(data_driven_bkg_list):
            arrays = ut.ReadFile(data_driven_path+"{}/{}_fake.root".format(chan, bkg), tree=chan, selections = selections)
            if bkg != "Data":
                arrays["weight"] = arrays["weight"] * (-1)
            data_driven_bkg_sb += np.sum(arrays.query("H_mass<122 | H_mass>128")["weight"])
            data_driven_bkg += np.sum(arrays["weight"])
            fake_bkg_hist, fake_yield_h = ut.AddHist(arrays, "H_mass", fake_bkg_hist)
        fake_sf = fake_bkg / data_driven_bkg
        fake_bkg_hist.Scale(fake_sf)
        data_driven_bkg_sb = data_driven_bkg_sb * fake_sf
        
        mc_sf = data_yield_h / (data_driven_bkg_sb + true_bkg_sb)
        fake_bkg_hist.Scale(mc_sf)
        true_bkg_hist.Scale(mc_sf)
        
        bkg_hist = true_bkg_hist.Clone()
        bkg_hist.Add(fake_bkg_hist)
            
        bkg_hist.Write(name)
        
        true_bkg_hist.Scale(1./true_bkg_hist.Integral())
        true_bkg_hist.Write("true_"+name)
        fake_bkg_hist.Scale(1./fake_bkg_hist.Integral())
        fake_bkg_hist.Write("fake_"+name)
        
        name = "fake_bkg_zg_enriched_{}_cat{}".format(chan, cat)
        fake_bkg_hist = ROOT.TH1D(name,name, 320, 100, 180)
        fake_bkg_hist.Sumw2()
        for i, bkg in enumerate(data_driven_bkg_list):
            arrays = ut.ReadFile(data_driven_path+"{}/{}_fake.root".format(chan, bkg), tree=chan, selections = selections + ["gamma_relpt>0.4"])
            if bkg != "Data":
                arrays["weight"] = arrays["weight"] * (-1)
            fake_bkg_hist, fake_yield_h = ut.AddHist(arrays, "H_mass", fake_bkg_hist)
        fake_bkg_hist.Scale(1./fake_bkg_hist.Integral())
        fake_bkg_hist.Write(name)
        
        name = "fake_bkg_dy_enriched_{}_cat{}".format(chan, cat)
        fake_bkg_hist = ROOT.TH1D(name,name, 320, 100, 180)
        fake_bkg_hist.Sumw2()
        for i, bkg in enumerate(data_driven_bkg_list):
            arrays = ut.ReadFile(data_driven_path+"{}/{}_fake.root".format(chan, bkg), tree=chan, selections = selections + ["gamma_relpt<=0.4"])
            if bkg != "Data":
                arrays["weight"] = arrays["weight"] * (-1)
            fake_bkg_hist, fake_yield_h = ut.AddHist(arrays, "H_mass", fake_bkg_hist)
        fake_bkg_hist.Scale(1./fake_bkg_hist.Integral())
        fake_bkg_hist.Write(name)
        
        name = "true_bkg_zg_enriched_{}_cat{}".format(chan, cat)
        true_bkg_hist = ROOT.TH1D(name,name, 320, 100, 180)
        true_bkg_hist.Sumw2()
        for i, bkg in enumerate(true_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), selections = selections + ["gamma_relpt>0.4"])
            true_bkg_hist, true_yield_h = ut.AddHist(arrays, "H_mass", true_bkg_hist)
        true_bkg_hist.Scale(1./true_bkg_hist.Integral())
        true_bkg_hist.Write(name)
        
        name = "true_bkg_dy_enriched_{}_cat{}".format(chan, cat)
        true_bkg_hist = ROOT.TH1D(name,name, 320, 100, 180)
        true_bkg_hist.Sumw2()
        for i, bkg in enumerate(true_bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), selections = selections + ["gamma_relpt<=0.4"])
            true_bkg_hist, true_yield_h = ut.AddHist(arrays, "H_mass", true_bkg_hist)
        true_bkg_hist.Scale(1./true_bkg_hist.Integral())
        true_bkg_hist.Write(name)

print("We get background and signal template: {}".format(target_path))

