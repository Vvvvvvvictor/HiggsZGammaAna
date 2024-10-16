import utils as ut
import uproot
import ROOT
from pdb import set_trace

log_path = "/eos/home-j/jiehan/root/outputs/"
bkg_list = ["ZGToLLG", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
sig_list = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
data_list = ["Data"]
# channels = ["zero_jet", "one_jet", "two_jet", "VH_ttH"]
channels = ["zero_to_one_jet", "two_jet"]
com_sel = []
target_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"

f = ROOT.TFile(target_path, "RECREATE")

xmin, xmax, xbin = 80, 180, 100

for chan in channels:
    result = ut.GetBDTBinary(log_path+"significances/bin_binaries_1D_{}.txt".format(chan))
    num_bin = int(result[0])
    binaries = result[1:num_bin+2]

    for cat in range(num_bin):
        selections = com_sel + ["(bdt_score_t>{}) & (bdt_score_t<{})".format(binaries[cat], binaries[cat+1])]

        name = "bkg_{}_cat{}".format(chan, cat)
        bkg_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        bkg_hist.Sumw2()
        for i, bkg in enumerate(bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), tree = chan, selections = selections)
            bkg_hist, yield_h = ut.AddHist(arrays, "H_mass", bkg_hist)
        bkg_hist.Write(name)

        name = "sig_{}_cat{}".format(chan, cat)
        sig_hist = ROOT.TH1D(name,name, 2*xbin, xmin, xmax)
        sig_hist.Sumw2()
        for i, sig in enumerate(sig_list):
            arrays = ut.ReadFile(log_path+"{}/{}_M125.root".format(chan, sig), tree = chan, selections = selections)
            sig_hist, yield_h = ut.AddHist(arrays, "H_mass", sig_hist)
        sig_hist.Write(name)

        name = "data_full_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree = chan, selections = selections)
            data_hist, yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)

        name = "data_{}_cat{}".format(chan, cat)
        data_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        data_hist.Sumw2()
        selections = selections + ["(H_mass<120) | (H_mass>130)"]
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree = chan, selections = selections)
            data_hist, yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)

print("We get background and signal template: {}".format(target_path))

