import utils as ut
import uproot
import ROOT
from pdb import set_trace

log_path = "/eos/home-j/jiehan/root/outputs/test/"
bkg_list = ["ZGToLLG", "EWKZ2J", "DYJetsToLL"]# "ZG2JToG2L2J", "TT", "TTGJets", "TGJets", "WZ", "ZZ", "ttZJets", "ttWJets"]
sig_list = ["ggH", "VBF"]#, "WminusH", "WplusH", "ZH", "ttH"]
data_list = ["Data"]
# channels = ["zero_jet", "one_jet", "two_jet", "VH_ttH"]
channels = ["zero_to_one_jet", "two_jet"]
names = ["ggH", "VBF"]
com_sel = []
target_path = "/eos/user/j/jiehan/root/input_finalfit/templates/template_all.root"

f = ROOT.TFile(target_path, "RECREATE")

xmin, xmax, xbin = 95, 180, 340

for j, chan in enumerate(channels):
    with open(log_path+"significances/bin_boundaries_1D_{}.txt".format(chan)) as log_file:
        boundaries = log_file.readlines()[0].split(" ")
        num_bin = len(boundaries)-1
        boundaries = [float(i) for i in boundaries]

    for cat in range(num_bin):
        selections = com_sel + ["(bdt_score_t>{}) & (bdt_score_t<{})".format(boundaries[cat], boundaries[cat+1])]

        name = "bkg_all_{}{}".format(names[j], cat)
        bkg_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        bkg_hist.Sumw2()
        for i, bkg in enumerate(bkg_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, bkg), tree = chan, selections = selections)
            bkg_hist, yield_h = ut.AddHist(arrays, "H_mass", bkg_hist)
        bkg_hist.Write(name)

        name = "sig_all_{}{}".format(names[j], cat)
        sig_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        sig_hist.Sumw2()
        for i, sig in enumerate(sig_list):
            arrays = ut.ReadFile(log_path+"{}/{}_M125.root".format(chan, sig), tree = chan, selections = selections)
            sig_hist, yield_h = ut.AddHist(arrays, "H_mass", sig_hist)
        sig_hist.Write(name)

        name = "data_full_all_{}{}".format(names[j], cat)
        data_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        data_hist.Sumw2()
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree = chan, selections = selections)
            data_hist, yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)

        name = "data_all_{}{}".format(names[j], cat)
        data_hist = ROOT.TH1D(name,name, xbin, xmin, xmax)
        data_hist.Sumw2()
        selections = selections + ["(H_mass<120) | (H_mass>130)"]
        for i, data in enumerate(data_list):
            arrays = ut.ReadFile(log_path+"{}/{}.root".format(chan, data), tree = chan, selections = selections)
            data_hist, yield_h = ut.AddHist(arrays, "H_mass", data_hist)
        data_hist.Write(name)

print("We get background and signal template: {}".format(target_path))

