# from ROOT import *
from pdb import set_trace

# target_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/sig_template.root"

# f = TFile(target_path, "RECREATE")
# for i in range(4):
#     file_sig = TFile.Open(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/HZGamma_data_sig_Hm125_workspace_cat{i}.root")

#     inWS_sig = file_sig.Get('CMS_hzg_workspace')
#     sig = inWS_sig.data(f"ggh_125_13TeV_untag_cat{i}_hist")
    
#     # set_trace()
#     vars = sig.get()
#     x = vars.first()
#     nbins = x.getBins()
#     xlow = x.getMin()
#     xhigh = x.getMax()
#     sig_hist = TH1D(f"sig_zero_to_one_jet_cat{i}", f"sig_zero_to_one_jet_cat{i}", nbins, xlow, xhigh)
#     for i in range(1, nbins+1):
#         x.setVal(sig_hist.GetBinCenter(i))
#         sig_hist.SetBinContent(i, sig.weight(vars))
#     file_sig.Close()
#     sig_hist.Write(f"sig_zero_to_one_jet_cat{i}")

import uproot

data = uproot.open(" /eos/home-j/jiehan/root/outputs/zero_to_one_jet/signal.root")['zero_to_one_jet'].arrays(library='pd')
data["bdt_score"] = data["BDT_score"]
with uproot.recreate("/eos/home-j/jiehan/root/outputs/zero_to_one_jet/signal.root") as f:
    f["zero_to_one_jet"] = data
