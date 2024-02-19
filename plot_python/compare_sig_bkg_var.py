import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
mc_legend = ["SM ZG", "DYJets", "LWK ZG", "TT", "Diboson"]
sig_legend = []
channel = "test"
var = "delta_phi_zgjj"
bins = 50
x_range = (0., 6.3)
ratio = 500
x_title = "#Delta#phi_{Z#gammajj}"
y_title = "Events/{:.3f}".format((x_range[1]-x_range[0])/bins)
sub_y_title = "Sigx{}/Bkg".format(ratio)
selections = []

# Dataset list
sig_file_list = [
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/ggH.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/VBF.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/WminusH.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/WplusH.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/ZH.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/ttH.root"
]
mc_file_list = [
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/ZGToLLG.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/DYJetsToLL.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/LLAJJ.root",
    "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/TT.root",
    ["/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/WW.root", "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/WZ.root", "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/two_jet/ZZ.root"]
]

# sig_file_list = [
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ggH/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/VBF/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/WminusH/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/WplusH/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ZH/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ttH/2017.root"
# ]
# mc_file_list = [
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ZGToLLG/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/DYJetsToLL/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/LLAJJ/2017.root",
#     "/eos/home-j/jiehan/root/2017/skimmed_ntuples/TT/2017.root",
#     ["/eos/home-j/jiehan/root/2017/skimmed_ntuples/WW/2017.root", "/eos/home-j/jiehan/root/2017/skimmed_ntuples/WZ/2017.root", "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ZZ/2017.root"]
# ]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "bkg sig Comparison")
# pads = plot.TwoPadSplit(0.29, 0.005, 0.005)
pad = plot.OnePad()[0]

print("============================")
print("Finish setting picture style")
print("============================")

# get hists
sig_yields = []
sig_yield = 0
for sig in sig_file_list:
    arrays = pic.read_file(sig, var, channel, selections)
    if "sig_hist" in globals():
        sig_hist, yields, yield_sum = pic.get_hist(arrays, var, 1, "sig", bins, x_range, sig_hist)
    else:
        sig_hist, yields, yield_sum = pic.get_hist(arrays, var, 1, "sig", bins, x_range)
    sig_yields.append(yields)
sig_yield = yield_sum
# print(sig_yield)
sig_hist.Scale(1/sig_yield)
    

h_stack = ROOT.THStack()
mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
mc_yields, mc_hists = [], []

for i, bkg in enumerate(mc_file_list):
    if isinstance(bkg, list):
        for sub_bkg in bkg:
            arrays = pic.read_file(sub_bkg, var, channel, selections)
            if file_hist in globals():
                file_hist, yields, yield_sum = pic.get_hist(arrays, var, 1, "mc_{}".format(i), bins, x_range, file_hist)
            else:
                file_hist, yields, yield_sum = pic.get_hist(arrays, var, 1, "mc_{}".format(i), bins, x_range)
    else: 
        arrays = pic.read_file(bkg, var, channel, selections)
        file_hist, yields, yield_sum = pic.get_hist(arrays, var, 1, "mc_{}".format(i), bins, x_range)
    mc_hists.append(file_hist)
    mc_yields.append(yields)
mc_yield = sum(mc_yields)
# print(mc_yield)
for i, hist in enumerate(mc_hists):
    hist.Scale(1./mc_yield)
    mc_hist.Add(hist)
    plot.Set(hist, LineWidth=0, FillColor=ROOT.TColor.GetColorDark(i+2))
    h_stack.Add(hist)

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

rp = pic.get_ratio_hist(sig_hist, mc_hist, (0,1))
rp_err = pic.get_ratio_hist(mc_hist, mc_hist, (0,1))

c1.Update()
# pads[0].cd()
# plot.Set(pads[0], Logy=1)
# plot.Set(pad, Logy=0)
h_stack.Draw("hist")
# plot.Set(h_stack.GetXaxis(), LabelSize=0)
plot.Set(h_stack.GetXaxis(), Title=x_title)
plot.Set(h_stack.GetYaxis(), Title=y_title)
sig_hist.Draw("E same")
plot.Set(sig_hist, MarkerStyle=8)

# pads[0].Update()
plot.Set(h_stack, Maximum=1.5*pad.GetFrame().GetY2())

legend = plot.PositionedLegend(0.6, 0.10, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.023)
legend.AddEntry("sig", "sig({:.2f})".format(sum(sig_yields)), "lep")
for i in range(len(mc_file_list)):
    legend.AddEntry("mc_{}".format(i), mc_legend[i]+"({:.2f})".format(mc_yields[i]), "f")
legend.Draw()

print("========================")
print("Finish drawing upper pad")
print("========================")

# pads[1].cd()
# rp_err.Draw("E2 same")
# rp.Draw("E same")
# plot.Set(rp_err, FillColorAlpha=(ROOT.kGray, 1))
# line = ROOT.TLine()
# plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
# plot.DrawHorizontalLine(pads[1], line, 1)
# plot.Set(rp, MarkerStyle=8)
# rp_err.GetXaxis().SetTitle(x_title)
# rp_err.GetYaxis().SetTitle(sub_y_title)
# pads[1].Update()

# print("========================")
# print("Finish drawing lower pad")
# print("========================")

plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55);
plot.DrawCMSLogo(c1, "41.5 fb^{-1} (13 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45);

print("========================")
print("Finish adding CMS word")
print("========================")

c1.Draw()
c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/sig_bkg_{}_{}.png".format(channel, var))
