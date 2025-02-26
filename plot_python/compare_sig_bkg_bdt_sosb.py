import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
ratio = 1
mc_legend = ["SM ZG", "DYJets", "EWK Z+Jets", "EWK ZG", "TT", "TTG+Jets", "TTVJets", "Diboson"]#"DYJets", "EWK Z+Jets", "EWK ZG", "TT", "TTG+Jets", "TTVJets", "Diboson", "SM ZG", "EWK ZG", "fake photon"]
sig_legend = ["sig", "ggH", "VBF"]
# path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/"
path = "/eos/home-j/jiehan/root/outputs/test/"
channel = "two_jet"
tree = "two_jet"
var = "bdt_score_t"
bins = 100
x_range = (0, 1)
blind_range = (120, 130)
x_title = "BDT score"
y_title = "Events/{:.2f}".format((x_range[1]-x_range[0])/bins)
sub_y_title = "S/#sqrt{B}"
selections = ["H_mass>120", "H_mass<130"]

# Dataset list
sig_file_list = [
    ["sig.root"],
    "ggH_M125.root",
    "VBF_M125.root"
]
mc_file_list = [
    "ZGToLLG.root",
    "DYJetsToLL.root",
    "EWKZ2J.root",
    "ZG2JToG2L2J.root",
    "TT.root",
    "TTGJets.root",
    ["ttWJets.root", "ttZJets.root"],
    ["WW.root", "WZ.root", "ZZ.root"]
    # "data_driven_bkg_v3.root"
]
data_file_list = [
    "data.root"
]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "bkg sig Comparison")
pads = plot.TwoPadSplit(0.29, 0.005, 0.005)

print("============================")
print("Finish setting picture style")
print("============================")

# get hists
# data_sb_yields, mc_sb_yields = 0, 0
# for i, bkg in enumerate(mc_file_list):
#     if isinstance(bkg, list):
#         for sub_bkg in bkg:
#             arrays = pic.read_root_file(path+channel+"/"+sub_bkg, var, tree, selections)
#             if file_hist in globals():
#                 file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range, file_hist)
#             else:
#                 file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
#     else: 
#         arrays = pic.read_root_file(path+channel+"/"+bkg, var, tree, selections)
#         file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
#     mc_sb_yields = mc_sb_yields + yields

# for data in data_file_list:
#     if isinstance(data, list):
#         for sub_data in data:
#             arrays = pic.read_root_file(path+channel+"/"+sub_data, var, tree, selections)
#             if file_hist in globals():
#                 file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range, file_hist)
#             else:
#                 file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
#     else: 
#         arrays = pic.read_root_file(path+channel+"/"+data, var, tree, selections)
#         file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
#     data_sb_yields = data_sb_yields + yields

# sb_ratio = data_sb_yields / mc_sb_yields

sb_ratio=1.

# selections+=["H_mass<130", "H_mass>120"]

sig_yields, sig_hist_list = [], []
for i, sig in enumerate(sig_file_list):
    if isinstance(sig, list):
        for sub_sig in sig:
            arrays, _ = pic.read_root_file(path+channel+"/"+sub_sig, var, tree, selections)
            if "sig_hist" in globals():
                sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range, sig_hist)
            else:
                sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range)
    else:
        arrays, _ = pic.read_root_file(path+channel+"/"+sig, var, tree, selections)
        sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range)
    plot.Set(sig_hist, LineWidth=2, LineStyle=i+1)
    sig_hist_list.append(sig_hist)
    sig_yields.append(yields)

h_stack = ROOT.THStack()
mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
mc_yields, mc_hist_list = [], []

for i, bkg in enumerate(mc_file_list):
    if isinstance(bkg, list):
        for sub_bkg in bkg:
            arrays, _ = pic.read_root_file(path+channel+"/"+sub_bkg, var, tree, selections)
            if file_hist in globals():
                file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range, file_hist)
            else:
                file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range)
    else: 
        arrays, _ = pic.read_root_file(path+channel+"/"+bkg, var, tree, selections)
        file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range)
    mc_hist.Add(file_hist)
    mc_hist_list.append(file_hist)
    mc_yields.append(yields)
for i, file_hist in enumerate(mc_hist_list):
    plot.Set(file_hist, LineWidth=0, FillColor=ROOT.TColor.GetColorDark(i+2))
    file_hist.Scale(1./sum(mc_yields))
    h_stack.Add(file_hist)

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

rp = pic.get_S_over_sqrtB(sig_hist_list[0], mc_hist, ratio, (0.,1.0005))

line = ROOT.TLine()
plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)

c1.Update()
pads[0].cd()
h_stack.Draw("hist")
plot.Set(h_stack.GetXaxis(), LabelSize=0)
plot.Set(h_stack.GetYaxis(), Title=y_title)
for sig_hist in sig_hist_list:
    sig_hist.Scale(1./sig_yields[0])
    sig_hist.Draw("hist same")
plot.Set(sig_hist, MarkerStyle=1, MarkerSize=3)

plot.Set(h_stack, Maximum=1.3*pads[0].GetFrame().GetY2())
plot.Set(h_stack, Minimum=0.00005)
# plot.Set(h_stack, Maximum=0.045)

# x1, x2 = 0.61, 0.8
# line.DrawLine(x1, 0, x1, 0.5*pads[0].GetFrame().GetY2())
# line.DrawLine(x2, 0, x2, 0.5*pads[0].GetFrame().GetY2())

plot.Set(pads[0], Logy=1)
pads[0].Update()

legend = plot.PositionedLegend(0.75, 0.15, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.025, TextFont=62, FillStyle=0)
for i in range(len(sig_hist_list)):
    legend.AddEntry("sig_{}".format(i), sig_legend[i]+"({:4.2f})".format(sig_yields[i]), "l")
for i in range(len(mc_file_list)):
    legend.AddEntry("mc_{}".format(i), mc_legend[i]+"({:4.2f})".format(mc_yields[i]), "f")
legend.Draw()

print("========================")
print("Finish drawing upper pad")
print("========================")

pads[1].cd()
rp.Draw("E")
boundaries = [0.07, 0.3, 0.61, 0.9]
for i in range(4):
    line.DrawLine(boundaries[i], 0, boundaries[i], 1)
# plot.Set(pads[1], Logy=1)
# line.DrawLine(0.29, 0, 0.29, 0.26)
# line.DrawLine(x1, 0, x1, 0.2)
# line.DrawLine(x2, 0, x2, 0.2)
# plot.DrawHorizontalLine(pads[1], line, 1)
# plot.Set(rp, MarkerStyle=7, MarkerSize=2)
rp.GetXaxis().SetTitle(x_title)
rp.GetYaxis().SetTitle(sub_y_title)
pads[1].Update()

print("========================")
print("Finish drawing lower pad")
print("========================")

plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55);
plot.DrawCMSLogo(c1, "137 fb^{-1} (13 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45);

print("========================")
print("Finish adding CMS word")
print("========================")

c1.Draw()
c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/sig_bkg_{}_{}.png".format(channel, var))
