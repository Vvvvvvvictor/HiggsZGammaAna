import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
ratio = 500
mc_legend = ["SM ZG", "DYJets", "LWK ZG", "TT", "Diboson"]
sig_legend = ["sig", "ggH", "VBF"]
path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/InputData/outputs/"
channel = "two_jet"
tree = "test"
var = "bdt_score_t"
bins = 100
x_range = (0, 1)
blind_range = (122, 128)
x_title = "BDT score"
y_title = "Events/{:.2f}".format((x_range[1]-x_range[0])/bins)
sub_y_title = "S/#sqrt{B^{MC}_{scale}}"
selections = []

# Dataset list
sig_file_list = [
    ["ggH.root", "VBF.root", "WminusH.root", "WplusH.root", "ZH.root", "ttH.root"],
    "ggH.root",
    "VBF.root"
]
mc_file_list = [
    "ZGToLLG.root",
    "DYJetsToLL.root",
    "LLAJJ.root",
    "TT.root",
    ["WW.root", "WZ.root", "ZZ.root"]
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
data_sb_yields, mc_sb_yields = 0, 0
for i, bkg in enumerate(mc_file_list):
    if isinstance(bkg, list):
        for sub_bkg in bkg:
            arrays = pic.read_file(path+channel+"/"+sub_bkg, var, tree, selections)
            if file_hist in globals():
                file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range, file_hist)
            else:
                file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
    else: 
        arrays = pic.read_file(path+channel+"/"+bkg, var, tree, selections)
        file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
    mc_sb_yields = mc_sb_yields + yields

for data in data_file_list:
    if isinstance(data, list):
        for sub_data in data:
            arrays = pic.read_file(path+channel+"/"+sub_data, var, tree, selections)
            if file_hist in globals():
                file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range, file_hist)
            else:
                file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
    else: 
        arrays = pic.read_file(path+channel+"/"+data, var, tree, selections)
        file_hist, _, yields = pic.get_hist_sb(arrays, var, 1, "mc_{}".format(i), bins, x_range, blind_range)
    data_sb_yields = data_sb_yields + yields

sb_ratio = data_sb_yields / mc_sb_yields

selections+=["H_mass<128", "H_mass>122"]

sig_yields, sig_hist_list = [], []
for i, sig in enumerate(sig_file_list):
    if isinstance(sig, list):
        for sub_sig in sig:
            arrays = pic.read_file(path+channel+"/"+sub_sig, var, tree, selections)
            if "sig_hist" in globals():
                sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range, sig_hist)
            else:
                sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range)
    else:
        arrays = pic.read_file(path+channel+"/"+sig, var, tree, selections)
        sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range)
    plot.Set(sig_hist, LineWidth=2, LineStyle=i+1)
    sig_hist_list.append(sig_hist)
    sig_yields.append(yields)

h_stack = ROOT.THStack()
mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
mc_yields = []

for i, bkg in enumerate(mc_file_list):
    if isinstance(bkg, list):
        for sub_bkg in bkg:
            arrays = pic.read_file(path+channel+"/"+sub_bkg, var, tree, selections)
            if file_hist in globals():
                file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range, file_hist)
            else:
                file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range)
    else: 
        arrays = pic.read_file(path+channel+"/"+bkg, var, tree, selections)
        file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range)
    mc_hist.Add(file_hist)
    mc_yields.append(yields*sb_ratio)
    plot.Set(file_hist, LineWidth=0, FillColor=ROOT.TColor.GetColorDark(i+2))
    h_stack.Add(file_hist)

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

rp = pic.get_S_over_sqrtB(sig_hist_list[0], mc_hist, ratio, (0.,0.26))

c1.Update()
pads[0].cd()
# plot.Set(pads[0], Logy=1)
h_stack.Draw("hist")
plot.Set(h_stack.GetXaxis(), LabelSize=0)
plot.Set(h_stack.GetYaxis(), Title=y_title)
for sig_hist in sig_hist_list:
    sig_hist.Draw("hist same")
plot.Set(sig_hist, MarkerStyle=1, MarkerSize=3)

pads[0].Update()
# plot.Set(h_stack, Maximum=1.3*pads[0].GetFrame().GetY2())

legend = plot.PositionedLegend(0.65, 0.10, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.023)
for i in range(len(sig_hist_list)):
    legend.AddEntry("sig_{}".format(i), sig_legend[i]+"x{:.0f}({:.2f})".format(ratio, sig_yields[i]), "l")
for i in range(len(mc_file_list)):
    legend.AddEntry("mc_{}".format(i), mc_legend[i]+"({:.2f})".format(mc_yields[i]), "f")
legend.Draw()

print("========================")
print("Finish drawing upper pad")
print("========================")

pads[1].cd()
rp.Draw("E")
# plot.Set(pads[1], Logy=1)
line = ROOT.TLine()
plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
line.DrawLine(0.29, 0, 0.29, 0.26)
line.DrawLine(0.57, 0, 0.57, 0.26)
line.DrawLine(0.73, 0, 0.73, 0.26)
# plot.DrawHorizontalLine(pads[1], line, 1)
# plot.Set(rp, MarkerStyle=7, MarkerSize=2)
rp.GetXaxis().SetTitle(x_title)
rp.GetYaxis().SetTitle(sub_y_title)
pads[1].Update()

print("========================")
print("Finish drawing lower pad")
print("========================")

plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55);
plot.DrawCMSLogo(c1, "41.5 fb^{-1} (13 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45);

print("========================")
print("Finish adding CMS word")
print("========================")

c1.Draw()
c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/sig_bkg_{}_{}.png".format(channel, var))
