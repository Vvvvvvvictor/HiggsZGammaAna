import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
mc_legend = [
    "SM ZG", "DYJets", "EWK Z(G)", "TT(G)", "Diboson", "WG", "ttVJets", "TGJets"
    ]
sig_legend = ["Signal"]
path = "/eos/home-j/jiehan/root/outputs/"
channel = "two_jet"
tree = "test"
var = "Z_pt"
bins = 100
x_range = (0, 50)
ratio = 1
x_title = "pt of ll"
y_title = "Events/{:.2f}".format((x_range[1]-x_range[0])/bins)
sub_y_title = "Sig/Bkg"
selections = []

# Dataset list
sig_file_list = [
    ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
]
mc_file_list = [
    "ZGToLLG",
    "DYJetsToLL",
    ["ZG2JToG2L2J", "EWKZ2J"],
    ["TT", "TTGJets"],
    ["WW", "WZ", "ZZ"],
    "WGToLNuG",
    ["ttWJets", "ttZJets"],
    "TGJets"
]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "bkg sig Comparison")
pads = plot.TwoPadSplit(0.29, 0.005, 0.005)

print("============================")
print("Finish setting picture style")
print("============================")

sb_ratio=1.

# get hists
sig_yields, sig_hist_list = [], []
for i, sig in enumerate(sig_file_list):
    if isinstance(sig, list):
        for sub_sig in sig:
            arrays = pic.read_file(path+channel+"/"+sub_sig+".root", var, tree, selections)
            if "sig_hist" in globals():
                sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range, sig_hist)
            else:
                sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range)
    else:
        arrays = pic.read_file(path+channel+"/"+sig+".root", var, tree, selections)
        sig_hist, _, yields = pic.get_hist(arrays, var, ratio, "sig_{}".format(i), bins, x_range)
    plot.Set(sig_hist, LineWidth=2, LineStyle=i+1)
    sig_hist_list.append(sig_hist)
    sig_yields.append(yields)
print(sum(sig_yields))
sig_hist.Scale(1./sum(sig_yields))

h_stack = ROOT.THStack()
mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
mc_yields, mc_hist_list = [], []

for i, bkg in enumerate(mc_file_list):
    if isinstance(bkg, list):
        for sub_bkg in bkg:
            arrays = pic.read_file(path+channel+"/"+sub_bkg+".root", var, tree, selections)
            if file_hist in globals():
                file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range, file_hist)
            else:
                file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range)
    else: 
        arrays = pic.read_file(path+channel+"/"+bkg+".root", var, tree, selections)
        file_hist, _, yields = pic.get_hist(arrays, var, sb_ratio, "mc_{}".format(i), bins, x_range)
    mc_hist.Add(file_hist)
    mc_hist_list.append(file_hist)
    mc_yields.append(yields)
print(sum(mc_yields))
for i, file_hist in enumerate(mc_hist_list):
    plot.Set(file_hist, LineWidth=0, FillColor=ROOT.TColor.GetColorDark(i+2))
    file_hist.Scale(1./sum(mc_yields))
    h_stack.Add(file_hist)

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

rp = pic.get_ratio_hist(sig_hist, mc_hist)
# rp_err = pic.get_ratio_hist(mc_hist, mc_hist, (0,4.5))
# rp = pic.get_S_over_sqrtB(sig_hist, mc_hist, ratio, (0.,0.))

c1.Update()
pads[0].cd()
plot.Set(pads[0], Logy=1)
h_stack.Draw("hist")
plot.Set(h_stack.GetXaxis(), LabelSize=0)
plot.Set(h_stack.GetYaxis(), Title=y_title)
sig_hist.Draw("E same")
plot.Set(sig_hist, MarkerStyle=8)

pads[0].Update()
# plot.Set(h_stack, Maximum=1.2*pads[0].GetFrame().GetY2())

legend = plot.PositionedLegend(0.6, 0.10, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.023)
# legend.AddEntry("sig", "sigx{:.0f}({:.2f})".format(ratio, sum(sig_yields)), "lep")
legend.AddEntry("sig", sig_legend[0]+"({:.2f})".format(sig_yields[0]), "lep")
for i in range(len(mc_file_list)):
    legend.AddEntry("mc_{}".format(i), mc_legend[i]+"({:.2f})".format(mc_yields[i]), "f")
legend.Draw()

print("========================")
print("Finish drawing upper pad")
print("========================")

pads[1].cd()
# plot.Set(pads[1], Logy=1)
# rp_err.Draw("E2 same")
rp.Draw("E same")
# plot.Set(rp_err, FillColorAlpha=(ROOT.kGray, 1))
line = ROOT.TLine()
plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
plot.DrawHorizontalLine(pads[1], line, 0)
plot.Set(rp, MarkerStyle=8)
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
