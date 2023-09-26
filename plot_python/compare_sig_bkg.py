import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
mc_legend = ["SM ZG", "DYJets", "LWK ZG", "TT", "TG"]
sig_legend = []
channel = "VH_ttH"
var = "additional_lepton_1_pt"
bins = 40
x_range = (10, 50)
ratio = 500
x_title = "P_{T}^{l_{3}}(GeV/c)"
y_title = "Events/{:.2f}(GeV/c)".format((x_range[1]-x_range[0])/bins)
sub_y_title = "Sig/Bkg"
selections = [ "gamma_mvaID_WPL==1"]

# Dataset list
sig_file_list = [
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/WminusH/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/WplusH/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ZH/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ttH/2017.root"
]
mc_file_list = [
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ZGToLLG/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/DYJetsToLL/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/LLAJJ/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/TT/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/TGJets/2017.root"
]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "bkg sig Comparison")
pads = plot.TwoPadSplit(0.29, 0.005, 0.005)

print("============================")
print("Finish setting picture style")
print("============================")

# get hists
sig_yields = []
for sig in sig_file_list:
    arrays = pic.read_file(sig, var, channel, selections)
    if "sig_hist" in globals():
        sig_hist, yields = pic.get_hist(arrays, var, ratio, "sig", bins, x_range, sig_hist)
    else:
        sig_hist, yields = pic.get_hist(arrays, var, ratio, "sig", bins, x_range)
    sig_yields.append(yields)
    

h_stack = ROOT.THStack()
mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
mc_yields = []

for i, bkg in enumerate(mc_file_list):
    arrays = pic.read_file(bkg, var, channel, selections)
    file_hist, yields = pic.get_hist(arrays, var, 1, "mc_{}".format(i), bins, x_range)
    mc_hist.Add(file_hist)
    mc_yields.append(yields)
    # mc_hist.Sumw2()
    plot.Set(file_hist, LineWidth=0, FillColor=ROOT.TColor.GetColorDark(i+2))
    h_stack.Add(file_hist)

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

rp = pic.get_ratio_hist(sig_hist, mc_hist, (0,2.4))
rp_err = pic.get_ratio_hist(mc_hist, mc_hist, (0,2.4))

c1.Update()
pads[0].cd()
plot.Set(pads[0], Logy=1)
h_stack.Draw("hist")
plot.Set(h_stack.GetXaxis(), LabelSize=0)
plot.Set(h_stack.GetYaxis(), Title=y_title)
sig_hist.Draw("E same")
plot.Set(sig_hist, MarkerStyle=8)

pads[0].Update()
# plot.Set(h_stack, Maximum=1.3*pads[0].GetFrame().GetY2())

legend = plot.PositionedLegend(0.6, 0.10, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.023)
legend.AddEntry("sig", "sigx{:.0f}({:.2f})".format(ratio, sum(sig_yields)), "lep")
for i in range(len(mc_file_list)):
    legend.AddEntry("mc_{}".format(i), mc_legend[i]+"({:.2f})".format(mc_yields[i]), "f")
legend.Draw()

print("========================")
print("Finish drawing upper pad")
print("========================")

pads[1].cd()
rp_err.Draw("E2 same")
rp.Draw("E same")
plot.Set(rp_err, FillColorAlpha=(ROOT.kGray, 1))
line = ROOT.TLine()
plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
plot.DrawHorizontalLine(pads[1], line, 1)
plot.Set(rp, MarkerStyle=8)
rp_err.GetXaxis().SetTitle(x_title)
rp_err.GetYaxis().SetTitle(sub_y_title)
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
c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/test.png")
