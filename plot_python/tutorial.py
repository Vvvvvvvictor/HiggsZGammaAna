import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
mc_legend = ["SM ZG", "DYJets", "LWK ZG", "TTG", "TG"]
data_legend = []
var = "Z_mass"
bins = 40
x_range = (80, 100)
x_title = "M_{ll}(GeV)"
y_title = "Events/GeV"
sub_y_title = "Data/MC"
selections = ["(H_mass<120) | (H_mass>130)", "gamma_mvaID_WP80==1", "(Z_mass>80) & (Z_mass<100)"]

# Dataset list
data_file_list = [
    "/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/data/2017.root"
]
mc_file_list = [
    "/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZGToLLG/2017.root",
    "/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/DYJetsToLL/2017.root",
    "/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/ZG2JToG2L2J/2017.root",
    "/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TTGJets/2017.root",
    "/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/TGJets/2017.root"
]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "MC Data Comparison")
pads = plot.TwoPadSplit(0.29, 0.005, 0.005)

print("============================")
print("Finish setting picture style")
print("============================")

# get hists
data_yields = []
for data in data_file_list:
    arrays = pic.read_file(data, "inclusive", selections)
    if "data_hist" in globals():
        data_hist, yields = pic.get_hist(arrays, var, "data", bins, x_range, data_hist)
    else:
        data_hist, yields = pic.get_hist(arrays, var, "data", bins, x_range)
    data_yields.append(yields)
    

h_stack = ROOT.THStack()
mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
mc_yields = []

for i, mc in enumerate(mc_file_list):
    arrays = pic.read_file(mc, "inclusive", selections)
    file_hist, yields = pic.get_hist(arrays, var, "mc_{}".format(i), bins, x_range)
    mc_hist.Add(file_hist)
    mc_yields.append(yields)
    # mc_hist.Sumw2()
    plot.Set(file_hist, LineWidth=0, FillColor=i+2)
    h_stack.Add(file_hist)

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

rp = pic.get_ratio_hist(data_hist, mc_hist)
rp_err = pic.get_ratio_hist(mc_hist, mc_hist)

c1.Update()
pads[0].cd()
h_stack.Draw("hist")
plot.Set(h_stack.GetXaxis(), LabelSize=0)
plot.Set(h_stack.GetYaxis(), Title=y_title)
data_hist.Draw("E same")

pads[0].Update()
plot.Set(h_stack, Maximum=1.2*pads[0].GetFrame().GetY2())

legend = plot.PositionedLegend(0.6, 0.10, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.023)
legend.AddEntry("data", "Data({:.0f})".format(sum(data_yields)), "lep")
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
rp_err.GetXaxis().SetTitle(x_title)
rp_err.GetYaxis().SetTitle(sub_y_title)
pads[1].Update()

print("========================")
print("Finish drawing lower pad")
print("========================")

plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55);
plot.DrawCMSLogo(c1, "5.67 fb^{-1} (13.6 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45);

print("========================")
print("Finish adding CMS word")
print("========================")

c1.Draw()
c1.SaveAs("test.pdf")
