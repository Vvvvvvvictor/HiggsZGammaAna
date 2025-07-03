import pic_template as pic
import plotting as plot
import ROOT
import os

print("=====================================================================")
print("Starting to plot background MC shapes in different BDT score bins")
print("=====================================================================")

# Basic set of picture's content
channel = "zero_to_one_jet"  # Channel name, e.g., "zero_to_one_jet", "one_to_two_jets", etc.
tree = "zero_to_one_jet"
var = "H_mass"  # 绘制 H_mass 分布
bins = 85
x_range = (95, 180)
x_title = "m_{ll#gamma} [GeV]"
y_title = "Events/{:.1f} GeV".format((x_range[1]-x_range[0])/bins)
path = "/eos/user/j/jiehan/root/outputs/test/"

# BDT score bins from 0.2 to 0.3 with step 0.02
bdt_bins = []
for i in range(4):  # 0.2, 0.22, 0.24, 0.26, 0.28, 0.30
    bdt_start = 0.20 + i * 0.05
    bdt_end = bdt_start + 0.05
    bdt_bins.append((bdt_start, bdt_end))

# Background file list from plot_cats_hmass_dis_mod.py
bkg_file_list = ["DYJetsToLL.root"] #"ZGToLLG.root", "DYJetsToLL.root", "EWKZ2J.root", "ZG2JToG2L2J.root", "TT.root", "TTGJets.root", "TGJets.root", "WWG.root", "WZG.root", "ZZG.root", "ttZJets.root", "ttWJets.root"]
bkg_labels = ["Z+#gamma", "Z+Fake #gamma", "VBSZ+#gamma"] #, "ZG2J", "t#bar{t}", "t#bar{t}#gamma", "t#gamma", "WW#gamma", "WZ#gamma", "ZZ#gamma", "t#bar{t}Z", "t#bar{t}W"]

# Color scheme
colors = [ROOT.kBlue, ROOT.kGreen+2, ROOT.kCyan, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange, ROOT.kYellow+2, ROOT.kViolet, ROOT.kPink, ROOT.kGray, ROOT.kTeal, ROOT.kSpring]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "Background MC BDT bins Comparison", 1200, 800)
pad = plot.OnePad()

print("============================")
print("Finish setting picture style")
print("============================")

# Get histograms for each BDT bin
bdt_hist_list = []
bdt_legend_list = []

for i, (bdt_start, bdt_end) in enumerate(bdt_bins):
    print("Processing BDT bin [{:.2f}, {:.2f}]".format(bdt_start, bdt_end))
    
    # Selection for this BDT bin
    bdt_selection = ["bdt_score_t >= {:.2f}".format(bdt_start), "bdt_score_t < {:.2f}".format(bdt_end)]
    
    # Create combined histogram for this BDT bin
    combined_hist = ROOT.TH1D("bkg_combined_{}".format(i), "bkg_combined_{}".format(i), bins, x_range[0], x_range[1])
    total_yields = 0
    
    for j, bkg_file in enumerate(bkg_file_list):
        arrays, _ = pic.read_root_file(path + channel + "/" + bkg_file, var, tree, bdt_selection)
        if len(arrays[var]) > 0:  # Check if there are events in this bin
            file_hist, _, yields = pic.get_hist(arrays, var, 1, "bkg_{}_{}".format(i, j), bins, x_range)
            combined_hist.Add(file_hist)
            total_yields += yields
    
    # Normalize histogram
    combined_hist.Scale(1.0 / total_yields)
    
    # Set histogram style
    plot.Set(combined_hist, LineWidth=3, LineColor=colors[i % len(colors)], MarkerColor=colors[i % len(colors)])
    
    bdt_hist_list.append(combined_hist)
    bdt_legend_list.append("BDT [{:.2f}, {:.2f}]".format(bdt_start, bdt_end))

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

# Draw histograms
c1.Update()
plot.Set(pad[0], Logy=0)

# Find maximum for y-axis scaling
max_val = 0
for hist in bdt_hist_list:
    if hist.GetMaximum() > max_val:
        max_val = hist.GetMaximum()

# Draw all histograms
for i, hist in enumerate(bdt_hist_list):
    if i == 0:
        hist.Draw("E1")
        plot.Set(hist.GetXaxis(), Title=x_title)
        plot.Set(hist.GetYaxis(), Title="Normalized " + y_title)
        hist.SetMaximum(max_val * 1.2)
    else:
        hist.Draw("E1 same")

total_hist = bdt_hist_list[0].Clone("total_background")
for hist in bdt_hist_list[1:]:
    total_hist.Add(hist)
total_hist.Scale(1.0 / total_hist.Integral())
# Draw total histogram
total_hist.SetLineColor(ROOT.kBlack)
total_hist.SetLineWidth(3)
total_hist.Draw("HIST E1 same")

# # Set X-axis range to 110-140 GeV
# for hist in bdt_hist_list:
#     hist.GetXaxis().SetRangeUser(110, 140)
# total_hist.GetXaxis().SetRangeUser(110, 140)

pad[0].Update()

# Create legend
legend = plot.PositionedLegend(0.65, 0.15, 3, 0.015)
plot.Set(legend, NColumns=2, TextSize=0.03)
for i, legend_text in enumerate(bdt_legend_list):
    legend.AddEntry(bdt_hist_list[i], legend_text, "pe")
legend.AddEntry(total_hist, "Average Background", "lpe")
legend.Draw()

print("========================")
print("Finish drawing")
print("========================")

# Add CMS labels
plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55)
plot.DrawCMSLogo(c1, "137.2 fb^{-1} (13 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45)

print("========================")
print("Finish adding CMS word")
print("========================")

c1.Draw()
if not os.path.exists("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/bump_check/"):
    os.makedirs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/bump_check/")
c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/bump_check/bdt_bins_bkg_{}.png".format(channel))

print("=====================================================================")
print("Plot saved to: /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/bump_check/bdt_bins_bkg_{}.png".format(channel))
print("=====================================================================")
