import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
ratio = 500
bkg_legend = ["bkg"]
sig_legend = ["sig"]
fold_legend = ["1^{st}", "2^{nd}", "3^{rd}", "4^{th}"]
path = "/eos/user/j/jiehan/root/outputs/"
channel = "two_jet"
tree = "two_jet"
var = "bdt_score_t"
bins = 25
x_range = (0, 1)
x_title = "Transformed BDT score"
y_title = "Events/{:.2f}".format((x_range[1]-x_range[0])/bins)
selections = []

# Dataset list
sig_file_list = [
    ["ggH_M125.root", "VBF_M125.root", "WminusH_M125.root", "WplusH_M125.root", "ZH_M125.root", "ttH_M125.root"]
]
bkg_file_list = [
    ["ZGToLLG.root", "DYJetsToLL.root", "EWKZ2J.root", "ZG2JToG2L2J.root", "TT.root", "TTGJets.root", "TGJets.root", "WWG.root", "WZG.root", "ZZG.root", "ttZJets.root", "ttWJets.root"]
]
data_file_list = [
    "data.root"
]

# Set initial style
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", "bkg sig Comparison")
pad = plot.OnePad()

print("============================")
print("Finish setting picture style")
print("============================")

# get hists
sig_hist_list, bkg_hist_list = [], []
for i in range(8):
    path = "/eos/user/j/jiehan/root/outputs/{}/".format("test" if i < 4 else "val")
    print("Reading from {}".format(path))
    sig_yields,  bkg_yields = [], []
    fold_selection = ["event%314159%4 == {}".format(i%4)]
    hist_selections = selections + fold_selection

    bkg_hist = ROOT.TH1D("bkg_mc_{}".format(i),"bkg_mc_{}".format(i),bins,x_range[0],x_range[1])
    for j, bkg in enumerate(bkg_file_list):
        if isinstance(bkg, list):
            for k, sub_bkg in enumerate(bkg):
                arrays, _ = pic.read_root_file(path+channel+"/"+sub_bkg, var, tree, hist_selections)
                if k != 0:
                    file_hist, _, yields = pic.get_hist(arrays, var, 1, "subbkg_{}".format(j), bins, x_range, file_hist)
                else:
                    file_hist, _, yields = pic.get_hist(arrays, var, 1, "subbkg_{}".format(j), bins, x_range)
        else: 
            arrays, _ = pic.read_file(path+channel+"/"+bkg, var, tree, hist_selections)
            file_hist, _, yields = pic.get_hist(arrays, var, 1, "subbkg_{}".format(j), bins, x_range)
        bkg_yields.append(yields)
        bkg_hist.Add(file_hist)
    bkg_hist.Scale(1./sum(bkg_yields))
    plot.Set(bkg_hist, LineWidth=3, MarkerSize=1, MarkerStyle=8, MarkerColor=ROOT.TColor.GetColorDark(i+2))
    bkg_hist_list.append(bkg_hist)

    # sig_hist = ROOT.TH1D("sig_{}".format(i),"sig_{}".format(i),bins,x_range[0],x_range[1])
    # for j, sig in enumerate(sig_file_list):
    #     if isinstance(sig, list):
    #         for k, sub_sig in enumerate(sig):
    #             arrays, _ = pic.read_root_file(path+channel+"/"+sub_sig, var, tree, hist_selections)
    #             if k != 0:
    #                 file_hist, _, yields = pic.get_hist(arrays, var, 1, "subsig_{}".format(j), bins, x_range, file_hist)
    #             else:
    #                 file_hist, _, yields = pic.get_hist(arrays, var, 1, "subsig_{}".format(j), bins, x_range)
    #     else:
    #         arrays, _ = pic.read_file(path+channel+"/"+sig, var, 1, tree, hist_selections)
    #         file_hist, _, yields = pic.get_hist(arrays, var, ratio, "subsig_{}".format(j), bins, x_range)
    #     sig_yields.append(yields)
    #     sig_hist.Add(file_hist)
    # sig_hist.Scale(1./sum(sig_yields))
    # plot.Set(sig_hist, LineWidth=3, MarkerSize=1, MarkerStyle=8, MarkerColor=ROOT.TColor.GetColorDark(i+2))
    # sig_hist_list.append(sig_hist)
    

print("==============================")
print("Finish reading skimmed ntuples")
print("==============================")

c1.Update()
plot.Set(pad[0], Logy=1)
for bkg_hist in bkg_hist_list:
    bkg_hist.Draw("E0 E1 X0 same")
    plot.Set(bkg_hist.GetXaxis(), Title=x_title)
    plot.Set(bkg_hist.GetYaxis(), Title=y_title)
# for sig_hist in sig_hist_list:
#     sig_hist.Draw("E0 E1 X0 same")
#     plot.Set(sig_hist.GetXaxis(), Title=x_title)
#     plot.Set(sig_hist.GetYaxis(), Title=y_title)
# plot.Set(sig_hist, MarkerStyle=1, MarkerSize=3)

pad[0].Update()
# plot.Set(h_stack, Maximum=1.3*pads[0].GetFrame().GetY2())

legend = plot.PositionedLegend(0.65, 0.10, 3, 0.015)
plot.Set(legend, NColumns=3, TextSize=0.03)
for i in range(8):
    legend.AddEntry("bkg_mc_{}".format(i), "{} fold: {}".format(fold_legend[i%4], "test" if i < 4 else "val"), "p")
# for i in range(8):
    # legend.AddEntry("sig_{}".format(i), "{} fold: {}".format(fold_legend[i%4], "test" if i < 4 else "val"), "p")
# for i in range(len(sig_hist_list)):
#     legend.AddEntry("sig_{}".format(i), sig_legend[i]+"x{:.0f}({:.2f})".format(ratio, sig_yields[i]), "l")
# for i in range(len(bkg_file_list)):
#     legend.AddEntry("bkg_{}".format(i), bkg_legend[i]+"({:.2f})".format(bkg_yields[i]), "f")
legend.Draw()

print("========================")
print("Finish drawing")
print("========================")

# pads[1].cd()
# rp.Draw("E")
# # plot.Set(pads[1], Logy=1)
# line = ROOT.TLine()
# plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
# line.DrawLine(0.29, 0, 0.29, 0.26)
# line.DrawLine(0.57, 0, 0.57, 0.26)
# line.DrawLine(0.73, 0, 0.73, 0.26)
# # plot.DrawHorizontalLine(pads[1], line, 1)
# # plot.Set(rp, MarkerStyle=7, MarkerSize=2)
# rp.GetXaxis().SetTitle(x_title)
# rp.GetYaxis().SetTitle(sub_y_title)
# pads[1].Update()

# print("========================")
# print("Finish drawing lower pad")
# print("========================")

plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55);
plot.DrawCMSLogo(c1, "137.2 fb^{-1} (13 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45);

print("========================")
print("Finish adding CMS word")
print("========================")

c1.Draw()
c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/4fold_{}_bkg.png".format(channel))
