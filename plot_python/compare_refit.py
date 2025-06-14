import pic_template as pic
import plotting as plot
import ROOT

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

legend = ["refit", "original"]
samples = ["DYJetsToLL"]
path = "/eos/user/j/jiehan/root/skimmed_ntuples_run3/"
years = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]
channel = "inclusive"
selections = ["z_ee == 1"]

settings = {
        "H_mass": {"range": (100, 180, 80), "x_title": r"$m_{ll\gamma}$ [GeV]"},
        "Z_mass": {"range": (80, 100, 80), "x_title": r"$m_{ll}$ [GeV]"},
        }

for var, config in settings.items():
    bins = config["range"][2]
    x_range = (config["range"][0], config["range"][1])
    x_title = config["x_title"]
    plot.ModTDRStyle()
    c1 = ROOT.TCanvas("c1", "refit vs original")
    pad = plot.OnePad()

    print("============================")
    print("Finish setting picture style")
    print("============================")

    original_hist = ROOT.TH1D("original", "original", bins, x_range[0], x_range[1])
    refit_hist = ROOT.TH1D("refit", "refit", bins, x_range[0], x_range[1])
    original_yield, refit_yield = 0, 0
    for i in samples:
        for j in years:
            print("Reading from {}/{}/{}.root".format(path, i, j))
            arrays, _ = pic.read_root_file("{}/{}/{}.root".format(path, i, j), [var, "{}_refit".format(var)], channel, selections)
            original_sub_hist, original_sub_yield, original_sub_inte = pic.get_hist(arrays, var, 1, "original_{}_{}".format(i, j), bins, x_range)
            refit_sub_hist, refit_sub_yield, refit_sub_inte = pic.get_hist(arrays, "{}_refit".format(var), 1, "refit_{}_{}".format(i, j), bins, x_range)
            original_hist.Add(original_sub_hist)
            refit_hist.Add(refit_sub_hist)
            original_yield += original_sub_yield
            refit_yield += refit_sub_yield
    print("Original yield: {}, Refitted yield: {}".format(original_yield, refit_yield))
    plot.Set(original_hist, LineWidth=3, MarkerSize=1, MarkerStyle=8, MarkerColor=ROOT.TColor.GetColorDark(2))
    plot.Set(refit_hist, LineWidth=3, MarkerSize=1, MarkerStyle=8, MarkerColor=ROOT.TColor.GetColorDark(3))

    print("==============================")
    print("Finish reading skimmed ntuples")
    print("==============================")

    c1.Update()
    refit_hist.Draw("E0 same")
    original_hist.Draw("E0 same")

    legend = plot.PositionedLegend(0.3, 0.10, 3, 0.015)
    plot.Set(legend, TextSize=0.03)
    legend.AddEntry(original_hist, "Original", "l")
    legend.AddEntry(refit_hist, "Refitted", "l")
    legend.Draw()

    print("========================")
    print("Finish drawing")
    print("========================")

    plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55)
    plot.DrawCMSLogo(c1, "62.32 fb^{-1} (13.6 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45)

    print("========================")
    print("Finish drawing CMS logo")
    print("========================")

    c1.Draw()
    c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/refit_vs_original_{}.pdf".format(var))