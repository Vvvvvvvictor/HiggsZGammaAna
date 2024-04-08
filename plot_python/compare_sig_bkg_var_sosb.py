import pic_template as pic
import plotting as plot
import pandas as pd
import ROOT
from pdb import set_trace

COLORS = ["Red", "Blue"]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

var_names = [
    "H_mass", 
    "delta_eta_jj",
    "delta_phi_jj",
    "delta_phi_zgjj",
    "gamma_eta",
    "gamma_mvaID",
    "H_ptt",
    "jet_1_pt",
    "jet_2_pt",
    "jet1G_deltaR",
    "jet2G_deltaR",
    "jet_1_btagDeepFlavB",
    "jet_2_btagDeepFlavB",
    "l1g_deltaR",
    "l2g_deltaR",
    "lep_cos_theta",
    "lep_phi",
    "photon_zeppenfeld",
    "pt_balance",
    "Z_cos_theta",
    "Z_lead_lepton_eta",
    "Z_sublead_lepton_eta",
    "H_relpt",
    "HZ_deltaRap"
    ]
x_ranges = [
    (100, 180),
    (0, 9),
    (0, 3.2),
    (0, 3.2),
    (-2.5, 2.5),
    (0.14, 1),
    (0, 160),
    (30, 330),
    (30, 150),
    (0.4, 6.4),
    (0.4, 6.4),
    (0, 0.1),
    (0, 0.15),
    (0.4, 4.8),
    (0.4, 3.4),
    (-1, 1),
    (-3.2, 3.2),
    (0, 5),
    (0, 1),
    (-1, 1),
    (-2.5, 2.5),
    (-2.5, 2.5),
    (0, 3),
    (-0.7, 0.7)
    ]
x_titles = [
    "m_{ll#gamma}(GeV/c^{2})", 
    "#Delta#eta_{jj}",
    "#Delta#phi_{jj}",
    "#Delta#phi_{z#gamma,jj}",
    "#eta_{#gamma}",
    "photon MVA",
    "pT_{t}^{Z#gamma}",
    "pT_{j1}(GeV/c)",
    "pT_{j2}(GeV/c)",
    "#DeltaR(#gamma,j1)",
    "#DeltaR(#gamma,j2)",
    "j1 btag",
    "j2 btag",
    "max(#DeltaR(l,#gamma))",
    "min(#DeltaR(l,#gamma))",
    "#cos#theta",
    "#phi",
    "Zeppenfeld #gamma",
    "system balance",
    "#cos#Theta",
    "#eta_{l1}",
    "#eta_{l2}",
    "pT_{ll#gamma}*c/m_{ll#gamma}",
    "#Deltay(ll,ll#gamma)"
    ]
var_df = pd.DataFrame({'var_name': var_names, 'x_range': x_ranges, 'x_title': x_titles})
    
# Dataset list
'''
sig_file_list = [
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ggH/2017.root",
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/VBF/2017.root",
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
    ["/eos/home-j/jiehan/root/2017/skimmed_ntuples/WW/2017.root", "/eos/home-j/jiehan/root/2017/skimmed_ntuples/WZ/2017.root", "/eos/home-j/jiehan/root/2017/skimmed_ntuples/ZZ/2017.root"]
]
data_file_list = [
    "/eos/home-j/jiehan/root/2017/skimmed_ntuples/data/2017.root"
]
'''

sig_file_list = [
    "/eos/user/j/jiehan/root/outputs/two_jet/ggH.root",
    "/eos/user/j/jiehan/root/outputs/two_jet/VBF.root",
    "/eos/user/j/jiehan/root/outputs/two_jet/WminusH.root",
    "/eos/user/j/jiehan/root/outputs/two_jet/WplusH.root",
    "/eos/user/j/jiehan/root/outputs/two_jet/ZH.root",
    "/eos/user/j/jiehan/root/outputs/two_jet/ttH.root"
]

data_file_list = [
    "/eos/user/j/jiehan/root/outputs/two_jet/data.root"
]

mc_file_list = [
    "/eos/user/j/jiehan/root/outputs/two_jet/ZGToLLG.root",
    "/eos/user/j/jiehan/root/outputs/two_jet/DYJetsToLL.root",
    [
        "/eos/user/j/jiehan/root/outputs/two_jet/EWKZ2J.root",
        "/eos/user/j/jiehan/root/outputs/two_jet/ZG2JToG2L2J.root"
    ],
    [
        "/eos/user/j/jiehan/root/outputs/two_jet/TT.root",
        "/eos/user/j/jiehan/root/outputs/two_jet/TTGJets.root"
    ],
    "/eos/user/j/jiehan/root/outputs/two_jet/TGJets.root",
    [
        "/eos/user/j/jiehan/root/outputs/two_jet/ttWJets.root",
        "/eos/user/j/jiehan/root/outputs/two_jet/ttZJets.root"
    ],
    "/eos/user/j/jiehan/root/outputs/two_jet/WGToLNuG.root",
    [
        "/eos/user/j/jiehan/root/outputs/two_jet/WW.root",
        "/eos/user/j/jiehan/root/outputs/two_jet/WZ.root",
        "/eos/user/j/jiehan/root/outputs/two_jet/ZZ.root"
    ]
]

# Basic set of picture's content
mc_legend = ["SM ZG", "DYJets", "EWK Z(G)", "TT(G)", "TGJets", "ttVJets", "WGamma", "Diboson"]
# mc_legend = []
data_legend = ["Data"]
sig_legend = ["Signal"]
channel = "test"
bins = 100
ratio = 1
sub_y_title = "S/#sqrt{B}"
selections = ["H_mass>100", "H_mass<180"]

# Convert root file to DataFrame
sig_arr = pd.DataFrame()
for i, sig in enumerate(sig_file_list):
    array, var_not_exist = pic.read_file(sig, var_df["var_name"], channel, selections)
    sig_arr = pd.concat([array, sig_arr], ignore_index=True)
    if (i == 0):
        var_df = var_df[~var_df['var_name'].isin(var_not_exist)]

mc_arr_list = []
for i, bkg in enumerate(mc_file_list):
    arrays = pd.DataFrame()
    if isinstance(bkg, list):
        for sub_bkg in bkg:
            array, var_not_exist = pic.read_file(sub_bkg, var_df["var_name"], channel, selections)
            arrays = pd.concat([arrays, array], ignore_index=True)
    else: 
        arrays, var_not_exist = pic.read_file(bkg, var_df["var_name"], channel, selections)
    mc_arr_list.append(arrays)

for var in var_df['var_name']:
    # var = "gamma_ptRelErr"
    target_row = var_df.loc[var_df['var_name'] == var]
    x_range = target_row['x_range'].values[0]
    x_title = target_row['x_title'].values[0]
    print(x_range, x_title)
    y_title = "Events/{:.3f}".format((x_range[1]-x_range[0])/bins)
    selections = ["{}>{}".format(var, x_range[0]), "{}<{}".format(var, x_range[1])]

    # Set initial style
    plot.ModTDRStyle()
    c1 = ROOT.TCanvas("c1", "bkg sig Comparison")
    pads = plot.TwoPadSplit(0.29, 0.005, 0.005)

    print("============================")
    print("Finish setting picture style")
    print("============================")

    maximum = 0.

    # get arrays
    sig_yields= []
    sig_hist, yields, _ = pic.get_hist(sig_arr, var, ratio, "sig", bins, x_range, selections=selections)
    sig_yields.append(yields)
    sig_hist.Scale(1./sum(sig_yields))
    # sig_hist.Scale(1./yields)
    # plot.Set(sig_hist, LineColor=2)
    sig_max = sig_hist.GetBinContent(sig_hist.GetMaximumBin())
    if sig_max > maximum:
        maximum = sig_max

    h_stack = ROOT.THStack()
    mc_hist = ROOT.TH1D("","",bins,x_range[0],x_range[1])
    mc_yields, mc_hist_list = [], []

    for i, arrays in enumerate(mc_arr_list):
        file_hist, yields, _ = pic.get_hist(arrays, var, 1, "mc_{}".format(i), bins, x_range, selections=selections)
        mc_hist_list.append(file_hist)
        mc_yields.append(yields)
    for i, file_hist in enumerate(mc_hist_list):
        file_hist.Scale(1./sum(mc_yields))
        mc_hist.Add(file_hist)
        plot.Set(file_hist, LineWidth=0, FillColor=ROOT.TColor.GetColorDark(i+2))
        h_stack.Add(file_hist)

    mc_max = mc_hist.GetBinContent(mc_hist.GetMaximumBin())
    if mc_max > maximum:
        maximum = mc_max

    # data_hists = ROOT.TH1D("","",bins,x_range[0],x_range[1])
    # data_yields = []

    # for i, data in enumerate(data_file_list):
    #     if isinstance(data, list):
    #         for sub_data in data:
    #             arrays = pic.read_file(sub_bkg, var, channel, selections)
    #             if data_hist in globals():
    #                 data_hist, yields, _ = pic.get_hist(arrays, var, 1, "data_{}".format(i), bins, x_range, data_hist)
    #             else:
    #                 data_hist, yields, _  = pic.get_hist(arrays, var, 1, "data_{}".format(i), bins, x_range)
    #     else: 
    #         arrays = pic.read_file(data, var, channel, selections)
    #         data_hist, yields, _ = pic.get_hist(arrays, var, 1, "data_{}".format(i), bins, x_range)
    #     data_hists.Add(data_hist)
    #     data_yields.append(yields)
    # data_hists.Scale(1./sum(data_yields))
    # plot.Set(data_hists, LineWidth=3)

    # data_max = data_hists.GetBinContent(data_hists.GetMaximumBin())
    # if data_max > maximum:
    #     maximum = data_max

    print("==============================")
    print("Finish reading skimmed ntuples")
    print("==============================")

    rp = pic.get_S_over_sqrtB(sig_hist, mc_hist, ratio)
    # rp = pic.get_S_over_sqrtB(sig_hist, data_hists, ratio, (0.,0.05))

    c1.Update()
    pads[0].cd()
    plot.Set(pads[0], Logy=0)
    h_stack.Draw("hist")
    plot.Set(h_stack.GetXaxis(), LabelSize=0)
    plot.Set(h_stack.GetYaxis(), Title=y_title)
    sig_hist.Draw("E same")
    plot.Set(sig_hist, MarkerStyle=8)
    # data_hists.Draw("E same")

    pads[0].Update()
    plot.Set(h_stack, Maximum=1.3*maximum)

    legend = plot.PositionedLegend(0.75, 0.15, 3, 0.005)
    plot.Set(legend, NColumns=3, TextSize=0.025, FillStyle=0)
    legend.AddEntry("sig", "sig({:4.2f})".format(sum(sig_yields)), "lep")
    for i in range(len(mc_file_list)):
        legend.AddEntry("mc_{}".format(i), mc_legend[i]+"({:.2f})".format(mc_yields[i]), "f")
    # for i in range(len(data_file_list)):
    #     legend.AddEntry("data_{}".format(i), data_legend[i]+"({:.2f})".format(data_yields[i]), "lep")
    legend.Draw()

    print("========================")
    print("Finish drawing upper pad")
    print("========================")

    pads[1].cd()
    rp.Draw("E same")
    line = ROOT.TLine()
    # plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
    # plot.DrawHorizontalLine(pads[1], line, 1)
    plot.Set(rp, MarkerStyle=8)
    rp.GetXaxis().SetTitle(x_title)
    rp.GetYaxis().SetTitle(sub_y_title)
    pads[1].Update()

    print("========================")
    print("Finish drawing lower pad")
    print("========================")

    plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55);
    plot.DrawCMSLogo(c1, "138 fb^{-1} (13 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45);

    print("========================")
    print("Finish adding CMS word")
    print("========================")

    c1.Draw()
    c1.SaveAs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/sig_bkg_{}_{}.png".format(channel, var))

    del sig_hist, file_hist, 
