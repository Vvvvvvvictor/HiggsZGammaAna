import pic_template as pic
import plotting as plot
import ROOT
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plot 4-fold comparison with average curves')
parser.add_argument('--mode', choices=['test', 'val'], default='test', help='Choose test or val data')
parser.add_argument('--type', choices=['background', 'signal', 'data'], default='background', help='Choose plot type')
parser.add_argument('--channel', default='two_jet', help='Channel to plot')
parser.add_argument('--variable', choices=['bdt_score', 'H_mass'], default='bdt_score', help='Variable to plot')
parser.add_argument('--score_type', default='score_t', help='type of score to plot (only for BDT)')
args = parser.parse_args()

COLORS = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange, ROOT.kViolet, ROOT.kCyan, ROOT.kMagenta, ROOT.kYellow+2]

print("=====================================================================")
print("Starting to run the tutorial of plotting in Higgs to Z Gamma analysis")
print("=====================================================================")

# Basic set of picture's content
bkg_legend = ["bkg"]
sig_legend = ["sig"]
fold_legend = ["1^{st}", "2^{nd}", "3^{rd}", "4^{th}"]
path = "/eos/user/j/jiehan/root/outputs/"
channel = args.channel
tree = channel

# Set variable and plotting parameters based on choice
if args.variable == 'H_mass':
    var = "H_mass"
    bins = 85 if args.type == 'background' else 60
    x_range = (95, 180) if args.type == 'background' else (110, 135)
    x_title = "$m_{\ell\ell\gamma}$ [GeV/c^2]"
    plot_title_suffix = "H_mass"
    
    # Read bin boundaries for H_mass
    bin_boundaries_file = f"/eos/user/j/jiehan/root/outputs/test/significances/bin_boundaries_1D_{channel}.txt"
    try:
        with open(bin_boundaries_file, 'r') as f:
            boundaries_line = f.readline().strip()
            boundaries = [float(x) for x in boundaries_line.split()]
        print(f"Read bin boundaries: {boundaries}")
        num_bins = len(boundaries) - 1
        print(f"Number of bins: {num_bins}")
    except FileNotFoundError:
        print(f"Warning: Could not find bin boundaries file {bin_boundaries_file}")
        boundaries = None
        num_bins = 1
else:  # bdt_score
    score_type = args.score_type
    var = f"bdt_{score_type}"
    bins = 56
    x_range = (-0.06, 1.06)
    x_title = "Transformed BDT score" if score_type == "score_t" else "Raw BDT score"
    plot_title_suffix = f"bdt_{score_type}"
    boundaries = None
    num_bins = 1

y_title = "Events/{:.2f}".format((x_range[1]-x_range[0])/bins)
sub_y_title = "Ratio to Aver."
selections = []

# Dataset list
sig_file_list = [
    ["ggH_M125.root", "VBF_M125.root"] #, "WminusH_M125.root", "WplusH_M125.root", "ZH_M125.root", "ttH_M125.root"]
]
bkg_file_list = [
    ["ZGToLLG.root", "DYJetsToLL.root", "EWKZ2J.root"] #, "ZG2JToG2L2J.root", "TT.root", "TTGJets.root", "TGJets.root", "WWG.root", "WZG.root", "ZZG.root", "ttZJets.root", "ttWJets.root"]
]
data_file_list = [
    "data.root"
]

# Select which file list to use
if args.type == 'background':
    file_list = bkg_file_list
    plot_title = "Background"
    legend_list = bkg_legend
elif args.type == 'signal':
    file_list = sig_file_list
    plot_title = "Signal"
    legend_list = sig_legend
else:  # data
    file_list = data_file_list
    plot_title = "Data"
    legend_list = ["data"]

# Define color function
def get_color(i):
    return COLORS[i % len(COLORS)]

def create_histograms(bin_selection, bin_suffix):
    """Create histograms for 4 folds with given bin selection"""
    hist_list = []
    for i in range(4):
        path_full = "/eos/user/j/jiehan/root/outputs/{}/".format(args.mode)
        print("Reading fold {} from {}{}".format(i+1, path_full, " (no bin selection)" if not bin_selection else ""))
        yields = []
        fold_selection = ["event%314159%4 == {}".format(i)]
        hist_selections = selections + fold_selection + bin_selection

        hist = ROOT.TH1D("hist{}_{}".format(bin_suffix, i), "hist{}_{}".format(bin_suffix, i), bins, x_range[0], x_range[1])
        
        for j, file_item in enumerate(file_list):
            if isinstance(file_item, list):
                for k, sub_file in enumerate(file_item):
                    arrays, _ = pic.read_root_file(path_full+channel+"/"+sub_file, var, tree, hist_selections)
                    if k != 0:
                        file_hist, _, file_yields = pic.get_hist(arrays, var, 1, "subfile_{}".format(j), bins, x_range, file_hist)
                    else:
                        file_hist, _, file_yields = pic.get_hist(arrays, var, 1, "subfile_{}".format(j), bins, x_range)
            else: 
                arrays, _ = pic.read_root_file(path_full+channel+"/"+sub_file, var, tree, hist_selections)
                file_hist, _, file_yields = pic.get_hist(arrays, var, 1, "subfile_{}".format(j), bins, x_range)
            yields.append(file_yields)
            hist.Add(file_hist)
        
        # # Normalize
        # if sum(yields) > 0:
        #     hist.Scale(1./sum(yields))
        
        plot.Set(hist, LineWidth=2, MarkerSize=0.25, MarkerStyle=20, MarkerColor=get_color(i), LineColor=get_color(i))
        hist_list.append(hist)
    
    return hist_list

def calculate_average_and_chi2(hist_list, hist_suffix):
    """Calculate average histogram and chi2 values"""
    avg_hist = ROOT.TH1D("avg_hist{}".format(hist_suffix), "avg_hist{}".format(hist_suffix), bins, x_range[0], x_range[1])
    for hist_bin_idx in range(1, bins + 1):
        sum_val = 0
        sum_err2 = 0
        for hist in hist_list:
            sum_val += hist.GetBinContent(hist_bin_idx)
            sum_err2 += hist.GetBinError(hist_bin_idx)**2
        avg_hist.SetBinContent(hist_bin_idx, sum_val / 4.0)
        avg_hist.SetBinError(hist_bin_idx, (sum_err2**0.5) / 4.0)

    plot.Set(avg_hist, LineWidth=3, MarkerSize=0.5, MarkerStyle=21, MarkerColor=ROOT.kBlack, LineColor=ROOT.kBlack)

    # Calculate chi2 for each fold against the average
    chi2_values = []
    for hist in hist_list:
        chi2 = hist.Chi2Test(avg_hist, "CHI2")
        chi2_values.append(chi2)
    
    return avg_hist, chi2_values

def draw_and_save_plots(hist_list, avg_hist, chi2_values, pads, hist_suffix, c1, output_suffix):
    """Draw main and ratio plots, add CMS logo and save plot"""
    # Draw main pad
    pads[0].cd()
    
    # Draw all 4-fold histograms
    for i, hist in enumerate(hist_list):
        hist.Draw("E0 E1 X0 same" if i > 0 else "E0 E1 X0")
        plot.Set(hist.GetXaxis(), Title="", LabelSize=0)
        hist.SetMaximum(hist.GetMaximum() * 1.6)
        plot.Set(hist.GetYaxis(), Title=y_title)

    # Draw average histogram
    avg_hist.Draw("E0 E1 X0 same")
    plot.Set(avg_hist.GetXaxis(), Title="", LabelSize=0)
    plot.Set(avg_hist.GetYaxis(), Title=y_title)

    # Set y-axis minimum to 0
    for hist in hist_list:
        hist.SetMinimum(0)
    avg_hist.SetMinimum(0)
    
    # Update the pad before creating legend
    pads[0].Update()

    # Create legend
    legend = plot.PositionedLegend(0.35, 0.2, 3, 0.015)
    plot.Set(legend, NColumns=1, TextSize=0.04)
    for i in range(4):
        hist_name = "hist{}_{}".format(hist_suffix, i)
        legend.AddEntry(hist_name, "{} fold (#chi^{{2}}: {:.2f})".format(fold_legend[i], chi2_values[i]), "lep")
    avg_hist_name = "avg_hist{}".format(hist_suffix)
    legend.AddEntry(avg_hist_name, "Average", "lep")
    total_chi2 = sum(chi2_values)
    legend.AddEntry("", "Total #chi^{{2}}: {:.2f}".format(total_chi2), "")
    legend.Draw()

    pads[0].Update()

    # Draw ratio plot
    pads[1].cd()
    ratio_hists = []
    for i, hist in enumerate(hist_list):
        ratio_hist = hist.Clone("ratio{}_{}".format(hist_suffix, i))
        ratio_hist.Divide(avg_hist)
        plot.Set(ratio_hist, LineWidth=2, MarkerSize=0.25, MarkerStyle=20, MarkerColor=get_color(i), LineColor=get_color(i))
        ratio_hist.Draw("E0 E1 X0 same" if i > 0 else "E0 E1 X0")
        ratio_hists.append(ratio_hist)

    # Set axis labels for ratio plot
    if ratio_hists:
        plot.Set(ratio_hists[0].GetXaxis(), Title=x_title, LabelSize=0.04)
        plot.Set(ratio_hists[0].GetYaxis(), Title=sub_y_title)
        ratio_hists[0].GetYaxis().SetRangeUser(0, 2.75)
        for ratio_hist in ratio_hists:
            ratio_hist.SetMinimum(0)

    # Draw y=1 horizontal line
    line = ROOT.TLine()
    plot.Set(line, LineStyle=2, LineWidth=2, LineColor=ROOT.kRed)
    line.DrawLine(x_range[0], 1, x_range[1], 1)
    pads[1].Update()
    
    # Add CMS logo and save plot
    plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0., -0.02, 1.2, cmsTextSize=0.55)
    plot.DrawCMSLogo(c1, "137.61 fb^{-1} (13 TeV) & 62.32 fb^{-1} (13.6 TeV)", "", 3, 0., -0.02, 1.2, cmsTextSize=0.45)

    c1.Draw()
    output_filename = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{}/4fold_{}_{}_{}_{}.png".format(plot_title_suffix, channel, args.mode, args.type, output_suffix)

    # Check and create output directory if needed
    output_dir = os.path.dirname(output_filename)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Created directory: {}".format(output_dir))

    c1.SaveAs(output_filename)
    print("Saved plot to: {}".format(output_filename))

print("============================")
print("Starting plotting process")
print("============================")

# Loop over bins (for H_mass) or just one iteration (for BDT)
for bin_idx in range(num_bins):
    if boundaries is not None:
        # Create bin selection for H_mass
        bin_min = boundaries[bin_idx]
        bin_max = boundaries[bin_idx + 1]
        bin_selection = [f"bdt_score_t >= {bin_min} & bdt_score_t < {bin_max}"]
        print(f"Processing bin {bin_idx + 1}/{num_bins}: BDT score [{bin_min:.2f}, {bin_max:.2f})")
        bin_suffix = f"_bin{bin_idx + 1}"
    else:
        bin_selection = []
        bin_suffix = ""
        print("Processing single plot (no binning)")
    
    # Set initial style for each plot
    plot.ModTDRStyle()
    c1 = ROOT.TCanvas("c1_{}".format(bin_idx), "{} Comparison{}".format(plot_title, bin_suffix))
    pads = plot.TwoPadSplit(0.3, 0.005, 0.005)

    # Create histograms and calculate statistics
    hist_list = create_histograms(bin_selection, bin_suffix)
    avg_hist, chi2_values = calculate_average_and_chi2(hist_list, bin_suffix)
    
    # Draw plots and save
    draw_and_save_plots(hist_list, avg_hist, chi2_values, pads, bin_suffix, c1, bin_suffix.lstrip("_") if bin_suffix else "single")
    
    # Clean up for next iteration
    c1.Delete()
    
print("========================")
print("Finished processing all bins")
print("========================")

# If we have multiple bins, create an additional plot without bin selection
if num_bins > 1:
    print("============================")
    print("Creating additional plot without bin selection")
    print("============================")
    
    # Set initial style for the combined plot
    plot.ModTDRStyle()
    c1 = ROOT.TCanvas("c1_combined", "{} Comparison (All Bins Combined)".format(plot_title))
    pads = plot.TwoPadSplit(0.3, 0.005, 0.005)

    # Create histograms and calculate statistics
    hist_list = create_histograms([], "_combined")
    avg_hist, chi2_values = calculate_average_and_chi2(hist_list, "_combined")
    
    # Draw plots and save
    draw_and_save_plots(hist_list, avg_hist, chi2_values, pads, "_combined", c1, "combined")
    
    # Clean up
    c1.Delete()

print("========================")
print("All plots completed successfully")
print("========================")
