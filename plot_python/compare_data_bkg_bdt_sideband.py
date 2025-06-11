import pic_template as pic
import plotting as plot
import ROOT
import argparse

ROOT.gROOT.SetBatch(True) # Run in batch mode (no GUI)

# Argument parser
parser = argparse.ArgumentParser(description="Compare data and MC BDT score in sideband for a given channel.")
parser.add_argument("--channel", type=str, default="zero_to_one_jet", 
                    help="Analysis channel (e.g., zero_to_one_jet, two_jet, etc.)")
args = parser.parse_args()

# Use the channel from arguments
channel = args.channel

print("=====================================================================")
print(f"Comparing Data vs MC in sideband for channel: {channel}")
print("=====================================================================")

# --- Configuration ---
folder = "test"  # Can be "val" or "test"
path = f"/eos/home-j/jiehan/root/outputs/{folder}/"
tree = channel  # Tree name is assumed to be the same as channel name
var = "bdt_score_t" # Variable to plot
bins = 100  # Number of bins
x_range = (0.0, 1.0)  # Range of the variable
x_title = "Transformed BDT score"
bin_width = (x_range[1] - x_range[0]) / bins
y_title = f"Events / {bin_width:.3f}"
sub_y_title = "Data / MC"

# Sideband selections: Example: H_mass < 120 GeV or H_mass > 130 GeV
# Ensure 'H_mass' is the correct variable name in your TTrees.
sideband_selection = ["(H_mass < 120 | H_mass > 130)"]

# MC files and legends - adjust as needed for your analysis
mc_file_list = [
    "ZGToLLG.root",
    "DYJetsToLL.root", # Assuming this is a single component now
    "EWKZ2J.root",
    # Add more MC files or groups:
    # "TT.root",
    # ["TTGJets.root", "TGJets.root"], # Example of a group
]
mc_legend = [
    "Z+#gamma", 
    "Z+Jets (DY)", 
    "EWK Z+#gamma",
    # Add corresponding legends:
    # "t#bar{t}",
    # "t#gamma/t#bar{t}#gamma",
]
# Colors for MC components
bkg_color = ["#3f90da", "#ffa90e", "#92dadd", "#e76300", "#bd1f01", "#94a4a2", "#832db6", "#33a02c", "#fb9a99"]


data_file_list = [
    "Data.root"
]

# Path for bin boundaries (optional)
boundaries_path_val = f"/eos/home-j/jiehan/root/outputs/test/significances/bin_boundaries_1D_{channel}.txt"

# --- Plotting Setup ---
plot.ModTDRStyle()
c1 = ROOT.TCanvas("c1", f"Data_MC_Sideband_{channel}", 800, 800)
pads = plot.TwoPadSplit(0.30, 0.005, 0.005) # Main pad and ratio pad

# --- Process MC Samples ---
h_stack_mc = ROOT.THStack("h_stack_mc", f";{x_title};{y_title}")
mc_hist_list = []
mc_total_hist = ROOT.TH1D("mc_total_hist", "Total MC", bins, x_range[0], x_range[1])
mc_total_hist.SetDirectory(0) # Avoid issues with ROOT's global directory

print("Processing MC files...")
for i, bkg_source in enumerate(mc_file_list):
    # This histogram will store the sum for the current source (single file or list of files)
    hist_for_stack_component = ROOT.TH1D(f"mc_stack_component_{i}", f"MC Component {i}", bins, x_range[0], x_range[1])
    hist_for_stack_component.SetDirectory(0)

    current_files_to_process = []
    if isinstance(bkg_source, list):
        current_files_to_process.extend(bkg_source)
    else:
        current_files_to_process.append(bkg_source)

    for mc_file_name in current_files_to_process:
        file_path = f"{path}{channel}/{mc_file_name}"
        print(f"  Reading MC file: {file_path} for component {i}")
        arrays, _ = pic.read_root_file(file_path, var, tree, sideband_selection)
        
        if arrays is None or (hasattr(arrays, 'size') and arrays.size == 0) or \
           (isinstance(arrays, (list, tuple)) and not arrays):
            print(f"    No data in {mc_file_name} for current selection, skipping.")
            continue
        
        # Create a temporary histogram for this file's arrays
        # Name must be unique if pic.get_hist registers it globally, or use None and it creates an anonymous one
        temp_mc_file_hist, _, _ = pic.get_hist(arrays, var, 1.0, f"temp_mc_{i}_{mc_file_name.replace('.root','')}", bins, x_range)
        hist_for_stack_component.Add(temp_mc_file_hist)

    if i < len(bkg_color):
        plot.Set(hist_for_stack_component, LineWidth=0, FillColor=ROOT.TColor.GetColor(bkg_color[i]), FillStyle=1001)
    else: # Fallback color if not enough defined
        plot.Set(hist_for_stack_component, LineWidth=0, FillColor=ROOT.kGray + i - len(bkg_color), FillStyle=1001)
        
    mc_hist_list.append(hist_for_stack_component)
    print(f"  Processed MC source {i} ({bkg_source}), Integral: {hist_for_stack_component.Integral()}")

# Add to stack and total sum AFTER all files for that component are processed
for hist in mc_hist_list:
    mc_total_hist.Add(hist)

# Add individual MC hists to stack in reverse order for correct drawing
for hist in reversed(mc_hist_list):
    h_stack_mc.Add(hist.Clone()) # Clone to be safe as THStack might take ownership

print(f"Total MC Integral (for ratio): {mc_total_hist.Integral()}")

# --- Process Data Sample ---
data_hist = ROOT.TH1D("data_hist", "Data", bins, x_range[0], x_range[1])
data_hist.SetDirectory(0)

print("Processing Data files...")
for data_file_name in data_file_list:
    file_path = f"{path}{channel}/{data_file_name}"
    print(f"  Reading data file: {file_path}")
    arrays, _ = pic.read_root_file(file_path, var, tree, sideband_selection)

    if arrays is None or (hasattr(arrays, 'size') and arrays.size == 0) or \
       (isinstance(arrays, (list, tuple)) and not arrays):
        print(f"    No data in {data_file_name} for current selection, skipping.")
        continue

    temp_data_file_hist, _, _ = pic.get_hist(arrays, var, 1.0, f"temp_data_{data_file_name.replace('.root','')}", bins, x_range)
    data_hist.Add(temp_data_file_hist)

plot.Set(data_hist, MarkerStyle=20, MarkerSize=1.0, LineColor=ROOT.kBlack) # Standard data style
print(f"Total Data Integral: {data_hist.Integral()}")

# --- Create Ratio Plot ---
ratio_hist = None
# Define base font sizes from ratio plot settings
base_font_title_size = 0.12
base_font_label_size = 0.11

if mc_total_hist.Integral() > 1e-6: # Avoid division by zero or tiny numbers
    ratio_hist = data_hist.Clone("ratio_hist")
    ratio_hist.Divide(mc_total_hist)
    plot.Set(ratio_hist, MarkerStyle=20, MarkerSize=1.0, LineColor=ROOT.kBlack)
    ratio_hist.GetYaxis().SetTitle(sub_y_title)
    ratio_hist.GetYaxis().SetNdivisions(505) # Standard for ratio plots
    # ratio_hist.GetYaxis().SetTitleSize(base_font_title_size)
    # ratio_hist.GetYaxis().SetLabelSize(base_font_label_size)
    # ratio_hist.GetYaxis().SetTitleOffset(pads[1].GetLeftMargin()*0.25 / base_font_title_size ) # Adjust offset based on margin
    ratio_hist.GetXaxis().SetTitle(x_title)
    # ratio_hist.GetXaxis().SetTitleSize(base_font_title_size)
    # ratio_hist.GetXaxis().SetLabelSize(base_font_label_size)
    ratio_hist.GetXaxis().SetTitleOffset(0.9) # Standard offset for X title on ratio
    ratio_hist.SetMinimum(0.0) # Typical Y-axis range for ratio
    ratio_hist.SetMaximum(2.0)
else:
    print("Warning: MC total histogram has zero or very small integral. Cannot create ratio plot.")

# --- Drawing ---
# Upper Pad (Main Plot)
pads[0].cd()
pads[0].SetLogy(True) # Set log scale for Y-axis

h_stack_mc.Draw("HIST") # Draw stacked MC

# Calculate font sizes for the main plot (pads[0]) to match apparent size of ratio plot (pads[1])
split_point_fraction = 0.30 # This is the fraction of the canvas height for the ratio pad (pads[1])
main_pad_height_fraction = 1.0 - split_point_fraction
ratio_pad_height_fraction = split_point_fraction

# Scale factor: (height_ratio_pad / height_main_pad)
font_scale_factor = ratio_pad_height_fraction / main_pad_height_fraction

scaled_title_size_main_pad_y = base_font_title_size * font_scale_factor
scaled_label_size_main_pad_y = base_font_label_size * font_scale_factor

# Apply scaled font sizes and Y-axis title offset to the main plot (h_stack_mc)
if h_stack_mc.GetXaxis(): # X-axis of the main plot
    h_stack_mc.GetXaxis().SetLabelSize(0)  # Hide X-axis labels on main plot
    h_stack_mc.GetXaxis().SetTitleSize(0)  # Hide X-axis title on main plot

if h_stack_mc.GetYaxis(): # Y-axis of the main plot
    h_stack_mc.GetYaxis().SetTitle(y_title) # Ensure Y title is set on the stack
    h_stack_mc.GetYaxis().SetTitleSize(scaled_title_size_main_pad_y)
    h_stack_mc.GetYaxis().SetLabelSize(scaled_label_size_main_pad_y)
    # Adjust Y-axis title offset for the main pad
    # Increase the multiplier (e.g., from 0.25) to create more space.
    # A factor of ~0.4 to 0.5 in the numerator part of the formula often works well.
    # Target offset value around 1.2-1.5.
    # NEW_FACTOR = TargetOffset * scaled_title_size / LeftMargin
    # e.g. NEW_FACTOR = 1.3 * 0.0514 / 0.16 = 0.417
    y_title_offset_main_pad = pads[0].GetLeftMargin() * 0.42 / scaled_title_size_main_pad_y
    h_stack_mc.GetYaxis().SetTitleOffset(y_title_offset_main_pad)


data_hist.Draw("E1 SAME") # Draw data with error bars on top

# Adjust Y-axis for stack and data after drawing to get maxima
pads[0].Update() # Ensure GetMaximum works correctly
max_val_mc = h_stack_mc.GetMaximum()
max_val_data = 0
for bin_idx in range(1, data_hist.GetNbinsX() + 1):
    if data_hist.GetBinContent(bin_idx) > 0: # Only consider bins with data
      max_val_data = max(max_val_data, data_hist.GetBinContent(bin_idx) + data_hist.GetBinError(bin_idx))

overall_max = max(max_val_mc, max_val_data)

if pads[0].GetLogy():
    h_stack_mc.SetMaximum(overall_max * 100) # Factor for log scale
    h_stack_mc.SetMinimum(max(0.01, mc_total_hist.GetMinimum(0)/100.0) if mc_total_hist.Integral() > 0 else 0.01) # Avoid zero or negative minimum for log
else:
    h_stack_mc.SetMaximum(overall_max * 1.3)
    h_stack_mc.SetMinimum(0)

# Legend
legend_x_start = 0.25
legend_y_start = 0.90 - (len(mc_legend) + 1) * 0.045 # Adjust Y start based on number of entries
legend = plot.PositionedLegend(legend_x_start, 0.25, 2, 0.04) # x, y_bottom_of_legend, n_columns, y_spacing
plot.Set(legend, TextSize=0.030, FillStyle=0, BorderSize=0)
legend.AddEntry(data_hist, f"Data ({data_hist.Integral():.1f})", "pe")
for i, mc_h in enumerate(mc_hist_list):
    if i < len(mc_legend):
        legend.AddEntry(mc_h, f"{mc_legend[i]} ({mc_h.Integral():.1f})", "f")
    else:
        legend.AddEntry(mc_h, f"MC Comp {i} ({mc_h.Integral():.1f})", "f") # Fallback legend
legend.Draw()

pads[0].Update()
pads[0].RedrawAxis()

# Lower Pad (Ratio Plot)
pads[1].cd()
pads[1].SetGridy(True)

# Scale factor: (height_ratio_pad / height_main_pad)
font_scale_factor = ratio_pad_height_fraction / main_pad_height_fraction

scaled_title_size_main_pad_y = base_font_title_size * font_scale_factor
scaled_label_size_main_pad_y = base_font_label_size * font_scale_factor

if ratio_hist:
    ratio_hist.Draw("E1")
    # Draw a line at y=1 for reference
    line_at_one = ROOT.TLine(x_range[0], 1.0, x_range[1], 1.0)
    line_at_one.SetLineStyle(2)
    line_at_one.SetLineColor(ROOT.kRed)
    line_at_one.Draw("SAME")
    ratio_hist.GetYaxis().SetTitleSize(scaled_title_size_main_pad_y)
    ratio_hist.GetYaxis().SetLabelSize(scaled_label_size_main_pad_y)
    ratio_hist.GetYaxis().SetTitleOffset(pads[1].GetLeftMargin() * 0.42 / scaled_title_size_main_pad_y) # Adjust offset based on margin
    ratio_hist.GetXaxis().SetTitleSize(scaled_title_size_main_pad_y)
    ratio_hist.GetXaxis().SetLabelSize(scaled_label_size_main_pad_y)

pads[1].Update()
pads[1].RedrawAxis()

# --- Final Touches ---
c1.cd() # Go back to canvas level for logos
plot.DrawCMSLogo(c1, "CMS", "Preliminary", 0, 0.1, -0.02, 1.2, cmsTextSize=0.55)
# Adjust lumi text as needed for your specific dataset combination
plot.DrawCMSLogo(c1, "Run 2 & Run 3 (13-13.6 TeV)", "", 3, 0.1, -0.02, 1.2, cmsTextSize=0.40) 

c1.Update()
output_filename = f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/data_mc_sideband_{channel}_{var}.png"
c1.SaveAs(output_filename)
c1.SaveAs(output_filename.replace(".png",".pdf")) # Also save as PDF
print(f"Plot saved as {output_filename} and .pdf")

print("========================")
print("Script finished.")
print("========================")

