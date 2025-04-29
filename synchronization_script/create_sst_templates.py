import os
import sys
import uproot
import ROOT
import numpy as np
import awkward as ak
from argparse import ArgumentParser
import logging

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=logging.INFO)

def get_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description='Generate templates directly from ROOT files with BDT categories')
    parser.add_argument('-i', '--inputFolder', default='/eos/home-j/jiehan/root/VBF_output/',
                        help='Path to the input folder')
    parser.add_argument('-o', '--outputFolder', default='/eos/user/j/jiehan/root/input_finalfit/templates/',
                        help='Path to the output folder')
    return parser.parse_args()

def add_to_hist(file_path, hist, var_name="llphoton_refit_m", 
                bdt_var="BDT_score_2j", bdt_min=None, bdt_max=None, 
                selection=None):
    """
    Add data from a ROOT file to a histogram with possible BDT score filtering.
    
    Args:
        file_path: Path to the ROOT file
        hist: ROOT histogram to fill
        var_name: Variable to plot in histogram
        bdt_var: BDT score variable name for filtering
        bdt_min: Minimum BDT score (inclusive)
        bdt_max: Maximum BDT score (exclusive)
        selection: Additional selection criteria
    """
    try:
        with uproot.open(file_path) as infile:
            if "outtree" not in infile:
                logging.warning(f"Tree 'outtree' not found in {file_path}")
                return hist, 0
                
            # Read all variables to ensure we get everything needed
            arrays = infile["outtree"].arrays(library="ak")
            
            # Skip if no events
            if len(arrays) == 0:
                return hist, 0
            
            # Apply BDT score filtering if requested
            if bdt_var in arrays.fields:
                # Extract BDT scores
                bdt_scores = arrays[bdt_var]
                
                # Apply BDT range selection
                mask = np.ones(len(arrays), dtype=bool)
                if bdt_min is not None and bdt_min != float('-inf'):
                    mask = mask & (bdt_scores >= bdt_min)
                if bdt_max is not None and bdt_max != float('inf'):
                    mask = mask & (bdt_scores < bdt_max)
                
                arrays = arrays[mask]
                
                # Skip if no events after BDT selection
                if len(arrays) == 0:
                    return hist, 0
            else:
                logging.warning(f"BDT variable '{bdt_var}' not found in {file_path}")
                
            # Rename the mass variable for consistency with template format
            # In the input files it's called "llphoton_refit_m", in the template it's "H_mass"
            if var_name == "llphoton_refit_m" and var_name in arrays.fields:
                # Make a copy with the new name for easier histogram filling
                arrays["H_mass"] = arrays[var_name]
                var_name = "H_mass"
            
            # Apply mass window selection if requested
            if selection == "blind" and var_name in arrays.fields:
                mass_mask = (arrays[var_name] < 120) | (arrays[var_name] > 130)
                arrays = arrays[mass_mask]
                
            # Skip if no events after selection
            if len(arrays) == 0:
                return hist, 0
                
            # Get variable data and weights
            var_data = arrays[var_name].to_numpy()
            weights = arrays["weight_corr"].to_numpy() if "weight_corr" in arrays.fields else np.ones(len(var_data))
            
            # Fill histogram
            for i in range(len(var_data)):
                hist.Fill(var_data[i], weights[i])
            
            return hist, np.sum(weights)
            
    except Exception as e:
        logging.error(f"Error processing {file_path}: {e}")
        return hist, 0

def create_templates(input_folder, output_folder):
    """Create Run2 and Run3 templates from ROOT files with BDT categories."""
    # Ensure output directory exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Define run groups
    run2_years = ["2016APV", "2016", "2017", "2018"]
    run3_years = ["2022", "2022EE", "2023", "2023BPix"]
    all_years = run2_years + run3_years
    
    # Sample groups
    data_samples = ["data"]
    bkg_samples = ["DY0", "EWK", "SM1"]  # SM1 is assumed to be ZG1 based on file naming
    signal_samples = ["GGF", "VBF", "ttH", "WH", "ZH"]
    
    # BDT score boundaries and categories
    # Boundaries are given as: inf, 0.518, 0.402, 0.083, -inf
    # We reverse them to go from lowest to highest for easier processing
    bdt_boundaries = [-float('inf'), 0.083, 0.286, 0.489, float('inf')]
    category_names = ["VBF0", "VBF1", "VBF2", "VBF3"]
    x_ranges = {
        "VBF0": [105, 170],
        "VBF1": [100, 165],
        "VBF2": [95, 161],
        "VBF3": [95, 160]
    }
    # category_names_run3 = ["VBF01", "VBF11", "VBF21", "VBF31"]
    
    # Histogram parameters
    xmin, xmax, xbin = 95, 170, 300
    channel = "two_jet"
    
    # Process each run separately
    for run_name, run_years in [
        # ("run2", run2_years),
        # ("run3", run3_years),
        ("all", all_years)
    ]:
        logging.info(f"Processing {run_name}...")

        # Special handling for bkg with dy mixing
        sstFile = ROOT.TFile("/eos/user/j/jiehan/root/input_finalfit/sst_vbf_hist_nodrop.root", "READ")
        hists = []
        name_map = {
            "vbf4": "VBF0",
            "vbf3": "VBF1",
            "vbf2": "VBF2",
            "vbf1": "VBF3",
        }
        for name in name_map.keys():
            xmin, xmax = x_ranges[name_map[name]]
            xbin = int((xmax - xmin) * 4)
            new_hist = ROOT.TH1D(f"bkg_all_{name_map[name]}", f"bkg_all_{name_map[name]}", 4*(xmax-xmin), xmin, xmax)
            new_hist.Sumw2()
            hist = sstFile.Get(f"{name}_dy_sm")
            for i in range(1, hist.GetNbinsX() + 1):
                x = hist.GetBinCenter(i)
                if x >= xmin and x <= xmax:
                    content = hist.GetBinContent(i)
                    error = hist.GetBinError(i)
                    new_bin = new_hist.FindBin(x)
                    new_hist.SetBinContent(new_bin, content)
                    new_hist.SetBinError(new_bin, error)
            new_hist.SetDirectory(0)
            hists.append(new_hist)
        sstFile.Close()

        output_file = ROOT.TFile(os.path.join(output_folder, f"template_{run_name}.root"), "RECREATE")
        
        # Process each BDT category
        for cat_idx in range(len(category_names)):
            cat_name = category_names[cat_idx] #if run_name == "run2" else category_names_run3[cat_idx]
            xmin, xmax = x_ranges[cat_name]
            xbin = int((xmax - xmin) * 4)
            bdt_min = bdt_boundaries[cat_idx]
            bdt_max = bdt_boundaries[cat_idx + 1]
            
            logging.info(f"Processing category {cat_name} with BDT range [{bdt_min}, {bdt_max})")
            
            # # Create background histogram

            for hist in hists:
                hist.Write()

            # hist_name = f"bkg_{run_name}_{cat_name}"
            # bkg_hist = ROOT.TH1D(hist_name, hist_name, xbin, xmin, xmax)
            # bkg_hist.Sumw2()
            
            # # Process background samples
            # for bkg in bkg_samples:
            #     for year in run_years:
            #         bkg_file = os.path.join(input_folder, f"{bkg}_{year}_output.root")
            #         if os.path.exists(bkg_file):
            #             logging.info(f"Adding {bkg_file} to background for {cat_name}")
            #             bkg_hist, yield_h = add_to_hist(bkg_file, bkg_hist, 
            #                                           bdt_min=bdt_min, bdt_max=bdt_max)
            #         else:
            #             logging.warning(f"Background file not found: {bkg_file}")
            
            # bkg_hist.Write(hist_name)

            # Create signal histogram
            hist_name = f"sig_{run_name}_{cat_name}"
            sig_hist = ROOT.TH1D(hist_name, hist_name, xbin, xmin, xmax)
            sig_hist.Sumw2()
            
            # Process signal samples
            for sig in signal_samples:
                for year in run_years:
                    sig_file = os.path.join(input_folder, f"{sig}_{year}_output.root")
                    if os.path.exists(sig_file):
                        logging.info(f"Adding {sig_file} to signal for {cat_name}")
                        sig_hist, yield_h = add_to_hist(sig_file, sig_hist, 
                                                      bdt_min=bdt_min, bdt_max=bdt_max)
                    else:
                        logging.warning(f"Signal file not found: {sig_file}")
            
            sig_hist.Write(hist_name)
            
            # Create full data histogram
            hist_name = f"data_full_{run_name}_{cat_name}"
            data_full_hist = ROOT.TH1D(hist_name, hist_name, xbin, xmin, xmax)
            data_full_hist.Sumw2()
            
            # Process data samples
            for data in data_samples:
                for year in run_years:
                    data_file = os.path.join(input_folder, f"{data}_{year}_output.root")
                    if os.path.exists(data_file):
                        logging.info(f"Adding {data_file} to full data for {cat_name}")
                        data_full_hist, yield_h = add_to_hist(data_file, data_full_hist, 
                                                           bdt_min=bdt_min, bdt_max=bdt_max)
                    else:
                        logging.warning(f"Data file not found: {data_file}")
            
            data_full_hist.Write(hist_name)
            
            # Create blinded data histogram
            hist_name = f"data_{run_name}_{cat_name}"
            data_hist = ROOT.TH1D(hist_name, hist_name, xbin, xmin, xmax)
            data_hist.Sumw2()
            
            # Process data samples with blinding
            for data in data_samples:
                for year in run_years:
                    data_file = os.path.join(input_folder, f"{data}_{year}_output.root")
                    if os.path.exists(data_file):
                        logging.info(f"Adding {data_file} to blinded data for {cat_name}")
                        data_hist, yield_h = add_to_hist(data_file, data_hist, 
                                                      bdt_min=bdt_min, bdt_max=bdt_max,
                                                      selection="blind")
                    else:
                        logging.warning(f"Data file not found: {data_file}")
            
            data_hist.Write(hist_name)
        output_file.Close()
    
    logging.info(f"Templates saved to: {output_folder}")
    return True

if __name__ == "__main__":
    args = get_args()
    create_templates(args.inputFolder, args.outputFolder)
    logging.info("Template generation complete")
