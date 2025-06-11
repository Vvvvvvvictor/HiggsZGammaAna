import uproot
import pandas as pd
import numpy as np # Add numpy import
import os
import logging
import argparse # Import argparse here

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

filesPath = "/eos/project-h/htozg-dy-privatemc/rzou/bdt/VBF_output"
dataName = "data"
bkgName = ["DY0", "EWK", "SM"]
sigName = ["GGF", "VBF", "WH", "ZH", "ttH"]
boundaries = [-np.inf, 0.083, 0.286, 0.489, np.inf]

def get_yields(file_path, name, boundaries, is_sb=False):
    """
    Get yields from the ROOT file.
    """
    fileList = [file for file in os.listdir(file_path) if name in file]
    if len(fileList) == 0:
        logging.warning(f"No files found for {name} in {file_path}")
        return None
    yields = 0
    for file in fileList:
        data = uproot.open(os.path.join(file_path, file))["outtree"].arrays(["weight_corr", "BDT_score_2j", "llphoton_refit_m"], library="pd")
        data = data.query("BDT_score_2j > @boundaries[0] and BDT_score_2j < @boundaries[1]")
        if is_sb:
            data = data.query("llphoton_refit_m < 120 or llphoton_refit_m > 130")
        else:
            data = data.query("llphoton_refit_m > 120 and llphoton_refit_m < 130")
        yields += data["weight_corr"].sum()
    return yields

for i in range(len(boundaries) - 1):
    logging.info(f"Processing boundary {i + 1}: {boundaries[i]} to {boundaries[i + 1]}")
    sigYield, bkgYield, bkgSBYield, dataSBYield = 0, 0, 0, 0
    for sig in sigName:
        sigYield += get_yields(filesPath, sig, boundaries[i:i + 2])
    for bkg in bkgName:
        bkgYield += get_yields(filesPath, bkg, boundaries[i:i + 2])
        bkgSBYield += get_yields(filesPath, bkg, boundaries[i:i + 2], is_sb=True)
    dataSBYield = get_yields(filesPath, dataName, boundaries[i:i + 2], is_sb=True)
    logging.info(f"Signal yield: {sigYield}, Background yield: {bkgYield}, Background SB yield: {bkgSBYield}, Data SB yield: {dataSBYield}")
    logging.info(f"Category {i + 1}: sig/sqrt(bkg) = {sigYield / np.sqrt(bkgYield)}, sig/sqrt(bkgCorr) = {sigYield / (np.sqrt(bkgYield / bkgSBYield * dataSBYield))}") 