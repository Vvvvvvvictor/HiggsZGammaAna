#!/usr/bin/env python
import argparse
from higgs_dna.utils.logger_utils import setup_logger
import pandas as pd
import uproot
import awkward as ak
import numpy as np
import os
import json
from importlib import resources

parser = argparse.ArgumentParser(
    description="Simple utility script to convert one parquet file into one ROOT file."
)
parser.add_argument("--source", type=str, help="Path to input file.")
parser.add_argument("--target", type=str, help="Path to desired output file.")
parser.add_argument("--type", type=str, help="Type of dataset (data or mc).")
parser.add_argument(
    "--log", dest="log", type=str, default="INFO", help="Logger info level"
)
parser.add_argument("--process", type=str, default="", help="Production mode.")
parser.add_argument(
    "--notag",
    dest="notag",
    action="store_true",
    default=False,
    help="create NOTAG dataset as well.",
)
parser.add_argument(
    "--do_syst",
    dest="do_syst",
    action="store_true",
    default=False,
    help="create branches for systematic variations",
)
parser.add_argument(
    "--cats",
    type=str,
    dest="cats_dict",
    default="",
    help="Dictionary containing category selections.",
)
parser.add_argument(
    "--vars",
    type=str,
    dest="vars_dict",
    default="",
    help="Dictionary containing variations.",
)
args = parser.parse_args()
source_path = args.source
target_path = args.target
type = args.type
notag = True if (type == "mc" and args.notag == True) else False
process = args.process if (args.process != "") else "data"

logger = setup_logger(level=args.log)

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
) + "/higgs_dna/"

# Dictionary for renaming variables in ROOT tree output for final fits
rename_dict = {
    "H_mass": "CMS_hzg_mass"
}

df_dict = {}
outfiles = {
    "ch": target_path.replace(
        "merged.root", "output_cHToGG_M125_13TeV_amcatnloFXFX_pythia8.root"
    ),
    "ggh": target_path.replace(
        "merged.root", "output_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8.root"
    ),
    "vbf": target_path.replace(
        "merged.root", "output_VBFHToGG_M125_13TeV_amcatnlo_pythia8.root"
    ),
    "vh": target_path.replace(
        "merged.root", "output_VHToGG_M125_13TeV_amcatnlo_pythia8.root"
    ),
    "tth": target_path.replace(
        "merged.root", "output_TTHToGG_M125_13TeV_amcatnlo_pythia8.root"
    ),
    "data": target_path.replace("merged.root", "allData_2017.root"),
}

# Loading category informations (used for naming of files to read/write)
if args.cats_dict != "":
    with open(BASEDIR + args.cats_dict) as pf:
    # with resources.open_text("higgs_dna", args.cats_dict) as pf:
        cat_dict = json.load(pf)
    for cat in cat_dict:
        logger.debug(f"Found category: {cat}")
else:
    logger.info(
        "You provided an invalid dictionary containing categories information, have a look at your version of prepare_output_file.py"
    )
    logger.info(
        "An inclusive NOTAG category is used as default"
    )
    cat_dict = {"NOTAG": {"cat_filter": [("pt", ">", -1.0)]}}

# Loading variation informations (used for naming of files to read/write)
# Active object systematics, weight systematics are just different sets of weights contained in the nominal file
if args.vars_dict != "":
    # with resources.open_text("higgs_dna", args.vars_dict) as pf:
    with open(BASEDIR + args.vars_dict) as pf:
        variation_dict = json.load(pf)
    for var in variation_dict:
        logger.debug(f"Found variation: {var}")
else:
    if args.do_syst:
        raise Exception(
            "You provided an invalid dictionary containing systematic variations information, have a look at your version of prepare_output_file.py"
        )

if args.do_syst:
    # object systematics come from a different file (you are supposed to have merged .parquet with the merge_parquet.py script)
    for var in variation_dict:
        df_dict[var] = {}
        for cat in cat_dict:
            var_path = source_path.replace(
                "merged.parquet", f"{variation_dict[var]}/{cat}_merged.parquet"
            )
            logger.info(
                f"Starting conversion of one parquet file to ROOT. Attempting to read file {var_path} for category: {cat}."
            )

            eve = ak.from_parquet(var_path)

            logger.info("Successfully read from parquet file with awkward.")

            dict = {}
            for i in eve.fields:
                i_re = rename_dict[i] if i in rename_dict else i
                dict[i_re] = eve[i]
                
            df_dict[var][cat] = dict

            logger.debug(
                f"Successfully created dict from awkward arrays for {var} variation for category: {cat}."
            )
            
    logger.info(f"Attempting to write dict to ROOT file {target_path}.")
else:
    for cat in cat_dict:
        var_path = source_path.replace("merged.parquet", f"{cat}_merged.parquet")
        logger.info(
            f"Starting conversion of one parquet file to ROOT. Attempting to read file {var_path}."
        )

        eve = ak.from_parquet(var_path)

        logger.info("Successfully read from parquet file wit awkward.")

        dict = {}
        for i in eve.fields:
            i_re = rename_dict[i] if i in rename_dict else i
            dict[i_re] = eve[i]
            
        df_dict[cat] = dict

        logger.debug(
            f"Successfully created dict from awkward arrays without variation for category: {cat}."
        )

cat_postfix = {"ggh": "GG2H", "vbf": "VBF", "tth": "TTH", "vh": "VH"}

# For MC: {inputTreeDir}/{production-mode}_{mass}_{sqrts}_{category}
# For data: {inputTreeDir}/Data_{sqrts}_{category}
labels = {}
names = {}
if type == "mc":
    for cat in cat_dict:
        names[
            cat
        ] = f"LLGammaTree/{process}_125_13TeV_{cat}"  # _"+cat_postfix[process]
        labels[cat] = []
    name_notag = "LLGammaTree/" + process + "_125_13TeV_NOTAG"
    # flashggFinalFit needs to have each systematic variation in a different branch
    if args.do_syst:
        for var in variation_dict:
            for cat in cat_dict:
                # for object systematics we have different files storing the variated collections with the nominal weights
                syst_ = var
                logger.info("found syst: %s for category: %s" % (syst_, cat))
                labels[cat].append(
                    [
                        "LLGammaTree/" + process + f"_125_13TeV_{cat}_" + syst_,
                        "weight",
                        syst_,
                        cat,
                    ]
                )
    else:
        for cat in cat_dict:
            syst_ = ""
            labels[cat].append(
                [
                    "LLGammaTree/" + process + f"_125_13TeV_{cat}_" + syst_,
                    "weight",
                    syst_,
                    cat,
                ]
            )

else:
    for cat in cat_dict:
        labels[cat] = []
        labels[cat].append([f"LLGammaTree/Data_13TeV_{cat}", cat])
        names[cat] = f"LLGammaTree/Data_13TeV_{cat}"
        # name = "LLGammaTree/Data_13TeV_WStest"

# Now we want to write the dictionary to a root file, since object systematics don't come from
# the nominal file we have to separate again the treatment of them from the object ones
with uproot.recreate(outfiles[process]) as file:
    logger.debug(outfiles[process])
    # Final fit want a separate tree for each category and variation,
    # the naming of the branches are quite rigid:
    # For MC: {inputTreeDir}/{production-mode}_{mass}_{sqrts}_{category}_{syst}
    # For data: {inputTreeDir}/Data_{sqrts}_{category}
    for cat in cat_dict:
        logger.debug(f"writing category: {cat}")

        if args.do_syst:
            # check that the category actually contains something, otherwise the slattening step will make the script crash,
            # an improvement (not sure if needed) may be to also write an empty TTree to not confuse FinalFit
            if len(df_dict["NOMINAL"][cat]["weight"]):
                for branch in df_dict["NOMINAL"][cat]:
                    # here I had to add a flattening step to help uproot with the type of the awkward arrays,
                    # if you don't flatten (event if you don't have a nested field) you end up having a type like (len_of_array) * ?type, which make uproot very mad apparently
                    df_dict["NOMINAL"][cat][branch] = ak.flatten(df_dict["NOMINAL"][cat][branch], axis=0)
                file[names[cat]] = df_dict["NOMINAL"][cat]
                if notag:
                    file[name_notag] = df_dict["NOMINAL"][cat]  # this is wrong, to be fixed
                for syst_name, weight, syst_, c in labels[cat]:
                    # Skip "NOMINAL" as information included in nominal tree
                    if syst_ == "NOMINAL":
                        continue
                    logger.debug(syst_name, weight, syst_, c)
                    # If the name is not in the variation dictionary it is assumed to be a weight systematic
                    if syst_ not in variation_dict:
                        red_dict = {
                            new_key: df_dict["NOMINAL"][cat][key]
                            for key, new_key in (
                                ["CMS_hgg_mass", "CMS_hgg_mass"],
                                [weight, "weight"],
                            )
                        }
                        logger.info(f"Adding {syst_name}01sigma to out tree...")
                        file[syst_name + "01sigma"] = red_dict
                    else:
                        red_dict = {
                            new_key: ak.flatten(df_dict[syst_][cat][key], 0)
                            for key, new_key in (
                                ["CMS_hgg_mass", "CMS_hgg_mass"],
                                [weight, "weight"],
                            )
                        }
                        logger.info(f"Adding {syst_name}01sigma to out tree...")
                        file[syst_name + "01sigma"] = red_dict
            else:
                logger.info(f"no events survived category selection for cat: {cat}")

        else:
            # if there are no syst there is no df_dict["NOMINAL"] entry in the dict
            if len(df_dict[cat][[*df_dict[cat]][0]]):
                # same as before
                for branch in df_dict[cat]:
                    df_dict[cat][branch] = ak.flatten(df_dict[cat][branch], axis=0)
                file[names[cat]] = df_dict[cat]
                if notag:
                    file[name_notag] = df_dict[cat]  # this is wrong, to be fixed
            else:
                logger.info(f"no events survived category selection for cat: {cat}")

    logger.info(
        f"Successfully converted parquet file to ROOT file for process {process}."
    )
