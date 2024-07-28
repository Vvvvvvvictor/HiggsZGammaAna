import uproot
import pandas as pd

for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
    standard = pd.read_parquet(f"/eos/user/j/jiehan/parquet/nanov9/bkgmc/DYJetsToLL_{year}/merged_nominal.parquet")["weight_central"][0]
    for type in ["DYJets_01J", "DYJets_2J"]:
        if type == "DYJets_01J":
            tree = "zero_to_one_jet"
        elif type == "DYJets_2J":
            tree = "two_jet"
        outputPath = f"/eos/home-j/jiehan/data_for_norm_float_v1/{type}/DYJets_{tree}_{year}.root"
        inputPath = f"/eos/home-j/jiehan/data_for_norm_float_v1/DYJets/DYJets_{tree}_{year}.root"
        print("Processing", inputPath)
        with uproot.recreate(outputPath) as outputFile:
            with uproot.open(inputPath) as inputFile:
                df = inputFile[tree].arrays(library="pd")
                print(df["weight"][0], standard)
                df["weight"] = abs(standard) * df["weight"] / abs(df.loc[0, "weight"])
                print(df["weight"][0])
                outputFile[tree] = df
            