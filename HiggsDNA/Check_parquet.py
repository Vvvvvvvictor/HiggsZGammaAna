import pandas as pd

try:
    df = pd.read_parquet("/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/run2/NanoV9/Bkg_MC/DYJetsToLL_2016postVFP/job_1/output_job_1_nominal.parquet")
    print(df.head())  # Print first few rows
except Exception as e:
    print(f"Error reading Parquet file: {e}")
