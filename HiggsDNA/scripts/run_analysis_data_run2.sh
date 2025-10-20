outdir="/eos/home-j/jiehan/parquet/nanov9/data"
# outdir="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/Data"

python scripts/run_analysis.py --config "metadata/zgamma_data_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --batch_system "condor" --unretire_jobs 