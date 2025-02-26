outdir="/eos/home-p/pelai/HZgamma/Parquet/NanoV9/run2/Data"

python scripts/run_analysis.py --config "metadata/zgamma_data_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --batch_system "condor" --unretire_jobs #local, condor