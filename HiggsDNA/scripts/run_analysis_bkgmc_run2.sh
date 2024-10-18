outdir="/eos/home-j/jiehan/parquet/nanov9/data_for_norm_v2"

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs #--batch_system "local" #--short 
