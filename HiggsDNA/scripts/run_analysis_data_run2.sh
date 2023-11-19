outdir="/eos/home-j/jiehan/parquet/2017/nanov9/data/"

python scripts/run_analysis.py --config "metadata/zgamma_data_run2.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs --batch_system "local" #--short
