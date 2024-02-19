outdir="/eos/home-j/jiehan/parquet/nanov9/signal"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs #--batch_system "local" #--short
