outdir="/eos/home-j/jiehan/parquet/nanov9/cutflow"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short
