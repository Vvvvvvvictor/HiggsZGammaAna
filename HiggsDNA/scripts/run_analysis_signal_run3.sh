outdir="/eos/home-j/jiehan/parquet/cutflow_ggf"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run3.json" --log-level "DEBUG" --n_cores 15 --output_dir $outdir --unretire_jobs --batch_system "local" --with_skimmed # --short
