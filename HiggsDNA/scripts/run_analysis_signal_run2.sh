outdir="/eos/home-j/jiehan/parquet/cutflow_ggf_var"

# python scripts/run_analysis.py --config "metadata/zgamma_signal_run2_store.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short
python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short
