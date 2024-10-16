outdir="/eos/user/z/zewang/HZGamma_data/run2UL_CR/signal"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "local" # --short
