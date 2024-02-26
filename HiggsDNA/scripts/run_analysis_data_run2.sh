outdir="/eos/user/m/mingtao/workspace/zgamma/parquet/data_run2"

python scripts/run_analysis.py --config "metadata/zgamma_data_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --batch_system "local" --unretire_jobs --merge_outputs #--short 
