outdir="/eos/home-j/jiehan/parquet/2017/nanov9/sig/"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs --batch_system "local" #--short

#python scripts/convert_parquet_to_root.py --source $outdir/ggH_M125_2017/merged_nominal.parquet --target $outdir/merged.root --type mc --log DEBUG --process ggh --notag
