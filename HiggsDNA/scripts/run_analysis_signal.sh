outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/signal2"

python scripts/run_analysis.py --config "metadata/zgamma_signal.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs --short

#python scripts/convert_parquet_to_root.py --source $outdir/ggH_M125_2017/merged_nominal.parquet --target $outdir/merged.root --type mc --log DEBUG --process ggh --notag
