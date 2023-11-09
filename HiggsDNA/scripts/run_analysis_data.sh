outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/data"

python scripts/run_analysis.py --config "metadata/zgamma_data.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --short
