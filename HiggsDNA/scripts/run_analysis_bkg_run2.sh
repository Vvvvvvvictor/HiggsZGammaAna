outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_run2_test"

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run2.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs --short
