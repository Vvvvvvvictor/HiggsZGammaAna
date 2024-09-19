outdir="/eos/home-j/jiehan/parquet/nanov12/bkgmc_jet25"

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short
