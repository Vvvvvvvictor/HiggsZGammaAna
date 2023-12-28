outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_run2_2017_CR2"

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs #--batch_system "local" #--short
