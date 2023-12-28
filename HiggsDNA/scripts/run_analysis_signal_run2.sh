outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/signal_run2_2017_CR"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs --batch_system "local" #--short
