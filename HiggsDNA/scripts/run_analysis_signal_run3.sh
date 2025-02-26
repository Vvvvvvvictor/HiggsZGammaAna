outdir="/eos/home-p/pelai/HZgamma/Parquet/NanoV12/run3/Sig_MC_WO_Systematic"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short
