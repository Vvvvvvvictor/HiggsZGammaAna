outdir="/eos/home-p/pelai/HZgamma/Parquet/NanoV9/run2/Sig_MC/WI_Systematic"
# outdir="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/run2/NanoV9/Sig_MC"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2_store.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short

# local, condor 
