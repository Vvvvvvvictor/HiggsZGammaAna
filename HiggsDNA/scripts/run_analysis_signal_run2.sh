outdir="/eos/home-j/jiehan/parquet/nanov9/signal"
# outdir="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/Sig_MC"

python scripts/run_analysis.py --config "metadata/zgamma_signal_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" --with_skimmed #--short 
# python scripts/run_analysis.py --config "metadata/zgamma_signal_run2_ext.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" --with_skimmed #--short #
