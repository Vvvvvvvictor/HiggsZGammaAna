outdir="/eos/home-j/jiehan/parquet/nanov9/background"
# outdir="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/Bkg_MC"

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" --with_skimmed #--short #--batch_system "local" condor
