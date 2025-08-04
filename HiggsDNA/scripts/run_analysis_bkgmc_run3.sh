# outdir="/eos/home-j/jiehan/parquet/test_bkg"
outdir="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/Bkg_MC"

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "local" --short --with_skimmed #--batch_system "local" condor