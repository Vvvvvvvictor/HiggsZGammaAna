outdir="private/HZGamma/HiggsZGammaAna/HiggsDNA/data_run2_DoubleMuon_Run2018C_CR"

python scripts/run_analysis.py --config "metadata/zgamma_data_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --batch_system "condor" --unretire_jobs #--short 
