# outdir="/eos/home-j/jiehan/parquet/nanov12/signal"
outdir="/eos/home-p/pelai/HZgamma/parquet_DNA/run3/signal"

rm -fr /eos/home-p/pelai/HZgamma/parquet_DNA/run3/signal/ggH_M125_2023postBPix
rm -fr /eos/home-p/pelai/HZgamma/parquet_DNA/run3/signal/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/zgamma_signal_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "local" #--short #--with_skimmed #condor #local
# python scripts/run_analysis.py --config "metadata/zgamma_signal_run3_ext.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--with_skimmed #--short #
