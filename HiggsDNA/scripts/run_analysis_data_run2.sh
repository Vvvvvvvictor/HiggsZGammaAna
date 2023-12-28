outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/data_run2_2017_nanov2_CR_noPhoID_noEleVeto_all"

python scripts/run_analysis.py --config "metadata/zgamma_data_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "local" #--short
