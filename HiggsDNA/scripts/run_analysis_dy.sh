outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_run2_2017_DY_madgraph_jiehan_test"

python scripts/run_analysis.py --config "metadata/zgamma_dy.json" --log-level "DEBUG" --n_cores 5 --output_dir $outdir --unretire_jobs --batch_system "local" --short
