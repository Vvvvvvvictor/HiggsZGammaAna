<<<<<<< HEAD
outdir="/afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_2016postVFP_CR"
=======
outdir="/eos/home-j/jiehan/parquet/nanov9/data_for_norm_v1"
>>>>>>> 0a1443ddd30e7de7324e36ae4a34dacbebf59ffc

python scripts/run_analysis.py --config "metadata/zgamma_bkgmc_run2.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs #--batch_system "local" #--short 
