rm -rf /eos/user/s/shsong/HiggsDNA/testsys/
python scripts/run_analysis.py --config "metadata/zgamma_signal_run2_syst_splitskimmed.json"  --log-level "DEBUG"  --sample_list "ttH_M125" --output_dir "/eos/user/s/shsong/HiggsDNA/testsys/" --short --batch_system "local"
