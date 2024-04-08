rm -rf eos_logs/data/Data_Double* eos_logs/data/Data_Single*
rm -rf /eos/user/j/jiehan/parquet/nanov9/data/Data_Double* /eos/user/j/jiehan/parquet/nanov9/data/Data_Single*
rm -rf /eos/user/j/jiehan/parquet/nanov9/data/analysis_manager.pkl 
bash scripts/run_analysis_data_run2.sh
