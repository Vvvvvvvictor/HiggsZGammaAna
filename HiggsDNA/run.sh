# rm -rf eos_logs/data/Data_Double* eos_logs/data/Data_Single*
# rm -rf /eos/user/j/jiehan/parquet/nanov9/data/Data_Double* /eos/user/j/jiehan/parquet/nanov9/data/Data_Single*
rm -rf /eos/user/j/jiehan/parquet/nanov9/data/analysis_manager.pkl 
rm -rf eos_logs/data/Data_DoubleMuon*
rm -rf /eos/user/j/jiehan/parquet/nanov9/data/Data_DoubleMuon*
bash scripts/run_analysis_data_run2.sh

# rm -rf eos_logs/signal/ggH*
# rm -rf /eos/user/j/jiehan/parquet/nanov9/signal/ggH*
# rm -rf /eos/user/j/jiehan/parquet/nanov9/signal/analysis_manager.pkl 
# bash scripts/run_analysis_signal_run2.sh
