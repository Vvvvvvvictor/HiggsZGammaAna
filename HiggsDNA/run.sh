# rm -rf eos_logs/data/Data_Double* eos_logs/data/Data_Single*
# rm -rf /eos/user/j/jiehan/parquet/nanov9/data/Data_Double* /eos/user/j/jiehan/parquet/nanov9/data/Data_Single*
rm -rf /eos/user/j/jiehan/parquet/nanov9/data/analysis_manager.pkl 
rm -rf eos_logs/data/Data_2018
rm -rf /eos/user/j/jiehan/parquet/nanov9/data/Data_2018
bash scripts/run_analysis_data_run2.sh

# rm -rf eos_logs/signal/ggH*
# rm -rf /eos/user/j/jiehan/parquet/nanov9/signal/ggH*
# rm -rf /eos/user/j/jiehan/parquet/nanov9/signal/analysis_manager.pkl 
# bash scripts/run_analysis_signal_run2.sh

# rm -rf eos_logs/bkgmc/DYJetsToLL_2017
# rm -rf /eos/user/j/jiehan/parquet/nanov9/bkgmc/DYJetsToLL_2017
# rm -rf /eos/user/j/jiehan/parquet/nanov9/bkgmc/analysis_manager.pkl 
# bash scripts/run_analysis_bkgmc_run2.sh
