baseDir=`pwd`

cd $baseDir/Trees2WS/
#!/bin/bash

# 设置输入目录
inputDir="/eos/user/j/jiehan/root/output_cor_syst"
outputWSDir="/eos/user/j/jiehan/root/ws_cor_syst"

# 遍历文件夹
# for proc in "$inputDir"/*; do
#     for year in "$proc"/*; do
#         for file in "$year"/*.root; do
#             if [[ -f "$file" ]]; then
#                 # 获取年份和文件名
#                 proc_name=$(basename "$proc")
#                 year_num=$(basename "$year")
#                 # 生成命令
#                 cmd="python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile $file --productionMode $proc_name --year $year_num --outputWSDir $outputWSDir --doSystematics" 
#                 echo $cmd
#                 $cmd
#             fi
#         done
#     done
# done
# python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile "/eos/user/j/jiehan/root/output_cor_syst/ggH/2017/output.root" --productionMode ggH --year 2017 --outputWSDir $outputWSDir --doSystematics

# python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile /eos/user/j/jiehan/root/output_cor_syst --productionMode ggH --year 2017

cd $baseDir/Signal/
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode fTest
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode calcPhotonSyst
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode getDiagProc
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode signalFit --groupSignalFitJobsByCat
python3 RunSignalScripts.py --inputConfig config_2017.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--useDCB"
for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 lep VH ZH ttHh ttHl;do
    python3 RunPlotter.py --years 2017 --ext 2017 --cats $cat --groupSignalFitByProc
done
# python3 RunPlotter.py --years 2017 --ext 2017