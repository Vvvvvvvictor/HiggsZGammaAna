#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

input="/eos/user/j/jiehan/parquet/nanov9/mc_cor/"
target="/eos/home-j/jiehan/root/mc_cor_syst/"
# target="./"

# years=(2016preVFP 2016postVFP 2017 2018 2022preEE 2022postEE 2023preBPix 2023postBPix)
years=(2016preVFP 2016postVFP 2017 2018)

# 函数定义：执行命令并处理错误
execute_command() {
    local cmd="$1"
    local max_retries=10  # 最大重试次数
    local attempt=1

    while [ $attempt -le $max_retries ]; do
        echo "Attempt $attempt: $cmd"
        $cmd && break  # 如果命令成功执行，则跳出循环
        echo "Command $cmd failed. Retrying..."
        ((attempt++))
    done

    if [ $attempt -gt $max_retries ]; then
        echo "Error: Maximum retries reached. Command failed: $cmd"
    fi
}

# 函数定义：处理样本数据
process_sample() {
    local sample="$1"
    local type="$2"
    local corr="$3"
    
    for year in "${years[@]}"; do
        command="python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py "
        if [ "$type" = "signal" ]; then
            # command+="-i ${input}${type}/${sample}_${year}/merged_nominal.parquet "
            command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
        else
            # command+="-i ${input}${type}/${sample}_${year}/merged_nominal.parquet "
            command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
        fi
        if [ "$type" = "Data" ]; then
            command+="-o ${target}data/${year}.root"
        else
            command+="-o ${target}${sample}_${corr}/${year}.root"
        fi
        
        # 使用函数执行命令
        execute_command "$command" &
        pid_list+=($!)
    done

    # 等待所有后台任务完成
    for pid in "${pid_list[@]}"; do
        wait $pid
    done

    echo "Sample $sample completed successfully."
}

# 处理 signal 样本

samples=(ggH VBF WplusH WminusH ZH ttH) #ggH_M125 VBF_M125 WplusH_M125 WminusH_M125 ZH_M125 ttH_M125 ggH_M120 VBFH_M120 WplusH_M120 WminusH_M120 ZH_M120 ttH_M120 ggH_M130 VBFH_M130 WplusH_M130 WminusH_M130 ZH_M130 ttH_M130 ggH_mix VBF_mix)
type="signal"
for sample in "${samples[@]}"; do
    mkdir -p "$target${sample}_nominal"
    # 存储后台任务的进程ID列表
    pid_list=()

    # 调用函数处理样本数据
    process_sample "$sample" "$type" "nominal"
    for sf in "up" "down"; do
        for corr in "fnuf" "material" "scale" "smear" "JER" "JES" "MET_JES" "MET_Unclustered" "Muon_pt"; do
            mkdir -p "$target${sample}_${corr}_${sf}"
            # 存储后台任务的进程ID列表
            pid_list=()

            # 调用函数处理样本数据
            process_sample "$sample" "$type" "${corr}_${sf}"
        done
    done
done

# # 处理 bkgmc 样本

# # samples=(ZGToLLG DYJetsToLL WGToLNuG ZG2JToG2L2J EWKZ2J TT TTGJets TGJets ttWJets ttZJets WW WZ ZZ DYGto2LG_10to50 DYGto2LG_50to100)
# samples=(DYJetsToLL ZG2JToG2L2J)
# type="bkgmc"
# for sample in "${samples[@]}"; do
#     mkdir -p "$target$sample"
#     # 存储后台任务的进程ID列表
#     pid_list=()

#     # 调用函数处理样本数据
#     process_sample "$sample" "$type"
# done

# # 处理 data 样本

# samples=(Data)
# type="data"
# for sample in "${samples[@]}"; do
#     mkdir -p "$target$sample"
#     # 存储后台任务的进程ID列表
#     pid_list=()

#     # 调用函数处理样本数据
#     process_sample "$sample" "$type"
# done

# Use fake photon background estimation with data-driven

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_med/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_fake/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_true/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_med/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/apply_weight.py

# ######################
# Non prompt MC sample
# ######################

echo "==============FINISHED==========="
