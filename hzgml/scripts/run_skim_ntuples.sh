#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

# input="/eos/home-j/jiehan/parquet/nanov9/"
# target="/eos/home-j/jiehan/root/skimmed_ntuples/"
input="/eos/home-j/jiehan/parquet/nanov9/data_for_norm_v1/"
target="/eos/home-j/jiehan/data_for_norm_v1/"


# years=(2016preVFP 2016postVFP 2017 2018)
years=(2018)

<<<<<<< HEAD
done


=======
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

    for year in "${years[@]}"; do
        command="python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py "
        if [ "$type" = "signal" ]; then
            command+="-i ${input}${type}/${sample}_${year}/merged_nominal.parquet "
        else
            # command+="-i ${input}${type}/${sample}_${year}/merged_nominal.parquet "
            command+="-i ${input}/${sample}_${year}/merged_nominal.parquet "
        fi
        if [ "$type" = "Data" ]; then
            command+="-o ${target}data/${year}.root"
        else
            command+="-o ${target}${sample}/${year}.root"
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

# # 处理 signal 样本
# 
# samples=(ggH VBF WplusH WminusH ZH ttH)
# type="signal"
# for sample in "${samples[@]}"; do
#     mkdir -p "$target$sample"
#     # 存储后台任务的进程ID列表
#     pid_list=()
# 
#     # 调用函数处理样本数据
#     process_sample "$sample" "$type"
# done
# 
# # 处理 data 样本
# 
# samples=(Data)
# type="data"
# for sample in "${samples[@]}"; do
#     mkdir -p "$target$sample"
#     # 存储后台任务的进程ID列表
#     pid_list=()
<<<<<<< HEAD
# 
=======

>>>>>>> main
#     # 调用函数处理样本数据
#     process_sample "$sample" "$type"
# done

# # 处理 bkgmc 样本

<<<<<<< HEAD
samples=(ZG2JToG2L2J) # ZGToLLG DYJetsToLL WGToLNuG ZG2JToG2L2J EWKZ2J TT TTGJets TGJets ttWJets ttZJets WW WZ ZZ)
# samples=(ZGToLLG)
type="bkgmc"
for sample in "${samples[@]}"; do
    mkdir -p "$target$sample"
    # 存储后台任务的进程ID列表
    pid_list=()
=======
# samples=(ZGToLLG DYJetsToLL WGToLNuG ZG2JToG2L2J EWKZ2J TT TTGJets TGJets ttWJets ttZJets WW WZ ZZ)
# # samples=(ZGToLLG)
# type="bkgmc"
# for sample in "${samples[@]}"; do
#     mkdir -p "$target$sample"
#     # 存储后台任务的进程ID列表
#     pid_list=()
>>>>>>> main

#     # 调用函数处理样本数据
#     process_sample "$sample" "$type"
# done

# Use fake photon background estimation with data-driven

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_med/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_fake/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_true/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_med/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/apply_weight.py

# ######################
# Non prompt MC sample
# ######################
>>>>>>> 0a1443ddd30e7de7324e36ae4a34dacbebf59ffc

echo "==============FINISHED==========="
