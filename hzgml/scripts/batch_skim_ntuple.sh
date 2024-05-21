#!/bin/bash

# 定义输入文件夹路径
input_dir="/eos/home-j/jiehan/parquet/nanov9/data_for_norm_v1"

# 定义目标输出文件夹路径
target_dir="/eos/home-j/jiehan/data_for_norm_float_v1"

# 定义脚本路径
script_path="/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py"

# 定义log、err、out文件夹路径
log_dir="/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/eos_log"

# 遍历文件夹中的所有 parquet 文件
find "$input_dir" -type f -name '*merged_nominal.parquet' | while read -r input_file; do
    # 解析文件路径中的 sample 和 year 信息
    # 假设文件路径格式为 /eos/home-j/jiehan/parquet/nanov9/data_for_norm_v1/{sample}/{year}/merged_nominal.parquet
    echo $input_file
    sample=$(echo "$input_file" | awk -F'/' '{print $(NF-2)}')
    year=$(echo "$input_file" | awk -F'/' '{print $(NF-1)}')

    # 构建输出文件路径
    output_file="${target_dir}/${sample}_${year}.root"

    # 确保输出目录存在
    mkdir -p "$(dirname "$output_file")"

    # 创建HTCondor submission脚本
    submission_script="${log_dir}/${sample}_${year}.sub"

    cat <<EOL > "$submission_script"
Universe               = Vanilla
Executable             = /bin/bash
Arguments              = /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/eos_log/run.sh $input_file $output_file

Error                  = ${log_dir}/err.\$(Cluster)-\$(Process)
Output                 = ${log_dir}/out.\$(Cluster)-\$(Process)
Log                    = ${log_dir}/log.\$(Cluster)-\$(Process)

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
+JobFlavour = "workday"

Queue 1
EOL

    # 创建run.sh脚本
    run_script="${log_dir}/run.sh"
    cat <<'EOF' > "$run_script"
#!/bin/bash
input_file=$1
output_file=$2

source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_102/ROOT/6.26.04/x86_64-centos9-gcc11-opt/bin/thisroot.sh
source /eos/user/${USER::1}/$USER/hzgmlenv/bin/activate

export PATH="`pwd`:${PATH}"
export PYTHONPATH="`pwd`:${PYTHONPATH}"
export THEANO_FLAGS="gcc.cxxflags='-march=core2'"

export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"

export PATH="`pwd`/hzgml:$PATH"
export PYTHONPATH="`pwd`/hzgml:$PYTHONPATH"

export PYTHONPATH="/eos/user/${USER::1}/$USER/hzgmlenv/lib/python3.9/site-packages/:$PYTHONPATH"

python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i "$input_file" -o "$output_file"
EOF
    chmod +x "$run_script"

    # 提交作业
    condor_submit "$submission_script"
done
