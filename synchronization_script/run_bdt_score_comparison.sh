#!/bin/bash

# BDT分数比较脚本
# 比较两个ROOT文件中的BDT分数

echo "开始BDT分数比较..."

# 设置Python环境（如果需要）
source /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/setup_env.sh

# 文件路径
FILE1="/eos/project-h/htozg-dy-privatemc/rzou/bdt/Output_ggF_rui_commonparam/SM1_2018_output.root"
FILE2="/eos/user/j/jiehan/root/outputs/bdt_comparison_20250818_153334/zero_to_one_jet/ZGToLLG_2018_zero_to_one_jet_bdt_comparison.root"

# 输出目录
OUTPUT_DIR="/eos/home-j/jiehan/root/outputs/bdt_score_comparison_$(date +%Y%m%d_%H%M%S)"

# 创建输出目录
mkdir -p $OUTPUT_DIR

echo "文件1: $FILE1"
echo "文件2: $FILE2"
echo "输出目录: $OUTPUT_DIR"

# 运行比较脚本
python3 compare_bdt_scores.py \
    --file1 "$FILE1" \
    --tree1 "outtree" \
    --score1 "BDT_score" \
    --event1 "event" \
    --file2 "$FILE2" \
    --tree2 "tree" \
    --score2 "external_bdt_score" \
    --event2 "event" \
    --output-dir "$OUTPUT_DIR" \
    #--max-events 10000  # 取消注释这行来限制事件数（用于测试）

echo "比较完成！结果保存在: $OUTPUT_DIR"
echo "查看结果文件:"
echo "  - 比较图: $OUTPUT_DIR/bdt_score_comparison.png"
echo "  - 2D密度图: $OUTPUT_DIR/bdt_score_comparison_2d.png"
echo "  - 匹配数据: $OUTPUT_DIR/matched_scores.csv"
echo "  - 统计摘要: $OUTPUT_DIR/comparison_summary.txt"
echo "  - ROOT文件: $OUTPUT_DIR/matched_scores.root"
