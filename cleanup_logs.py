#!/usr/bin/env python3
"""
清理 HiggsDNA/eos_logs/background/*_201*/job_* 目录中的多余日志文件
只保留每种类型(.log, .err, .out)的最新文件
"""

import os
import glob
import sys
from collections import defaultdict
from pathlib import Path

def get_file_groups(job_dir):
    """
    将同一个 job 目录下的日志文件按类型分组
    返回: {extension: [file_paths]}
    """
    file_groups = defaultdict(list)
    
    # 查找所有 .log, .err, .out 文件
    for ext in ['log', 'err', 'out']:
        pattern = os.path.join(job_dir, f"*.{ext}")
        files = glob.glob(pattern)
        if files:
            file_groups[ext] = files
    
    return file_groups

def get_latest_file(file_list):
    """
    从文件列表中找到修改时间最新的文件
    """
    if not file_list:
        return None
    
    latest_file = None
    latest_time = 0
    
    for file_path in file_list:
        try:
            mtime = os.path.getmtime(file_path)
            if mtime > latest_time:
                latest_time = mtime
                latest_file = file_path
        except OSError:
            continue
    
    return latest_file

def cleanup_job_directory(job_dir, dry_run=True):
    """
    清理单个 job 目录中的多余日志文件
    """
    if not os.path.isdir(job_dir):
        return
    
    file_groups = get_file_groups(job_dir)
    removed_count = 0
    
    print(f"\n处理目录: {job_dir}")
    
    for ext, files in file_groups.items():
        if len(files) <= 1:
            print(f"  {ext} 文件: {len(files)} 个，无需清理")
            continue
        
        # 找到最新的文件
        latest_file = get_latest_file(files)
        
        if latest_file is None:
            continue
        
        # 需要删除的文件
        files_to_remove = [f for f in files if f != latest_file]
        
        print(f"  {ext} 文件: {len(files)} 个，保留最新: {os.path.basename(latest_file)}")
        
        for file_to_remove in files_to_remove:
            print(f"    {'[DRY RUN] ' if dry_run else ''}删除: {os.path.basename(file_to_remove)}")
            if not dry_run:
                try:
                    os.remove(file_to_remove)
                    removed_count += 1
                except OSError as e:
                    print(f"    错误: 无法删除 {file_to_remove}: {e}")
            else:
                removed_count += 1
    
    return removed_count

def main():
    """
    主函数
    """
    # 获取基础路径
    base_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/eos_logs/signal/"
    # base_path = "/eos/user/j/jiehan/eos_logs/background/"
    
    if not os.path.exists(base_path):
        print(f"错误: 路径不存在: {base_path}")
        return 1
    
    # 检查命令行参数
    dry_run = True
    if len(sys.argv) > 1 and sys.argv[1] == "--execute":
        dry_run = False
        print("执行模式: 将实际删除文件")
    else:
        print("测试模式: 只显示将要删除的文件，使用 --execute 参数来实际执行删除")
    
    total_removed = 0
    processed_dirs = 0
    
    # 查找所有匹配的目录模式 *_20* (包含 2016, 2017, 2018 等)
    pattern = os.path.join(base_path, "*_20*")
    sample_dirs = glob.glob(pattern)
    
    if not sample_dirs:
        print(f"未找到匹配模式 *_20* 的目录在: {base_path}")
        return 1
    
    print(f"找到 {len(sample_dirs)} 个样本目录")
    
    for sample_dir in sorted(sample_dirs):
        sample_name = os.path.basename(sample_dir)
        
        # 查找所有 job_* 子目录
        job_pattern = os.path.join(sample_dir, "job_*")
        job_dirs = glob.glob(job_pattern)
        
        if not job_dirs:
            print(f"\n{sample_name}: 未找到 job_* 目录")
            continue
        
        print(f"\n{sample_name}: 找到 {len(job_dirs)} 个 job 目录")
        
        # 处理每个 job 目录
        for job_dir in sorted(job_dirs):
            removed = cleanup_job_directory(job_dir, dry_run)
            total_removed += removed
            processed_dirs += 1
    
    print(f"\n总结:")
    print(f"  处理的目录数: {processed_dirs}")
    print(f"  {'预计删除' if dry_run else '已删除'}的文件数: {total_removed}")
    
    if dry_run and total_removed > 0:
        print(f"\n如要实际执行删除，请运行: python {sys.argv[0]} --execute")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
