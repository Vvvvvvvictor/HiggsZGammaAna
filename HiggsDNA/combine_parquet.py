import os
import pandas as pd
from tqdm import tqdm
import psutil

def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss

def merge_parquet_files(folder_path, output_file):
    # 获取文件夹中所有子文件夹的路径
    subfolders = [f.path for f in os.scandir(folder_path) if f.is_dir()]
    
    # 获取第一个子文件夹中的第一个Parquet文件的schema
    first_parquet_file = next((f.path for f in os.scandir(subfolders[0]) if f.is_file() and f.name.endswith('.parquet')), None)
    if not first_parquet_file:
        print("No Parquet files found in the folder")
        return
    
    # 合并每个子文件夹中的 Parquet 文件
    for subfolder in tqdm(subfolders, desc="Combining parquet files:", bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
        parquet_files = [f.path for f in os.scandir(subfolder) if f.is_file() and f.name.endswith('.parquet')]
        # 用于存储所有 Parquet 文件的 DataFrame
        dfs = []
        for parquet_file in parquet_files:
            # 逐个读取 Parquet 文件并将其转换为 DataFrame
            df = pd.read_parquet(parquet_file)
            dfs.append(df)
            
            # 输出程序占用的内存空间
            print("Memory usage: {}".format(sizeof_fmt(get_memory_usage())))
            
            # 清理内存
            del df
        
        # 合并所有 DataFrame
        combined_df = pd.concat(dfs, ignore_index=True)
        
        # 将合并后的 DataFrame 写入到输出文件中
        combined_df.to_parquet(output_file, index=False)
        
        # 清理内存
        del combined_df

# 文件夹路径和输出文件名
folder_path = "/eos/home-j/jiehan/parquet/nanov9/data/Data_2018"
output_file = "/eos/home-j/jiehan/parquet/nanov9/data/Data_2018/merged_nominal.parquet"

# 执行合并
merge_parquet_files(folder_path, output_file)
