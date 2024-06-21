import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# 指定文件夹路径
folder_path = '/eos/user/j/jiehan/parquet/nanov9/data/Data_EGamma_2018A_2018/'
output_file = os.path.join(folder_path, 'merged.parquet')

# 查找所有符合条件的文件
parquet_files = []
for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.startswith('output_job_') and file.endswith('_nominal.parquet'):
            parquet_files.append(os.path.join(root, file))

# 检查是否找到了符合条件的文件
if not parquet_files:
    print("No matching parquet files found.")
else:
    # 初始化写入器
    writer = None

    # 逐个读取并追加到合并的 Parquet 文件
    for file in parquet_files:
        print(f"Processing {file}...")
        df = pd.read_parquet(file)
        table = pa.Table.from_pandas(df)

        if writer is None:
            # 对于第一个文件，创建新的 Parquet 文件
            writer = pq.ParquetWriter(output_file, table.schema)
        writer.write_table(table)
    
    if writer:
        writer.close()
    
    print(f"Merged file saved to {output_file}")
