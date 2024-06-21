import pandas as pd
import pyarrow.parquet as pq

# 指定合并后的 Parquet 文件路径
merged_file = '/eos/user/j/jiehan/parquet/nanov9/data/Data_EGamma_2018A_2018/merged.parquet'

# 使用 pyarrow 读取 Parquet 文件的分段
def read_parquet_in_chunks(file_path, chunk_size=1000):
    parquet_file = pq.ParquetFile(file_path)
    total_row_groups = parquet_file.num_row_groups

    for rg in range(total_row_groups):
        row_group = parquet_file.read_row_group(rg)
        df = row_group.to_pandas()

        # 分块处理每个row group的DataFrame
        for start in range(0, len(df), chunk_size):
            end = min(start + chunk_size, len(df))
            yield df.iloc[start:end]

cut_types = ["zgammas", "zgammas_ele", "zgammas_mu"]#, "zgammas_w", "zgammas_ele_w", "zgammas_mu_w"]
cut_selections = ["all", "N_lep_sel", "trig_cut", "lep_pt_cut", "has_g_cand", "os_cut", "has_z_cand", "sel_h_1", "sel_h_2", "sel_h_3"]
cut_flow = {cut_type: {cut_selection: 0 for cut_selection in cut_selections} for cut_type in cut_types}
num_events = 0
output_file = '/eos/user/j/jiehan/run_lumi_eve_no/Data_EGamma_2018A.txt'
# with open(output_file, 'w') as f:
#     for chunk in read_parquet_in_chunks(merged_file, chunk_size=10000):
#         # 在这里处理每个块的数据
#         num_events += len(chunk)
#         print(f"Processed {num_events} events")
#         for index, row in chunk.iterrows():
#             num_events += 1
#             if num_events % 100000 == 0:
#                 print(f"Processed {num_events} events")
#             f.write(f"{row['run']}\t{row['luminosityBlock']}\t{row['event']}\n")
for chunk in read_parquet_in_chunks(merged_file, chunk_size=50000):
    # 在这里处理每个块的数据
    data = chunk[(chunk['zgammas_ele_os_cut']==0) & (chunk['zgammas_ele_has_g_cand']==1)]
    # data = chunk
    if not data.empty:
        print(f"run: {data['run']}, lumi: {data['luminosityBlock']}, event: {data['event']}, lead lepton charge: {data['Z_lead_lepton_charge']}, sublead lepton charge: {data['Z_sublead_lepton_charge']}")
    num_events += len(chunk)
    # print(f"Processed {num_events} events")
    # for cut_type in cut_types:
    #     for cut_selection in cut_selections:
    #         if cut_selection == "all":
    #             cut_flow[cut_type][cut_selection] += len(chunk)
    #         else:
    #             cut_flow[cut_type][cut_selection] += len(chunk[(chunk[f'{cut_type}_{cut_selection}'])])        
print(f"Total number of events: {num_events}")
# for i in cut_flow:
#     print(f"\nCut type: {i}")
#     for j in cut_flow[i]:
#         print(f"{j:20s}\t & {cut_flow[i][j]}")
    
