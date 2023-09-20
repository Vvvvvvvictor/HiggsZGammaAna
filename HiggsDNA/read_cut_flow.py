import os
import numpy as np
import pandas as pd
import re

# log_path = '/eos/home-j/jiehan/parquet/2022/mva_based/data/Data_E_2022/'
# log_path = '/afs/cern.ch/user/j/jiehan/private/HiggsDNA/cutflow/DYJetsToLL_amc_2017/'
# log_path = '/afs/cern.ch/user/j/jiehan/private/HiggsDNA/cutflow/ggH_M125_2017/'
log_path = '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/eos_logs/data/Data_E_2022/'
cutflow = {'zgammas':np.array(np.zeros(10)), 'zgammas_ele':np.array(np.zeros(10)), 'zgammas_mu':np.array(np.zeros(10)), 'zgammas_w':np.array(np.zeros(10)), 'zgammas_ele_w':np.array(np.zeros(10)), 'zgammas_mu_w':np.array(np.zeros(10))}
cutflow_type = ['zgammas','zgammas_ele','zgammas_mu','zgammas_w','zgammas_ele_w','zgammas_mu_w']
cut_type = ['all', 'N_lep_sel','trig_cut','lead_lep_pt_cut','sub_lep_pt_cut','has_g_cand','has_z_cand','sel_h_1','sel_h_2','sel_h_3']

# data = pd.read_parquet(log_path+'merged_nominal.parquet')
# weight = data['weight_central'][1]/data['weight_central_initial'][1]
# # print(data)
# print(weight)
# del data

# weight = 0.00002678185575243
weight = 1

temp = np.zeros(10)
cuti = 0

for log_dir in os.listdir(log_path):
    log_dir = log_path + log_dir
    if os.path.isdir(log_dir):
        log_dir = log_dir + '/'
        for log_file in os.listdir(log_dir):
            # if '.log' in log_file:
            if '.out' in log_file:
                log_file = log_dir + log_file
                with open(log_file, 'r') as f:
                    lines = f.read().split("DEBUG")
                    if len(lines) < 10:
                        continue;
                    print(log_file)
                    for line in lines:
                        line = line.replace("\n", " ")
                        for cutflow_t in cutflow_type:
                            if cutflow_t + "," in line:
                                if cut_type[cuti] + "," in line:
                                    # print(line)
                                    # split_line = re.split(":|,", line)
                                    split_line = line.split(" ")
                                    for index, item in enumerate(split_line):
                                        if item == "yields":
                                            indice = index
                                    yields = None
                                    step = 1
                                    num_cand = "0"
                                    while yields is None:
                                        try: 
                                            num_cand = split_line[indice+step]
                                        except:
                                            print("No more in this block!")
                                        if "." in num_cand:
                                            yields = float(num_cand)
                                        else:
                                            step+=1
                                    if 'w' in cutflow_t:
                                        temp[cuti] += yields*weight
                                    else:
                                        temp[cuti] += yields
                                    if cuti == 9:
                                        cutflow[cutflow_t] += temp
                                        temp = np.zeros(10)
                                        cuti = 0
                                        continue
                                    cuti += 1
                                    continue
                    f.close()
for i in cutflow:
    print(i)
    for j in cutflow[i]:
        if 'w' in i:
            print('%.3f'%j)
        else:
            print(int(j))
    print(' ')