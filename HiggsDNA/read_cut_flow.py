import os
import numpy as np
import pandas as pd
import re

eos_path = '/eos/home-j/jiehan/parquet/nanov9/'
log_path = '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/eos_logs/'

# dataset_type = 'bkgmc'
# dataset_names = ["DYJetsToLL", "ZGToLLG"]
# dataset_years = ["2016", "2017", "2018"]

# dataset_type = 'signal'
# dataset_names = ["ggH_M125", "VBFH_M125", "ZH_M125", "ttH_M125"] # "WplusH_M125", "WminusH_M125",
# dataset_years = ["2018"]#, "2017", "2018"]

dataset_type = 'data'
dataset_names = ["Data"]
dataset_years = ["2016", "2017", "2018"]

cutflow_type = ['zgammas','zgammas_ele','zgammas_mu','zgammas_w','zgammas_ele_w','zgammas_mu_w']
cut_type = ['all', 'N_lep_sel','trig_cut','lep_pt_cut','has_g_cand', 'os_cut', 'has_z_cand','sel_h_1','sel_h_2','sel_h_3']
cut_num = len(cut_type)
cutflow = {'zgammas':np.array(np.zeros(cut_num)), 'zgammas_ele':np.array(np.zeros(cut_num)), 'zgammas_mu':np.array(np.zeros(cut_num)), 'zgammas_w':np.array(np.zeros(cut_num)), 'zgammas_ele_w':np.array(np.zeros(cut_num)), 'zgammas_mu_w':np.array(np.zeros(cut_num))}

temp = np.zeros(cut_num)
cuti = 0
for dataset in dataset_names:
    for year in dataset_years:
        # try:
        #     data = pd.read_parquet("{}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
        # except:
        #     continue
        # weight = data['weight_central'][1]/data['weight_central_initial'][1]
        weight = 1
        # print(data)
        print(weight)
        # del data
        for log_dir in os.listdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
            log_dir = "{}{}/{}_{}/{}".format(log_path, dataset_type, dataset, year, log_dir)
            if os.path.isdir(log_dir):
                for log_file in os.listdir(log_dir):
                    # if '.log' in log_file:
                    if '.out' in log_file:
                        f = open("{}/{}".format(log_dir, log_file), 'r')
                        lines = f.read().split("DEBUG")
                        if len(lines) < cut_num:
                            continue;
                        print("reading: {}/{}".format(log_dir, log_file))
                        for line in lines:
                            line = line.replace("\n", " ")
                            for cutflow_t in cutflow_type:
                                if cutflow_t + "," in line:
                                    # print(line)
                                    if cut_type[cuti] + "," in line:
                                        # print(line)
                                        # split_line = re.split(":|,", line)
                                        split_line = line.split(":")
                                        for index, item in enumerate(split_line):
                                            if "yields" in item:
                                                # print(item)
                                                indice = index
                                        yields = None
                                        step = 1
                                        num_cand = "0"
                                        while yields is None:
                                            try: 
                                                num_cand = split_line[indice+step]
                                                # print(num_cand)
                                            except:
                                                print("No more in this block!")
                                                break
                                            if "." in num_cand:
                                                yields = float(num_cand)
                                            else:
                                                step+=1
                                        if 'w' in cutflow_t:
                                            temp[cuti] += yields*weight
                                        else:
                                            temp[cuti] += yields
                                        if cuti == cut_num-1:
                                            # print(cutflow[cutflow_t])
                                            cutflow[cutflow_t] += temp
                                            temp = np.zeros(cut_num)
                                            cuti = 0
                                            continue
                                        cuti += 1
                                        continue
                        f.close()
                        break
            # break
for i in cutflow:
    print(i)
    for j in cutflow[i]:
        if 'w' in i:
            print('%.3f'%j)
        else:
            print(int(j))
    print(' ')