import os
import numpy as np
import pandas as pd
import re
from pdb import set_trace
import time

start_time = time.time()

eos_path = '/eos/home-j/jiehan/parquet/nanov9/'
log_path = '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/eos_logs/'

# dataset_type = 'bkgmc'
# dataset_names = ["DYJetsToLL"]
# dataset_years = ["2016preVFP", "2016postVFP", "2017", "2018"]

# dataset_type = 'signal'
# dataset_names = ["ggH_M125"] #, "VBFH_M125", "ZH_M125", "ttH_M125"] # "WplusH_M125", "WminusH_M125",
# dataset_years = ["2017"]#, "2017", "2018"]

dataset_type = 'data'
dataset_names = ["Data"]
dataset_years = ["2016preVFP", "2016postVFP", "2017", "2018"]

cutflow_type = ['zgammas','zgammas_ele','zgammas_mu','zgammas_w','zgammas_ele_w','zgammas_mu_w']
type_num = len(cutflow_type)
cut_type = ['all', 'N_lep_sel','trig_cut','lep_pt_cut','has_g_cand', 'os_cut', 'has_z_cand','sel_h_1','sel_h_2','sel_h_3']
cut_num = len(cut_type)
cutflow = {'zgammas':np.array(np.zeros(cut_num)), 'zgammas_ele':np.array(np.zeros(cut_num)), 'zgammas_mu':np.array(np.zeros(cut_num)), 'zgammas_w':np.array(np.zeros(cut_num)), 'zgammas_ele_w':np.array(np.zeros(cut_num)), 'zgammas_mu_w':np.array(np.zeros(cut_num))}

for dataset in dataset_names:
    for year in dataset_years:
        try:
            data = pd.read_parquet("{}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
            if 'weight_central_initial' in data.keys():
                weight = data['weight_central'][1]/data['weight_central_initial'][1]
            else:
                weight = 1
            del data
        except:
            weight = 1
        print(weight)
        
        for log_dir in os.listdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
            log_dir = "{}{}/{}_{}/{}".format(log_path, dataset_type, dataset, year, log_dir)
            temp = np.zeros(cut_num)
            typei, cuti = 0, 0
            if not os.path.isdir(log_dir):
                continue
            for log_file in os.listdir(log_dir):
                # if '.log' in log_file:
                if '.out' not in log_file:
                    continue
                f = open("{}/{}".format(log_dir, log_file), 'r')
                lines = f.read().split("DEBUG")
                if len(lines) < (cut_num*6+20):
                    continue;
                print("reading: {}/{}".format(log_dir, log_file))
                for line in lines:
                    line = line.replace("\n", " ")
                    if "zgammas" not in line:
                        continue
                    if cuti == cut_num:
                        cutflow[cutflow_type[typei]] += temp
                        temp = np.zeros(cut_num)
                        cuti = 0
                        typei += 1
                        if typei == type_num:
                            break
                    yields = 0
                    match = re.search(r'{}.*?{}.*?yields\s*:\s*.*?(\d+\.\d+)'.format(cutflow_type[typei], cut_type[cuti]), line)
                    if match:
                        yields = float(match.group(1))
                    if 'w' in cutflow_type[typei]:
                        temp[cuti] += yields*weight
                    else:
                        temp[cuti] += yields
                    cuti += 1
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
end_time = time.time()
run_time = end_time - start_time
print(f"Processing time: {run_time:.2f} s")
