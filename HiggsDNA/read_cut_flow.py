import os
import numpy as np
import pandas as pd
import re
from pdb import set_trace
import time

start_time = time.time()

eos_path = '/eos/home-j/jiehan/parquet/nanov9/'
# log_path = '/eos/user/j/jiehan/eos_logs/'
log_path = '/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/eos_logs/'

# dataset_type = 'bkgmc'
# dataset_names = ["DYJetsToLL"]
# dataset_years = ["2016preVFP", "2016postVFP", "2017", "2018"]

# dataset_type = 'signal'
# dataset_names = ["ggH_M125"] #, "VBFH_M125", "ZH_M125", "ttH_M125"] # "WplusH_M125", "WminusH_M125",
# dataset_years = ["2017"]#, "2017", "2018"]

dataset_type = 'data'
dataset_names = ["Data_DoubleMuon"] # "Data_SingleMuon", "Data_DoubleMuon", "Data_SingleElectron", "Data_DoubleEG"
dataset_years = ["2017"] #2016preVFP", "2016postVFP", "2017", "2018"]

cutflow_type = ['zgammas','zgammas_ele','zgammas_mu','zgammas_w','zgammas_ele_w','zgammas_mu_w']
type_num = len(cutflow_type)
cut_type = ['all', 'N_lep_sel','trig_cut','lep_pt_cut','has_g_cand', 'os_cut', 'has_z_cand','sel_h_1','sel_h_2','sel_h_3']
cut_name = {
    'zgammas':['No selction', r'$N_{l}\geq2$', r'e, ee trigger || $\mu$, $\mu\mu$ trigger', r'lepton trigger $p_{T}$ cut', r'$N_{\gamma}\geq1$', 'Opposite charge', r'80GeV<m_{ll}<100GeV', r'p_{T}^{\gamma}>15./110.', r'm_{ll}+m_{ll\gamma}<185GeV', r'100GeV<m_{ll\gamma}<180GeV'], 
    'zgammas_ele':['No selction', r'$N_{e}\geq2$', 'e, ee trigger', r'e trigger $p_{T}$ cut', r'$N_{\gamma}\geq1$', 'Opposite charge', r'80GeV<$m_{ee}$<100GeV', r'$p_{T}^{\gamma}$>15./110.', r'$m_{ee}+m_{ee\gamma}$<185GeV', r'100GeV<$m_{ee\gamma}$<180GeV'],
    'zgammas_mu':['No selction', r'$N_{\mu}\geq2$', r'$\mu$, $\mu\mu$ trigger', r'$\mu$ trigger $p_{T}$ cut', r'$N_{\gamma}\geq1$', 'Opposite charge', r'80GeV<$m_{\mu\mu}$<100GeV', r'$p_{T}^{\gamma}$>15./110.', r'$m_{\mu\mu}+m_{\mu\mu\gamma}$<185GeV', r'100GeV<$m_{\mu\mu\gamma}$<180GeV'], 
    'zgammas_w':['No selction', r'$N_{l}\geq2$', r'e, ee trigger || $\mu$, $\mu\mu$ trigger', r'lepton trigger $p_{T}$ cut', r'$N_{\gamma}\geq1$', 'Opposite charge', r'80GeV<m_{ll}<100GeV', r'p_{T}^{\gamma}>15./110.', r'm_{ll}+m_{ll\gamma}<185GeV', r'100GeV<m_{ll\gamma}<180GeV'],  
    'zgammas_ele_w':['No selction', r'$N_{e}\geq2$', 'e, ee trigger', r'e trigger $p_{T}$ cut', r'$N_{\gamma}\geq1$', 'Opposite charge', r'80GeV<$m_{ee}$<100GeV', r'$p_{T}^{\gamma}$>15./110.', r'$m_{ee}+m_{ee\gamma}$<185GeV', r'100GeV<$m_{ee\gamma}$<180GeV'],
    'zgammas_mu_w':['No selction', r'$N_{\mu}\geq2$', r'$\mu$, $\mu\mu$ trigger', r'$\mu$ trigger $p_{T}$ cut', r'$N_{\gamma}\geq1$', 'Opposite charge', r'80GeV<$m_{\mu\mu}$<100GeV', r'$p_{T}^{\gamma}$>15./110.', r'$m_{\mu\mu}+m_{\mu\mu\gamma}$<185GeV', r'100GeV<$m_{\mu\mu\gamma}$<180GeV']
}
cut_num = len(cut_type)
cutflow = {'zgammas':np.array(np.zeros(cut_num)), 'zgammas_ele':np.array(np.zeros(cut_num)), 'zgammas_mu':np.array(np.zeros(cut_num)), 'zgammas_w':np.array(np.zeros(cut_num)), 'zgammas_ele_w':np.array(np.zeros(cut_num)), 'zgammas_mu_w':np.array(np.zeros(cut_num))}

electron_yields_dict, muon_yields_dict = {}, {}
for dataset in dataset_names:
    for year in dataset_years:
        weight = 1
        # try:
        #     data = pd.read_parquet("{}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
        #     print("{}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
        #     if 'weight_central_initial' in data.keys():
        #         weight = data['weight_central'][1]/data['weight_central_initial'][1]
        #     else:
        #         weight = 1
        #         print("No weight exists, set it as 1.")
        #     del data
        # except:
        #     print("No weight exists, set it as 1.")
        print(weight)
        
        for log_dir in os.listdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
            flag = 1
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
                if len(lines) < (cut_num*3+20):
                    continue;
                print("reading: {}/{}".format(log_dir, log_file))
                for line in lines:
                    line = line.replace("\n", " ")
                    # print(line)
                    if "!!!!charge of lead lepton" in line and "[]" not in line:
                        print(line)
                    match = re.match(r".*?cut type : SelectedElectron,.*? cut.\s*:\s*.*?(.+),.*? yields\s*:\s*.*?(\d+)", line)
                    if match:
                        cut = match.group(1)
                        yields = int(match.group(2))
                        if cut in electron_yields_dict:
                            electron_yields_dict[cut] += yields
                        else:
                            electron_yields_dict[cut] = yields
                    # set_trace()
                    match = re.match(r".*?cut type : SelectedMuon,.*? cut.\s*:\s*.*?(.+),.*? yields\s*:\s*.*?(\d+)", line)
                    if match:
                        cut = match.group(1)
                        yields = int(match.group(2))
                        if cut not in muon_yields_dict:
                            muon_yields_dict[cut] = yields
                        else:
                            muon_yields_dict[cut] += yields
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
                        flag = 0
                        yields = float(match.group(1))
                    else:
                        continue
                    if 'w' in cutflow_type[typei]:
                        temp[cuti] += yields*weight
                    else:
                        temp[cuti] += yields
                    cuti += 1
                f.close()
                break
            if flag:
                print(log_dir)
            # for i in range(cut_num):
            #     if cutflow['zgammas_mu'][i] + cutflow['zgammas_ele'][i] != cutflow['zgammas'][i]:
            #         print(log_dir)
            # break
    lepton_flag = 0

for name in cutflow:
    print(name)
    cuts = cut_name[name]
    for i, yields in enumerate(cutflow[name]):
        cut = cuts[i]
        # if i == 1 and 'w' not in name:
        #     if 'mu' in name:
        #         # print(muon_yields_dict)
        #         print('\n'.join([f'{cut:>40} {yields:.0f}' for cut, yields in muon_yields_dict.items()]))
        #     if "ele" in name:
        #         # print(electron_yields_dict)
        #         print('\n'.join([f'{cut:>40} {yields:.0f}' for cut, yields in electron_yields_dict.items()]))
        if 'w' in name:
            print(f'{cut:>40} {yields:.3f}')
        else:
            print(f'{cut:>40} {yields:.0f}')
    print(' ')
end_time = time.time()
run_time = end_time - start_time
print(f"Processing time: {run_time:.2f} s")
