import os
import numpy as np
import pandas as pd
import re
from pdb import set_trace
import time

start_time = time.time()


# Specify the output file path
# output_file = "syn_list/cut_yields_output_Data_run2.txt"
output_file = "syn_list/cut_yields_output_Data_run3.txt"

eos_path = '/eos/home-p/pelai/HZgamma/Parquet/NanoV12/run3/'
log_path = '/afs/cern.ch/work/p/pelai/HZgamma/Output/HiggsDNA/Syn/eos_logs/run3/'
# Reading log files from eos is too slow

dataset_type = 'Data'
dataset_names = ["Data"]
# dataset_years = ["2016preVFP", "2016postVFP", "2017", "2018"] #"2016preVFP", "2016postVFP", "2017", "2018", "2023preBPix"
dataset_years = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"] #"2022preEE", "2022postEE", "2023preBPix", "2023postBPix"

# eos_path = '/eos/home-p/pelai/HZgamma/Parquet/NanoV9/run2/Sig_MC/'
# log_path = '/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/eos_logs/Sig_MC/'

# dataset_type = 'signal'
# dataset_type = 'WI_Systematic'
# dataset_names = ["ggH_M125"] #"ggH", "VBFH", "ZH", "ttH", "WplusH", "WminusH" "ggH_M125", "VBFH_M125", "ZH_M125", "ttH_M125", "WplusH_M125", "WminusH_M125"
# dataset_years = ["2017"]#"2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"

# dataset_type = 'bkgmc'
# dataset_names = ["DYJetsToLL"] # "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J" "ZGToLLG" "Data_SingleMuon", "Data_DoubleMuon", "Data_SingleElectron", "Data_DoubleEG"
# dataset_years = ["2017"] #"2016preVFP", "2016postVFP", "2017", "2018"]

cutflow_type = ['zgammas','zgammas_ele','zgammas_mu','zgammas_w','zgammas_ele_w','zgammas_mu_w']
type_num = len(cutflow_type)
cut_type = ['all', 'N_lep_sel','trig_cut','lep_pt_cut','has_g_cand', 'has_z_cand','sel_h_1','sel_h_2','sel_h_3', 'event']
cut_name = {
    'zgammas':['Initial Events', r'$N_{l} \geq 2$', 'e, ee trigger || $\mu\mu$, $\mu$ triggers ', r'lepton trigger pT cut', r'$N_{\gamma} \geq 1$', r'$80 \text{ GeV} < m_{ll} < 100 \text{ GeV}$', r'$p_{T}^{\gamma}/m_{ll\gamma} > 15/110$', r'$m_{ll} + m_{ll\gamma} > 185 \text{ GeV}$', r'$100 \text{ GeV} < m_{ll\gamma} < 180 \text{ GeV}$', 'Event Filtering', 'Total Baseline Events'],
    'zgammas_ele':['Initial Events', r'$N_{\textrm{e}} \geq 2$', 'ee, e triggers ', r'ee, e trigger pT cut', r'$N_{\gamma} \geq 1$', r'$80 \text{ GeV} < m_{\textrm{ee}} < 100 \text{ GeV}$', r'$p_{T}^{\gamma}/m_{\textrm{ee}\gamma} > 15/110$', r'$m_{\textrm{ee}} + m_{\textrm{ee}\gamma} > 185 \text{ GeV}$', r'$100 \text{ GeV} < m_{\textrm{ee}\gamma} < 180 \text{ GeV}$', 'Event Filtering', 'Total Baseline Events'],
    'zgammas_mu':['Initial Events', r'$N_{\mu} \geq 2$', '$\mu\mu$, $\mu$ triggers ', r'$\mu\mu$, $\mu$ trigger pT cut', r'$N_{\gamma} \geq 1$', r'$80 \text{ GeV} < m_{\mu\mu} < 100 \text{ GeV}$', r'$p_{T}^{\gamma}/m_{\mu\mu\gamma} > 15/110$', r'$m_{\mu\mu} + m_{\mu\mu\gamma} > 185 \text{ GeV}$', r'$100 \text{ GeV} < m_{\mu\mu\gamma} < 180 \text{ GeV}$', 'Event Filtering', 'Total Baseline Events'],
    'zgammas_w':['Initial Events', r'$N_{l} \geq 2$', 'e, ee trigger || $\mu\mu$, $\mu$ triggers ', r'lepton trigger pT cut', r'$N_{\gamma} \geq 1$', r'$80 \text{ GeV} < m_{ll} < 100 \text{ GeV}$', r'$p_{T}^{\gamma}/m_{ll\gamma} > 15/110$', r'$m_{ll} + m_{ll\gamma} > 185 \text{ GeV}$', r'$100 \text{ GeV} < m_{ll\gamma} < 180 \text{ GeV}$', 'Event Filtering', 'Total Baseline Events'],
    'zgammas_ele_w':['Initial Events', r'$N_{\textrm{e}} \geq 2$', 'ee, e triggers ', r'ee, e trigger pT cut', r'$N_{\gamma} \geq 1$', r'$80 \text{ GeV} < m_{\textrm{ee}} < 100 \text{ GeV}$', r'$p_{T}^{\gamma}/m_{\textrm{ee}\gamma} > 15/110$', r'$m_{\textrm{ee}} + m_{\textrm{ee}\gamma} > 185 \text{ GeV}$', r'$100 \text{ GeV} < m_{\textrm{ee}\gamma} < 180 \text{ GeV}$', 'Event Filtering', 'Total Baseline Events'],
    'zgammas_mu_w':['Initial Events', r'$N_{\mu} \geq 2$', '$\mu\mu$, $\mu$ triggers ', r'$\mu\mu$, $\mu$ trigger pT cut', r'$N_{\gamma} \geq 1$', r'$80 \text{ GeV} < m_{\mu\mu} < 100 \text{ GeV}$', r'$p_{T}^{\gamma}/m_{\mu\mu\gamma} > 15/110$', r'$m_{\mu\mu} + m_{\mu\mu\gamma} > 185 \text{ GeV}$', r'$100 \text{ GeV} < m_{\mu\mu\gamma} < 180 \text{ GeV}$', 'Event Filtering', 'Total Baseline Events'],
}
cut_num = len(cut_type)
cutflow = {'zgammas':np.array(np.zeros(cut_num)), 'zgammas_ele':np.array(np.zeros(cut_num)), 'zgammas_mu':np.array(np.zeros(cut_num)), 'zgammas_w':np.array(np.zeros(cut_num)), 'zgammas_ele_w':np.array(np.zeros(cut_num)), 'zgammas_mu_w':np.array(np.zeros(cut_num))}

# To renew the output file
with open(output_file, 'w') as f:
    # Output BEFORE replacement
    f.write("📌📌📌📌📌📌📌📌📌📌📌📌📌")

electron_yields_dict, muon_yields_dict = {}, {}
for dataset in dataset_names:
    for year in dataset_years:
        cutflow = {'zgammas':np.array(np.zeros(cut_num)), 'zgammas_ele':np.array(np.zeros(cut_num)), 'zgammas_mu':np.array(np.zeros(cut_num)), 'zgammas_w':np.array(np.zeros(cut_num)), 'zgammas_ele_w':np.array(np.zeros(cut_num)), 'zgammas_mu_w':np.array(np.zeros(cut_num))}
        weight = 1
        try:
            print("reading: {}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
            data = pd.read_parquet("{}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
            print(data["weight_central"].to_numpy().astype('float64'), data["weight_central_no_lumi"].to_numpy().astype('float64'))
            print("{}{}/{}_{}/merged_nominal.parquet".format(eos_path, dataset_type, dataset, year))
            if 'weight_central_initial' in data.keys():
                weight = data['weight_central'].to_numpy().astype('float64')[1]/data['weight_central_initial'].to_numpy().astype('float64')[1]
            else:
                weight = 1
                print("No weight exists, set it as 1.")
            del data
        except:
            print("No weight exists, set it as 1.")
        print(weight)

        # Storage dictionaries
        yields_dict = {}

        if not os.path.isdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
            continue      
        for log_dir in os.listdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
        # for i, log_dir in enumerate(os.listdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year))):
            # if i == 10:
                # break
            flag = 1
            log_dir = "{}{}/{}_{}/{}".format(log_path, dataset_type, dataset, year, log_dir)
            # print(log_dir)
            temp = np.zeros(cut_num)
            if not os.path.isdir(log_dir):
                continue
            # Step 1: Collect all .out files
            out_files = []
            for log_file in os.listdir(log_dir):
                if '.out' not in log_file:  # Skip files without .out
                    continue
                out_files.append(log_file)

            # Step 2: If there are .out files, find the one with the largest number
            if out_files:
                # Extract the numerical part (before '.out') and find the max
                largest_file = max(out_files, key=lambda x: float(x.split('.out')[0]))
                
                # Step 3: Process only the largest file
                typei, cuti = 0, 0
                f = open("{}/{}".format(log_dir, largest_file), 'r')

                lines = f.read().split("DEBUG")
                if len(lines) < (cut_num*3+40):
                    continue
                # print("reading: {}/{}".format(log_dir, log_file), end = ' ')
                
                # Define regex pattern
                pattern = re.compile(
                    r"\[Tagger\]\s*:.*?cut type\s*:\s*(\w+),\s*cut\s*:\s*"  # Captures cut type
                    r"([\w\s|<>.=!'\-]+?),?\s*yields\s*:?\s*([\d.]+)"  # Captures cut condition & yields, allowing `'`
                )

                for line in lines:
                    # line = line.replace("\n", " ")
                    line = re.sub(r"tagger\.py:190", "", line)
                    line = re.sub(r"tagger\.py:224", "", line)
                    line = re.sub(r"analysis\.py:\d+\s*", "", line)
                    
                    # Normalize spaces
                    line = " ".join(line.split())  # Removes extra spaces while keeping words intact

                    if "[Tagger]" in line: 
                        if "WARNING" in line or "[[INFO]]" in line:
                            continue
                        
                        # print(f"🔍 Checking line: {repr(line)}")  # Show exact formatting
                        
                        match = pattern.match(line)

                        if match:
                            # Extract values from the match
                            cut_type = match.group(1).strip()  # Normalize spaces & case
                            cut = match.group(2).strip().strip()
                            yields = float(match.group(3))

                            # print(f"🔍 Checking storage → Cut Type: {cut_type}, Cut: {cut}, Yields: {yields:.3f}")
                            
                            # Ensure `cut_type` exists in yields_dict
                            if cut_type not in yields_dict:
                                # print(f"📌 Creating new entry for cut type: {repr(cut_type)}")
                                yields_dict[cut_type] = {}  # Initialize as an empty dictionary

                            # Ensure `cut` exists and accumulates yields
                            if cut in yields_dict[cut_type]:
                                # print(f"🔄 Updating existing cut: {cut} (Previous: {yields_dict[cut_type][cut]}, Adding: {yields})")
                                yields_dict[cut_type][cut] += yields  # ✅ Accumulate the value
                            else:
                                # print(f"🆕 Adding new cut: {cut} with initial value {yields}")
                                yields_dict[cut_type][cut] = yields  # ✅ Initialize the first time
                            
                        # else:
                            # 🔍 Debugging output for missing matches
                            # print(f"❌ No match for line: {repr(line)}")
            f.close()

        # Output BEFORE replacement
        print("\n📌 Final Stored Data:")
        for cut_type, cuts in yields_dict.items():
            print(f"Cut Type: {cut_type}")
            for cut, yield_value in cuts.items():
                print(f"   - {cut:30} → {yield_value:.0f}")

        # Output AFTER replacement
        print("\n📌 Final Stored Data (After Replacement):")
        for cut_type_key, cuts in yields_dict.items():
            print(f"Cut Type: {cut_type_key}")
            if cut_type_key in cut_name:  # Check if cut_type_key exists in cut_name
                for i, (cut, yield_value) in enumerate(cuts.items()):
                    descriptive_name = cut_name[cut_type_key][i]  # Use mapped name
                    print(f"     {descriptive_name:30} & {yield_value:.0f} \\\\")
            else:
                # If not in cut_name, use the original cut names
                for cut, yield_value in cuts.items():
                    print(f"     {cut:30} & {yield_value:.0f} \\\\")


        # Open the file in write mode and redirect output
        with open(output_file, 'a') as f:
            # Output BEFORE replacement
            f.write(f"\n📌 Years {year}\n")
            
            f.write("\n📌 Final Stored Data (Before Replacement):\n")
            for cut_type_key, cuts in yields_dict.items():
                f.write(f"Cut Type: {cut_type_key}\n")
                for cut, yield_value in cuts.items():
                    f.write(f"   - {cut:30} → {yield_value:.0f}\n")

            # Output AFTER replacement
            f.write("\n📌 Final Stored Data (After Replacement):\n")
            for cut_type_key, cuts in yields_dict.items():
                f.write(f"Cut Type: {cut_type_key}\n")
                if cut_type_key in cut_name:  # Check if cut_type_key exists in cut_name
                    for i, (cut, yield_value) in enumerate(cuts.items()):
                        descriptive_name = cut_name[cut_type_key][i]  # Use mapped name
                        f.write(f"     {descriptive_name:30} & {yield_value:.0f} \\\\ \n")
                else:
                    # If not in cut_name, use the original cut names
                    for cut, yield_value in cuts.items():
                        f.write(f"     {cut:30} & {yield_value:.0f} \\\\ \n")

        print(f"Output has been saved to {output_file}")

        # print("\n🔹 Cutflow Summary\n")
        # for cut_type, cuts in yields_dict.items():
        #     print(f"📌 Cut Type: {cut_type}")
        #     for cut, yield_value in cuts.items():
        #         print(f"   - {cut:30} → {yield_value:,.0f}")  # Format with comma for readability
        #         print("-" * 50)  # Separator for readability

        # for name in cutflow:
        #     print(name)
        #     cuts = cut_name[name]
        #     for i, yields in enumerate(cutflow[name]):
        #         cut = cuts[i]
        #         if i == 1 and 'w' not in name:
        #             if 'mu' in name:
        #                 # print("muon_yields_dict keys:", muon_yields_dict.keys())
        #                 # print("Expected cut names:", [cuts[i] for i in range(len(cutflow[name]))])

        #                 print(muon_yields_dict)
        #                 print('\n'.join([f'{cut:>40} {yields:.0f}' for cut, yields in muon_yields_dict.items()]))
        #             if "ele" in name:
        #                 print(electron_yields_dict)
        #                 print('\n'.join([f'{cut:>40} {yields:.0f}' for cut, yields in electron_yields_dict.items()]))
        #         if 'w' in name:
        #             print(f'{cut:>40} {yields:.3f}')
        #         else:
        #             print(f'{cut:>40} {yields:.0f}')
        #     print(' ')
# print(cutflow['zgammas_ele'] + cutflow['zgammas_mu'])
end_time = time.time()
run_time = end_time - start_time
print(f"Processing time: {run_time:.2f} s")

# # events_sum = np.zeros(2)
# for dataset in dataset_names:
#     for year in dataset_years:
#         if not os.path.isdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
#             continue      
#         for log_dir in os.listdir("{}{}/{}_{}".format(log_path, dataset_type, dataset, year)):
#             flag = 1
#             log_dir = "{}{}/{}_{}/{}".format(log_path, dataset_type, dataset, year, log_dir)
#             if not os.path.isdir(log_dir):
#                 continue
#             for log_file in os.listdir(log_dir): 
#                 if '.out' not in log_file:
#                     continue
#                 f = open("{}/{}".format(log_dir, log_file), 'r')
#                 lines = f.read().split("DEBUG")
#                 if len(lines) < 10:
#                     continue;
#                 print("reading: {}/{}".format(log_dir, log_file))
#                 for line in lines:
#                     line = line.replace("\n", " ")
#                     # matched = re.match(r".*?Duplicated_samples:  selected samples:.*?(\d+).*?\((\d+)\).*?", line)
#                     # if matched:
#                     #     events_sum[0] += int(matched.group(1))
#                     #     temp = int(matched.group(2))
#                     #     print(line)
#                     matched = re.match(r".*?2022C*.?EGamma/NANOAOD/.*?", line)
#                     if matched:
#                         print(log_file)
#                         # events_sum[1] += temp
#                         break
#                 f.close()
# # print(events_sum)