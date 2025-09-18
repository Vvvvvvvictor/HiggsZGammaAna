#!/usr/bin/env python3
import os

def get_input_base(year):
    if any(run2_year in year for run2_year in ["2016", "2017", "2018"]):
        return "/eos/user/j/jiehan/parquet/nanov9/"
    return "/eos/user/j/jiehan/parquet/nanov12/"

def get_systs(year):
    if any(run2_year in year for run2_year in ["2016", "2017", "2018"]):
        return ["Photon_scale", "Photon_smear", "Electron_scale", "Electron_smear", "JER", "JES", "MET_JES", "MET_Unclustered", "Muon_pt_smear"]
    return ["Photon_scale", "Photon_smear", "Electron_scale", "Electron_smear", "JER", "JES", "MET_JES", "MET_Unclustered", "Muon_pt_scale", "Muon_pt_smear"]

def generate_condor_submission():
    # --- Configuration ---
    target_base = "/eos/home-j/jiehan/root/skimmed_ntuples/"
    script_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py"
    
    years = ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"] #["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]

    signal_samples = ["ggH_M125", "VBF_M125", "WplusH_M125", "WminusH_M125", "ZH_M125", "ttH_M125"]
    bkg_samples = [
        # "ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J", "EWKZ2J",
        # "DYGto2LG_10to50", "DYGto2LG_50to100", "DYGto2LG_10to100"
    ]
    data_samples = ["Data"]

    all_commands = []

    # 1. Process nominal signal samples
    for sample in signal_samples:
        for year in years:
            input_base = get_input_base(year)
            input_file = f"{input_base}signal/{sample}_{year}/merged_nominal.parquet"
            if not os.path.exists(input_file):
                print(f"Warning: Input file not found, skipping: {input_file}")
                continue
            output_dir = f"{target_base}{sample}/"
            os.makedirs(output_dir, exist_ok=True)
            output_file = f"{output_dir}{year}.root"
            command = f"python {script_path} -i {input_file} -o {output_file}"
            all_commands.append(command)

    # 2. Process systematic variations for signal samples
    for sample in signal_samples:
        for year in years:
            systs = get_systs(year)
            input_base = get_input_base(year)
            for syst in systs:
                for uod in ["up", "down"]:
                    corr = f"{syst}_{uod}"
                    input_file = f"{input_base}signal/{sample}_{year}/merged_{corr}.parquet"
                    if not os.path.exists(input_file):
                        print(f"Warning: Input file not found, skipping: {input_file}")
                        continue
                    output_dir = f"{target_base}{sample}_{syst}_{uod}/"
                    os.makedirs(output_dir, exist_ok=True)
                    output_file = f"{output_dir}{year}.root"
                    command = f"python {script_path} -i {input_file} -o {output_file}"
                    all_commands.append(command)

    # 3. Process nominal background samples
    for sample in bkg_samples:
        for year in years:
            input_base = get_input_base(year)
            input_file = f"{input_base}background/{sample}_{year}/merged_nominal.parquet"
            if not os.path.exists(input_file):
                print(f"Warning: Input file not found, skipping: {input_file}")
                continue
            output_dir = f"{target_base}{sample}/"
            os.makedirs(output_dir, exist_ok=True)
            output_file = f"{output_dir}{year}.root"
            command = f"python {script_path} -i {input_file} -o {output_file}"
            all_commands.append(command)
            
    # 4. Process data samples
    for sample in data_samples:
        for year in years:
            input_base = get_input_base(year)
            # Data samples might have a different naming convention in the path
            input_file = f"{input_base}data/{sample}_{year}/merged_nominal.parquet"
            if not os.path.exists(input_file):
                print(f"Warning: Input file not found, skipping: {input_file}")
                continue
            output_dir = f"{target_base}Data/"
            os.makedirs(output_dir, exist_ok=True)
            output_file = f"{output_dir}{year}.root"
            command = f"python {script_path} -i {input_file} -o {output_file}"
            all_commands.append(command)

    # --- Setup output directories and file paths ---
    current_dir = os.getcwd()
    base_log_dir = os.path.join(current_dir, "eos_log")
    condor_log_dir = os.path.join(base_log_dir, "skimmed_ntuple_log")
    os.makedirs(condor_log_dir, exist_ok=True)

    commands_file_path = os.path.join(base_log_dir, "condor_commands.txt")
    sub_file_path = os.path.join(base_log_dir, "submit.sub")

    # --- Write commands to a file ---
    with open(commands_file_path, "w") as f:
        for cmd in all_commands:
            f.write(f'{cmd}\n')

    # --- Create Condor submission file ---
    sub_file_content = f"""
executable              = /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/run_condor_job.sh
arguments               = "$(command)"
output                  = {condor_log_dir}/job_$(ClusterId)_$(ProcId).out
error                   = {condor_log_dir}/job_$(ClusterId)_$(ProcId).err
log                     = {condor_log_dir}/job_$(ClusterId).log
+JobFlavour             = "longlunch"

# Use this to avoid transferring large root files
should_transfer_files   = NO

queue command from {commands_file_path}
"""
    
    with open(sub_file_path, "w") as f:
        f.write(sub_file_content)

    print(f"Generated {len(all_commands)} commands in '{commands_file_path}'")
    print(f"Generated Condor submission file '{sub_file_path}'")
    print(f"To submit, run: condor_submit {sub_file_path}")

if __name__ == "__main__":
    generate_condor_submission()
