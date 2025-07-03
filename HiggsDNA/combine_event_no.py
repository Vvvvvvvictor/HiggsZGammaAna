import os
import re

for dataset in ["ggH_M125_2023postBPix"]: #"Data_SingleElectron_2017", "Data_SingleMuon_2017", "Data_DoubleEG_2017", "Data_DoubleMuon_2017"
    # Define folder path
    # folder_path = f"/eos/user/j/jiehan/eos_logs/data/{dataset}"
    # folder_path = f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/eos_logs/data/{dataset}"
    folder_path = f"/eos/user/j/jiehan/parquet/cutflow_ggf/{dataset}"

    for event_type in ["electron"]:
        # Define output file path
        output_file = f"/eos/user/j/jiehan/run_lumi_eve_no/{dataset}_{event_type}.txt"

        # Define regular expression patterns
        start_pattern = f"!!!start check events tag({event_type})!!!"
        end_pattern = f"!!!end check events tag({event_type})!!!"

        # Open the output file in write mode
        with open(output_file, 'w') as out_file:

            # Traverse all subfolders in the folder
            for root, dirs, files in os.walk(folder_path):
                writed = False
                for file_name in files:
                    # print(file_name)
                    if ".log" not in file_name:
                    # if ".out" not in file_name:
                        continue
                    # Concatenate the file path
                    file_path = os.path.join(root, file_name)

                    print(f"Reading {file_path}...")
                    
                    # Read the file content
                    with open(file_path, 'r') as in_file:
                        lines = in_file.read()

                        # Use regular expression to extract content between start and end patterns
                        match = re.search(f"{re.escape(start_pattern)}(.*?){re.escape(end_pattern)}", lines, re.DOTALL)
                        if match:
                            # Extracted content is in the first group of the matching object
                            extracted_content = match.group(1).strip()
                            
                            # Write the extracted content to file
                            out_file.write(extracted_content)
                            if len(extracted_content) > 0:
                                out_file.write("\n")
                            writed = True
                            break
                            # print(f"Extracted content has been written to {output_file}.")
                if not writed:
                    print(f"{root} No matching content found.")

        print(f"Merging {dataset}_{event_type} completed!")
