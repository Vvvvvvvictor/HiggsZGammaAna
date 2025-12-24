import os
import re
import argparse
from glob import glob
from collections import defaultdict
import pandas as pd

def parse_yields_from_file(file_path):
    """
    Parses a single log file by linearizing its content and then extracting 
    the 'yields' for each cut of 'SelectedJet'.
    Returns a dictionary of {cut_name: yield_value}.
    """
    yields = {}
    try:
        with open(file_path, 'r') as f:
            # Read the whole file and replace newlines and multiple spaces with a single space.
            content = f.read()
            linearized_content = ' '.join(content.split())
            
            # This regex is now much simpler and more robust thanks to the linearization.
            regex = re.compile(r"syst variation : nominal, cut type : SelectedJet, cut : (.*?), yields : (\d+)")
            matches = regex.finditer(linearized_content)
            
            for match in matches:
                # The cut name is cleaned by stripping any leading/trailing whitespace.
                cut_name = match.group(1).strip()
                yield_value = int(match.group(2))
                yields[cut_name] = yield_value
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
    return yields

def main():
    parser = argparse.ArgumentParser(description="Extract and sum yields from log files for a specific year to generate a cut flow.")
    parser.add_argument("year", default="2018", type=str, help="The year to process (e.g., 2018).")
    args = parser.parse_args()

    year = args.year
    # Using glob to find all matching log files based on the user-provided path structure.
    path_pattern = f'HiggsDNA/nanov12/signal_test/ggH_M125_{year}/job_1/*.log'
    # path_pattern = f'HiggsDNA/nanov12/data_test/Data_{year}/job_*/*.log'
    log_files = glob(path_pattern)

    if not log_files:
        print(f"No log files found for the year {year} with pattern: {path_pattern}")
        return

    total_yields = defaultdict(int)
    cut_order = None
    
    for file_path in log_files:
        yields = parse_yields_from_file(file_path)
        
        # Establish the cut order from the first file that has yields.
        if cut_order is None and yields:
            cut_order = list(yields.keys())

        for cut_name, yield_value in yields.items():
            total_yields[cut_name] += yield_value

    if cut_order is None:
        print("No 'SelectedJet' yields found in any log files.")
        return

    # Prepare data for pandas DataFrame, ensuring the established order is respected.
    data = {'Cut': [], 'Total Yield': []}
    for cut_name in cut_order:
        if cut_name in total_yields:
            data['Cut'].append(cut_name)
            data['Total Yield'].append(total_yields[cut_name])

    df = pd.DataFrame(data)

    print(f"\n--- SelectedJet Cut Flow for {year} ---")
    print(df.to_string(index=False))
    print("---------------------------------------\n")

if __name__ == "__main__":
    main()
