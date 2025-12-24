import uproot
import numpy as np
import pandas as pd

var_map = {
    "event": "event",
    "n_jets": "n_jets",
    "nbdfm": "n_b_jets",
    "jet_1_pt": "jet_1_pt",
    "jet_1_eta": "jet_1_eta",
    "jet_2_pt": "jet_2_pt",
    "jet_2_eta": "jet_2_eta",
    "j2_phi": "jet_2_phi",
    "j3_pt": "jet_3_pt",
    "j3_eta": "jet_3_eta",
    "j3_phi": "jet_3_phi",
    "met": "MET_pt",
    "nel": "n_electrons",
    "nmu": "n_muons",
}

# match the events with `event` variable in `zero_to_one_jet` tree

# file1: /eos/user/j/jiehan/root/skimmed_ntuples_rui_new/Data/2023preBPix.root
# file2: /eos/user/j/jiehan/root/skimmed_ntuples/Data/2023preBPix.root
# var_map = {var in file1: var in file2}

# store the following log in a text file with good formatting(\t) in /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/synchronization_script/log/

# for matched events in zero_to_one_jet tree in file1
# print out var. value in one row in file1
# then the value from file2
# the difference value in the third row
# a blank row between two event

# After finishing matched events

# for non-matched events
# find it in `inclusive` tree in file2
# print out var. value in one row in file1
# then the value from file2
# the difference value in the third row

file1_path = "/eos/user/j/jiehan/root/skimmed_ntuples_rui_new/Data/2023postBPix.root"
file2_path = "/eos/user/j/jiehan/root/skimmed_ntuples/Data/2023postBPix.root"
log_file_path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/synchronization_script/log/event_comparison_2023postBPix.log"

def get_df(file_path, tree_name, branches):
    with uproot.open(file_path) as f:
        tree = f[tree_name]
        df = tree.arrays(branches, library="pd")
    return df

def format_event_info(df, event_id, columns):
    event_data = df[df["event"] == event_id]
    if event_data.empty:
        return "Event not found"
    
    info_str = f"Event: {event_id}\n"
    for col in columns:
        info_str += f"{col}: {event_data[col].iloc[0]:.4f}\t"
    return info_str.strip()

def format_diff_info(df1, df2, event_id, var_map):
    event1 = df1[df1["event"] == event_id]
    event2 = df2[df2["event"] == event_id]
    if event1.empty or event2.empty:
        return "Event not found in one of the files"

    diff_str = "Differences:\n"
    has_diff = False
    for var1, var2 in var_map.items():
        val1 = event1[var1].iloc[0]
        val2 = event2[var2].iloc[0]
        diff = val1 - val2
        if abs(diff) > 1e-6: # Use a threshold for float comparison
            has_diff = True
            diff_str += f"{var1}-{var2}: {diff:.4f}\t"
    
    if not has_diff:
        return "No significant differences."
    return diff_str.strip()


def main():
    branches1 = list(var_map.keys())
    branches2 = list(var_map.values())

    df1_two_jet = get_df(file1_path, "two_jet", branches=branches1)
    df2_two_jet = get_df(file2_path, "two_jet", branches=branches2)

    events1 = set(df1_two_jet["event"])
    events2 = set(df2_two_jet["event"])

    # The user wants to find events that are in file2 but not in file1.
    unmatched_events = sorted(list(events2.difference(events1)))

    with open(log_file_path, "w") as log_file:
        log_file.write("="*50 + "\n")
        log_file.write("Processing events in file2's two_jet but not in file1's\n")
        log_file.write("="*50 + "\n\n")

        if unmatched_events:
            # For these events, find them in file1's zero_to_one_jet tree
            df1_fallback = get_df(file1_path, "zero_to_one_jet", branches=branches1)
            events1_fallback = set(df1_fallback["event"])

            branches2_extra = branches2 + ['run', 'luminosityBlock', "Z_lead_lepton_pt", "Z_lead_lepton_eta", "Z_lead_lepton_phi", "Z_sublead_lepton_pt", "Z_sublead_lepton_eta", "Z_sublead_lepton_phi"]
            df2_two_jet_extra = get_df(file2_path, "two_jet", branches=branches2_extra)

            for event_id in unmatched_events:
                if event_id in events1_fallback:
                    log_file.write(f"--- Event: {event_id} (in file2:two_jet, in file1:zero_to_one_jet) ---\n")

                    # File 2 (two_jet)
                    row2_str = "File2 (two_jet):\t"
                    event2_data = df2_two_jet[df2_two_jet["event"] == event_id]
                    for var in branches2:
                        row2_str += f"{var}: {event2_data[var].iloc[0]:.4f}\t"
                    log_file.write(row2_str + "\n")

                    # File 1 (zero_to_one_jet)
                    row1_str = "File1 (zero_to_one_jet):\t"
                    event1_data = df1_fallback[df1_fallback["event"] == event_id]
                    for var in branches1:
                        row1_str += f"{var}: {event1_data[var].iloc[0]:.4f}\t"
                    log_file.write(row1_str + "\n")
                    
                    # extra info for file2
                    event2_extra_data = df2_two_jet_extra[df2_two_jet_extra["event"] == event_id]
                    extra_info = f"File2 extra info:\t"
                    extra_info += f"run: {event2_extra_data['run'].iloc[0]}\t"
                    extra_info += f"luminosityBlock: {event2_extra_data['luminosityBlock'].iloc[0]}\t"
                    extra_info += f"event: {event2_extra_data['event'].iloc[0]}\n"
                    extra_info += f"Z_lead_lepton_pt: {event2_extra_data['Z_lead_lepton_pt'].iloc[0]:.4f}\t"
                    extra_info += f"Z_lead_lepton_eta: {event2_extra_data['Z_lead_lepton_eta'].iloc[0]:.4f}\t"
                    extra_info += f"Z_lead_lepton_phi: {event2_extra_data['Z_lead_lepton_phi'].iloc[0]:.4f}\n"
                    extra_info += f"Z_sublead_lepton_pt: {event2_extra_data['Z_sublead_lepton_pt'].iloc[0]:.4f}\t"
                    extra_info += f"Z_sublead_lepton_eta: {event2_extra_data['Z_sublead_lepton_eta'].iloc[0]:.4f}\t"
                    extra_info += f"Z_sublead_lepton_phi: {event2_extra_data['Z_sublead_lepton_phi'].iloc[0]:.4f}\n"
                    log_file.write(extra_info)

                    # Differences
                    diff_row_str = "Diff:\t"
                    has_diff = False
                    for var1, var2 in var_map.items():
                        val1 = event1_data[var1].iloc[0]
                        val2 = event2_data[var2].iloc[0]
                        diff = val1 - val2
                        diff_row_str += f"{var1}-{var2}: {diff:.4f}\t"
                        if abs(diff) > 1e-6:
                            has_diff = True
                    
                    if has_diff:
                        log_file.write(diff_row_str + "\n")
                    else:
                        log_file.write("Diff:\tNo significant differences.\n")

                    log_file.write("\n")
                else:
                    log_file.write(f"--- Event: {event_id} (in file2:two_jet, NOT found in file1:zero_to_one_jet) ---\n\n")


    print(f"Comparison finished. Log saved to {log_file_path}")

if __name__ == "__main__":
    main()