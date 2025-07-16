import uproot
import pandas as pd
import os

def filter_root_file(input_path, output_path):
    """
    Filters events in all trees of a ROOT file based on 'weight_corr' > 0.1
    and saves the result to a new file.

    Args:
        input_path (str): Path to the input ROOT file.
        output_path (str): Path to save the filtered ROOT file.
    """
    try:
        with uproot.open(input_path) as original_file:
            with uproot.recreate(output_path) as new_file:
                print(f"Opened input file: {input_path}")
                print(f"Will write to output file: {output_path}")

                tree_names = [key.split(';')[0] for key in original_file.keys() if isinstance(original_file[key], uproot.TTree)]
                
                if not tree_names:
                    print("No trees found in the input file.")
                    return

                print(f"Found trees: {tree_names}")

                for tree_name in tree_names:
                    print(f"Processing tree: '{tree_name}'")
                    tree = original_file[tree_name]
                    
                    # Check if 'weight_corr' branch exists
                    if 'weight_corr' not in tree.keys():
                        print(f"  Branch 'weight_corr' not found in tree '{tree_name}'. Copying tree as is.")
                        new_file[tree_name] = tree.arrays()
                        continue

                    # Read all branches into a pandas DataFrame
                    df = tree.arrays(library="pd")
                    
                    initial_events = len(df)
                    print(f"  Initial number of events: {initial_events}")

                    # Apply the filter
                    filtered_df = df[df['weight_corr'] <= 0.1]
                    
                    final_events = len(filtered_df)
                    print(f"  Number of events after filtering: {final_events}")
                    print(f"  Events removed: {initial_events - final_events}")

                    # Write the filtered DataFrame to the new ROOT file
                    new_file[tree_name] = filtered_df
                
                print(f"\nSuccessfully processed all trees. Filtered file saved to: {output_path}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Define input and output paths
    input_file = "/eos/home-j/jiehan/root/skimmed_ntuples_rui_new/ggH_M125/2022postEE.root"
    
    # Create an output directory if it doesn't exist
    output_dir = "filtered_root_files"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_file = os.path.join(output_dir, "2022postEE.root")

    # Run the filtering function
    filter_root_file(input_file, output_file)
