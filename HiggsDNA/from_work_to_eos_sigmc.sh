#!/bin/bash

# Define base paths
SRC_DIR="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/run2/NanoV9/Sig_MC"
DEST_DIR="/eos/home-p/pelai/HZgamma/Parquet/NanoV9/run2/Sig_MC/WI_Systematic"

# Define the base dataset names
# DATASETS=("ZH_M125" "WplusH_M125" "WminusH_M125" "VBF_M125" "ggH_M125" "ttH_M125") # All Processes
DATASETS=("WplusH_M125" "WminusH_M125" "VBFH_M125")

# Define the year-based suffixes
ERAS=("2018" "2017" "2016postVFP" "2016preVFP")

# Loop through each dataset
for DATASET in "${DATASETS[@]}"; do
    for ERA in "${ERAS[@]}"; do
        # Construct full directory names
        DIR="${DATASET}_${ERA}"
        SRC_PATH="$SRC_DIR/$DIR"
        DEST_PATH="$DEST_DIR/$DIR"

        # Ensure the destination directory exists
        if [ ! -d "$DEST_PATH" ]; then
            echo "Creating directory: $DEST_PATH"
            mkdir -p "$DEST_PATH"
        fi

        # (Optional) Copy files from source to destination
        if [ -d "$SRC_PATH" ]; then
            echo "Copying contents from $SRC_PATH to $DEST_PATH"
            cp -r "$SRC_PATH"/* "$DEST_PATH/"
        else
            echo "Warning: Source directory $SRC_PATH does not exist, skipping copy."
        fi
    done
done

echo "All directories processed!"
