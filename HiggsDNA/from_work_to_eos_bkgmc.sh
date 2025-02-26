#!/bin/bash

# Define base paths
SRC_DIR="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/HiggsDNA/Parquet/run2/NanoV9/Bkg_MC"
DEST_DIR="/eos/home-p/pelai/HZgamma/Parquet/NanoV9/run2/Bkg_MC"

# Define the base dataset names
# DATASETS=("EWKZ2J" "ZG2JToG2L2J" "ZGToLLG" "DYJetsToLL")
# DATASETS=("TT" "TTGJets" "TGJets" "WW" "ZZ" "WZ")
DATASETS=("ttWJets" "ttZJets")


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
            cp -r "$SRC_PATH"/* "$DEST_PATH"
        else
            echo "Warning: Source directory $SRC_PATH does not exist, skipping copy."
        fi
    done
done

echo "All directories processed!"
