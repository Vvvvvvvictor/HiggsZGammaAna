#!/bin/bash

#!/bin/bash

# Script to run BDT model comparison for multiple samples and regions
# Usage: ./run_bdt_comparison.sh

# Activate environment
echo "Activating hzgml environment..."
# source /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/setup.sh || echo "Failed to activate environment, continuing..."

# Set output directory with timestamp
OUTPUT_DIR="/eos/user/j/jiehan/root/outputs/bdt_comparison_$(date +%Y%m%d_%H%M%S)"
mkdir -p $OUTPUT_DIR

echo "Output directory: $OUTPUT_DIR"
echo "Starting BDT model comparison..."

# Define samples and regions to process
SAMPLES=("ZGToLLG") #"VBF_M125" "ggH_M125" "ZGToLLG" "EWKZ2J" "DYJetsToLL"
REGIONS=("zero_to_one_jet" "two_jet") # 
YEARS=("2018") #"2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBpix"

# Create region subdirectories
for region in "${REGIONS[@]}"; do
    mkdir -p $OUTPUT_DIR/$region
done

# Process each combination
for region in "${REGIONS[@]}"; do
    for year in "${YEARS[@]}"; do
        echo "=== Processing region: $region, year: $year ==="
        
        cd /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/synchronization_script
        
        # Process all samples for this region and year
        python3 compare_two_bdt_models.py \
            --region $region \
            --sample ${SAMPLES[@]} \
            --year $year \
            --output-dir $OUTPUT_DIR/$region/ \
            --input-dir /eos/user/j/jiehan/root/skimmed_ntuples_rui_new/ \
            --your-model-dir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/models \
            --external-model-dir /eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/ \
            --config-dir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data
        
        if [ $? -eq 0 ]; then
            echo "✓ Completed region: $region, year: $year"
        else
            echo "✗ Failed region: $region, year: $year"
        fi
        
        echo ""
    done
done

echo "=== BDT Comparison Summary ==="
echo "Output directory: $OUTPUT_DIR"
echo "Generated files:"
find $OUTPUT_DIR -name "*.png" -o -name "*.root" | sort
echo ""
echo "BDT model comparison completed!"
