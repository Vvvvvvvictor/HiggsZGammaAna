#!/bin/bash

# Script to apply external models to all skimmed ntuples produced by generate_condor_sub.py
# This script handles signal samples (including nominal-only and systematic variations) and background/data samples
# Following the naming conventions from apply_bdt_*.py scripts

# Configuration
SCRIPT_DIR="/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/synchronization_script"
EXTERNAL_MODEL_SCRIPT="${SCRIPT_DIR}/apply_external_model.py"
CONFIG_FILE="${SCRIPT_DIR}/external_model_config.json"
EXTERNAL_MODEL_FOLDER="/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/"

# Input folder - where generate_condor_sub.py outputs are stored
INPUT_FOLDER="/eos/home-j/jiehan/root/skimmed_ntuples/"

# Output folders - matching apply_bdt_*.py conventions
SIGNAL_OUTPUT_FOLDER="/eos/home-j/jiehan/root/fitting_signal"
BACKGROUND_OUTPUT_FOLDER="/eos/home-j/jiehan/root/fitting_bkg"

# Regions to process
REGIONS=("zero_to_one_jet" "two_jet")

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Check if required files exist
check_requirements() {
    log_info "Checking requirements..."
    
    if [[ ! -f "$EXTERNAL_MODEL_SCRIPT" ]]; then
        log_error "External model script not found: $EXTERNAL_MODEL_SCRIPT"
        exit 1
    fi
    
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "Config file not found: $CONFIG_FILE"
        exit 1
    fi
    
    if [[ ! -d "$EXTERNAL_MODEL_FOLDER" ]]; then
        log_error "External model folder not found: $EXTERNAL_MODEL_FOLDER"
        exit 1
    fi
    
    if [[ ! -d "$INPUT_FOLDER" ]]; then
        log_error "Input folder not found: $INPUT_FOLDER"
        exit 1
    fi
    
    log_success "All requirements satisfied"
}

# Function to run external model application
run_external_model() {
    local process_type=$1
    local output_folder=$2
    local region=$3
    
    log_info "Processing $process_type samples for region $region"
    log_info "Input: $INPUT_FOLDER"
    log_info "Output: $output_folder"
    
    # Create output directory
    mkdir -p "$output_folder"
    
    # Run the external model application
    python "$EXTERNAL_MODEL_SCRIPT" \
        --config "$CONFIG_FILE" \
        --inputFolder "$INPUT_FOLDER" \
        --externalModelFolder "$EXTERNAL_MODEL_FOLDER" \
        --outputFolder "$output_folder" \
        --process-type "$process_type" \
        --region "$region"
    
    local exit_code=$?
    
    if [[ $exit_code -eq 0 ]]; then
        log_success "Successfully processed $process_type samples for region $region"
    else
        log_error "Failed to process $process_type samples for region $region (exit code: $exit_code)"
        return $exit_code
    fi
}

# Function to run with time measurement
run_with_timing() {
    local start_time=$(date +%s)
    "$@"
    local exit_code=$?
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [[ $exit_code -eq 0 ]]; then
        log_success "Completed in ${duration}s"
    else
        log_error "Failed after ${duration}s"
    fi
    
    return $exit_code
}

# Function to check sample counts
check_sample_counts() {
    local folder=$1
    local label=$2
    
    if [[ -d "$folder" ]]; then
        local file_count=$(find "$folder" -name "*.root" | wc -l)
        local dir_count=$(find "$folder" -mindepth 1 -maxdepth 1 -type d | wc -l)
        log_info "$label: $dir_count sample directories, $file_count ROOT files"
        
        # Show some sample directories
        if [[ $dir_count -gt 0 ]]; then
            log_info "Sample directories in $label:"
            find "$folder" -mindepth 1 -maxdepth 1 -type d | head -5 | while read dir; do
                echo "  - $(basename "$dir")"
            done
            if [[ $dir_count -gt 5 ]]; then
                echo "  ... and $((dir_count - 5)) more"
            fi
        fi
    else
        log_warn "$label folder not found: $folder"
    fi
}

# Main execution
main() {
    log_info "Starting external model application batch processing"
    log_info "Script: $EXTERNAL_MODEL_SCRIPT"
    log_info "Config: $CONFIG_FILE"
    log_info "External models: $EXTERNAL_MODEL_FOLDER"
    log_info "Input: $INPUT_FOLDER"
    
    # Check requirements
    check_requirements
    
    # Check input sample counts
    log_info "=== Input Sample Overview ==="
    check_sample_counts "$INPUT_FOLDER" "Input samples"
    
    # Track overall success
    local overall_success=true
    local total_start_time=$(date +%s)
    
    # Process each region
    for region in "${REGIONS[@]}"; do
        log_info "=== Processing region: $region ==="
        
        # Skip signal processing if background-only mode
        if [[ "$PROCESS_BACKGROUND_ONLY" != true ]]; then
            log_info "--- Processing signal samples for $region ---"
            if ! run_with_timing run_external_model "signal" "$SIGNAL_OUTPUT_FOLDER" "$region"; then
                log_error "Signal processing failed for region $region"
                overall_success=false
            fi
        fi
        
        # Skip background processing if signal-only mode
        if [[ "$PROCESS_SIGNAL_ONLY" != true ]]; then
            log_info "--- Processing background samples for $region ---"
            if ! run_with_timing run_external_model "background" "$BACKGROUND_OUTPUT_FOLDER" "$region"; then
                log_error "Background processing failed for region $region"
                overall_success=false
            fi
        fi
        
        log_info "=== Completed region: $region ==="
        echo
    done
    
    # Summary
    local total_end_time=$(date +%s)
    local total_duration=$((total_end_time - total_start_time))
    
    echo "======================================"
    if [[ "$overall_success" == true ]]; then
        log_success "All external model applications completed successfully!"
        log_success "Total time: ${total_duration}s"
    else
        log_error "Some external model applications failed!"
        log_error "Total time: ${total_duration}s"
        exit 1
    fi
    echo "======================================"
    
    # Print output summary
    log_info "Output folders:"
    if [[ "$PROCESS_BACKGROUND_ONLY" != true ]]; then
        log_info "  Signal: $SIGNAL_OUTPUT_FOLDER"
        check_sample_counts "$SIGNAL_OUTPUT_FOLDER" "Signal output"
    fi
    
    if [[ "$PROCESS_SIGNAL_ONLY" != true ]]; then
        log_info "  Background: $BACKGROUND_OUTPUT_FOLDER"
        check_sample_counts "$BACKGROUND_OUTPUT_FOLDER" "Background output"
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --signal-only)
            PROCESS_SIGNAL_ONLY=true
            shift
            ;;
        --background-only)
            PROCESS_BACKGROUND_ONLY=true
            shift
            ;;
        --region)
            SINGLE_REGION="$2"
            shift 2
            ;;
        --input-folder)
            INPUT_FOLDER="$2"
            shift 2
            ;;
        --signal-output)
            SIGNAL_OUTPUT_FOLDER="$2"
            shift 2
            ;;
        --background-output)
            BACKGROUND_OUTPUT_FOLDER="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --signal-only                Process only signal samples"
            echo "  --background-only            Process only background samples"
            echo "  --region REGION              Process only specified region (zero_to_one_jet|two_jet)"
            echo "  --input-folder PATH          Override input folder path"
            echo "  --signal-output PATH         Override signal output folder path"
            echo "  --background-output PATH     Override background output folder path"
            echo "  --dry-run                    Show what would be done without executing"
            echo "  --help, -h                   Show this help message"
            echo ""
            echo "Default paths:"
            echo "  Input: $INPUT_FOLDER"
            echo "  Signal output: $SIGNAL_OUTPUT_FOLDER"
            echo "  Background output: $BACKGROUND_OUTPUT_FOLDER"
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate mutually exclusive options
if [[ "$PROCESS_SIGNAL_ONLY" == true && "$PROCESS_BACKGROUND_ONLY" == true ]]; then
    log_error "Cannot use both --signal-only and --background-only options"
    exit 1
fi

# Override regions if single region specified
if [[ -n "$SINGLE_REGION" ]]; then
    if [[ "$SINGLE_REGION" != "zero_to_one_jet" && "$SINGLE_REGION" != "two_jet" ]]; then
        log_error "Invalid region: $SINGLE_REGION. Must be 'zero_to_one_jet' or 'two_jet'"
        exit 1
    fi
    REGIONS=("$SINGLE_REGION")
fi

# Show configuration if dry run
if [[ "$DRY_RUN" == true ]]; then
    log_info "DRY RUN MODE - Commands that would be executed:"
    log_info "Input folder: $INPUT_FOLDER"
    
    for region in "${REGIONS[@]}"; do
        if [[ "$PROCESS_BACKGROUND_ONLY" != true ]]; then
            echo "python $EXTERNAL_MODEL_SCRIPT --process-type signal --region $region --inputFolder $INPUT_FOLDER --outputFolder $SIGNAL_OUTPUT_FOLDER ..."
        fi
        if [[ "$PROCESS_SIGNAL_ONLY" != true ]]; then
            echo "python $EXTERNAL_MODEL_SCRIPT --process-type background --region $region --inputFolder $INPUT_FOLDER --outputFolder $BACKGROUND_OUTPUT_FOLDER ..."
        fi
    done
    
    log_info "Expected outputs:"
    if [[ "$PROCESS_BACKGROUND_ONLY" != true ]]; then
        log_info "  Signal: $SIGNAL_OUTPUT_FOLDER"
    fi
    if [[ "$PROCESS_SIGNAL_ONLY" != true ]]; then
        log_info "  Background: $BACKGROUND_OUTPUT_FOLDER"
    fi
    
    exit 0
fi

# Run main function
main