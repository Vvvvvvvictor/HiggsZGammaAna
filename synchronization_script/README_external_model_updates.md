# External Model Application - Updated Scripts

## Overview
Updated `apply_external_model.py` and `run_external_model_batch.sh` to process all skimmed ntuples produced by `generate_condor_sub.py` and organize outputs according to `apply_bdt_*.py` conventions.

## Key Changes

### 1. Input Processing
- **Source**: Now processes files directly from `/eos/home-j/jiehan/root/skimmed_ntuples/` (output of `generate_condor_sub.py`)
- **Structure**: Handles the complete sample structure including:
  - Signal samples: `ggH_M125`, `VBF_M125`, `WplusH_M125`, `WminusH_M125`, `ZH_M125`, `ttH_M125`
  - Nominal-only samples: `*_M120`, `*_M130`, `*_up`, `*_down`, `*_mu`
  - Systematic variations: `{sample}_{syst}_{up/down}` for M125 samples
  - Background/Data: `Data`

### 2. Special Naming Rules
- **`*_up`/`*_down` samples** → **`*_Tune_up`/`*_Tune_down`** (systematic: Tune)
- **`*_mu` samples** → **`*_Mmu`** (renaming convention)

### 3. Output Organization
Follows the same conventions as `apply_bdt_*.py` scripts:

#### Signal Output (`/eos/home-j/jiehan/root/fitting_signal`)
```
{sample}_{year}/
  └── output_{sample}.root
      └── DiphotonTree/
          ├── {proc}_ele_13TeV_{region}[{syst}]
          └── {proc}_mu_13TeV_{region}[{syst}]
```

#### Background Output (`/eos/home-j/jiehan/root/fitting_bkg`)
```
{proc}_{year}/
  └── output_{proc}_{year}.root
      └── DiphotonTree/
          └── Data_13TeV_{region}
```

### 4. External Model Integration
- **Variable mapping**: Uses `external_model_config.json` for variable name translation
- **Cross-validation**: Applies 4-fold CV using event ID modulo
- **Model types**: 
  - `zero_to_one_jet` → `ggf` models
  - `two_jet` → `vbf` models
- **Output**: Adds `external_bdt_score` column to all processed events

### 5. Systematic Processing
- **M125 samples**: Process all systematic variations (`{syst}_{up/down}`)
- **Other samples**: Nominal only
- **Naming**: Follows `apply_bdt_sig_corr.py` convention for systematic naming

## Usage

### Batch Processing (Recommended)
```bash
# Process both signal and background for all regions
./run_external_model_batch.sh

# Process only signals for a specific region
./run_external_model_batch.sh --signal-only --region zero_to_one_jet

# Process only backgrounds
./run_external_model_batch.sh --background-only

# Dry run to see what would be executed
./run_external_model_batch.sh --dry-run
```

### Direct Python Script
```bash
# Process signal samples for zero_to_one_jet region
python apply_external_model.py \
  --process-type signal \
  --region zero_to_one_jet \
  --inputFolder /eos/home-j/jiehan/root/skimmed_ntuples/ \
  --outputFolder /eos/home-j/jiehan/root/fitting_signal

# Process background samples for two_jet region
python apply_external_model.py \
  --process-type background \
  --region two_jet \
  --inputFolder /eos/home-j/jiehan/root/skimmed_ntuples/ \
  --outputFolder /eos/home-j/jiehan/root/fitting_bkg
```

## File Structure

### Input Files (from `generate_condor_sub.py`)
```
/eos/home-j/jiehan/root/skimmed_ntuples/
├── ggH_M125/{year}.root
├── ggH_M120/{year}.root
├── ggH_up/{year}.root              # → ggH_Tune_up
├── ggH_down/{year}.root            # → ggH_Tune_down
├── ggH_mu/{year}.root              # → ggH_Mmu
├── ggH_M125_Photon_scale_up/{year}.root
├── ggH_M125_Photon_scale_down/{year}.root
├── ...
└── Data/{year}.root
```

### Output Files (matching `apply_bdt_*.py`)
```
/eos/home-j/jiehan/root/fitting_signal/
├── ggH_125_{year}/output_ggH_125.root
├── ggH_120_{year}/output_ggH_120.root
├── ggH_Tune_up_{year}/output_ggH_Tune_up.root
├── ggH_Mmu_{year}/output_ggH_Mmu.root
└── ...

/eos/home-j/jiehan/root/fitting_bkg/
├── Data_{year}/output_Data_{year}.root
└── ...
```

## Configuration
- **External models**: `/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/`
- **Config file**: `external_model_config.json`
- **Variable mapping**: Defined in config file for name translation between internal and external formats

## Dependencies
- `uproot`: ROOT file I/O
- `pandas`: Data manipulation
- `xgboost`: External model inference
- `numpy`: Numerical operations

## Error Handling
- Missing input files: Warning logged, processing continues
- Missing variables: Filled with 0.0, warning logged
- Model loading failures: Error with descriptive message
- Cross-validation: Handles missing event column gracefully