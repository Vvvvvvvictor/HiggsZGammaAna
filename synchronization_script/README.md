# BDT Model Comparison Scripts

This directory contains scripts to compare two different BDT models:
1. Your models (from hzgml/models)
2. External models (from /eos/project/h/htozg-dy-privatemc/rzou/bdt/XGB_scores/)

## Files

### Main Scripts
- `compare_two_bdt_models.py`: Main comparison script
- `run_bdt_comparison.sh`: Batch processing script
- `external_model_config.json`: Configuration for external model variables

### Usage

#### Single Sample Comparison
```bash
# Activate environment first
hzgml

# Run comparison for a single sample
python3 compare_two_bdt_models.py \
    --region two_jet \
    --sample VBF_M125 \
    --year 2018 \
    --output-dir ./output/
```

#### Batch Processing
```bash
# Run comparison for multiple samples and regions
./run_bdt_comparison.sh
```

### Parameters

#### compare_two_bdt_models.py
- `--region`: Region to process (`two_jet`, `zero_to_one_jet`, `all_jet`)
- `--sample`: Sample(s) to process (can specify multiple)
- `--year`: Year to process (default: 2018)
- `--input-dir`: Directory containing input ROOT files
- `--output-dir`: Output directory for results
- `--your-model-dir`: Directory containing your trained models
- `--external-model-dir`: Directory containing external models
- `--config-dir`: Directory containing configuration files

### Output Files

For each sample and region, the script generates:
1. **2D histogram plot**: `bdt_comparison_2d_{sample}_{region}.png`
2. **Scatter plot**: `bdt_comparison_scatter_{sample}_{region}.png`
3. **ROOT file**: `bdt_comparison_{region}.root` containing:
   - `your_bdt_score`: Scores from your model
   - `external_bdt_score`: Scores from external model
   - `weight`: Event weights
   - `H_mass`: Higgs mass (if available)

### Model Details

#### Your Models
- Format: XGBoost models in H5 or JSON format
- Location: `/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/models`
- Transformers: WeightedQuantileTransformer for input features and score transformation
- Files: `BDT_{region}_{fold}.h5` and `BDT_tsf_{region}_{fold}.pkl`

#### External Models
- Format: XGBoost models in JSON format
- Location: `/eos/project/h/htozg-dy-privatemc/rzou/bdt/XGB_scores/`
- Files: `model_xgb_{fold}_{model_type}.json`
- Model types: `vbf` for two_jet region, `ggf` for other regions

### Variable Mapping

The script automatically maps between different variable naming conventions:
- Your model variables are defined in training configuration files
- External model variables are mapped through `external_model_config.json`
- Branch mapping handles ROOT file variable names

### Example Results

The comparison plots show:
- **Correlation coefficient** between the two BDT scores
- **2D distribution** of score pairs
- **Scatter plot** with reference diagonal line (y=x)

High correlation indicates good agreement between models.

### Troubleshooting

1. **Environment issues**: Make sure to run `hzgml` to activate the correct environment
2. **Missing files**: Check that input ROOT files and model files exist
3. **Variable mismatches**: Review variable mappings in configuration files
4. **Memory issues**: Process fewer samples at once for large datasets

### Configuration Files Required

1. `/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/training_config_BDT.json`
2. `/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/data/apply_config_BDT.json`
3. `./external_model_config.json`

Make sure these files are properly configured with the correct variable lists and mappings.
