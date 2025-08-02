# Kinematic DNN Weight Implementation

## Overview
This modification adds kinematic DNN reweighting functionality to the `skim_ntuples.py` script. The DNN model calculates event weights based on kinematic variables to improve analysis accuracy.

## Added Features

### 1. DNN Model Implementation (`KinR3Weighter` class)
- Complete Python implementation of the C++ DNN model from `kinr3_weighter.hpp`
- 3-layer neural network with ELU and sigmoid activation functions
- Proper input scaling using pre-trained mean and standard deviation values

### 2. New Variables
- `photon_mht_dphi`: Delta phi between photon and MHT
- `kin_weight`: DNN-calculated kinematic weight

### 3. Input Variables (in order)
1. `gamma_mvaID`: Photon MVA ID score
2. `photon_mht_dphi`: Delta phi between photon and MHT (calculated)
3. `n_jets`: Number of jets
4. `Z_pt`: Z boson transverse momentum
5. `gamma_pt`: Photon transverse momentum  
6. `H_pt`: Higgs boson transverse momentum
7. `MHT_pt`: Missing HT magnitude
8. `HT`: Hadronic transverse energy

## Usage

The DNN weight calculation is automatically applied in the `decorate()` function:

```python
data['photon_mht_dphi'] = data.apply(lambda x: compute_photon_mht_dphi(x), axis=1)
data['kin_weight'] = data.apply(lambda x: compute_kin_weight(x), axis=1)
```

## Error Handling
- Returns default weight of 1.0 if any input is NaN or infinite
- Returns default weight of 1.0 if computation fails
- Graceful handling of edge cases

## Files Modified
1. `hzgml/scripts/skim_ntuples.py`: Main implementation
2. Variable list updated to include `MHT_pt`, `MHT_phi`, and `HT`

## Testing
The implementation has been tested with:
- Normal kinematic ranges
- Edge cases (NaN values, large delta phi)
- Multiple input combinations

## Output
The `kin_weight` variable will be available in all output trees:
- inclusive
- zero_jet
- one_jet
- two_jet
- zero_to_one_jet
- VH_ttH
- VH
- ZH
- ttH_had
- ttH_lep
