#!/usr/bin/env python
"""
Test script for the kinematic DNN weight functionality
"""

import sys
import os
import pandas as pd
import numpy as np

# Add the hzgml scripts directory to path
sys.path.append('/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts')

# Import our modified skim_ntuples functions
from skim_ntuples import compute_photon_mht_dphi, compute_kin_weight, kin_weighter, true_delta_phi

def test_dnn_model():
    """Test the DNN model with some example inputs"""
    print("Testing DNN model...")
    
    # Test inputs (gamma_mvaID, photon_mht_dphi, n_jets, Z_pt, gamma_pt, H_pt, MHT_pt, HT)
    test_inputs = [
        [0.8, 1.5, 2.0, 45.0, 35.0, 80.0, 25.0, 200.0],
        [0.9, 0.5, 1.0, 50.0, 40.0, 90.0, 20.0, 150.0],
        [0.7, 2.0, 0.0, 55.0, 30.0, 85.0, 30.0, 100.0]
    ]
    
    for i, inputs in enumerate(test_inputs):
        weight = kin_weighter.evaluate(inputs)
        print(f"Test {i+1}: inputs = {inputs}")
        print(f"         output weight = {weight:.6f}")
        print()

def test_helper_functions():
    """Test the helper functions"""
    print("Testing helper functions...")
    
    # Create a mock data row
    class MockRow:
        def __init__(self):
            self.gamma_phi = 1.5
            self.MHT_phi = 2.0
            self.gamma_mvaID = 0.8
            self.n_jets = 2
            self.Z_pt = 45.0
            self.gamma_pt = 35.0
            self.H_pt = 80.0
            self.MHT_pt = 25.0
            self.HT = 200.0
    
    row = MockRow()
    
    # Test photon_mht_dphi calculation
    dphi = compute_photon_mht_dphi(row)
    print(f"photon_mht_dphi = {dphi:.6f}")
    
    # Test kin_weight calculation
    weight = compute_kin_weight(row)
    print(f"kin_weight = {weight:.6f}")
    print()

def test_edge_cases():
    """Test edge cases"""
    print("Testing edge cases...")
    
    class MockRow:
        def __init__(self, **kwargs):
            self.gamma_phi = kwargs.get('gamma_phi', 1.5)
            self.MHT_phi = kwargs.get('MHT_phi', 2.0)
            self.gamma_mvaID = kwargs.get('gamma_mvaID', 0.8)
            self.n_jets = kwargs.get('n_jets', 2)
            self.Z_pt = kwargs.get('Z_pt', 45.0)
            self.gamma_pt = kwargs.get('gamma_pt', 35.0)
            self.H_pt = kwargs.get('H_pt', 80.0)
            self.MHT_pt = kwargs.get('MHT_pt', 25.0)
            self.HT = kwargs.get('HT', 200.0)
    
    # Test with NaN values
    row_nan = MockRow(gamma_mvaID=float('nan'))
    weight_nan = compute_kin_weight(row_nan)
    print(f"Weight with NaN input: {weight_nan}")
    
    # Test with very large delta phi
    row_large_dphi = MockRow(gamma_phi=0.0, MHT_phi=3.5)
    weight_large_dphi = compute_kin_weight(row_large_dphi)
    print(f"Weight with large dphi: {weight_large_dphi:.6f}")
    print()

if __name__ == '__main__':
    print("Running kinematic DNN weight tests...")
    print("=" * 50)
    
    test_dnn_model()
    test_helper_functions()
    test_edge_cases()
    
    print("All tests completed!")
