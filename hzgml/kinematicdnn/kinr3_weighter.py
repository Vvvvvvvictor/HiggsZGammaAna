#!/usr/bin/env python
"""
Python wrapper for C++ kinematic DNN reweighting model.
This module provides a clean interface for the kinematic DNN model
by compiling and loading the C++ implementation directly.

Author: Jie Han
Date: August 2025
"""

import os
import ctypes
import tempfile
import subprocess
import numpy as np
from typing import List


class KinR3Weighter:
    """
    Wrapper class for the C++ kinematic DNN model.
    Compiles and loads the C++ implementation for efficient evaluation.
    """
    
    def __init__(self, compile_timeout=30):
        """
        Initialize the C++ DNN model by compiling and loading it.
        
        Args:
            compile_timeout: Timeout for C++ compilation in seconds (default: 30)
        """
        self._lib = None
        self._dnn_ptr = None
        self._compile_timeout = compile_timeout
        
        self._compile_cpp_model()
    
    def _compile_cpp_model(self):
        """
        Compile the C++ DNN model into a shared library and load it.
        """
        print("Starting C++ compilation process...")
        
        # Get the path to the C++ header file
        current_dir = os.path.dirname(os.path.abspath(__file__))
        cpp_header_path = os.path.join(current_dir, "kinr3_weighter.hpp")
        
        print(f"Looking for C++ header at: {cpp_header_path}")
        
        if not os.path.exists(cpp_header_path):
            raise FileNotFoundError(f"C++ header file not found: {cpp_header_path}")
        
        print("C++ header file found, creating temporary files...")
        
        # Create a temporary directory for compilation
        temp_dir = tempfile.mkdtemp()
        cpp_source = os.path.join(temp_dir, "kinr3_weighter_wrapper.cpp")
        shared_lib = os.path.join(temp_dir, "libkinr3_weighter.so")
        
        print(f"Temporary directory: {temp_dir}")
        
        # Create a C++ source file with a C interface
        cpp_code = f'''
#include "{cpp_header_path}"
#include <vector>

extern "C" {{
    kinr3_weighter* create_dnn() {{
        return new kinr3_weighter();
    }}
    
    void delete_dnn(kinr3_weighter* dnn) {{
        delete dnn;
    }}
    
    float evaluate_dnn(kinr3_weighter* dnn, float* input) {{
        std::vector<float> input_vec(input, input + 8);
        return dnn->evaluate(input_vec);
    }}
}}
'''
        
        print("Writing C++ wrapper source...")
        # Write the C++ source file
        with open(cpp_source, 'w') as f:
            f.write(cpp_code)
        
        print("Starting compilation...")
        # Compile the shared library
        compile_cmd = [
            'g++', '-shared', '-fPIC', '-O3', '-std=c++11', 
            '-o', shared_lib, cpp_source, '-lm'
        ]
        
        print(f"Compile command: {' '.join(compile_cmd)}")
        
        result = subprocess.run(compile_cmd, capture_output=True, text=True, timeout=self._compile_timeout)
        if result.returncode != 0:
            raise RuntimeError(f"Compilation failed: {result.stderr}")
        
        print("Compilation successful, loading shared library...")
        
        # Load the shared library
        self._lib = ctypes.CDLL(shared_lib)
        
        print("Defining function signatures...")
        # Define function signatures
        self._lib.create_dnn.restype = ctypes.c_void_p
        self._lib.delete_dnn.argtypes = [ctypes.c_void_p]
        self._lib.evaluate_dnn.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_float)]
        self._lib.evaluate_dnn.restype = ctypes.c_float
        
        print("Creating DNN instance...")
        # Create DNN instance
        self._dnn_ptr = self._lib.create_dnn()
        
        print("Successfully compiled and loaded C++ kinematic DNN model")
    
    def evaluate(self, input_vec: List[float]) -> float:
        """
        Evaluate the DNN model using C++ implementation.
        
        Args:
            input_vec: Input vector with 8 elements:
                [gamma_mvaID, photon_mht_dphi, n_jets, Z_pt, gamma_pt, H_pt, MHT_pt, HT]
                
        Returns:
            DNN output (kinematic weight)
        """
        if len(input_vec) != 8:
            raise ValueError(f"Input vector must have 8 elements, got {len(input_vec)}")
            
        # Check for invalid inputs
        for i, val in enumerate(input_vec):
            if np.isnan(val) or np.isinf(val):
                raise ValueError(f"Invalid input at position {i}: {val}")
        
        # Convert to ctypes array
        input_ctypes = (ctypes.c_float * 8)(*input_vec)
        # Evaluate DNN
        result = self._lib.evaluate_dnn(self._dnn_ptr, input_ctypes)
        return float(result)/(1-float(result))  # Convert to weight format
    
    def __del__(self):
        """Cleanup when object is destroyed."""
        if hasattr(self, '_dnn_ptr') and hasattr(self, '_lib') and self._dnn_ptr is not None:
            try:
                self._lib.delete_dnn(self._dnn_ptr)
            except:
                pass


# Global DNN instance (similar to photon_id_shape_systematics.py)
_weighter_instance = None

def get_weighter_instance():
    """Get or create the global kinematic weighter instance."""
    global _weighter_instance
    if _weighter_instance is None:
        _weighter_instance = KinR3Weighter()
    return _weighter_instance


def create_kin_weighter() -> KinR3Weighter:
    """
    Factory function to create a kinematic weighter
    
    Returns:
        Initialized KinR3Weighter instance
    """
    return KinR3Weighter()


# For backward compatibility
def get_default_weighter() -> KinR3Weighter:
    """Get a default weighter instance"""
    return get_weighter_instance()


if __name__ == "__main__":
    # Test the model
    print("Testing kinematic DNN weighter...")
    
    # Create weighter
    weighter = KinR3Weighter()
    print("Weighter created successfully!")
    
    # Test with some example inputs
    test_inputs = [
        [0.8, 1.5, 2.0, 45.0, 35.0, 80.0, 25.0, 200.0],
        [0.9, 0.5, 1.0, 50.0, 40.0, 90.0, 20.0, 150.0],
        [0.7, -1.1415926535, 0.0, 55.0, 30.0, 85.0, 30.0, 100.0]
    ]
    
    for i, inputs in enumerate(test_inputs):
        try:
            weight = weighter.evaluate(inputs)
            print(f"Test {i+1}: inputs = {inputs}")
            print(f"         output weight = {weight:.6f}")
        except Exception as e:
            print(f"Test {i+1} failed: {e}")
    
    print("Testing completed!")
