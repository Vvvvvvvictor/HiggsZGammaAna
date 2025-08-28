import numpy as np
import awkward as ak
import ctypes
import os
import tempfile
import subprocess
from pathlib import Path

from higgs_dna.utils.logger_utils import simple_logger
from higgs_dna.utils import awkward_utils, misc_utils

logger = simple_logger(__name__)


class PhotonIDShapeDNN:
    """
    Wrapper class for the C++ DNN model used in photon ID shape systematic calculations.
    """
    
    def __init__(self):
        self._compiled_lib = None
        self._compile_cpp_model()
    
    def _compile_cpp_model(self):
        """
        Compile the C++ DNN model into a shared library and load it.
        """
        try:
            # Get the path to the C++ header file
            cpp_header_path = misc_utils.expand_path("higgs_dna/systematics/data/rw_mmp_r3.hpp")
            
            # Create a temporary directory for compilation
            temp_dir = tempfile.mkdtemp()
            cpp_source = os.path.join(temp_dir, "rw_mmp_r3.cpp")
            shared_lib = os.path.join(temp_dir, "librw_mmp_r3.so")
            
            # Create a C++ source file with a C interface
            cpp_code = f'''
#include "{cpp_header_path}"
#include <vector>

extern "C" {{
    rw_mmp_r3* create_dnn() {{
        return new rw_mmp_r3();
    }}
    
    void delete_dnn(rw_mmp_r3* dnn) {{
        delete dnn;
    }}
    
    float evaluate_dnn(rw_mmp_r3* dnn, float* input) {{
        std::vector<float> input_vec(input, input + 4);
        return dnn->evaluate(input_vec);
    }}
}}
'''
            
            # Write the C++ source file
            with open(cpp_source, 'w') as f:
                f.write(cpp_code)
            
            # Compile the shared library
            compile_cmd = [
                'g++', '-shared', '-fPIC', '-O3', '-o', shared_lib, cpp_source, '-lm'
            ]
            
            result = subprocess.run(compile_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"Compilation failed: {result.stderr}")
            
            # Load the shared library
            self._lib = ctypes.CDLL(shared_lib)
            
            # Define function signatures
            self._lib.create_dnn.restype = ctypes.c_void_p
            self._lib.delete_dnn.argtypes = [ctypes.c_void_p]
            self._lib.evaluate_dnn.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_float)]
            self._lib.evaluate_dnn.restype = ctypes.c_float
            
            # Create DNN instance
            self._dnn_ptr = self._lib.create_dnn()
            
            logger.info("Successfully compiled and loaded C++ DNN model")
            
        except Exception as e:
            logger.error(f"Failed to compile C++ DNN model: {e}")
            # Fallback to pure Python implementation
            self._use_fallback = True
            
    def evaluate(self, input_array):
        """
        Evaluate the DNN model on input data.
        
        Args:
            input_array: numpy array with shape (N, 4) where N is number of photons
                        and the 4 features are [pt, abs_eta, idmva, relpterr]
        
        Returns:
            numpy array with DNN outputs
        """
        if hasattr(self, '_use_fallback') and self._use_fallback:
            return self._python_fallback(input_array)
            
        try:
            n_photons = input_array.shape[0]
            results = np.zeros(n_photons, dtype=np.float32)
            
            for i in range(n_photons):
                # Convert to ctypes array
                input_ctypes = (ctypes.c_float * 4)(*input_array[i])
                # Evaluate DNN
                results[i] = self._lib.evaluate_dnn(self._dnn_ptr, input_ctypes)
            
            return results
            
        except Exception as e:
            logger.warning(f"C++ evaluation failed, using Python fallback: {e}")
            return self._python_fallback(input_array)
    
    def _python_fallback(self, input_array):
        """
        Pure Python fallback implementation (simplified approximation).
        This should ideally implement the same network architecture.
        """
        logger.warning("Using simplified Python fallback for DNN evaluation")
        
        # This is a simplified fallback - in practice you might want to 
        # implement the full network or provide a different solution
        n_photons = input_array.shape[0]
        
        # Simple linear combination as placeholder
        # You might want to implement the actual network here
        weights = np.array([0.1, 0.2, 0.3, 0.4])
        bias = 0.5
        
        linear_output = np.dot(input_array, weights) + bias
        # Apply sigmoid activation
        sigmoid_output = 1.0 / (1.0 + np.exp(-linear_output))
        
        return sigmoid_output
    
    def __del__(self):
        """Cleanup when object is destroyed."""
        if hasattr(self, '_dnn_ptr') and hasattr(self, '_lib'):
            try:
                self._lib.delete_dnn(self._dnn_ptr)
            except:
                pass


# Global DNN instance
_dnn_instance = None

def get_dnn_instance():
    """Get or create the global DNN instance."""
    global _dnn_instance
    if _dnn_instance is None:
        _dnn_instance = PhotonIDShapeDNN()
    return _dnn_instance


def photon_id_shape_sf(events, year, central_only, input_collection="Photon"):
    """
    Calculate photon ID shape scale factors using DNN model.
    
    Args:
        events: awkward array of events
        year: data-taking year (string)
        central_only: if True, only return central SF (no variations)
        input_collection: name of photon collection (default: "Photon")
    
    Returns:
        dict with variations: {"central": sf_values, "up": sf_up, "down": sf_down}
    """
    
    required_fields = [
        (input_collection, "pt"),
        (input_collection, "eta"), 
        (input_collection, "mvaID"),
        (input_collection, "ptErrRel")  # relative pt error
    ]
    
    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.error(f"Missing required fields for photon ID shape SF: {missing_fields}")
        # Return dummy values
        photons = events[input_collection]
        n_photons = ak.num(photons)
        dummy_sf = ak.ones_like(photons.pt)
        variations = {"central": dummy_sf}
        if not central_only:
            variations["up"] = dummy_sf
            variations["down"] = dummy_sf
        return variations
    
    photons = events[input_collection]
    
    # Flatten photons for processing
    n_photons = ak.num(photons)
    photons_flattened = ak.flatten(photons)

    if ak.sum(n_photons) == 0:
        if not central_only:
            return {"central": ak.zeros_like(photons.pt), "up": ak.zeros_like(photons.pt), "down": ak.zeros_like(photons.pt)}
        else:
            return {"central": ak.zeros_like(photons.pt)}

    # Prepare input features for DNN
    pt = ak.to_numpy(photons_flattened.pt)
    abs_eta = np.abs(ak.to_numpy(photons_flattened.eta))
    idmva = ak.to_numpy(photons_flattened.mvaID)
    relpterr = ak.to_numpy(photons_flattened.ptErrRel)
    
    # Apply reasonable clipping to avoid extreme values
    pt = np.clip(pt, 15.0, 1000.0)
    abs_eta = np.clip(abs_eta, 0.0, 2.5)
    idmva = np.clip(idmva, -1.0, 1.0)
    relpterr = np.clip(relpterr, 0.0, 1.0)
    
    # Stack input features
    dnn_input = np.column_stack([pt, abs_eta, idmva, relpterr])
    
    # Get DNN instance and evaluate
    dnn = get_dnn_instance()
    dnn_output = dnn.evaluate(dnn_input)
    
    # Convert DNN output to scale factor
    # sf = dnn_output / (1.0 - dnn_output)
    # Avoid division by zero
    epsilon = 1e-7
    dnn_output = np.clip(dnn_output, epsilon, 1.0 - epsilon)
    logger.debug(f"DNN output shape: {dnn_output.shape}")
    print(f"Input features (first 10): {dnn_input[:10]}")
    print(f"DNN output: {dnn_output[:10]}")
    sf_values = dnn_output / (1.0 - dnn_output)
    
    # Apply reasonable SF limits
    sf_values = np.clip(sf_values, 0.5, 2.0)
    
    # Unflatten back to original structure
    sf_central = ak.unflatten(sf_values, n_photons)
    
    variations = {"central": sf_central}
    
    if not central_only:
        # Create systematic variations
        # For simplicity, use ±5% variation
        variation_size = 0.05
        
        sf_up = sf_values * (1.0 + variation_size)
        sf_down = sf_values * (1.0 - variation_size)
        
        # Apply limits to variations as well
        sf_up = np.clip(sf_up, 0.5, 2.0)
        sf_down = np.clip(sf_down, 0.5, 2.0)
        
        variations["up"] = ak.unflatten(sf_up, n_photons)
        variations["down"] = ak.unflatten(sf_down, n_photons)
    
    logger.debug(f"Calculated photon ID shape SF for {len(sf_values)} photons")
    logger.debug(f"SF range: [{np.min(sf_values):.3f}, {np.max(sf_values):.3f}]")
    
    return variations
