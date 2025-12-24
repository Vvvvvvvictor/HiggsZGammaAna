#!/usr/bin/env python

import numpy as np
import pandas as pd
import math
from tqdm import tqdm

# Z constraint refit parameters for electrons and muons
ELECTRON_LINESHAPE_PARAMS = {
    "meanCB": 90.8919,
    "sigmaCB": 4.00007,
    "alphaCB": 1.1981,
    "nCB": 3.25604,
    "meanGauss1": 96.4278,
    "meanGauss2": 91.1649,
    "meanGauss3": 91.1513,
    "sigmaGauss1": 6.17509,
    "sigmaGauss2": 0.856305,
    "sigmaGauss3": 1.77148,
    "f1": 0.86449,
    "f2": 0.514371,
    "f3": 0.648225
}

MUON_LINESHAPE_PARAMS = {
    "meanCB": 90.8919,
    "sigmaCB": 4.00007,
    "alphaCB": 1.1981,
    "nCB": 3.25604,
    "meanGauss1": 96.4278,
    "meanGauss2": 91.1649,
    "meanGauss3": 91.1513,
    "sigmaGauss1": 6.17509,
    "sigmaGauss2": 0.856305,
    "sigmaGauss3": 1.77148,
    "f1": 0.86449,
    "f2": 0.514371,
    "f3": 0.648225
}

def get_energy_resolution_em(corrected_energy, eta):
    """Calculate electromagnetic energy resolution."""
    if abs(eta) < 1.48:
        C = 0.35 / 100.
        S = 5.51 / 100.
        N = 98. / 1000.
    else:
        C = 0.
        S = 12.8 / 100.
        N = 440. / 1000.

    result = math.sqrt(C * C * corrected_energy * corrected_energy + 
                      S * S * corrected_energy + N * N)
    return result

def get_photon_pt_error(photon_pt, photon_eta):
    """Calculate photon pt error based on energy resolution.
    This matches the pterr function in KinZfitter.cpp"""
    # Convert pt to energy (assuming massless photon)
    theta = 2 * math.atan(math.exp(-photon_eta))
    photon_energy = photon_pt / math.sin(theta)
    
    # Calculate energy error
    energy_err = get_energy_resolution_em(photon_energy, photon_eta)
    
    return energy_err

def crystal_ball_pdf(x, mean, sigma, alpha, n):
    """Crystal Ball PDF implementation."""
    t = (x - mean) / sigma
    if alpha < 0:
        t = -t
        alpha = -alpha
    
    abs_alpha = abs(alpha)
    if t >= -abs_alpha:
        return math.exp(-0.5 * t * t)
    else:
        A = math.pow(n / abs_alpha, n) * math.exp(-0.5 * abs_alpha * abs_alpha)
        B = n / abs_alpha - abs_alpha
        return A * math.pow(B - t, -n)

def gaussian_pdf(x, mean, sigma):
    """Gaussian PDF implementation."""
    return math.exp(-0.5 * ((x - mean) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))

def z_lineshape_pdf(mz, params):
    """Z boson lineshape PDF using Crystal Ball + 3 Gaussians."""
    cb = crystal_ball_pdf(mz, params["meanCB"], params["sigmaCB"], 
                         params["alphaCB"], params["nCB"])
    
    gauss1 = gaussian_pdf(mz, params["meanGauss1"], params["sigmaGauss1"])
    gauss2 = gaussian_pdf(mz, params["meanGauss2"], params["sigmaGauss2"])
    gauss3 = gaussian_pdf(mz, params["meanGauss3"], params["sigmaGauss3"])
    
    # Combine components with fractions
    cb_plus_gauss = params["f1"] * cb + (1 - params["f1"]) * gauss1
    cb_plus_gauss_plus_gauss = params["f2"] * cb_plus_gauss + (1 - params["f2"]) * gauss2
    full_pdf = params["f3"] * cb_plus_gauss_plus_gauss + (1 - params["f3"]) * gauss3
    
    return full_pdf

def calculate_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2, 
                           gamma_pt=None, gamma_eta=None, gamma_phi=None):
    """Calculate invariant mass of two leptons (and optionally a photon)."""
    # Convert to 4-vectors
    theta1 = 2 * math.atan(math.exp(-eta1))
    theta2 = 2 * math.atan(math.exp(-eta2))
    
    e1 = math.sqrt(pt1**2 / math.sin(theta1)**2 + m1**2)
    e2 = math.sqrt(pt2**2 / math.sin(theta2)**2 + m2**2)
    
    px1 = pt1 * math.cos(phi1)
    py1 = pt1 * math.sin(phi1)
    pz1 = pt1 / math.tan(theta1)
    
    px2 = pt2 * math.cos(phi2)
    py2 = pt2 * math.sin(phi2)
    pz2 = pt2 / math.tan(theta2)
    
    if gamma_pt is not None and gamma_eta is not None and gamma_phi is not None:
        # Include FSR photon
        gamma_theta = 2 * math.atan(math.exp(-gamma_eta))
        gamma_e = gamma_pt / math.sin(gamma_theta)
        gamma_px = gamma_pt * math.cos(gamma_phi)
        gamma_py = gamma_pt * math.sin(gamma_phi)
        gamma_pz = gamma_pt / math.tan(gamma_theta)
        
        total_e = e1 + e2 + gamma_e
        total_px = px1 + px2 + gamma_px
        total_py = py1 + py2 + gamma_py
        total_pz = pz1 + pz2 + gamma_pz
    else:
        total_e = e1 + e2
        total_px = px1 + px2
        total_py = py1 + py2
        total_pz = pz1 + pz2
    
    mass_squared = total_e**2 - total_px**2 - total_py**2 - total_pz**2
    return math.sqrt(max(0, mass_squared))

def negative_log_likelihood(pt_vals, pt_reco, pt_errors, eta, phi, mass, 
                          gamma_pt_fsr_reco=None, gamma_eta_fsr=None, gamma_phi_fsr=None,
                          gamma_pt_fsr_err=None, num_fsr=0, lineshape_params=None):
    """Calculate negative log likelihood for Z constraint fit.
    This matches the likelihood in KinZfitter.cpp which includes:
    - Gaussian constraints on lepton pT measurements
    - Gaussian constraint on FSR photon pT measurement (if present, also fitted)
    - Z lineshape constraint
    
    Args:
        pt_vals: List of pT values to fit [pt1, pt2] or [pt1, pt2, pt_fsr] depending on num_fsr
        pt_reco: [pt1_reco, pt2_reco]
        pt_errors: [pt1_err, pt2_err]
        eta, phi, mass: Lepton kinematics
        gamma_pt_fsr_reco: Reconstructed FSR photon pT (or list if multiple)
        gamma_eta_fsr: FSR photon eta (or list if multiple)
        gamma_phi_fsr: FSR photon phi (or list if multiple)
        gamma_pt_fsr_err: FSR photon pT error (or list if multiple)
        num_fsr: Number of FSR photons (0, 1, or 2)
        lineshape_params: Z lineshape parameters
    """
    pt1 = pt_vals[0]
    pt2 = pt_vals[1]
    
    pt1_reco, pt2_reco = pt_reco
    pt1_err, pt2_err = pt_errors
    eta1, eta2 = eta
    phi1, phi2 = phi
    m1, m2 = mass
    
    # Gaussian constraints on pT measurements for leptons
    gauss1 = -0.5 * ((pt1 - pt1_reco) / pt1_err) ** 2
    gauss2 = -0.5 * ((pt2 - pt2_reco) / pt2_err) ** 2
    
    # Gaussian constraints for FSR photons (matching KinZfitter.cpp)
    gauss_fsr = 0.0
    gamma_pt_fsr_fit = []
    
    if num_fsr >= 1:
        # Process all FSR photons
        for i in range(num_fsr):
            pt_fsr = pt_vals[2 + i]
            
            # Handle single vs multiple FSR photons
            if num_fsr == 1:
                pt_fsr_reco = gamma_pt_fsr_reco
                pt_fsr_err = gamma_pt_fsr_err
                eta_fsr = gamma_eta_fsr
                phi_fsr = gamma_phi_fsr
            else:
                pt_fsr_reco = gamma_pt_fsr_reco[i]
                pt_fsr_err = gamma_pt_fsr_err[i]
                eta_fsr = gamma_eta_fsr[i]
                phi_fsr = gamma_phi_fsr[i]
            
            # Add Gaussian constraint for this FSR photon
            gauss_fsr += -0.5 * ((pt_fsr - pt_fsr_reco) / pt_fsr_err) ** 2
            
            # Store fitted FSR photon parameters
            gamma_pt_fsr_fit.append([pt_fsr, eta_fsr, phi_fsr])
    
    # Calculate Z mass with current pT values (including FSR if present)
    if num_fsr == 0:
        mz = calculate_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2)
    elif num_fsr == 1:
        pt_fsr, eta_fsr, phi_fsr = gamma_pt_fsr_fit[0]
        mz = calculate_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2,
                                    pt_fsr, eta_fsr, phi_fsr)
    else:  # num_fsr >= 2
        # For multiple FSR photons, calculate mass with all of them
        # First add both leptons
        theta1 = 2 * math.atan(math.exp(-eta1))
        theta2 = 2 * math.atan(math.exp(-eta2))
        
        e1 = math.sqrt(pt1**2 / math.sin(theta1)**2 + m1**2)
        e2 = math.sqrt(pt2**2 / math.sin(theta2)**2 + m2**2)
        
        px1 = pt1 * math.cos(phi1)
        py1 = pt1 * math.sin(phi1)
        pz1 = pt1 / math.tan(theta1)
        
        px2 = pt2 * math.cos(phi2)
        py2 = pt2 * math.sin(phi2)
        pz2 = pt2 / math.tan(theta2)
        
        total_e = e1 + e2
        total_px = px1 + px2
        total_py = py1 + py2
        total_pz = pz1 + pz2
        
        # Add all FSR photons
        for fsr_params in gamma_pt_fsr_fit:
            pt_fsr, eta_fsr, phi_fsr = fsr_params
            theta_fsr = 2 * math.atan(math.exp(-eta_fsr))
            e_fsr = pt_fsr / math.sin(theta_fsr)
            px_fsr = pt_fsr * math.cos(phi_fsr)
            py_fsr = pt_fsr * math.sin(phi_fsr)
            pz_fsr = pt_fsr / math.tan(theta_fsr)
            
            total_e += e_fsr
            total_px += px_fsr
            total_py += py_fsr
            total_pz += pz_fsr
        
        mass_squared = total_e**2 - total_px**2 - total_py**2 - total_pz**2
        mz = math.sqrt(max(0, mass_squared))
    
    # Z lineshape constraint
    z_constraint = math.log(max(1e-10, z_lineshape_pdf(mz, lineshape_params)))
    
    return -(gauss1 + gauss2 + gauss_fsr + z_constraint)

def minimize_simple(func, x0, bounds=None):
    """Simple grid search minimization for multi-dimensional optimization.
    Supports 2 or more parameters (pt1, pt2, [pt_fsr1], [pt_fsr2], ...)
    """
    best_val = float('inf')
    best_x = x0[:]
    n_params = len(x0)
    
    # Simple grid search around initial values
    # For more than 2 parameters, use coarser grid to avoid combinatorial explosion
    if n_params <= 2:
        ranges = [range(-5, 6) for _ in range(n_params)]
        step_size = 0.5
    elif n_params <= 4:
        ranges = [range(-3, 4) for _ in range(n_params)]
        step_size = 1.0
    else:
        # For many FSR photons, use even coarser grid
        ranges = [range(-2, 3) for _ in range(n_params)]
        step_size = 1.5
    
    # Generate all combinations
    import itertools
    for deltas in itertools.product(*ranges):
        x_test = []
        valid = True
        
        for i, delta in enumerate(deltas):
            if bounds and i < len(bounds):
                x_val = max(bounds[i][0], min(bounds[i][1], x0[i] + delta * step_size))
            else:
                x_val = x0[i] + delta * step_size
            x_test.append(x_val)
        
        try:
            val = func(x_test)
            if val < best_val:
                best_val = val
                best_x = x_test
        except:
            continue
    
    return best_x, True

def perform_z_refit(row):
    """Perform Z constraint refit for a single event.
    This matches the KinZfitter.cpp implementation including:
    - FSR photon pT fitting (not fixed)
    - Proper FSR photon assignment to leptons
    - Multiple FSR photons support
    """
    try:
        # Extract lepton information
        l1_pt = row['Z_lead_lepton_pt']
        l1_eta = row['Z_lead_lepton_eta']
        l1_phi = row['Z_lead_lepton_phi']
        l1_mass = row['Z_lead_lepton_mass']
        l1_id = row['Z_lead_lepton_id']
        
        l2_pt = row['Z_sublead_lepton_pt']
        l2_eta = row['Z_sublead_lepton_eta']
        l2_phi = row['Z_sublead_lepton_phi']
        l2_mass = row['Z_sublead_lepton_mass']
        
        # Calculate pT errors
        if abs(l1_id) == 11:  # electrons
            l1_energy = l1_pt / math.sin(2 * math.atan(math.exp(-l1_eta)))
            l2_energy = l2_pt / math.sin(2 * math.atan(math.exp(-l2_eta)))
            l1_energy_err = row.get('Z_lead_lepton_ptE_error', l1_energy * 0.01)
            l2_energy_err = row.get('Z_sublead_lepton_ptE_error', l2_energy * 0.01)
            
            l1_pt_err = abs(l1_pt * l1_energy_err / l1_energy)
            l2_pt_err = abs(l2_pt * l2_energy_err / l2_energy)
            lineshape_params = ELECTRON_LINESHAPE_PARAMS
        else:  # muons
            l1_pt_err = abs(row.get('Z_lead_lepton_ptE_error', l1_pt * 0.01))
            l2_pt_err = abs(row.get('Z_sublead_lepton_ptE_error', l2_pt * 0.01))
            lineshape_params = MUON_LINESHAPE_PARAMS
        
        # Collect FSR photons information
        # Check for multiple FSR photons (can be 2 or more)
        fsr_photons = []
        gamma_pt = row['gamma_pt']
        gamma_eta = row['gamma_eta'] 
        gamma_phi = row['gamma_phi']
        
        # Try to find all FSR photons
        # First check for columns named gamma_fsr_pt, gamma_fsr2_pt, gamma_fsr3_pt, etc.
        # Or gamma_fsr_pt with array/list structure
        
        # Method 1: Single FSR photon column
        gamma_pt_fsr = row.get('gamma_fsr_pt', -1)
        if gamma_pt_fsr > 0:
            gamma_eta_fsr = row.get('gamma_fsr_eta', 0)
            gamma_phi_fsr = row.get('gamma_fsr_phi', 0)
            
            # Check if FSR photon is clean (not overlapping with main photon or leptons)
            dr_gamma = math.sqrt((gamma_eta_fsr - gamma_eta)**2 + (gamma_phi_fsr - gamma_phi)**2)
            dr_l1 = math.sqrt((gamma_eta_fsr - l1_eta)**2 + (gamma_phi_fsr - l1_phi)**2)
            dr_l2 = math.sqrt((gamma_eta_fsr - l2_eta)**2 + (gamma_phi_fsr - l2_phi)**2)
            
            if dr_gamma > 0.001 and dr_l1 > 0.001 and dr_l2 > 0.001:
                # Calculate FSR photon pT error (matching KinZfitter.cpp pterr function)
                gamma_pt_fsr_err = get_photon_pt_error(gamma_pt_fsr, gamma_eta_fsr)
                
                # Determine which lepton this FSR photon belongs to
                # Based on minimum dR (closest lepton)
                ifsr = 0 if dr_l1 < dr_l2 else 1
                
                fsr_photons.append({
                    'pt': gamma_pt_fsr,
                    'eta': gamma_eta_fsr,
                    'phi': gamma_phi_fsr,
                    'pt_err': gamma_pt_fsr_err,
                    'ifsr': ifsr  # 0 = belongs to lepton1, 1 = belongs to lepton2
                })
        
        # Method 2: Check for additional FSR photons (gamma_fsr2_pt, gamma_fsr3_pt, etc.)
        fsr_idx = 2
        while True:
            fsr_key = f'gamma_fsr{fsr_idx}_pt' if fsr_idx > 1 else 'gamma_fsr_pt'
            gamma_pt_fsrN = row.get(fsr_key, -1)
            
            if gamma_pt_fsrN <= 0:
                # Also try without number suffix for arrays
                if fsr_idx == 2:
                    # Try checking if gamma_fsr_pt is an array/list
                    try:
                        fsr_pt_val = row.get('gamma_fsr_pt', -1)
                        if isinstance(fsr_pt_val, (list, np.ndarray)) and len(fsr_pt_val) > 1:
                            # It's an array, process all elements
                            fsr_eta_arr = row.get('gamma_fsr_eta', [])
                            fsr_phi_arr = row.get('gamma_fsr_phi', [])
                            
                            for i in range(1, len(fsr_pt_val)):  # Start from 1 since 0 was processed
                                if fsr_pt_val[i] > 0:
                                    dr_gamma = math.sqrt((fsr_eta_arr[i] - gamma_eta)**2 + (fsr_phi_arr[i] - gamma_phi)**2)
                                    dr_l1 = math.sqrt((fsr_eta_arr[i] - l1_eta)**2 + (fsr_phi_arr[i] - l1_phi)**2)
                                    dr_l2 = math.sqrt((fsr_eta_arr[i] - l2_eta)**2 + (fsr_phi_arr[i] - l2_phi)**2)
                                    
                                    if dr_gamma > 0.001 and dr_l1 > 0.001 and dr_l2 > 0.001:
                                        gamma_pt_fsr_err = get_photon_pt_error(fsr_pt_val[i], fsr_eta_arr[i])
                                        ifsr = 0 if dr_l1 < dr_l2 else 1
                                        
                                        fsr_photons.append({
                                            'pt': fsr_pt_val[i],
                                            'eta': fsr_eta_arr[i],
                                            'phi': fsr_phi_arr[i],
                                            'pt_err': gamma_pt_fsr_err,
                                            'ifsr': ifsr
                                        })
                    except:
                        pass
                break
            
            # Process FSR photon N
            gamma_eta_fsrN = row.get(f'gamma_fsr{fsr_idx}_eta', 0)
            gamma_phi_fsrN = row.get(f'gamma_fsr{fsr_idx}_phi', 0)
            
            dr_gamma = math.sqrt((gamma_eta_fsrN - gamma_eta)**2 + (gamma_phi_fsrN - gamma_phi)**2)
            dr_l1 = math.sqrt((gamma_eta_fsrN - l1_eta)**2 + (gamma_phi_fsrN - l1_phi)**2)
            dr_l2 = math.sqrt((gamma_eta_fsrN - l2_eta)**2 + (gamma_phi_fsrN - l2_phi)**2)
            
            if dr_gamma > 0.001 and dr_l1 > 0.001 and dr_l2 > 0.001:
                gamma_pt_fsr_err = get_photon_pt_error(gamma_pt_fsrN, gamma_eta_fsrN)
                ifsr = 0 if dr_l1 < dr_l2 else 1
                
                fsr_photons.append({
                    'pt': gamma_pt_fsrN,
                    'eta': gamma_eta_fsrN,
                    'phi': gamma_phi_fsrN,
                    'pt_err': gamma_pt_fsr_err,
                    'ifsr': ifsr
                })
            
            fsr_idx += 1
            # Safety limit to prevent infinite loop
            if fsr_idx > 10:
                break
        
        num_fsr = len(fsr_photons)
        
        # Define objective function for minimization
        def objective(pt_vals):
            gamma_pt_fsr_reco = None
            gamma_eta_fsr_list = None
            gamma_phi_fsr_list = None
            gamma_pt_fsr_err_list = None
            
            if num_fsr >= 1:
                gamma_pt_fsr_reco = [fsr['pt'] for fsr in fsr_photons]
                gamma_eta_fsr_list = [fsr['eta'] for fsr in fsr_photons]
                gamma_phi_fsr_list = [fsr['phi'] for fsr in fsr_photons]
                gamma_pt_fsr_err_list = [fsr['pt_err'] for fsr in fsr_photons]
                
                if num_fsr == 1:
                    gamma_pt_fsr_reco = gamma_pt_fsr_reco[0]
                    gamma_eta_fsr_list = gamma_eta_fsr_list[0]
                    gamma_phi_fsr_list = gamma_phi_fsr_list[0]
                    gamma_pt_fsr_err_list = gamma_pt_fsr_err_list[0]
            
            return negative_log_likelihood(
                pt_vals, [l1_pt, l2_pt], [l1_pt_err, l2_pt_err],
                [l1_eta, l2_eta], [l1_phi, l2_phi], [l1_mass, l2_mass],
                gamma_pt_fsr_reco,
                gamma_eta_fsr_list,
                gamma_phi_fsr_list,
                gamma_pt_fsr_err_list,
                num_fsr, lineshape_params
            )
        
        # Initial values and bounds
        x0 = [l1_pt, l2_pt]
        bounds = [
            (max(5.0, l1_pt - 2*l1_pt_err), l1_pt + 2*l1_pt_err),
            (max(5.0, l2_pt - 2*l2_pt_err), l2_pt + 2*l2_pt_err)
        ]
        
        # Add FSR photon pT to fitting parameters
        for fsr in fsr_photons:
            x0.append(fsr['pt'])
            fsr_pt_min = max(2.0, fsr['pt'] - 3*fsr['pt_err'])
            fsr_pt_max = fsr['pt'] + 3*fsr['pt_err'] if fsr['pt'] >= 2 else fsr_pt_min
            bounds.append((fsr_pt_min, fsr_pt_max))
        
        # Perform minimization
        try:
            result, success = minimize_simple(objective, x0, bounds)
            if success:
                # Extract fitted values
                l1_pt_refit = result[0]
                l2_pt_refit = result[1]
                
                # Extract fitted FSR photon pT values
                fsr_pt_refit = []
                fsr_ifsr = []  # Store which lepton each FSR belongs to
                for i in range(num_fsr):
                    fsr_pt_refit.append(result[2 + i])
                    # Get ifsr from fsr_photons
                    fsr_ifsr.append(fsr_photons[i]['ifsr'])
                
                return (l1_pt_refit, l2_pt_refit, l1_pt_err, l2_pt_err, 
                        fsr_pt_refit, fsr_ifsr, num_fsr)
        except:
            pass
        
        # If optimization fails, return original values
        return (l1_pt, l2_pt, l1_pt_err, l2_pt_err, 
                [fsr['pt'] for fsr in fsr_photons], 
                [fsr['ifsr'] for fsr in fsr_photons], num_fsr)
        
    except Exception as e:
        # If any error occurs, return original values
        return (row['Z_lead_lepton_pt'], row['Z_sublead_lepton_pt'], 
                0.0, 0.0, [], [], 0)

def apply_z_refit(df):
    """Apply Z constraint refit to DataFrame.
    This matches the KinZfitter.cpp behavior including:
    - FSR photon pT refitting
    - Proper assignment of FSR photons to leptons (ifsr==0 -> lepton1, ifsr==1 -> lepton2)
    - Calculation of refit Z and H kinematics
    """
    print("Applying Z constraint refit...")
    
    # Initialize refit columns
    df['Z_lead_lepton_pt_refit'] = df['Z_lead_lepton_pt'].copy()
    df['Z_sublead_lepton_pt_refit'] = df['Z_sublead_lepton_pt'].copy()
    df['Z_lead_lepton_pt_refit_err'] = 0.0
    df['Z_sublead_lepton_pt_refit_err'] = 0.0
    df['num_fsr_photons'] = 0
    
    # Store FSR photon info as lists (to support multiple FSR photons)
    df['gamma_fsr_pt_refit_list'] = [[] for _ in range(len(df))]
    df['gamma_fsr_ifsr_list'] = [[] for _ in range(len(df))]
    
    # Also keep first FSR photon info for backward compatibility
    df['gamma_fsr_pt_refit'] = 0.0
    df['gamma_fsr_ifsr'] = -1  # -1 means no FSR, 0 = belongs to lepton1, 1 = belongs to lepton2
    
    # Apply refit with progress bar
    refit_results = [
        perform_z_refit(row)
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Z refit progress")
    ]

    if refit_results:
        (pt1_refit_col, pt2_refit_col, pt1_err_col, pt2_err_col,
         fsr_pt_refit_col, fsr_ifsr_col, num_fsr_col) = zip(*refit_results)
    else:
        pt1_refit_col = pt2_refit_col = pt1_err_col = pt2_err_col = ()
        fsr_pt_refit_col = fsr_ifsr_col = num_fsr_col = ()

    # Update DataFrame with refit results in vectorized form
    if pt1_refit_col:
        df['Z_lead_lepton_pt_refit'] = pt1_refit_col
        df['Z_sublead_lepton_pt_refit'] = pt2_refit_col
        df['Z_lead_lepton_pt_refit_err'] = pt1_err_col
        df['Z_sublead_lepton_pt_refit_err'] = pt2_err_col
        df['num_fsr_photons'] = num_fsr_col
        df['gamma_fsr_pt_refit_list'] = list(fsr_pt_refit_col)
        df['gamma_fsr_ifsr_list'] = list(fsr_ifsr_col)
        df['gamma_fsr_pt_refit'] = [pts[0] if n >= 1 else 0.0 for pts, n in zip(fsr_pt_refit_col, num_fsr_col)]
        df['gamma_fsr_ifsr'] = [ifsr[0] if n >= 1 else -1 for ifsr, n in zip(fsr_ifsr_col, num_fsr_col)]
    
    # Calculate refit Z and H kinematics
    print("Calculating refit kinematics...")
    z_pt_vals = []
    z_eta_vals = []
    z_phi_vals = []
    z_mass_vals = []
    h_pt_vals = []
    h_eta_vals = []
    h_phi_vals = []
    h_mass_vals = []

    for row in tqdm(df.itertuples(index=True), total=len(df), desc="Kinematics"):
        try:
            l1_pt_refit = row.Z_lead_lepton_pt_refit
            l2_pt_refit = row.Z_sublead_lepton_pt_refit
            num_fsr = int(row.num_fsr_photons)

            l1_eta = row.Z_lead_lepton_eta
            l1_phi = row.Z_lead_lepton_phi
            l1_mass = row.Z_lead_lepton_mass

            l2_eta = row.Z_sublead_lepton_eta
            l2_phi = row.Z_sublead_lepton_phi
            l2_mass = row.Z_sublead_lepton_mass

            gamma_pt = row.gamma_pt
            gamma_eta = row.gamma_eta
            gamma_phi = row.gamma_phi

            theta1 = 2 * math.atan(math.exp(-l1_eta))
            theta2 = 2 * math.atan(math.exp(-l2_eta))

            l1_px = l1_pt_refit * math.cos(l1_phi)
            l1_py = l1_pt_refit * math.sin(l1_phi)
            l1_pz = l1_pt_refit / math.tan(theta1)
            l1_e = math.sqrt(l1_pt_refit**2 / math.sin(theta1)**2 + l1_mass**2)

            l2_px = l2_pt_refit * math.cos(l2_phi)
            l2_py = l2_pt_refit * math.sin(l2_phi)
            l2_pz = l2_pt_refit / math.tan(theta2)
            l2_e = math.sqrt(l2_pt_refit**2 / math.sin(theta2)**2 + l2_mass**2)

            if num_fsr >= 1:
                fsr_pt_refit_list = row.gamma_fsr_pt_refit_list
                fsr_ifsr_list = row.gamma_fsr_ifsr_list

                for fsr_idx in range(num_fsr):
                    if fsr_idx == 0:
                        gamma_pt_fsr_orig = getattr(row, 'gamma_fsr_pt', 0)
                        gamma_eta_fsr = getattr(row, 'gamma_fsr_eta', 0)
                        gamma_phi_fsr = getattr(row, 'gamma_fsr_phi', 0)
                    else:
                        fsr_key = f'gamma_fsr{fsr_idx+1}_pt'
                        gamma_pt_fsr_orig = getattr(row, fsr_key, 0)
                        gamma_eta_fsr = getattr(row, f'gamma_fsr{fsr_idx+1}_eta', 0)
                        gamma_phi_fsr = getattr(row, f'gamma_fsr{fsr_idx+1}_phi', 0)

                        if gamma_pt_fsr_orig == 0:
                            try:
                                fsr_pt_arr = getattr(row, 'gamma_fsr_pt', [])
                                if isinstance(fsr_pt_arr, (list, np.ndarray)) and len(fsr_pt_arr) > fsr_idx:
                                    gamma_pt_fsr_orig = fsr_pt_arr[fsr_idx]
                                    fsr_eta_arr = getattr(row, 'gamma_fsr_eta', [])
                                    fsr_phi_arr = getattr(row, 'gamma_fsr_phi', [])
                                    gamma_eta_fsr = fsr_eta_arr[fsr_idx]
                                    gamma_phi_fsr = fsr_phi_arr[fsr_idx]
                            except Exception:
                                continue

                    if fsr_idx < len(fsr_pt_refit_list) and fsr_idx < len(fsr_ifsr_list):
                        gamma_fsr_pt_refit = fsr_pt_refit_list[fsr_idx]
                        gamma_fsr_ifsr = int(fsr_ifsr_list[fsr_idx])
                    else:
                        continue

                    if gamma_fsr_pt_refit > 0 and gamma_fsr_ifsr >= 0:
                        gamma_theta_fsr = 2 * math.atan(math.exp(-gamma_eta_fsr))

                        gamma_px_fsr = gamma_fsr_pt_refit * math.cos(gamma_phi_fsr)
                        gamma_py_fsr = gamma_fsr_pt_refit * math.sin(gamma_phi_fsr)
                        gamma_pz_fsr = gamma_fsr_pt_refit / math.tan(gamma_theta_fsr)
                        gamma_e_fsr = gamma_fsr_pt_refit / math.sin(gamma_theta_fsr)

                        if gamma_fsr_ifsr == 0:
                            l1_px += gamma_px_fsr
                            l1_py += gamma_py_fsr
                            l1_pz += gamma_pz_fsr
                            l1_e += gamma_e_fsr
                        elif gamma_fsr_ifsr == 1:
                            l2_px += gamma_px_fsr
                            l2_py += gamma_py_fsr
                            l2_pz += gamma_pz_fsr
                            l2_e += gamma_e_fsr

            z_px = l1_px + l2_px
            z_py = l1_py + l2_py
            z_pz = l1_pz + l2_pz
            z_e = l1_e + l2_e

            z_pt_refit = math.sqrt(z_px**2 + z_py**2)
            z_phi_refit = math.atan2(z_py, z_px)
            z_p = math.sqrt(z_px**2 + z_py**2 + z_pz**2)
            z_eta_refit = 0.5 * math.log((z_p + z_pz) / (z_p - z_pz)) if abs(z_p - z_pz) > 1e-10 else 0
            z_mass_refit = math.sqrt(max(0, z_e**2 - z_px**2 - z_py**2 - z_pz**2))

            gamma_theta = 2 * math.atan(math.exp(-gamma_eta))
            gamma_px = gamma_pt * math.cos(gamma_phi)
            gamma_py = gamma_pt * math.sin(gamma_phi)
            gamma_pz = gamma_pt / math.tan(gamma_theta)
            gamma_e = gamma_pt / math.sin(gamma_theta)

            h_px = z_px + gamma_px
            h_py = z_py + gamma_py
            h_pz = z_pz + gamma_pz
            h_e = z_e + gamma_e

            h_pt_refit = math.sqrt(h_px**2 + h_py**2)
            h_phi_refit = math.atan2(h_py, h_px)
            h_p = math.sqrt(h_px**2 + h_py**2 + h_pz**2)
            h_eta_refit = 0.5 * math.log((h_p + h_pz) / (h_p - h_pz)) if abs(h_p - h_pz) > 1e-10 else 0
            h_mass_refit = math.sqrt(max(0, h_e**2 - h_px**2 - h_py**2 - h_pz**2))

            z_pt_vals.append(z_pt_refit)
            z_eta_vals.append(z_eta_refit)
            z_phi_vals.append(z_phi_refit)
            z_mass_vals.append(z_mass_refit)
            h_pt_vals.append(h_pt_refit)
            h_eta_vals.append(h_eta_refit)
            h_phi_vals.append(h_phi_refit)
            h_mass_vals.append(h_mass_refit)

        except Exception:
            z_pt_vals.append(getattr(row, 'Z_pt', 0))
            z_eta_vals.append(getattr(row, 'Z_eta', 0))
            z_phi_vals.append(getattr(row, 'Z_phi', 0))
            z_mass_vals.append(getattr(row, 'Z_mass', 0))
            h_pt_vals.append(getattr(row, 'H_pt', 0))
            h_eta_vals.append(getattr(row, 'H_eta', 0))
            h_phi_vals.append(getattr(row, 'H_phi', 0))
            h_mass_vals.append(getattr(row, 'H_mass', 0))

    df['Z_pt_refit'] = z_pt_vals
    df['Z_eta_refit'] = z_eta_vals
    df['Z_phi_refit'] = z_phi_vals
    df['Z_mass_refit'] = z_mass_vals
    df['H_pt_refit'] = h_pt_vals
    df['H_eta_refit'] = h_eta_vals
    df['H_phi_refit'] = h_phi_vals
    df['H_mass_refit'] = h_mass_vals

    print(f"Z constraint refit completed for {len(df)} events")
    return df
