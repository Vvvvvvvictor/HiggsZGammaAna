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

def gamma_fsr_pt_error(gamma_fsr_pt, gamma_fsr_eta, gamma_fsr_theta):
    """Calculate FSR photon pt error."""
    gamma_fsr_energy = gamma_fsr_pt / math.sin(gamma_fsr_theta)
    perr = get_energy_resolution_em(gamma_fsr_energy, gamma_fsr_eta)
    gamma_fsr_p = gamma_fsr_pt / math.tan(gamma_fsr_theta)
    pterr = perr * gamma_fsr_pt / gamma_fsr_p
    return pterr

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
                          gamma_pt_fsr=None, gamma_eta_fsr=None, gamma_phi_fsr=None,
                          has_fsr=False, lineshape_params=None):
    """Calculate negative log likelihood for Z constraint fit."""
    pt1, pt2 = pt_vals
    pt1_reco, pt2_reco = pt_reco
    pt1_err, pt2_err = pt_errors
    eta1, eta2 = eta
    phi1, phi2 = phi
    m1, m2 = mass
    
    # Gaussian constraints on pT measurements
    gauss1 = -0.5 * ((pt1 - pt1_reco) / pt1_err) ** 2
    gauss2 = -0.5 * ((pt2 - pt2_reco) / pt2_err) ** 2
    
    # Calculate Z mass with current pT values
    if has_fsr and gamma_pt_fsr is not None:
        mz = calculate_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2,
                                    gamma_pt_fsr, gamma_eta_fsr, gamma_phi_fsr)
    else:
        mz = calculate_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2)
    
    # Z lineshape constraint
    z_constraint = math.log(max(1e-10, z_lineshape_pdf(mz, lineshape_params)))
    
    return -(gauss1 + gauss2 + z_constraint)

def minimize_simple(func, x0, bounds=None):
    """Simple grid search minimization."""
    best_val = float('inf')
    best_x = x0
    
    # Simple grid search around initial values
    for i in range(-5, 6):
        for j in range(-5, 6):
            if bounds:
                x1 = max(bounds[0][0], min(bounds[0][1], x0[0] + i * 0.5))
                x2 = max(bounds[1][0], min(bounds[1][1], x0[1] + j * 0.5))
            else:
                x1 = x0[0] + i * 0.5
                x2 = x0[1] + j * 0.5
            
            try:
                val = func([x1, x2])
                if val < best_val:
                    best_val = val
                    best_x = [x1, x2]
            except:
                continue
    
    return best_x, True

def perform_z_refit(row):
    """Perform Z constraint refit for a single event."""
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
        
        # Check for FSR photon
        has_fsr = False
        gamma_pt_fsr = row.get('gamma_fsr_pt', -1)
        gamma_eta_fsr = row.get('gamma_fsr_eta', 0)
        gamma_phi_fsr = row.get('gamma_fsr_phi', 0)
        
        if gamma_pt_fsr > 0:
            # Check if FSR photon is clean (not overlapping with main photon or leptons)
            gamma_pt = row['gamma_pt']
            gamma_eta = row['gamma_eta'] 
            gamma_phi = row['gamma_phi']
            
            dr_gamma = math.sqrt((gamma_eta_fsr - gamma_eta)**2 + (gamma_phi_fsr - gamma_phi)**2)
            dr_l1 = math.sqrt((gamma_eta_fsr - l1_eta)**2 + (gamma_phi_fsr - l1_phi)**2)
            dr_l2 = math.sqrt((gamma_eta_fsr - l2_eta)**2 + (gamma_phi_fsr - l2_phi)**2)
            
            if dr_gamma > 0.001 and dr_l1 > 0.001 and dr_l2 > 0.001:
                has_fsr = True
        
        # Define objective function for minimization
        def objective(pt_vals):
            return negative_log_likelihood(
                pt_vals, [l1_pt, l2_pt], [l1_pt_err, l2_pt_err],
                [l1_eta, l2_eta], [l1_phi, l2_phi], [l1_mass, l2_mass],
                gamma_pt_fsr if has_fsr else None,
                gamma_eta_fsr if has_fsr else None,
                gamma_phi_fsr if has_fsr else None,
                has_fsr, lineshape_params
            )
        
        # Initial values and bounds
        x0 = [l1_pt, l2_pt]
        bounds = [(max(5.0, l1_pt - 2*l1_pt_err), l1_pt + 2*l1_pt_err),
                 (max(5.0, l2_pt - 2*l2_pt_err), l2_pt + 2*l2_pt_err)]
        
        # Perform minimization
        try:
            result, success = minimize_simple(objective, x0, bounds)
            if success:
                return result[0], result[1], l1_pt_err, l2_pt_err
        except:
            pass
        
        # If optimization fails, return original values
        return l1_pt, l2_pt, l1_pt_err, l2_pt_err
        
    except Exception as e:
        # If any error occurs, return original values
        return row['Z_lead_lepton_pt'], row['Z_sublead_lepton_pt'], 0.0, 0.0

def apply_z_refit(df):
    """Apply Z constraint refit to DataFrame."""
    print("Applying Z constraint refit...")
    
    # Initialize refit columns
    df['Z_lead_lepton_pt_refit'] = df['Z_lead_lepton_pt'].copy()
    df['Z_sublead_lepton_pt_refit'] = df['Z_sublead_lepton_pt'].copy()
    df['Z_lead_lepton_pt_refit_err'] = 0.0
    df['Z_sublead_lepton_pt_refit_err'] = 0.0
    
    # Apply refit with progress bar
    refit_results = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Z refit progress"):
        result = perform_z_refit(row)
        refit_results.append(result)
    
    # Update DataFrame with refit results
    for i, (pt1_refit, pt2_refit, pt1_err, pt2_err) in enumerate(refit_results):
        df.iloc[i, df.columns.get_loc('Z_lead_lepton_pt_refit')] = pt1_refit
        df.iloc[i, df.columns.get_loc('Z_sublead_lepton_pt_refit')] = pt2_refit
        df.iloc[i, df.columns.get_loc('Z_lead_lepton_pt_refit_err')] = pt1_err
        df.iloc[i, df.columns.get_loc('Z_sublead_lepton_pt_refit_err')] = pt2_err
    
    # Calculate refit Z and H kinematics
    print("Calculating refit kinematics...")
    
    for idx, row in df.iterrows():
        try:
            # Create refitted lepton 4-vectors
            l1_pt_refit = row['Z_lead_lepton_pt_refit']
            l2_pt_refit = row['Z_sublead_lepton_pt_refit']
            
            l1_eta = row['Z_lead_lepton_eta']
            l1_phi = row['Z_lead_lepton_phi']
            l1_mass = row['Z_lead_lepton_mass']
            
            l2_eta = row['Z_sublead_lepton_eta']
            l2_phi = row['Z_sublead_lepton_phi']
            l2_mass = row['Z_sublead_lepton_mass']
            
            gamma_pt = row['gamma_pt']
            gamma_eta = row['gamma_eta']
            gamma_phi = row['gamma_phi']
            gamma_mass = row['gamma_mass']
            
            # Check for clean FSR photon
            has_fsr = False
            gamma_pt_fsr = row.get('gamma_fsr_pt', -1)
            gamma_eta_fsr = row.get('gamma_fsr_eta', 0)
            gamma_phi_fsr = row.get('gamma_fsr_phi', 0)
            
            if gamma_pt_fsr > 0:
                dr_gamma = math.sqrt((gamma_eta_fsr - gamma_eta)**2 + (gamma_phi_fsr - gamma_phi)**2)
                dr_l1 = math.sqrt((gamma_eta_fsr - l1_eta)**2 + (gamma_phi_fsr - l1_phi)**2)
                dr_l2 = math.sqrt((gamma_eta_fsr - l2_eta)**2 + (gamma_phi_fsr - l2_phi)**2)
                
                if dr_gamma > 0.001 and dr_l1 > 0.001 and dr_l2 > 0.001:
                    has_fsr = True
            
            # Calculate refit Z mass
            if has_fsr:
                z_mass_refit = calculate_invariant_mass(
                    l1_pt_refit, l1_eta, l1_phi, l1_mass,
                    l2_pt_refit, l2_eta, l2_phi, l2_mass,
                    gamma_pt_fsr, gamma_eta_fsr, gamma_phi_fsr
                )
            else:
                z_mass_refit = calculate_invariant_mass(
                    l1_pt_refit, l1_eta, l1_phi, l1_mass,
                    l2_pt_refit, l2_eta, l2_phi, l2_mass
                )
            
            # Calculate refit Z pt, eta, phi
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
            
            if has_fsr:
                gamma_theta_fsr = 2 * math.atan(math.exp(-gamma_eta_fsr))
                gamma_px_fsr = gamma_pt_fsr * math.cos(gamma_phi_fsr)
                gamma_py_fsr = gamma_pt_fsr * math.sin(gamma_phi_fsr)
                gamma_pz_fsr = gamma_pt_fsr / math.tan(gamma_theta_fsr)
                gamma_e_fsr = gamma_pt_fsr / math.sin(gamma_theta_fsr)
                
                z_px = l1_px + l2_px + gamma_px_fsr
                z_py = l1_py + l2_py + gamma_py_fsr
                z_pz = l1_pz + l2_pz + gamma_pz_fsr
                z_e = l1_e + l2_e + gamma_e_fsr
            else:
                z_px = l1_px + l2_px
                z_py = l1_py + l2_py
                z_pz = l1_pz + l2_pz
                z_e = l1_e + l2_e
            
            z_pt_refit = math.sqrt(z_px**2 + z_py**2)
            z_phi_refit = math.atan2(z_py, z_px)
            z_eta_refit = -math.log(math.tan(0.5 * math.atan2(z_pt_refit, z_pz)))
            
            # Calculate refit H kinematics
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
            h_eta_refit = -math.log(math.tan(0.5 * math.atan2(h_pt_refit, h_pz)))
            h_mass_refit = math.sqrt(max(0, h_e**2 - h_px**2 - h_py**2 - h_pz**2))
            
            # Store results
            df.loc[idx, 'Z_pt_refit'] = z_pt_refit
            df.loc[idx, 'Z_eta_refit'] = z_eta_refit
            df.loc[idx, 'Z_phi_refit'] = z_phi_refit
            df.loc[idx, 'Z_mass_refit'] = z_mass_refit
            df.loc[idx, 'H_pt_refit'] = h_pt_refit
            df.loc[idx, 'H_eta_refit'] = h_eta_refit
            df.loc[idx, 'H_phi_refit'] = h_phi_refit
            df.loc[idx, 'H_mass_refit'] = h_mass_refit
            
        except Exception as e:
            # If calculation fails, use original values
            df.loc[idx, 'Z_pt_refit'] = row['Z_pt']
            df.loc[idx, 'Z_eta_refit'] = row['Z_eta'] 
            df.loc[idx, 'Z_phi_refit'] = row['Z_phi']
            df.loc[idx, 'Z_mass_refit'] = row['Z_mass']
            df.loc[idx, 'H_pt_refit'] = row['H_pt']
            df.loc[idx, 'H_eta_refit'] = row['H_eta']
            df.loc[idx, 'H_phi_refit'] = row['H_phi']
            df.loc[idx, 'H_mass_refit'] = row['H_mass']
    
    print(f"Z constraint refit completed for {len(df)} events")
    return df
