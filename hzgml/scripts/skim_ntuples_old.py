#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.21.2019
#
#
#
#
import os
import math
from argparse import ArgumentParser
from ROOT import Math, TVector2, TVector3, TLorentzVector, gROOT, gSystem
import numpy as np
#import time
import pandas as pd
import uproot
# from root_pandas import *
from tqdm import tqdm
import warnings
from z_refit import apply_z_refit
import ctypes

# Import the kinematic DNN weighter module
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'kinematicdnn'))
from kinr3_weighter import create_kin_weighter

warnings.simplefilter(action='ignore', category=FutureWarning)

# Initialize the DNN model globally
kin_weighter = create_kin_weighter()

def getArgs():
    parser = ArgumentParser(description="Skim the input ntuples for Hmumu XGBoost analysis.")
    parser.add_argument('-i', '--input', action='store', default='inputs', help='Path to the input ntuple')
    parser.add_argument('-o', '--output', action='store', default='outputs', help='Path to the output ntuple')
    parser.add_argument('--chunksize', type=int, default=500000, help='size to process at a time') 
    return  parser.parse_args()

def true_delta_phi(delta_phi):
    if delta_phi > math.pi:
        return 2 * math.pi - delta_phi
    return delta_phi
    
def true_delta_phi_vectorized(delta_phi):
    """Vectorized version of true_delta_phi"""
    return np.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)

def compute_lorentz_vectors(data, suffix=''):
    """Create all LorentzVectors at once to avoid repeated creation"""
    vectors = {}
    
    # Create Z vectors
    Z_data = data[[f'Z_pt{suffix}', f'Z_eta{suffix}', f'Z_phi{suffix}', f'Z_mass{suffix}']].values
    vectors['Z'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in Z_data])
    
    # Create H vectors
    H_data = data[[f'H_pt{suffix}', f'H_eta{suffix}', f'H_phi{suffix}', f'H_mass{suffix}']].values
    vectors['H'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in H_data])
    
    # Create gamma vectors
    gamma_data = data[['gamma_pt', 'gamma_eta', 'gamma_phi', 'gamma_mass']].values
    vectors['gamma'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in gamma_data])
    
    # Create lepton vectors
    lead_lep_pt = f'Z_lead_lepton_pt{suffix}' if f'Z_lead_lepton_pt{suffix}' in data.columns else 'Z_lead_lepton_pt'
    sublead_lep_pt = f'Z_sublead_lepton_pt{suffix}' if f'Z_sublead_lepton_pt{suffix}' in data.columns else 'Z_sublead_lepton_pt'
    
    lead_lep_data = data[[lead_lep_pt, 'Z_lead_lepton_eta', 'Z_lead_lepton_phi', 'Z_lead_lepton_mass']].values
    sublead_lep_data = data[[sublead_lep_pt, 'Z_sublead_lepton_eta', 'Z_sublead_lepton_phi', 'Z_sublead_lepton_mass']].values
    
    vectors['lead_lepton'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in lead_lep_data])
    vectors['sublead_lepton'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in sublead_lep_data])
    
    # Create jet vectors when needed
    if 'jet_1_pt' in data.columns:
        jet1_data = data[['jet_1_pt', 'jet_1_eta', 'jet_1_phi', 'jet_1_mass']].values
        jet2_data = data[['jet_2_pt', 'jet_2_eta', 'jet_2_phi', 'jet_2_mass']].values
        vectors['jet_1'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in jet1_data])
        vectors['jet_2'] = np.array([Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(*row) for row in jet2_data])
    
    return vectors

def compute_relative_variables_vectorized(data, suffix=''):
    """Compute all relative variables using vectorized operations"""
    H_mass_col = f'H_mass{suffix}'
    Z_mass_col = f'Z_mass{suffix}'
    H_pt_col = f'H_pt{suffix}'
    Z_pt_col = f'Z_pt{suffix}'
    
    lead_lep_pt = f'Z_lead_lepton_pt{suffix}' if f'Z_lead_lepton_pt{suffix}' in data.columns else 'Z_lead_lepton_pt'
    sublead_lep_pt = f'Z_sublead_lepton_pt{suffix}' if f'Z_sublead_lepton_pt{suffix}' in data.columns else 'Z_sublead_lepton_pt'
    
    rel_vars = {}
    rel_vars[f'HZ_relM{suffix}'] = data[H_mass_col] / data[Z_mass_col]
    rel_vars[f'H_relpt{suffix}'] = data[H_pt_col] / data[H_mass_col]
    rel_vars[f'Z_relpt{suffix}'] = data[Z_pt_col] / data[H_mass_col]
    rel_vars[f'Z_lead_lepton_relpt{suffix}'] = data[lead_lep_pt] / data[H_mass_col]
    rel_vars[f'Z_sublead_lepton_relpt{suffix}'] = data[sublead_lep_pt] / data[H_mass_col]
    rel_vars[f'gamma_relpt{suffix}'] = data['gamma_pt'] / data[H_mass_col]
    rel_vars[f'jet_1_relpt{suffix}'] = data['jet_1_pt'] / data[H_mass_col]
    rel_vars[f'jet_2_relpt{suffix}'] = data['jet_2_pt'] / data[H_mass_col]
    rel_vars[f'MET_relpt{suffix}'] = data['MET_pt'] / data[H_mass_col]
    
    return rel_vars

def compute_delta_phi_variables_vectorized(data):
    """Compute delta phi variables using vectorized operations"""
    delta_phi_vars = {}
    
    # H, Z, lepton delta phi with gamma
    delta_phi_vars['H_deltaphi'] = true_delta_phi_vectorized(np.abs(data['H_phi'] - data['gamma_phi']))
    delta_phi_vars['Z_deltaphi'] = true_delta_phi_vectorized(np.abs(data['Z_phi'] - data['gamma_phi']))
    delta_phi_vars['Z_lead_lepton_deltaphi'] = true_delta_phi_vectorized(np.abs(data['Z_lead_lepton_phi'] - data['gamma_phi']))
    delta_phi_vars['Z_sublead_lepton_deltaphi'] = true_delta_phi_vectorized(np.abs(data['Z_sublead_lepton_phi'] - data['gamma_phi']))
    delta_phi_vars['MET_deltaphi'] = true_delta_phi_vectorized(np.abs(data['MET_phi'] - data['gamma_phi']))
    
    # Jet delta phi with gamma
    for i in range(1, 5):
        jet_phi_col = f'jet_{i}_phi'
        if jet_phi_col in data.columns:
            mask = data['n_jets'] >= i
            delta_phi_vars[f'jet_{i}_deltaphi'] = np.where(
                mask,
                true_delta_phi_vectorized(np.abs(data[jet_phi_col] - data['gamma_phi'])),
                -9999
            )
        else:
            delta_phi_vars[f'jet_{i}_deltaphi'] = np.full(len(data), -9999)
    
    # Additional lepton delta phi
    for i in [1, 2]:
        col = f'additional_lepton_{i}_phi'
        if col in data.columns:
            delta_phi_vars[f'additional_lepton_{i}_deltaphi'] = true_delta_phi_vectorized(np.abs(data[col] - data['gamma_phi']))
        else:
            delta_phi_vars[f'additional_lepton_{i}_deltaphi'] = np.full(len(data), -9999)
    
    return delta_phi_vars

def compute_simple_jet_variables_vectorized(data):
    """Compute simple jet-related variables using vectorized operations"""
    jet_vars = {}
    
    # Jet zeppenfeld variables
    mask_2jets = data['n_jets'] >= 2
    avg_jet_eta = (data['jet_1_eta'] + data['jet_2_eta']) / 2
    
    jet_vars['photon_zeppenfeld'] = np.where(
        mask_2jets,
        np.abs(data['gamma_eta'] - avg_jet_eta),
        -9999
    )
    
    jet_vars['H_zeppenfeld'] = np.where(
        mask_2jets,
        np.abs(data['H_eta'] - avg_jet_eta),
        -9999
    )
    
    jet_vars['H_zeppenfeld_refit'] = np.where(
        mask_2jets,
        np.abs(data['H_eta_refit'] - avg_jet_eta),
        -9999
    )
    
    # Delta eta and phi between jets
    jet_vars['delta_eta_jj'] = np.where(
        mask_2jets,
        np.abs(data['jet_1_eta'] - data['jet_2_eta']),
        -9999
    )
    
    jet_vars['delta_phi_jj'] = np.where(
        mask_2jets,
        true_delta_phi_vectorized(np.abs(data['jet_1_phi'] - data['jet_2_phi'])),
        -9999
    )
    
    return jet_vars

def compute_is_center_vectorized(data, suffix=''):
    """Vectorized version of is_center computation"""
    H_mass_col = f'H_mass{suffix}'
    return ((data[H_mass_col] >= 120) & (data[H_mass_col] <= 130)).astype(int)

def compute_complex_variables_optimized(data, vectors, suffix=''):
    """Compute complex variables that need LorentzVector operations with optimization"""
    complex_vars = {}
    n_events = len(data)
    
    # Get the appropriate vectors
    Z_vecs = vectors['Z']
    H_vecs = vectors['H']
    gamma_vecs = vectors['gamma']
    lead_lep_vecs = vectors['lead_lepton']
    sublead_lep_vecs = vectors['sublead_lepton']
    
    # Compute variables that need LorentzVector operations
    # Use list comprehension with enumerate for better performance
    
    # DeltaR calculations
    complex_vars[f'll_deltaR{suffix}'] = np.array([
        Math.VectorUtil.DeltaR(lead_lep_vecs[i], sublead_lep_vecs[i]) 
        for i in range(n_events)
    ])
    
    complex_vars[f'leadLG_deltaR{suffix}'] = np.array([
        Math.VectorUtil.DeltaR(lead_lep_vecs[i], gamma_vecs[i]) 
        for i in range(n_events)
    ])
    
    complex_vars[f'subleadLG_deltaR{suffix}'] = np.array([
        Math.VectorUtil.DeltaR(sublead_lep_vecs[i], gamma_vecs[i]) 
        for i in range(n_events)
    ])
    
    complex_vars[f'ZG_deltaR{suffix}'] = np.array([
        Math.VectorUtil.DeltaR(Z_vecs[i], gamma_vecs[i]) 
        for i in range(n_events)
    ])
    
    # H_ptt, H_al, H_bt calculations
    complex_vars[f'H_ptt{suffix}'] = np.array([
        abs(Z_vecs[i].Px() * gamma_vecs[i].Py() - gamma_vecs[i].Px() * Z_vecs[i].Py()) / (Z_vecs[i] - gamma_vecs[i]).Pt() * 2.0
        for i in range(n_events)
    ])
    
    complex_vars[f'H_al{suffix}'] = np.array([
        (Z_vecs[i].Pt() ** 2 - gamma_vecs[i].Pt() ** 2) / (Z_vecs[i] - gamma_vecs[i]).Pt()
        for i in range(n_events)
    ])
    
    # H_bt calculation
    h_bt_values = []
    for i in range(n_events):
        Z, gamma = Z_vecs[i], gamma_vecs[i]
        r = abs(Z.Pt()) / abs(gamma.Pt())
        dev = ((Z.Px() - r * gamma.Px()) ** 2 + (Z.Py() - r * gamma.Py()) ** 2) ** 0.5
        h_bt_values.append(r * abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / dev * 2.0)
    complex_vars[f'H_bt{suffix}'] = np.array(h_bt_values)
    
    # Rapidity calculations
    complex_vars[f'HZ_deltaRap{suffix}'] = np.array([
        H_vecs[i].Rapidity() - Z_vecs[i].Rapidity()
        for i in range(n_events)
    ])
    
    # ECM calculations (boost to H rest frame)
    g_ecm_values = []
    z_ecm_values = []
    z_rap_cm_values = []
    for i in range(n_events):
        H_beta = TLorentzVector(H_vecs[i].Px(), H_vecs[i].Py(), H_vecs[i].Pz(), H_vecs[i].E()).BoostVector()
        g_ecm_values.append(Math.VectorUtil.boost(gamma_vecs[i], -H_beta).E())
        z_boosted = Math.VectorUtil.boost(Z_vecs[i], -H_beta)
        z_ecm_values.append(z_boosted.E())
        z_rap_cm_values.append(z_boosted.Rapidity())
    
    complex_vars[f'G_ECM{suffix}'] = np.array(g_ecm_values)
    complex_vars[f'Z_ECM{suffix}'] = np.array(z_ecm_values)
    complex_vars[f'Z_rapCM{suffix}'] = np.array(z_rap_cm_values)
    
    return complex_vars

def compute_jet_complex_variables_optimized(data, vectors, suffix=''):
    """Compute jet-related complex variables with optimization"""
    jet_vars = {}
    n_events = len(data)
    
    if 'jet_1' not in vectors or 'jet_2' not in vectors:
        # Return default values if jets not available
        for var_name in ['system_pt', 'pt_balance', 'delta_phi_zgjj', 'delta_eta_zgjj']:
            jet_vars[f'{var_name}{suffix}'] = np.full(n_events, -9999.0)
        return jet_vars
    
    Z_vecs = vectors['Z']
    gamma_vecs = vectors['gamma']
    jet1_vecs = vectors['jet_1']
    jet2_vecs = vectors['jet_2']
    
    mask_2jets = data['n_jets'] >= 2
    
    # System pt and pt balance calculations
    system_pt_values = []
    pt_balance_values = []
    delta_phi_zgjj_values = []
    delta_eta_zgjj_values = []
    
    for i in range(n_events):
        if mask_2jets.iloc[i]:
            total = Z_vecs[i] + gamma_vecs[i] + jet1_vecs[i] + jet2_vecs[i]
            system_pt_values.append(total.Pt())
            
            # pt balance
            pt_sum = data['Z_pt'].iloc[i] + data['gamma_pt'].iloc[i] + data['jet_1_pt'].iloc[i] + data['jet_2_pt'].iloc[i]
            if suffix == '_refit':
                pt_sum = data['Z_pt_refit'].iloc[i] + data['gamma_pt'].iloc[i] + data['jet_1_pt'].iloc[i] + data['jet_2_pt'].iloc[i]
            pt_balance_values.append(total.Pt() / pt_sum)
            
            # ZG-jj angles
            zg = Z_vecs[i] + gamma_vecs[i]
            j1j2 = jet1_vecs[i] + jet2_vecs[i]
            delta_phi_zgjj_values.append(abs(zg.Phi() - j1j2.Phi()))
            delta_eta_zgjj_values.append(abs(zg.Eta() - j1j2.Eta()))
        else:
            system_pt_values.append(-9999)
            pt_balance_values.append(-9999)
            delta_phi_zgjj_values.append(-9999)
            delta_eta_zgjj_values.append(-9999)
    
    jet_vars[f'system_pt{suffix}'] = np.array(system_pt_values)
    jet_vars[f'pt_balance{suffix}'] = np.array(pt_balance_values)
    jet_vars[f'delta_phi_zgjj{suffix}'] = true_delta_phi_vectorized(np.array(delta_phi_zgjj_values))
    jet_vars[f'delta_eta_zgjj{suffix}'] = np.array(delta_eta_zgjj_values)
    
    return jet_vars

def compute_optimized_decorate(data):
    """Optimized version of the decorate function using vectorized operations"""
    if data.shape[0] == 0: 
        return data
    
    # Create all LorentzVectors at once
    print("Creating LorentzVectors...")
    vectors_orig = compute_lorentz_vectors(data, '')
    vectors_refit = compute_lorentz_vectors(data, '_refit')
    
    print("Computing relative variables...")
    # Compute relative variables using vectorized operations
    rel_vars_orig = compute_relative_variables_vectorized(data, '')
    rel_vars_refit = compute_relative_variables_vectorized(data, '_refit')
    
    # Add to dataframe
    for key, value in rel_vars_orig.items():
        data[key] = value
    for key, value in rel_vars_refit.items():
        data[key] = value
    
    print("Computing delta phi variables...")
    # Compute delta phi variables
    delta_phi_vars = compute_delta_phi_variables_vectorized(data)
    for key, value in delta_phi_vars.items():
        data[key] = value
    
    print("Computing simple jet variables...")
    # Compute simple jet variables
    jet_simple_vars = compute_simple_jet_variables_vectorized(data)
    for key, value in jet_simple_vars.items():
        data[key] = value
    
    print("Computing complex variables...")
    # Compute complex variables that need LorentzVector operations
    complex_vars_orig = compute_complex_variables_optimized(data, vectors_orig, '')
    complex_vars_refit = compute_complex_variables_optimized(data, vectors_refit, '_refit')
    
    for key, value in complex_vars_orig.items():
        data[key] = value
    for key, value in complex_vars_refit.items():
        data[key] = value
    
    print("Computing jet complex variables...")
    # Compute jet complex variables
    jet_complex_orig = compute_jet_complex_variables_optimized(data, vectors_orig, '')
    jet_complex_refit = compute_jet_complex_variables_optimized(data, vectors_refit, '_refit')
    
    for key, value in jet_complex_orig.items():
        data[key] = value
    for key, value in jet_complex_refit.items():
        data[key] = value
    
    print("Computing is_center variables...")
    # Compute is_center variables
    data['is_center'] = compute_is_center_vectorized(data, '')
    data['is_center_refit'] = compute_is_center_vectorized(data, '_refit')
    
    print("Computing remaining simple variables...")
    # Add simple calculations that don't need optimization
    data['weight'] = data.weight_central
    
    # Jet pair pt (simple calculation)
    mask_2jets = data['n_jets'] >= 2
    jet_pair_pt_values = []
    for i, has_2jets in enumerate(mask_2jets):
        if has_2jets:
            jet1 = vectors_orig['jet_1'][i]
            jet2 = vectors_orig['jet_2'][i]
            jet_pair_pt_values.append((jet1 + jet2).Pt())
        else:
            jet_pair_pt_values.append(-9999)
    data['jet_pair_pt'] = jet_pair_pt_values
    
    # Add remaining apply-based calculations only for those that are truly complex
    # and cannot be easily vectorized (like the very complex angular calculations)
    complex_angular_vars = [
        'gamma_ptRelErr', 'l_rapCM', 'l_cosProdAngle', 'Z_cosProdAngle',
        'mass_jj', 'max_two_jet_btag', 'jet_ptt',
        'Z_cos_theta', 'lep_cos_theta', 'lep_phi',
        'l1g_deltaR', 'l2g_deltaR', 'pt_balance_0j', 'pt_balance_1j',
        'photon_mht_dphi', 'kin_weight'
    ]
    
    # For now, keep the original apply-based approach for the most complex calculations
    # These could be further optimized in future iterations
    print("Computing remaining complex angular variables...")
    data['gamma_ptRelErr'] = data.apply(lambda x: compute_gamma_relEerror(x), axis=1)
    
    # Add refit versions with suffix
    for var in ['l_rapCM', 'l_cosProdAngle', 'Z_cosProdAngle', 'Z_cos_theta', 'lep_cos_theta', 'lep_phi']:
        data[var] = data.apply(lambda x: globals()[f'compute_{var.replace("_", "_").split("_")[0]}_{var.replace("_", "_").split("_")[1] if len(var.split("_")) > 1 else var}'](x), axis=1)
        if f'compute_{var}_refit' in globals():
            data[f'{var}_refit'] = data.apply(lambda x: globals()[f'compute_{var}_refit'](x), axis=1)
    
    return data

    if (x.H_mass >= 120 and x.H_mass <= 130): return 1
    else: return 0

def compute_is_center_refit(x):

    if (x.H_mass_refit >= 120 and x.H_mass_refit <= 130): return 1
    else: return 0

def compute_photon_mht_dphi(x):
    """Compute delta phi between photon and MHT"""
    return true_delta_phi(abs(x.MHT_phi - x.gamma_phi))

def compute_kin_weight(x):
    """Compute kinematic weight using DNN model"""
    # Input order: gamma_mvaID, photon_mht_dphi, n_jets, Z_pt, gamma_pt, H_pt, MHT_pt, HT
    try:
        inputs = [
            x.gamma_mvaID,
            x.photon_mht_dphi,
            float(x.n_jets),
            x.Z_pt,
            x.gamma_pt,
            x.H_pt,
            x.MHT_pt,
            x.HT
        ]
        
        # Check for any invalid inputs
        if any(math.isnan(val) or math.isinf(val) for val in inputs):
            return 1.0  # Default weight
            
        return kin_weighter.evaluate(inputs)
    except:
        return 1.0  # Default weight if computation fails

def compute_kin_weight_refit(x):
    """Compute kinematic weight using DNN model with refit variables"""
    # Input order: gamma_mvaID, photon_mht_dphi, n_jets, Z_pt_refit, gamma_pt, H_pt_refit, MHT_pt, HT
    try:
        inputs = [
            x.gamma_mvaID,
            x.photon_mht_dphi,
            float(x.n_jets),
            x.Z_pt_refit,
            x.gamma_pt,
            x.H_pt_refit,
            x.MHT_pt,
            x.HT
        ]
        
        # Check for any invalid inputs
        if any(math.isnan(val) or math.isinf(val) for val in inputs):
            return 1.0  # Default weight
            
        return kin_weighter.evaluate(inputs)
    except:
        return 1.0  # Default weight if computation fails

# others

def compute_Z_cosTheta(x):
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    M, mll = x.H_mass, x.Z_mass
    lZ = math.sqrt((H.Dot(Z) / M) ** 2 - mll ** 2)

    H_transverse_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    H_transverse_beta.SetZ(0)
    hH = Math.VectorUtil.boost(H, -H_transverse_beta)
    hPz, hE = hH.Pz(), hH.E()
    q = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz + hE) / 2, (hE + hPz) / 2)
    q = Math.VectorUtil.boost(q, H_transverse_beta)
    qbar = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz - hE) / 2, (hE - hPz) / 2)
    qbar = Math.VectorUtil.boost(qbar, H_transverse_beta)

    cosTheta = (qbar - q).Dot(Z)/(M * lZ)

    return cosTheta

def compute_Z_cosTheta_refit(x):
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    M, mll = x.H_mass_refit, x.Z_mass_refit
    lZ = math.sqrt((H.Dot(Z) / M) ** 2 - mll ** 2)

    H_transverse_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    H_transverse_beta.SetZ(0)
    hH = Math.VectorUtil.boost(H, -H_transverse_beta)
    hPz, hE = hH.Pz(), hH.E()
    q = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz + hE) / 2, (hE + hPz) / 2)
    q = Math.VectorUtil.boost(q, H_transverse_beta)
    qbar = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz - hE) / 2, (hE - hPz) / 2)
    qbar = Math.VectorUtil.boost(qbar, H_transverse_beta)

    cosTheta = (qbar - q).Dot(Z)/(M * lZ)

    return cosTheta

def compute_l_costheta(x):
    if (x.Z_lead_lepton_charge < 0 and x.Z_sublead_lepton_charge > 0): 
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    elif (x.Z_sublead_lepton_charge < 0 and x.Z_lead_lepton_charge > 0):
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    else: print('leptons have same sign')
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    Z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    l1_BZ = Math.VectorUtil.boost(l1, -Z_beta)
    gamma_BZ = Math.VectorUtil.boost(gamma, -Z_beta)

    ## photon and lepton
    costheta = - (gamma_BZ.Vect()).Dot(l1_BZ.Vect())/gamma_BZ.Vect().R()/l1_BZ.Vect().R()

    return costheta

def compute_l_costheta_refit(x):
    if (x.Z_lead_lepton_charge < 0 and x.Z_sublead_lepton_charge > 0): 
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    elif (x.Z_sublead_lepton_charge < 0 and x.Z_lead_lepton_charge > 0):
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    else: print('leptons have same sign')
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    Z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    l1_BZ = Math.VectorUtil.boost(l1, -Z_beta)
    gamma_BZ = Math.VectorUtil.boost(gamma, -Z_beta)

    ## photon and lepton
    costheta = - (gamma_BZ.Vect()).Dot(l1_BZ.Vect())/gamma_BZ.Vect().R()/l1_BZ.Vect().R()

    return costheta

def compute_l_phi(x):
    if (x.Z_lead_lepton_charge < 0 and x.Z_sublead_lepton_charge > 0): 
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    elif (x.Z_sublead_lepton_charge < 0 and x.Z_lead_lepton_charge > 0):
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    else: print('leptons have same sign')
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    
    M, mll = x.H_mass, x.Z_mass
    lZ = math.sqrt((H.Dot(Z) / M) ** 2 - mll ** 2)

    H_transverse_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    H_transverse_beta.SetZ(0)
    hH = Math.VectorUtil.boost(H, -H_transverse_beta)
    hPz, hE = hH.Pz(), hH.E()
    q = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz + hE) / 2, (hE + hPz) / 2)
    q = Math.VectorUtil.boost(q, H_transverse_beta)
    qbar = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz - hE) / 2, (hE - hPz) / 2)
    qbar = Math.VectorUtil.boost(qbar, H_transverse_beta)

    cosTheta = (qbar - q).Dot(Z)/(M * lZ)
    if (abs(cosTheta) > 1): 
        print("cosTheta = ", cosTheta)
        cosTheta = 1./cosTheta
    sinTheta = math.sqrt(1 - cosTheta ** 2)
    
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    l1_BH = Math.VectorUtil.boost(l1, -H_beta)
    l2_BH = Math.VectorUtil.boost(l2, -H_beta)
    Z_BH = Math.VectorUtil.boost(Z, -H_beta)
    N1_BH = l1_BH.Vect().Cross(l2_BH.Vect())
    
    q = Math.VectorUtil.boost(q, -H_beta)
    
    NSC_BH = (q.Vect().Cross(Z_BH.Vect()))
    tmpSgnPhi1 = - N1_BH.Dot(q.Vect())/q.Vect().R()/N1_BH.R()/sinTheta
    
    sgnPhi1 = 0.
    if (abs(tmpSgnPhi1)>0.): sgnPhi1 = tmpSgnPhi1/abs(tmpSgnPhi1)
    dot_BH1SC = - N1_BH.Dot(NSC_BH)/NSC_BH.R()/N1_BH.R()
    if (abs(dot_BH1SC)>=1.): 
      print("dot_BH1SC = ", dot_BH1SC)
      dot_BH1SC *= 1./abs(dot_BH1SC)
    Phi1 = sgnPhi1 * math.acos(dot_BH1SC)

    return Phi1

def compute_l_phi_refit(x):
    if (x.Z_lead_lepton_charge < 0 and x.Z_sublead_lepton_charge > 0): 
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    elif (x.Z_sublead_lepton_charge < 0 and x.Z_lead_lepton_charge > 0):
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    else: print('leptons have same sign')
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    
    M, mll = x.H_mass_refit, x.Z_mass_refit
    lZ = math.sqrt((H.Dot(Z) / M) ** 2 - mll ** 2)

    H_transverse_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    H_transverse_beta.SetZ(0)
    hH = Math.VectorUtil.boost(H, -H_transverse_beta)
    hPz, hE = hH.Pz(), hH.E()
    q = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz + hE) / 2, (hE + hPz) / 2)
    q = Math.VectorUtil.boost(q, H_transverse_beta)
    qbar = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz - hE) / 2, (hE - hPz) / 2)
    qbar = Math.VectorUtil.boost(qbar, H_transverse_beta)

    cosTheta = (qbar - q).Dot(Z)/(M * lZ)
    if (abs(cosTheta) > 1): 
        print("cosTheta = ", cosTheta)
        cosTheta = 1./cosTheta
    sinTheta = math.sqrt(1 - cosTheta ** 2)
    
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    l1_BH = Math.VectorUtil.boost(l1, -H_beta)
    l2_BH = Math.VectorUtil.boost(l2, -H_beta)
    Z_BH = Math.VectorUtil.boost(Z, -H_beta)
    N1_BH = l1_BH.Vect().Cross(l2_BH.Vect())
    
    q = Math.VectorUtil.boost(q, -H_beta)
    
    NSC_BH = (q.Vect().Cross(Z_BH.Vect()))
    tmpSgnPhi1 = - N1_BH.Dot(q.Vect())/q.Vect().R()/N1_BH.R()/sinTheta
    
    sgnPhi1 = 0.
    if (abs(tmpSgnPhi1)>0.): sgnPhi1 = tmpSgnPhi1/abs(tmpSgnPhi1)
    dot_BH1SC = - N1_BH.Dot(NSC_BH)/NSC_BH.R()/N1_BH.R()
    if (abs(dot_BH1SC)>=1.): 
      print("dot_BH1SC = ", dot_BH1SC)
      dot_BH1SC *= 1./abs(dot_BH1SC)
    Phi1 = sgnPhi1 * math.acos(dot_BH1SC)

    return Phi1

def compute_dR1lg(x):

    if (((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5) > ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5:
        return ((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5
    else:
        return ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5
    
def compute_dR2lg(x):

    if (((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5) < ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5:
        return ((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5
    else:
        return ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5

# mine

def compute_l_prodAngle(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    if (x.Z_lead_lepton_charge > 0):
        a = Math.VectorUtil.boost(l1, -z_beta).Vect()
    else :
        a = Math.VectorUtil.boost(l2, -z_beta).Vect()
    return a.Unit().Dot(z_beta.Unit())

def compute_l_prodAngle_refit(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    if (x.Z_lead_lepton_charge > 0):
        a = Math.VectorUtil.boost(l1, -z_beta).Vect()
    else :
        a = Math.VectorUtil.boost(l2, -z_beta).Vect()
    return a.Unit().Dot(z_beta.Unit())

def compute_Z_prodAngle(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    Z_BH = Math.VectorUtil.boost(Z, -H_beta).Vect()
    gamma_BH = Math.VectorUtil.boost(gamma, -H_beta).Vect()
    return Z_BH.Unit().Dot(H.Vect().Unit())

def compute_Z_prodAngle_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    Z_BH = Math.VectorUtil.boost(Z, -H_beta).Vect()
    gamma_BH = Math.VectorUtil.boost(gamma, -H_beta).Vect()
    return Z_BH.Unit().Dot(H.Vect().Unit())

def compute_gamma_relEerror(x):

    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return x.gamma_energyErr / gamma.E()

def compute_G_ECM(x):

    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    return Math.VectorUtil.boost(gamma, -H_beta).E()

def compute_G_ECM_refit(x):

    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    return Math.VectorUtil.boost(gamma, -H_beta).E()

def compute_Z_ECM(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    return Math.VectorUtil.boost(Z, -H_beta).E()

def compute_Z_ECM_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    return Math.VectorUtil.boost(Z, -H_beta).E()

def compute_l_rapCM(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    a = Math.VectorUtil.boost(l1, -z_beta)
    return a.Rapidity()

def compute_l_rapCM_refit(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    a = Math.VectorUtil.boost(l1, -z_beta)
    return a.Rapidity()

def compute_Z_rapCM(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    a = Math.VectorUtil.boost(Z, -H_beta)
    return a.Rapidity()

def compute_Z_rapCM_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    a = Math.VectorUtil.boost(Z, -H_beta)
    return a.Rapidity()

def compute_HZ_deltaRap(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    return H.Rapidity() - Z.Rapidity()

def compute_HZ_deltaRap_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt_refit, x.H_eta_refit, x.H_phi_refit, x.H_mass_refit)
    return H.Rapidity() - Z.Rapidity()

def compute_ll_deltaR(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    return Math.VectorUtil.DeltaR(l1, l2)

def compute_ll_deltaR_refit(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    return Math.VectorUtil.DeltaR(l1, l2)

def compute_leadLG_deltaR(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return Math.VectorUtil.DeltaR(l1, gamma)

def compute_leadLG_deltaR_refit(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt_refit, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return Math.VectorUtil.DeltaR(l1, gamma)

def compute_subleadLG_deltaR(x):

    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return Math.VectorUtil.DeltaR(l2, gamma)

def compute_subleadLG_deltaR_refit(x):

    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt_refit, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return Math.VectorUtil.DeltaR(l2, gamma)

def compute_ZG_deltaR(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return Math.VectorUtil.DeltaR(Z, gamma)

def compute_ZG_deltaR_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return Math.VectorUtil.DeltaR(Z, gamma)

def compute_H_ptt(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / (Z - gamma).Pt() * 2.0

def compute_H_ptt_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / (Z - gamma).Pt() * 2.0


def compute_H_al(x):
    
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return (Z.Pt() ** 2 - gamma.Pt() ** 2) / (Z - gamma).Pt()

def compute_H_al_refit(x):
    
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return (Z.Pt() ** 2 - gamma.Pt() ** 2) / (Z - gamma).Pt()

def compute_H_bt(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    r = abs(Z.Pt()) / abs(gamma.Pt())
    dev = ((Z.Px() - r * gamma.Px()) ** 2 + (Z.Py() - r * gamma.Py()) ** 2) ** 0.5

    return r * abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / dev * 2.0

def compute_H_bt_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    r = abs(Z.Pt()) / abs(gamma.Pt())
    dev = ((Z.Px() - r * gamma.Px()) ** 2 + (Z.Py() - r * gamma.Py()) ** 2) ** 0.5

    return r * abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / dev * 2.0

def compute_jet_ptt(x):

    if x.n_jets < 2:
        return -9999
    else:
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)

        return abs(Jets_1.Px() * Jets_2.Py() - Jets_2.Px() * Jets_1.Py()) / (Jets_1 - Jets_2).Pt() * 2.0

def compute_max_two_jet_btag(x):

    n_jets = x.n_jets
    jet_btag_array = []

    for i in range(1, min(n_jets+1, 5)):
        jet_name = f"jet_{i}"
        jet_btag_array.append(getattr(x, f"{jet_name}_btagDeepFlavB"))

    jet_btag_array = np.array(jet_btag_array)
    jet_btag_array_sorted = np.sort(jet_btag_array)[::-1]

    return np.sum(jet_btag_array_sorted[:2])


def compute_mass_jj(x):

    if x.n_jets < 2:
        return -9999
    else:
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return j1j2.M()

def compute_delta_eta_jj(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.jet_1_eta-x.jet_2_eta)

def compute_delta_phi_jj(x):

    if x.n_jets >= 2:
        return true_delta_phi(abs(x.jet_1_phi-x.jet_2_phi))
    return -9999

def compute_delta_phi_zg_jj(x):

    if x.n_jets >= 2:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        zg = Z+gamma
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return true_delta_phi(abs(zg.Phi()-j1j2.Phi()))
    return -9999

def compute_delta_phi_zg_jj_refit(x):

    if x.n_jets >= 2:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        zg = Z+gamma
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return true_delta_phi(abs(zg.Phi()-j1j2.Phi()))
    return -9999

def compute_delta_eta_zg_jj(x):

    if x.n_jets >= 2:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        zg = Z+gamma
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return abs(zg.Eta()-j1j2.Eta())
    return -9999

def compute_delta_eta_zg_jj_refit(x):

    if x.n_jets >= 2:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        zg = Z+gamma
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return abs(zg.Eta()-j1j2.Eta())
    return -9999

def compute_photon_zeppenfeld(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.gamma_eta-(x.jet_1_eta+x.jet_2_eta)/2)

def compute_H_zeppenfeld(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.H_eta-(x.jet_1_eta+x.jet_2_eta)/2)

def compute_H_zeppenfeld_refit(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.H_eta_refit-(x.jet_1_eta+x.jet_2_eta)/2)

def compute_pt_balance(x):

    if x.n_jets < 2:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        total = Z+gamma+Jet_1+Jet_2
        return total.Pt()/(x.Z_pt+x.gamma_pt+x.jet_1_pt+x.jet_2_pt)

def compute_pt_balance_refit(x):

    if x.n_jets < 2:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        total = Z+gamma+Jet_1+Jet_2
        return total.Pt()/(x.Z_pt_refit+x.gamma_pt+x.jet_1_pt+x.jet_2_pt)

def compute_system_pt(x):
    
    if x.n_jets < 2:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        total = Z+gamma+Jet_1+Jet_2
        return total.Pt()

def compute_system_pt_refit(x):
    
    if x.n_jets < 2:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        total = Z+gamma+Jet_1+Jet_2
        return total.Pt()

def compute_jet_pair_pt(x):
    
    if x.n_jets < 2:
        return -9999
    else:
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        Jet_pair = Jet_1+Jet_2
        return Jet_pair.Pt()

def compute_pt_balance_0j(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    total = Z+gamma
    return total.Pt()/(x.Z_pt+x.gamma_pt)

def compute_pt_balance_0j_refit(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    total = Z+gamma
    return total.Pt()/(x.Z_pt_refit+x.gamma_pt)

def compute_pt_balance_1j(x):

    if x.n_jets < 1:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        total = Z+gamma+Jets_1
        return total.Pt()/(x.Z_pt+x.gamma_pt+x.jet_1_pt)

def compute_pt_balance_1j_refit(x):

    if x.n_jets < 1:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt_refit, x.Z_eta_refit, x.Z_phi_refit, x.Z_mass_refit)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        total = Z+gamma+Jets_1
        return total.Pt()/(x.Z_pt_refit+x.gamma_pt+x.jet_1_pt)

def compute_Delta_Phi(x, var = "gamma_phi", min_jet=0):

    if min_jet:
        if x.n_jets < min_jet: return -9999
        return true_delta_phi(abs(getattr(x, "jet_%d_phi" % min_jet) - x.gamma_phi))
    else:
        return true_delta_phi(abs(getattr(x, var) - x.gamma_phi))

def compute_Delta_R(x, min_jet=0):

    if x.n_jets < min_jet: return -9999
    else: 
        jet = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(getattr(x, "jet_%d_pt" % min_jet), getattr(x, "jet_%d_eta" % min_jet), getattr(x, "jet_%d_phi" % min_jet), getattr(x, "jet_%d_mass" % min_jet))
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

        return Math.VectorUtil.DeltaR(jet, gamma)

# default

def compute_QG(x):

    if x.Jets_jetMultip >= 1 and (abs(x.Jets_Eta_Lead) > 2.1 or x.Jets_PT_Lead < 50):
        Jets_QGscore_Lead, Jets_QGflag_Lead = -1, -1
    else:
        Jets_QGscore_Lead = x.Jets_NTracks_Lead
        Jets_QGflag_Lead = np.heaviside(x.Jets_NTracks_Lead - 11, 0) + np.heaviside(-x.Jets_NTracks_Lead - 9999, -9999)

    if x.Jets_jetMultip >= 2 and (abs(x.Jets_Eta_Sub) > 2.1 or x.Jets_PT_Sub < 50):
        Jets_QGscore_Sub, Jets_QGflag_Sub = -1, -1
    else:
        Jets_QGscore_Sub = x.Jets_NTracks_Sub
        Jets_QGflag_Sub = np.heaviside(x.Jets_NTracks_Sub - 11, 0) + np.heaviside(-x.Jets_NTracks_Sub - 9999, -9999)

    return Jets_QGscore_Lead, Jets_QGflag_Lead, Jets_QGscore_Sub, Jets_QGflag_Sub

def preselect(data):

    #data.query('(Muons_Minv_MuMu_Paper >= 110) | (Event_Paper_Category >= 17)', inplace=True)
    #data.query('Event_Paper_Category > 0', inplace=True)

    return data

def decorate(data):

    if data.shape[0] == 0: return data

    data['jet_pair_pt'] = data.apply(lambda x: compute_jet_pair_pt(x), axis=1)
    data['system_pt'] = data.apply(lambda x: compute_system_pt(x), axis=1)
    data['HZ_relM'] = data.H_mass / data.Z_mass
    data['H_relpt'] = data.H_pt / data.H_mass
    data['Z_relpt'] = data.Z_pt / data.H_mass
    data['Z_lead_lepton_relpt'] = data.Z_lead_lepton_pt / data.H_mass
    data['Z_sublead_lepton_relpt'] = data.Z_sublead_lepton_pt / data.H_mass
    data['gamma_relpt'] = data.gamma_pt / data.H_mass
    data['jet_1_relpt'] = data.jet_1_pt / data.H_mass
    data['jet_2_relpt'] = data.jet_2_pt / data.H_mass
    data['MET_relpt'] = data.MET_pt / data.H_mass
    
    # Add refit relative variables  
    data['system_pt_refit'] = data.apply(lambda x: compute_system_pt_refit(x), axis=1)
    data['HZ_relM_refit'] = data.H_mass_refit / data.Z_mass_refit
    data['H_relpt_refit'] = data.H_pt_refit / data.H_mass_refit
    data['Z_relpt_refit'] = data.Z_pt_refit / data.H_mass_refit
    data['Z_lead_lepton_relpt_refit'] = data.Z_lead_lepton_pt_refit / data.H_mass_refit
    data['Z_sublead_lepton_relpt_refit'] = data.Z_sublead_lepton_pt_refit / data.H_mass_refit
    data['gamma_relpt_refit'] = data.gamma_pt / data.H_mass_refit
    data['jet_1_relpt_refit'] = data.jet_1_pt / data.H_mass_refit
    data['jet_2_relpt_refit'] = data.jet_2_pt / data.H_mass_refit
    data['MET_relpt_refit'] = data.MET_pt / data.H_mass_refit
    data['gamma_ptRelErr'] = data.apply(lambda x:compute_gamma_relEerror(x), axis=1)
    data['G_ECM'] = data.apply(lambda x:compute_G_ECM(x), axis=1)
    data['Z_ECM'] = data.apply(lambda x:compute_Z_ECM(x), axis=1)
    data['Z_rapCM'] = data.apply(lambda x:compute_Z_rapCM(x), axis=1)
    data['l_rapCM'] = data.apply(lambda x:compute_l_rapCM(x), axis=1)
    data['HZ_deltaRap'] = data.apply(lambda x:compute_HZ_deltaRap(x), axis=1)
    data['l_cosProdAngle'] = data.apply(lambda x:compute_l_prodAngle(x), axis=1)
    data['Z_cosProdAngle'] = data.apply(lambda x:compute_Z_prodAngle(x), axis=1)
    data['ll_deltaR'] = data.apply(lambda x:compute_ll_deltaR(x), axis=1)
    data['leadLG_deltaR'] = data.apply(lambda x:compute_leadLG_deltaR(x), axis=1)
    data['ZG_deltaR'] = data.apply(lambda x:compute_ZG_deltaR(x), axis=1)
    data['subleadLG_deltaR'] = data.apply(lambda x:compute_subleadLG_deltaR(x), axis=1)
    
    # Add refit versions
    data['G_ECM_refit'] = data.apply(lambda x:compute_G_ECM_refit(x), axis=1)
    data['Z_ECM_refit'] = data.apply(lambda x:compute_Z_ECM_refit(x), axis=1)
    data['Z_rapCM_refit'] = data.apply(lambda x:compute_Z_rapCM_refit(x), axis=1)
    data['l_rapCM_refit'] = data.apply(lambda x:compute_l_rapCM_refit(x), axis=1)
    data['HZ_deltaRap_refit'] = data.apply(lambda x:compute_HZ_deltaRap_refit(x), axis=1)
    data['l_cosProdAngle_refit'] = data.apply(lambda x:compute_l_prodAngle_refit(x), axis=1)
    data['Z_cosProdAngle_refit'] = data.apply(lambda x:compute_Z_prodAngle_refit(x), axis=1)
    data['ll_deltaR_refit'] = data.apply(lambda x:compute_ll_deltaR_refit(x), axis=1)
    data['leadLG_deltaR_refit'] = data.apply(lambda x:compute_leadLG_deltaR_refit(x), axis=1)
    data['ZG_deltaR_refit'] = data.apply(lambda x:compute_ZG_deltaR_refit(x), axis=1)
    data['subleadLG_deltaR_refit'] = data.apply(lambda x:compute_subleadLG_deltaR_refit(x), axis=1)
    data['H_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'H_phi'), axis=1)
    data['Z_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'Z_phi'), axis=1)
    data['Z_lead_lepton_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'Z_lead_lepton_phi'), axis=1)
    data['Z_sublead_lepton_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'Z_sublead_lepton_phi'), axis=1)
    for i in np.arange(1,5):
        data['jet_%d_deltaphi' %i] = data.apply(lambda x: compute_Delta_Phi(x, "jet", min_jet=i), axis=1)
        data['jet%dG_deltaR' %i] = data.apply(lambda x: compute_Delta_R(x, min_jet=i), axis=1)
    data['max_jet_deltaR'] = data[['jet1G_deltaR', 'jet2G_deltaR']].max(axis=1)
    data['min_jet_deltaR'] = data[['jet1G_deltaR', 'jet2G_deltaR']].min(axis=1)
    data['additional_lepton_1_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'additional_lepton_1_phi', min_jet=0), axis=1)
    data['additional_lepton_2_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'additional_lepton_2_phi', min_jet=0), axis=1) 
    data['MET_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'MET_phi'), axis=1)
    data['weight'] = data.weight_central
    data['mass_jj'] = data.apply(lambda x: compute_mass_jj(x), axis=1)
    data['max_two_jet_btag'] = data.apply(lambda x: compute_max_two_jet_btag(x), axis=1)
    data['jet_ptt'] = data.apply(lambda x: compute_jet_ptt(x), axis=1)
    data['H_ptt'] = data.apply(lambda x: compute_H_ptt(x), axis=1)
    data['H_al'] = data.apply(lambda x: compute_H_al(x), axis=1)
    data['H_bt'] = data.apply(lambda x: compute_H_bt(x), axis=1)
    data['Z_cos_theta'] = data.apply(lambda x:compute_Z_cosTheta(x), axis=1)
    data['lep_cos_theta'] = data.apply(lambda x: compute_l_costheta(x), axis=1)
    data['lep_phi'] = data.apply(lambda x: compute_l_phi(x), axis=1)
    
    # Add refit versions
    data['H_ptt_refit'] = data.apply(lambda x: compute_H_ptt_refit(x), axis=1)
    data['H_al_refit'] = data.apply(lambda x: compute_H_al_refit(x), axis=1)
    data['H_bt_refit'] = data.apply(lambda x: compute_H_bt_refit(x), axis=1)
    data['Z_cos_theta_refit'] = data.apply(lambda x:compute_Z_cosTheta_refit(x), axis=1)
    data['lep_cos_theta_refit'] = data.apply(lambda x: compute_l_costheta_refit(x), axis=1)
    data['lep_phi_refit'] = data.apply(lambda x: compute_l_phi_refit(x), axis=1)
    data['l1g_deltaR'] = data.apply(lambda x: compute_dR1lg(x), axis=1) 
    data['l2g_deltaR'] = data.apply(lambda x: compute_dR2lg(x), axis=1)
    data['delta_eta_jj'] = data.apply(lambda x: compute_delta_eta_jj(x), axis=1)
    data['delta_phi_jj'] = data.apply(lambda x: compute_delta_phi_jj(x), axis=1)
    data['delta_phi_zgjj'] = data.apply(lambda x: compute_delta_phi_zg_jj(x), axis=1)
    data['delta_eta_zgjj'] = data.apply(lambda x: compute_delta_eta_zg_jj(x), axis=1)
    data['photon_zeppenfeld'] = data.apply(lambda x: compute_photon_zeppenfeld(x), axis=1)
    data['H_zeppenfeld'] = data.apply(lambda x: compute_H_zeppenfeld(x), axis=1)
    data['pt_balance'] = data.apply(lambda x: compute_pt_balance(x), axis=1)
    data['pt_balance_0j'] = data.apply(lambda x: compute_pt_balance_0j(x), axis=1)
    data['pt_balance_1j'] = data.apply(lambda x: compute_pt_balance_1j(x), axis=1)
    data['is_center'] = data.apply(lambda x: compute_is_center(x), axis=1)
    
    # Add refit versions
    data['delta_phi_zgjj_refit'] = data.apply(lambda x: compute_delta_phi_zg_jj_refit(x), axis=1)
    data['delta_eta_zgjj_refit'] = data.apply(lambda x: compute_delta_eta_zg_jj_refit(x), axis=1)
    data['H_zeppenfeld_refit'] = data.apply(lambda x: compute_H_zeppenfeld_refit(x), axis=1)
    data['pt_balance_refit'] = data.apply(lambda x: compute_pt_balance_refit(x), axis=1)
    data['pt_balance_0j_refit'] = data.apply(lambda x: compute_pt_balance_0j_refit(x), axis=1)
    data['pt_balance_1j_refit'] = data.apply(lambda x: compute_pt_balance_1j_refit(x), axis=1)
    data['is_center_refit'] = data.apply(lambda x: compute_is_center_refit(x), axis=1)
    
    # Add kinematic DNN weighting
    data['photon_mht_dphi'] = data.apply(lambda x: compute_photon_mht_dphi(x), axis=1)
    data['kin_weight'] = data.apply(lambda x: compute_kin_weight(x), axis=1)
    data['kin_weight_refit'] = data.apply(lambda x: compute_kin_weight_refit(x), axis=1)
    #data[['Jets_QGscore_Lead', 'Jets_QGflag_Lead', 'Jets_QGscore_Sub', 'Jets_QGflag_Sub']] = data.apply(lambda x: compute_QG(x), axis=1, result_type='expand')

    #data.rename(columns={'Muons_Minv_MuMu_Paper': 'm_mumu', 'Muons_Minv_MuMu_VH': 'm_mumu_VH', 'EventInfo_EventNumber': 'eventNumber', 'Jets_jetMultip': 'n_j'}, inplace=True)
    #data.drop(['PassesttHSelection', 'PassesVHSelection', 'GlobalWeight', 'SampleOverlapWeight', 'EventWeight_MCCleaning_5'], axis=1, inplace=True)
    data = data.astype(float)
    data = data.astype({'is_center': int, 'is_center_refit': int, 'Z_lead_lepton_charge': int, 'Z_lead_lepton_id': int, 'Z_sublead_lepton_charge': int, 'Z_sublead_lepton_id': int, "n_jets": int, "n_b_jets": int, "n_leptons": int, "n_electrons": int, "n_muons": int, 'event': int})

    return data
    

def main():
    
    args = getArgs()

    variables = [
        'H_pt', 'H_eta', 'H_phi', 'H_mass',
        'Z_pt', 'Z_eta', 'Z_phi', 'Z_mass',
        'Z_lead_lepton_pt', 'Z_lead_lepton_eta', 'Z_lead_lepton_phi', 'Z_lead_lepton_mass',
        'Z_lead_lepton_charge', 'Z_lead_lepton_id', 'Z_lead_lepton_ptE_error',
        'Z_sublead_lepton_pt', 'Z_sublead_lepton_eta', 'Z_sublead_lepton_phi', 'Z_sublead_lepton_mass',
        'Z_sublead_lepton_charge', 'Z_sublead_lepton_id', 'Z_sublead_lepton_ptE_error',
        # Add refit variables
        'H_pt_refit', 'H_eta_refit', 'H_phi_refit', 'H_mass_refit',
        'Z_pt_refit', 'Z_eta_refit', 'Z_phi_refit', 'Z_mass_refit',
        'Z_lead_lepton_pt_refit', 'Z_sublead_lepton_pt_refit',
        'Z_lead_lepton_pt_refit_err', 'Z_sublead_lepton_pt_refit_err',
        'gamma_pt', 'gamma_eta', 'gamma_phi', 'gamma_mass',
        'gamma_mvaID', 'gamma_energyErr',
        "gamma_fsr_pt","gamma_fsr_eta","gamma_fsr_phi","gamma_fsr_mass","gamma_fsr_clean",
        'jet_1_pt', 'jet_1_eta', 'jet_1_phi', 'jet_1_mass', 'jet_1_btagDeepFlavB',
        'jet_2_pt', 'jet_2_eta', 'jet_2_phi', 'jet_2_mass', 'jet_2_btagDeepFlavB',
        'n_jets', 'n_leptons', 'n_electrons', 'n_muons', 'n_iso_photons', 'n_b_jets',
        'MET_pt', 'MET_phi', 
        'MHT_pt', 'MHT_phi', 'HT',  # Added for kinematic DNN
        'weight_central',
        'event'
    ]

    print("Now!!! Porcessing {:s} ......".format(args.input))

    if os.path.isfile(args.output): os.remove(args.output)

    initial_events = 0
    final_events = 0

#    for data in tqdm(read_root(args.input, key='DiMuonNtuple', columns=variables, chunksize=args.chunksize), desc='Processing %s' % args.input, bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):

    data = pd.read_parquet(args.input)
    # Remove columns with object dtype
    data = data.drop(columns=data.select_dtypes(include=['object']).columns)
    initial_events += data.shape[0]
    #data = preprocess(data)
    data = preselect(data) #TODO add cutflow
    # Apply Z constraint refit
    data = apply_z_refit(data)
    data = decorate(data)
    
    final_events += data.shape[0]
    data_zero_jet = data.query("n_jets == 0 & n_leptons == 2 & MET_pt < 90")
    data_one_jet = data.query("n_jets == 1 & n_leptons == 2 & MET_pt < 90")
    data_two_jet = data.query("n_jets >= 2 & n_leptons == 2 & n_b_jets == 0")
    data_zero_to_one_jet = data.query("n_jets <= 1 & n_leptons == 2 & MET_pt < 90")
    data_VH_ttH = data[(data.n_leptons > 2) & (data.n_b_jets == 0)]
    data_VH =  data.query("n_leptons >= 3 & n_b_jets == 0 & max_I_mini < 0.15 & H_relpt > 0.3 & MET_pt > 30 & Z_mass > 85 & Z_mass < 95")
    data_ZH = data.query("n_leptons == 2 & n_jets <= 1 & MET_pt > 90 & H_relpt > 0.4 & Z_mass > 85 & Z_mass < 95") 
    data_ttH_had = data.query("n_leptons == 2 & n_jets >= 5 & n_b_jets >= 1 & Z_mass > 85 & Z_mass < 95")
    data_ttH_lep = data.query("((n_leptons == 3 & n_jets >= 3 & n_b_jets >= 1) | (n_leptons >= 4 & n_jets >= 1 & n_b_jets >= 1)) & (max_I_mini < 0.1 & Z_mass > 85 & Z_mass < 95)")
    with uproot.recreate(args.output) as f:
        f['inclusive'] = data
        f['zero_jet'] = data_zero_jet
        f['one_jet'] = data_one_jet
        f['zero_to_one_jet'] = data_zero_to_one_jet
        f['two_jet'] = data_two_jet
        f['VH_ttH'] = data_VH_ttH
        f['VH'] = data_VH
        f['ZH'] = data_ZH
        f['ttH_had'] = data_ttH_had
        f['ttH_lep'] = data_ttH_lep
    # data.to_root(args.output, key='inclusive', mode='a', index=False)
    # data_zero_jet.to_root(args.output, key='zero_jet', mode='a', index=False)
    # data_one_jet.to_root(args.output, key='one_jet', mode='a', index=False)
    # data_zero_to_one_jet.to_root(args.output, key='zero_to_one_jet', mode='a', index=False)
    # data_two_jet.to_root(args.output, key='two_jet', mode='a', index=False)
    # data_VH_ttH.to_root(args.output, key='VH_ttH', mode='a', index=False)

    # meta_data = pd.DataFrame({'initial_events': [initial_events], 'final_events': [final_events]})
    # meta_data.to_root(args.output, key='MetaData', mode='a', index=False)

    print("Finished!!! Have gotten the skimmed data in {:s}".format(args.output))

if __name__ == '__main__':
    main()
