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
from kinr2_weighter import create_kin_weighter as create_kinr2_weighter
from kinr3_weighter import create_kin_weighter as create_kinr3_weighter

warnings.simplefilter(action='ignore', category=FutureWarning)

# Initialize the DNN models globally
kin_weighter_r2 = None
kin_weighter_r3 = None
kin_weighter = None

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
    
    # MET and MHT delta phi with gamma
    delta_phi_vars['MET_deltaphi'] = true_delta_phi_vectorized(np.abs(data['MET_phi'] - data['gamma_phi']))
    delta_phi_vars['photon_mht_deltaphi'] = true_delta_phi_vectorized(np.abs(data['MHT_phi'] - data['gamma_phi']))
    
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

def compute_additional_delta_phi_vectorized(data):
    """Vectorized computation of additional delta phi variables"""
    additional_vars = {}
    
    # Additional lepton delta phi (if columns exist)
    if 'additional_lepton_1_phi' in data.columns:
        additional_vars['additional_lepton_1_deltaphi'] = true_delta_phi_vectorized(
            np.abs(data['additional_lepton_1_phi'] - data['gamma_phi'])
        )
    else:
        additional_vars['additional_lepton_1_deltaphi'] = np.full(len(data), -9999.0)
        
    if 'additional_lepton_2_phi' in data.columns:
        additional_vars['additional_lepton_2_deltaphi'] = true_delta_phi_vectorized(
            np.abs(data['additional_lepton_2_phi'] - data['gamma_phi'])
        )
    else:
        additional_vars['additional_lepton_2_deltaphi'] = np.full(len(data), -9999.0)
    
    return additional_vars

def compute_jet_delta_phi_variables_vectorized(data):
    """Vectorized computation of jet delta phi and deltaR variables"""
    jet_vars = {}
    
    # Jet delta phi with gamma
    for i in range(1, 5):
        jet_phi_col = f'jet_{i}_phi'
        jet_eta_col = f'jet_{i}_eta'
        
        if jet_phi_col in data.columns:
            jet_vars[f'jet_{i}_deltaphi'] = true_delta_phi_vectorized(
                np.abs(data[jet_phi_col] - data['gamma_phi'])
            )
            # Jet-gamma deltaR
            if jet_eta_col in data.columns:
                delta_eta = data[jet_eta_col] - data['gamma_eta']
                delta_phi = jet_vars[f'jet_{i}_deltaphi']
                jet_vars[f'jet{i}G_deltaR'] = np.sqrt(delta_eta**2 + delta_phi**2)
            else:
                jet_vars[f'jet{i}G_deltaR'] = np.full(len(data), -9999.0)
        else:
            jet_vars[f'jet_{i}_deltaphi'] = np.full(len(data), -9999.0)
            jet_vars[f'jet{i}G_deltaR'] = np.full(len(data), -9999.0)
    
    # Max and min jet deltaR
    jet_deltaR_cols = [f'jet{i}G_deltaR' for i in [1, 2] if f'jet{i}G_deltaR' in jet_vars]
    if len(jet_deltaR_cols) >= 2:
        jet_deltaR_array = np.column_stack([jet_vars[col] for col in jet_deltaR_cols])
        jet_vars['max_jet_deltaR'] = np.max(jet_deltaR_array, axis=1)
        jet_vars['min_jet_deltaR'] = np.min(jet_deltaR_array, axis=1)
    else:
        jet_vars['max_jet_deltaR'] = np.full(len(data), -9999.0)
        jet_vars['min_jet_deltaR'] = np.full(len(data), -9999.0)
    
    return jet_vars

def compute_l1g_l2g_deltaR_vectorized(data):
    """Vectorized computation of l1g and l2g deltaR"""
    # Calculate delta phi
    Z_lead_lepton_deltaphi = true_delta_phi_vectorized(
        np.abs(data['Z_lead_lepton_phi'] - data['gamma_phi'])
    )
    Z_sublead_lepton_deltaphi = true_delta_phi_vectorized(
        np.abs(data['Z_sublead_lepton_phi'] - data['gamma_phi'])
    )
    
    # Calculate deltaR for both leptons
    lead_deltaR = np.sqrt(
        (data['Z_lead_lepton_eta'] - data['gamma_eta'])**2 + Z_lead_lepton_deltaphi**2
    )
    sublead_deltaR = np.sqrt(
        (data['Z_sublead_lepton_eta'] - data['gamma_eta'])**2 + Z_sublead_lepton_deltaphi**2
    )
    
    # l1g_deltaR: larger deltaR, l2g_deltaR: smaller deltaR
    l1g_deltaR = np.where(lead_deltaR > sublead_deltaR, lead_deltaR, sublead_deltaR)
    l2g_deltaR = np.where(lead_deltaR < sublead_deltaR, lead_deltaR, sublead_deltaR)
    
    return {'l1g_deltaR': l1g_deltaR, 'l2g_deltaR': l2g_deltaR}

def compute_l1g_l2g_deltaR_refit_vectorized(data):
    """Vectorized computation of l1g and l2g deltaR for refit versions"""
    # Calculate delta phi (using original phi since no refit version)
    Z_lead_lepton_deltaphi = true_delta_phi_vectorized(
        np.abs(data['Z_lead_lepton_phi'] - data['gamma_phi'])
    )
    Z_sublead_lepton_deltaphi = true_delta_phi_vectorized(
        np.abs(data['Z_sublead_lepton_phi'] - data['gamma_phi'])
    )
    
    # Calculate deltaR for both leptons (using original eta since no refit version)
    lead_deltaR = np.sqrt(
        (data['Z_lead_lepton_eta'] - data['gamma_eta'])**2 + Z_lead_lepton_deltaphi**2
    )
    sublead_deltaR = np.sqrt(
        (data['Z_sublead_lepton_eta'] - data['gamma_eta'])**2 + Z_sublead_lepton_deltaphi**2
    )
    
    # l1g_deltaR_refit: larger deltaR, l2g_deltaR_refit: smaller deltaR
    l1g_deltaR_refit = np.where(lead_deltaR > sublead_deltaR, lead_deltaR, sublead_deltaR)
    l2g_deltaR_refit = np.where(lead_deltaR < sublead_deltaR, lead_deltaR, sublead_deltaR)
    
    return {'l1g_deltaR_refit': l1g_deltaR_refit, 'l2g_deltaR_refit': l2g_deltaR_refit}

def compute_jet_pair_variables_vectorized(data, vectors):
    """Vectorized computation of jet pair variables"""
    jet_vars = {}
    n_events = len(data)
    
    if 'jet_1' not in vectors or 'jet_2' not in vectors:
        for var_name in ['delta_eta_jj', 'delta_phi_jj', 'mass_jj', 'jet_pair_pt']:
            jet_vars[var_name] = np.full(n_events, -9999.0)
        return jet_vars
    
    jet1_vecs = vectors['jet_1']
    jet2_vecs = vectors['jet_2']
    
    # Delta eta and delta phi between jets
    delta_eta_jj = []
    delta_phi_jj = []
    mass_jj = []
    jet_pair_pt = []
    
    mask_2jets = data['n_jets'] >= 2
    
    for i in range(n_events):
        if mask_2jets.iloc[i]:
            delta_eta_jj.append(abs(jet1_vecs[i].Eta() - jet2_vecs[i].Eta()))
            delta_phi_jj.append(true_delta_phi(abs(jet1_vecs[i].Phi() - jet2_vecs[i].Phi())))
            dijet = jet1_vecs[i] + jet2_vecs[i]
            mass_jj.append(dijet.M())
            jet_pair_pt.append(dijet.Pt())
        else:
            delta_eta_jj.append(-9999)
            delta_phi_jj.append(-9999)
            mass_jj.append(-9999)
            jet_pair_pt.append(-9999)
    
    jet_vars['delta_eta_jj'] = np.array(delta_eta_jj)
    jet_vars['delta_phi_jj'] = np.array(delta_phi_jj)
    jet_vars['mass_jj'] = np.array(mass_jj)
    jet_vars['jet_pair_pt'] = np.array(jet_pair_pt)
    
    return jet_vars

def compute_jet_btag_variables_vectorized(data):
    """Vectorized computation of jet btag variables"""
    jet_vars = {}
    n_events = len(data)
    
    # Max two jet btag
    max_btag_values = []
    jet_ptt_values = []
    
    mask_2jets = data['n_jets'] >= 2
    
    for i in range(n_events):
        if mask_2jets.iloc[i]:
            btag1 = data['jet_1_btagDeepFlavB'].iloc[i] if 'jet_1_btagDeepFlavB' in data.columns else 0
            btag2 = data['jet_2_btagDeepFlavB'].iloc[i] if 'jet_2_btagDeepFlavB' in data.columns else 0
            max_btag_values.append(max(btag1, btag2))
            
            # jet_ptt calculation
            jet1_pt = data['jet_1_pt'].iloc[i] if 'jet_1_pt' in data.columns else 0
            jet2_pt = data['jet_2_pt'].iloc[i] if 'jet_2_pt' in data.columns else 0
            jet1_phi = data['jet_1_phi'].iloc[i] if 'jet_1_phi' in data.columns else 0
            jet2_phi = data['jet_2_phi'].iloc[i] if 'jet_2_phi' in data.columns else 0
            
            if jet1_pt > 0 and jet2_pt > 0:
                delta_phi = true_delta_phi(abs(jet1_phi - jet2_phi))
                jet_ptt_values.append(abs(jet1_pt * np.sin(delta_phi)))
            else:
                jet_ptt_values.append(-9999)
        else:
            max_btag_values.append(-9999)
            jet_ptt_values.append(-9999)
    
    jet_vars['max_two_jet_btag'] = np.array(max_btag_values)
    jet_vars['jet_ptt'] = np.array(jet_ptt_values)
    
    return jet_vars

def compute_gamma_ptRelErr_vectorized(data):
    """Vectorized computation of gamma_ptRelErr"""
    return data['gamma_energyErr'] / data['gamma_pt']

def compute_zeppenfeld_variables_vectorized(data):
    """Vectorized computation of zeppenfeld variables"""
    zep_vars = {}
    n_events = len(data)
    
    # Photon and H zeppenfeld
    photon_zep = []
    h_zep = []
    h_zep_refit = []
    
    mask_2jets = data['n_jets'] >= 2
    
    for i in range(n_events):
        if mask_2jets.iloc[i]:
            jet1_eta = data['jet_1_eta'].iloc[i] if 'jet_1_eta' in data.columns else 0
            jet2_eta = data['jet_2_eta'].iloc[i] if 'jet_2_eta' in data.columns else 0
            avg_jet_eta = (jet1_eta + jet2_eta) / 2.0;
            
            photon_zep.append(data['gamma_eta'].iloc[i] - avg_jet_eta)
            h_zep.append(data['H_eta'].iloc[i] - avg_jet_eta)
            h_zep_refit.append(data['H_eta_refit'].iloc[i] - avg_jet_eta)
        else:
            photon_zep.append(-9999)
            h_zep.append(-9999)
            h_zep_refit.append(-9999)
    
    zep_vars['photon_zeppenfeld'] = np.array(photon_zep)
    zep_vars['H_zeppenfeld'] = np.array(h_zep)
    zep_vars['H_zeppenfeld_refit'] = np.array(h_zep_refit)
    
    return zep_vars

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

def compute_is_center(x):
    if (x.H_mass >= 120 and x.H_mass <= 130): return 1
    else: return 0

def compute_is_center_refit(x):
    if (x.H_mass_refit >= 120 and x.H_mass_refit <= 130): return 1
    else: return 0

def compute_photon_mht_deltaphi(x):
    """Compute delta phi between photon and MHT"""
    return true_delta_phi(abs(x.MHT_phi - x.gamma_phi))

def compute_photon_mht_deltaphi_refit(x):
    """Compute delta phi between photon and MHT for refit version"""
    return true_delta_phi(abs(x.MHT_phi - x.gamma_phi))  # MHT_phi doesn't have refit version

def compute_kin_weight(x):
    """Compute kinematic weight using DNN model"""
    # Input order: gamma_mvaID, photon_mht_deltaphi, n_jets, Z_pt, gamma_pt, H_pt, MHT_pt, HT
    try:
        inputs = [
            x.gamma_mvaID,
            x.photon_mht_deltaphi,
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
    except Exception as e:
        print(f"Error in compute_kin_weight: {e}")
        return 1.0  # Default weight if computation fails

def compute_kin_weight_refit(x):
    """Compute kinematic weight using DNN model with refit variables"""
    # Input order: gamma_mvaID, photon_mht_deltaphi, n_jets, Z_pt_refit, gamma_pt, H_pt_refit, MHT_pt, HT
    try:
        inputs = [
            x.gamma_mvaID,
            x.photon_mht_deltaphi,
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
    except Exception as e:
        print(f"Error in compute_kin_weight_refit: {e}")
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

def preselect(data):

    #data.query('(Muons_Minv_MuMu_Paper >= 110) | (Event_Paper_Category >= 17)', inplace=True)
    #data.query('Event_Paper_Category > 0', inplace=True)

    return data

def decorate(data):
    """Main decoration function with optimized vectorized operations"""
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
    if 'weight_central' in data.columns:
        data['weight'] = data.weight_central
    else:
        print("Warning: weight_central column not found, using default weight of 1.0")
        data['weight'] = np.ones(len(data))
    
    # Use vectorized gamma_ptRelErr computation
    data['gamma_ptRelErr'] = compute_gamma_ptRelErr_vectorized(data)
    
    # Add vectorized additional delta phi variables
    print("Computing additional delta phi variables...")
    additional_vars = compute_additional_delta_phi_vectorized(data)
    for key, value in additional_vars.items():
        data[key] = value
    
    # Add vectorized jet delta phi and deltaR variables
    print("Computing jet delta phi and deltaR variables...")
    jet_delta_vars = compute_jet_delta_phi_variables_vectorized(data)
    for key, value in jet_delta_vars.items():
        data[key] = value
    
    # Add vectorized l1g and l2g deltaR variables
    print("Computing l1g and l2g deltaR variables...")
    l1g_l2g_vars = compute_l1g_l2g_deltaR_vectorized(data)
    for key, value in l1g_l2g_vars.items():
        data[key] = value
    
    l1g_l2g_refit_vars = compute_l1g_l2g_deltaR_refit_vectorized(data)
    for key, value in l1g_l2g_refit_vars.items():
        data[key] = value
    
    # Add vectorized jet pair variables
    print("Computing jet pair variables...")
    jet_pair_vars = compute_jet_pair_variables_vectorized(data, vectors_orig)
    for key, value in jet_pair_vars.items():
        data[key] = value
    
    # Add vectorized jet btag variables
    print("Computing jet btag variables...")
    jet_btag_vars = compute_jet_btag_variables_vectorized(data)
    for key, value in jet_btag_vars.items():
        data[key] = value
    
    # Add vectorized zeppenfeld variables
    print("Computing zeppenfeld variables...")
    zep_vars = compute_zeppenfeld_variables_vectorized(data)
    for key, value in zep_vars.items():
        data[key] = value
    
    # Complex angular calculations that need apply (most complex physics calculations)
    print("Computing complex angular variables...")
    data['l_rapCM'] = data.apply(lambda x:compute_l_rapCM(x), axis=1)
    data['l_cosProdAngle'] = data.apply(lambda x:compute_l_prodAngle(x), axis=1)
    data['Z_cosProdAngle'] = data.apply(lambda x:compute_Z_prodAngle(x), axis=1)
    data['Z_cos_theta'] = data.apply(lambda x:compute_Z_cosTheta(x), axis=1)
    data['lep_cos_theta'] = data.apply(lambda x: compute_l_costheta(x), axis=1)
    data['lep_phi'] = data.apply(lambda x: compute_l_phi(x), axis=1)
    
    # Refit versions of angular calculations
    data['l_rapCM_refit'] = data.apply(lambda x:compute_l_rapCM_refit(x), axis=1)
    data['l_cosProdAngle_refit'] = data.apply(lambda x:compute_l_prodAngle_refit(x), axis=1)
    data['Z_cosProdAngle_refit'] = data.apply(lambda x:compute_Z_prodAngle_refit(x), axis=1)
    data['Z_cos_theta_refit'] = data.apply(lambda x:compute_Z_cosTheta_refit(x), axis=1)
    data['lep_cos_theta_refit'] = data.apply(lambda x: compute_l_costheta_refit(x), axis=1)
    data['lep_phi_refit'] = data.apply(lambda x: compute_l_phi_refit(x), axis=1)
    
    # Kinematic DNN weighting 
    print("Computing kinematic weighting...")
    data['kin_weight'] = data.apply(lambda x: compute_kin_weight(x), axis=1)
    data['kin_weight_refit'] = data.apply(lambda x: compute_kin_weight_refit(x), axis=1)

    data['pythia_weight'] = np.ones(data.shape[0]) / 0.96934
    
    # Type conversions
    data = data.astype(float)
    data = data.astype({'is_center': int, 'is_center_refit': int, 'Z_lead_lepton_charge': int, 'Z_lead_lepton_id': int, 'Z_sublead_lepton_charge': int, 'Z_sublead_lepton_id': int, "n_jets": int, "n_b_jets": int, "n_leptons": int, "n_electrons": int, "n_muons": int, 'event': int})

    return data
    

def main():
    
    args = getArgs()

    # Initialize weighters based on year from input file name
    global kin_weighter, kin_weighter_r2, kin_weighter_r3
    if any(year in args.input for year in ["2016", "2017", "2018"]):
        if kin_weighter_r2 is None:
            kin_weighter_r2 = create_kinr2_weighter()
        kin_weighter = kin_weighter_r2
        print("Using kinr2_weighter for 201X data.")
    else: # Assume 202X for others
        if kin_weighter_r3 is None:
            kin_weighter_r3 = create_kinr3_weighter()
        kin_weighter = kin_weighter_r3
        print("Using kinr3_weighter for 202X data.")

    variables = [
        'H_pt', 'H_eta', 'H_phi', 'H_mass',
        'Z_pt', 'Z_eta', 'Z_phi', 'Z_mass',
        'Z_lead_lepton_pt', 'Z_lead_lepton_eta', 'Z_lead_lepton_phi', 'Z_lead_lepton_mass',
        'Z_lead_lepton_charge', 'Z_lead_lepton_id', 'Z_lead_lepton_ptE_error',
        'Z_sublead_lepton_pt', 'Z_sublead_lepton_eta', 'Z_sublead_lepton_phi', 'Z_sublead_lepton_mass',
        'Z_sublead_lepton_charge', 'Z_sublead_lepton_id', 'Z_sublead_lepton_ptE_error',
        # Add refit variables - these will be created by apply_z_refit
        'H_pt_refit', 'H_eta_refit', 'H_phi_refit', 'H_mass_refit',
        'Z_pt_refit', 'Z_eta_refit', 'Z_phi_refit', 'Z_mass_refit',
        'Z_lead_lepton_pt_refit', 'Z_sublead_lepton_pt_refit',
        'Z_lead_lepton_pt_refit_err', 'Z_sublead_lepton_pt_refit_err',
        'gamma_pt', 'gamma_eta', 'gamma_phi', 'gamma_mass',
        'gamma_mvaID', 'gamma_energyErr',
        "gamma_fsr_pt","gamma_fsr_eta","gamma_fsr_phi","gamma_fsr_mass",  # Removed gamma_fsr_clean
        'jet_1_pt', 'jet_1_eta', 'jet_1_phi', 'jet_1_mass', 'jet_1_btagDeepFlavB',
        'jet_2_pt', 'jet_2_eta', 'jet_2_phi', 'jet_2_mass', 'jet_2_btagDeepFlavB',
        'n_jets', 'n_leptons', 'n_electrons', 'n_muons', 'n_b_jets',  # Removed n_iso_photons
        'MET_pt', 'MET_phi', 
        'MHT_pt', 'MHT_phi', 'HT',  # Added for kinematic DNN
        'weight_central',
        'event'
    ]

    print("Now!!! Porcessing {:s} ......".format(args.input))

    if os.path.isfile(args.output): os.remove(args.output)

    initial_events = 0
    final_events = 0

    data = pd.read_parquet(args.input)
    # print(f"Loaded data with shape: {data.shape}")
    # print(f"Available columns: {list(data.columns)}")
    
    # Remove columns with object dtype
    data = data.drop(columns=data.select_dtypes(include=['object']).columns)
    initial_events += data.shape[0]
    # data = preselect(data) #TODO add cutflow
    # Apply Z constraint refit
    data = apply_z_refit(data)
    data = decorate(data)

    weight_list = [
        "weight_hlt_sf_central", 
        "weight_pu_reweight_sf_central", 
        "weight_btag_deepjet_wp_sf_SelectedJet_central",
        "weight_electron_iso_sf_SelectedElectron_central",
        "weight_electron_reco_sf_SelectedElectron_central",
        "weight_electron_wplid_sf_SelectedElectron_central",
        "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central",
        "weight_muon_iso_sf_SelectedMuon_central",
        "weight_muon_looseid_sf_SelectedMuon_central",
        "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central",
        "weight_muon_reco_sf_SelectedMuon_central",
        "weight_photon_csev_sf_SelectedPhoton_central",
        "weight_photon_id_sf_SelectedPhoton_central",
        "weight_photon_id_shape_sf_SelectedPhoton_central"
    ]
    ggH_weight_list = [
        "weight_nnlo_sf_GenHzgHiggs_central"
    ]
    signal_weight_list = [
        "pythia_weight"
    ]
    background_weight_list = [
        "weight_photon_fake_photon_sf_SelectedPhoton_central",
        "kin_weight"
    ]

    # Check if required weight columns exist
    missing_weights = []
    for weight in ['weight'] + weight_list + ggH_weight_list + signal_weight_list + background_weight_list:
        if weight not in data.columns:
            missing_weights.append(weight)
    
    if missing_weights:
        print(f"Warning: Missing weight columns: {missing_weights}")
        # Create default weight if weight column is missing
        if 'weight' not in data.columns:
            data['weight'] = np.ones(len(data))
        # Create default values for missing weight columns
        for weight in missing_weights:
            if weight != 'weight':  # weight is already handled above
                data[weight] = np.ones(len(data))
    
    data["weight_corr"] = data["weight_central"]
    for weight in weight_list:
        if weight in data.columns:
            data["weight_corr"] *= data[weight]
    if "ggH" in args.input:
        for weight in ggH_weight_list:
            if weight in data.columns:
                data["weight_corr"] *= data[weight]
    if "signal" in args.input:
        for weight in signal_weight_list:
            if weight in data.columns:
                data["weight_corr"] *= data[weight]
    if "background" in args.input:
        for weight in background_weight_list:
            if weight in data.columns:
                data["weight_corr"] *= data[weight]

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

    print("Finished!!! Have gotten the skimmed data in {:s}".format(args.output))

if __name__ == '__main__':
    main()
