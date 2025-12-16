import awkward
import numpy
from correctionlib import _core
from higgs_dna.utils import misc_utils

# Base path placeholder
base_path = "higgs_dna/systematics/data/"

# -----------------------------------------------------------------------------------
# File Definitions
# -----------------------------------------------------------------------------------
SingleElectron_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_eltrig27_2016APV_efficiencies.json",
    "2016postVFP":  f"{base_path}2016postVFP_UL/hzg_eltrig27_2016_efficiencies.json",
    "2017":         f"{base_path}2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2018":         f"{base_path}2018_UL/hzg_eltrig32_2018_efficiencies.json",
    "2022preEE":    f"{base_path}2022preEE_UL/hzg_eltrig30_2022_efficiencies.json",
    "2022postEE":   f"{base_path}2022postEE_UL/hzg_eltrig30_2022EE_efficiencies.json",
    "2023preBPix":  f"{base_path}2023preBPix_UL/hzg_eltrig30_2023_efficiencies.json",
    "2023postBPix": f"{base_path}2023postBPix_UL/hzg_eltrig30_2023BPix_efficiencies.json",
}
DoubleElectron_HighLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_eltrig23_2016APV_efficiencies.json",
    "2016postVFP":  f"{base_path}2016postVFP_UL/hzg_eltrig23_2016_efficiencies.json",
    "2017":         f"{base_path}2017_UL/hzg_eltrig23_2017_efficiencies.json",
    "2018":         f"{base_path}2018_UL/hzg_eltrig23_2018_efficiencies.json",
    "2022preEE":    f"{base_path}2022preEE_UL/hzg_eltrig23_2022_efficiencies.json",
    "2022postEE":   f"{base_path}2022postEE_UL/hzg_eltrig23_2022EE_efficiencies.json",
    "2023preBPix":  f"{base_path}2023preBPix_UL/hzg_eltrig23_2023_efficiencies.json",
    "2023postBPix": f"{base_path}2023postBPix_UL/hzg_eltrig23_2023BPix_efficiencies.json",
}
DoubleElectron_LowLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_eltrig12_2016APV_efficiencies.json",
    "2016postVFP":  f"{base_path}2016postVFP_UL/hzg_eltrig12_2016_efficiencies.json",
    "2017":         f"{base_path}2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2018":         f"{base_path}2018_UL/hzg_eltrig12_2018_efficiencies.json",
    "2022preEE":    f"{base_path}2022preEE_UL/hzg_eltrig12_2022_efficiencies.json",
    "2022postEE":   f"{base_path}2022postEE_UL/hzg_eltrig12_2022EE_efficiencies.json",
    "2023preBPix":  f"{base_path}2023preBPix_UL/hzg_eltrig12_2023_efficiencies.json",
    "2023postBPix": f"{base_path}2023postBPix_UL/hzg_eltrig12_2023BPix_efficiencies.json",
}

# Hole files for 2023postBPix
SingleElectron_HLT_FILE_HOLE = { "2023postBPix" : f"{base_path}2023postBPix_UL/hzg_eltrig30_2023BPixHole_efficiencies.json" }
DoubleElectron_LowLeg_HLT_FILE_HOLE = { "2023postBPix" : f"{base_path}2023postBPix_UL/hzg_eltrig12_2023BPixHole_efficiencies.json" }
DoubleElectron_HighLeg_HLT_FILE_HOLE = { "2023postBPix" : f"{base_path}2023postBPix_UL/hzg_eltrig23_2023BPixHole_efficiencies.json" }

# Muon HLT files
SingleMuon_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_mutrig24_2016APV_efficiencies.json",
    "2016postVFP":  f"{base_path}2016postVFP_UL/hzg_mutrig24_2016_efficiencies.json",
    "2017":         f"{base_path}2017_UL/hzg_mutrig27_2017_efficiencies.json",
    "2018":         f"{base_path}2018_UL/hzg_mutrig24_2018_efficiencies.json",
    "2022preEE":    f"{base_path}2022preEE_UL/hzg_mutrig24_2022_efficiencies.json",
    "2022postEE":   f"{base_path}2022postEE_UL/hzg_mutrig24_2022EE_efficiencies.json",
    "2023preBPix":  f"{base_path}2023preBPix_UL/hzg_mutrig24_2023_efficiencies.json",
    "2023postBPix": f"{base_path}2023postBPix_UL/hzg_mutrig24_2023BPix_efficiencies.json",
}
DoubleMuon_HighLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_mutrig17_2016APV_efficiencies.json",
    "2016postVFP":  f"{base_path}2016postVFP_UL/hzg_mutrig17_2016_efficiencies.json",
    "2017":         f"{base_path}2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2018":         f"{base_path}2018_UL/hzg_mutrig17_2018_efficiencies.json",
    "2022preEE":    f"{base_path}2022preEE_UL/hzg_mutrig17_2022_efficiencies.json",
    "2022postEE":   f"{base_path}2022postEE_UL/hzg_mutrig17_2022EE_efficiencies.json",
    "2023preBPix":  f"{base_path}2023preBPix_UL/hzg_mutrig17_2023_efficiencies.json",
    "2023postBPix": f"{base_path}2023postBPix_UL/hzg_mutrig17_2023BPix_efficiencies.json",
}
DoubleMuon_LowLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_mutrig8_2016APV_efficiencies.json",
    "2016postVFP":  f"{base_path}2016postVFP_UL/hzg_mutrig8_2016_efficiencies.json",
    "2017":         f"{base_path}2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2018":         f"{base_path}2018_UL/hzg_mutrig8_2018_efficiencies.json",
    "2022preEE":    f"{base_path}2022preEE_UL/hzg_mutrig8_2022_efficiencies.json",
    "2022postEE":   f"{base_path}2022postEE_UL/hzg_mutrig8_2022EE_efficiencies.json",
    "2023preBPix":  f"{base_path}2023preBPix_UL/hzg_mutrig8_2023_efficiencies.json",
    "2023postBPix": f"{base_path}2023postBPix_UL/hzg_mutrig8_2023BPix_efficiencies.json",
}


def get_lepton_efficiencies(leptons, single_eval, high_eval, low_eval, 
                            is_electron, is_data, year, sigma=0.0,
                            single_eval_hole=None, high_eval_hole=None, low_eval_hole=None):
    """
    Retrieves efficiencies from correctionlib for a specific shift (sigma).
    Returns Awkward arrays of shape [N_events, N_leptons].
    """
    leptons_flat = awkward.flatten(leptons)
    lep_pt = awkward.to_numpy(leptons_flat.pt)
    lep_eta = awkward.to_numpy(leptons_flat.eta)
    if not is_electron:
        lep_eta = numpy.abs(lep_eta)

    eff_type = "effdata" if is_data else "effmc"
    syst_type = "systdata" if is_data else "systmc"

    # 1. Get Nominal Efficiencies
    eff_s = single_eval[eff_type].evalv(lep_pt, lep_eta)
    eff_h = high_eval[eff_type].evalv(lep_pt, lep_eta)
    eff_l = low_eval[eff_type].evalv(lep_pt, lep_eta)

    # 2. Apply Systematics (if sigma != 0)
    if sigma != 0.0:
        syst_s = single_eval[syst_type].evalv(lep_pt, lep_eta)
        syst_h = high_eval[syst_type].evalv(lep_pt, lep_eta)
        syst_l = low_eval[syst_type].evalv(lep_pt, lep_eta)
        
        eff_s = eff_s + (sigma * syst_s)
        eff_h = eff_h + (sigma * syst_h)
        eff_l = eff_l + (sigma * syst_l)

    # 3. Handle 2023postBPix Electron Holes
    if is_electron and year == "2023postBPix" and single_eval_hole is not None:
        lep_phi = awkward.to_numpy(leptons_flat.phi)
        hole_mask = (lep_eta > -1.566) & (lep_eta < 0) & (lep_phi > -1.2) & (lep_phi < -0.8)
        
        eff_s_hole = single_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_h_hole = high_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_l_hole = low_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        
        if sigma != 0.0:
            syst_s_hole = single_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_h_hole = high_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_l_hole = low_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            
            eff_s_hole = eff_s_hole + (sigma * syst_s_hole)
            eff_h_hole = eff_h_hole + (sigma * syst_h_hole)
            eff_l_hole = eff_l_hole + (sigma * syst_l_hole)

        eff_s = numpy.where(hole_mask, eff_s_hole, eff_s)
        eff_h = numpy.where(hole_mask, eff_h_hole, eff_h)
        eff_l = numpy.where(hole_mask, eff_l_hole, eff_l)

    # 4. Enforce Physical Hierarchy (Clipping)
    # Ensure 0 <= Prob(Single) <= Prob(High) <= Prob(Low) <= 1
    eff_s = numpy.clip(eff_s, 0.0, 1.0)
    eff_h = numpy.clip(eff_h, 0.0, 1.0)
    eff_l = numpy.clip(eff_l, 0.0, 1.0)
    
    eff_h = numpy.maximum(eff_h, eff_s)
    eff_l = numpy.maximum(eff_l, eff_h)

    # Restructure back to awkward arrays
    eff_s = awkward.unflatten(eff_s, awkward.num(leptons))
    eff_h = awkward.unflatten(eff_h, awkward.num(leptons))
    eff_l = awkward.unflatten(eff_l, awkward.num(leptons))

    return eff_s, eff_h, eff_l


def calculate_exclusive_probability(eff_s, eff_h, eff_l):
    """
    Calculates the total event probability using the Exclusive Probability Method.
    Supports an arbitrary number of leptons by dynamic padding.
    
    P_Event = 1 - (Fail_A + Fail_B)
    """
    
    # 1. Define Per-Lepton Probabilities for Exclusive States
    # State S (Single): eff_s
    # State H (High only): eff_h - eff_s
    # State L (Low only): eff_l - eff_h
    # State F (Fail Low): 1 - eff_l
    
    p_fail_high = 1.0 - eff_h  # Fails High (and Single)
    p_excl_H    = eff_h - eff_s # Passes High, Fails Single
    p_F         = 1.0 - eff_l  # Fails Low (and High and Single)

    # --- Scenario A: All leptons fail High Leg ---
    # P(A) = Product_i (1 - epsilon_High_i)
    # awkard.prod supports variable length arrays directly
    term_A = awkward.prod(p_fail_high, axis=1)

    # --- Scenario B: Exactly one exclusive High, others Fail Low ---
    # P(B) = Sum_i [ p_excl_H_i * Product_{j!=i} (p_F_j) ]
    
    # To handle arbitrary N vectorized, we pad to the max multiplicity in the batch
    # and loop over the columns (which is fast as max_N is small, e.g., 2, 3, 4).
    max_n = awkward.max(awkward.num(eff_s))
    if max_n is None or max_n == 0:
        # Handle case with empty event batch or no leptons
        return awkward.zeros_like(term_A)

    # Pad arrays: 
    # - p_excl_H padded with 0.0 (sum identity)
    # - p_F padded with 1.0 (product identity)
    def to_padded_numpy(arr, target_N, pad_val):
        padded = awkward.fill_none(awkward.pad_none(arr, target_N, axis=1), pad_val)
        return awkward.to_numpy(padded)

    ph_np = to_padded_numpy(p_excl_H, max_n, 0.0)
    pf_np = to_padded_numpy(p_F, max_n, 1.0)
    
    # Calculate term B by iterating columns to avoid 0/0 division issues
    term_B = numpy.zeros_like(ph_np[:, 0])
    
    # Create a reusable mask of all Trues
    # We will compute Product_{j!=i} using numpy.prod with a 'where' mask
    all_indices = numpy.ones((ph_np.shape[0], max_n), dtype=bool)
    
    for i in range(max_n):
        # Current lepton 'i' is the one Passing High (p_H)
        # All others must Fail Low (p_F)
        
        # Mask out column i, keep others
        all_indices[:, i] = False
        
        # Product of p_F for all other columns
        # where=False elements are replaced by 1.0 (identity) in prod
        prod_others = numpy.prod(pf_np, axis=1, where=all_indices)
        
        term_B += ph_np[:, i] * prod_others
        
        # Reset mask for next iteration
        all_indices[:, i] = True

    # Total Fail Probability
    p_fail_event = term_A + term_B
    
    # Total Pass Probability
    p_pass_event = numpy.clip(1.0 - p_fail_event, 0.0, 1.0)
    
    return awkward.Array(p_pass_event)


def get_combined_flavor_sf(events, year, central_only):
    """
    Main driver to calculate Scale Factors.
    """
    # 1. Load Evaluators
    ele_single_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleElectron_HLT_FILE[year]))
    ele_high_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_HighLeg_HLT_FILE[year]))
    ele_low_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_LowLeg_HLT_FILE[year]))
    
    ele_hole_evals = [None, None, None]
    if year == "2023postBPix":
        ele_hole_evals = [
            _core.CorrectionSet.from_file(misc_utils.expand_path(SingleElectron_HLT_FILE_HOLE[year])),
            _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_HighLeg_HLT_FILE_HOLE[year])),
            _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_LowLeg_HLT_FILE_HOLE[year]))
        ]

    muon_single_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleMuon_HLT_FILE[year]))
    muon_high_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_HighLeg_HLT_FILE[year]))
    muon_low_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_LowLeg_HLT_FILE[year]))

    electrons = events.Electron
    muons = events.Muon
    
    # Sort by pT
    electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
    muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]

    # 2. Define Variations to Calculate
    variations_to_run = [("central", 0.0)]
    if not central_only:
        variations_to_run += [("up", 1.0), ("down", -1.0)]

    results = {}

    for var_name, sigma in variations_to_run:
        # --- Electron Probability ---
        e_s, e_h, e_l = get_lepton_efficiencies(
            electrons, ele_single_eval, ele_high_eval, ele_low_eval,
            is_electron=True, is_data=True, year=year, sigma=sigma,
            single_eval_hole=ele_hole_evals[0], high_eval_hole=ele_hole_evals[1], low_eval_hole=ele_hole_evals[2]
        )
        p_ele_data = calculate_exclusive_probability(e_s, e_h, e_l)

        e_s_mc, e_h_mc, e_l_mc = get_lepton_efficiencies(
            electrons, ele_single_eval, ele_high_eval, ele_low_eval,
            is_electron=True, is_data=False, year=year, sigma=sigma,
            single_eval_hole=ele_hole_evals[0], high_eval_hole=ele_hole_evals[1], low_eval_hole=ele_hole_evals[2]
        )
        p_ele_mc = calculate_exclusive_probability(e_s_mc, e_h_mc, e_l_mc)

        # --- Muon Probability ---
        m_s, m_h, m_l = get_lepton_efficiencies(
            muons, muon_single_eval, muon_high_eval, muon_low_eval,
            is_electron=False, is_data=True, year=year, sigma=sigma
        )
        p_muon_data = calculate_exclusive_probability(m_s, m_h, m_l)

        m_s_mc, m_h_mc, m_l_mc = get_lepton_efficiencies(
            muons, muon_single_eval, muon_high_eval, muon_low_eval,
            is_electron=False, is_data=False, year=year, sigma=sigma
        )
        p_muon_mc = calculate_exclusive_probability(m_s_mc, m_h_mc, m_l_mc)

        # --- Combine Flavors (Union of Independent Flavors) ---
        # P_total = 1 - (1 - P_ele)(1 - P_muon)
        p_total_data = 1.0 - (1.0 - p_ele_data) * (1.0 - p_muon_data)
        p_total_mc   = 1.0 - (1.0 - p_ele_mc) * (1.0 - p_muon_mc)

        # --- Calculate Scale Factor ---
        # Safe division: convert to numpy for ufunc compatibility
        p_total_data_np = awkward.to_numpy(p_total_data)
        p_total_mc_np = awkward.to_numpy(p_total_mc)
        
        sf = numpy.divide(p_total_data_np, p_total_mc_np, out=numpy.ones_like(p_total_data_np), where=p_total_mc_np > 1e-6)
        
        # Clip reasonable range
        results[var_name] = awkward.Array(numpy.clip(sf, 0.0, 10.0))

    return results

def HLT_sf(events, year, central_only):
    """
    Entry point function.
    """
    return get_combined_flavor_sf(events, year, central_only)