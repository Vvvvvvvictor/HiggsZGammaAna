import awkward
import numpy
from correctionlib import _core

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.utils import systematic_from_bins
import logging
logger = logging.getLogger(__name__)

# Note: these are placeholders. You should update them with the actual paths to your correction files.
base_path = "higgs_dna/systematics/data/"

# Electron HLT files
SingleElectron_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_eltrig27_2016APV_efficiencies.json",
    "2016postVFP":      f"{base_path}2016postVFP_UL/hzg_eltrig27_2016_efficiencies.json",
    "2017":      f"{base_path}2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2018":      f"{base_path}2018_UL/hzg_eltrig32_2018_efficiencies.json",
    "2022preEE":      f"{base_path}2022preEE_UL/hzg_eltrig30_2022_efficiencies.json",
    "2022postEE":    f"{base_path}2022postEE_UL/hzg_eltrig30_2022EE_efficiencies.json",
    "2023preBPix":      f"{base_path}2023preBPix_UL/hzg_eltrig30_2023_efficiencies.json",
    "2023postBPix":  f"{base_path}2023postBPix_UL/hzg_eltrig30_2023BPix_efficiencies.json",
}
DoubleElectron_HighLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_eltrig23_2016APV_efficiencies.json",
    "2016postVFP":      f"{base_path}2016postVFP_UL/hzg_eltrig23_2016_efficiencies.json",
    "2017":      f"{base_path}2017_UL/hzg_eltrig23_2017_efficiencies.json",
    "2018":      f"{base_path}2018_UL/hzg_eltrig23_2018_efficiencies.json",
    "2022preEE":      f"{base_path}2022preEE_UL/hzg_eltrig23_2022_efficiencies.json",
    "2022postEE":    f"{base_path}2022postEE_UL/hzg_eltrig23_2022EE_efficiencies.json",
    "2023preBPix":      f"{base_path}2023preBPix_UL/hzg_eltrig23_2023_efficiencies.json",
    "2023postBPix":  f"{base_path}2023postBPix_UL/hzg_eltrig23_2023BPix_efficiencies.json",
}
DoubleElectron_LowLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_eltrig12_2016APV_efficiencies.json",
    "2016postVFP":      f"{base_path}2016postVFP_UL/hzg_eltrig12_2016_efficiencies.json",
    "2017":      f"{base_path}2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2018":      f"{base_path}2018_UL/hzg_eltrig12_2018_efficiencies.json",
    "2022preEE":      f"{base_path}2022preEE_UL/hzg_eltrig12_2022_efficiencies.json",
    "2022postEE":    f"{base_path}2022postEE_UL/hzg_eltrig12_2022EE_efficiencies.json",
    "2023preBPix":      f"{base_path}2023preBPix_UL/hzg_eltrig12_2023_efficiencies.json",
    "2023postBPix":  f"{base_path}2023postBPix_UL/hzg_eltrig12_2023BPix_efficiencies.json",
}

# Hole files for 2023postBPix
SingleElectron_HLT_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig30_2023BPixHole_efficiencies.json"
}
DoubleElectron_LowLeg_HLT_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig12_2023BPixHole_efficiencies.json"
}
DoubleElectron_HighLeg_HLT_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig23_2023BPixHole_efficiencies.json"
}


# Muon HLT files
SingleMuon_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_mutrig24_2016APV_efficiencies.json",
    "2016postVFP":      f"{base_path}2016postVFP_UL/hzg_mutrig24_2016_efficiencies.json",
    "2017":      f"{base_path}2017_UL/hzg_mutrig27_2017_efficiencies.json",
    "2018":      f"{base_path}2018_UL/hzg_mutrig24_2018_efficiencies.json",
    "2022preEE":      f"{base_path}2022preEE_UL/hzg_mutrig24_2022_efficiencies.json",
    "2022postEE":    f"{base_path}2022postEE_UL/hzg_mutrig24_2022EE_efficiencies.json",
    "2023preBPix":      f"{base_path}2023preBPix_UL/hzg_mutrig24_2023_efficiencies.json",
    "2023postBPix":  f"{base_path}2023postBPix_UL/hzg_mutrig24_2023BPix_efficiencies.json",
}
DoubleMuon_HighLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_mutrig17_2016APV_efficiencies.json",
    "2016postVFP":      f"{base_path}2016postVFP_UL/hzg_mutrig17_2016_efficiencies.json",
    "2017":      f"{base_path}2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2018":      f"{base_path}2018_UL/hzg_mutrig17_2018_efficiencies.json",
    "2022preEE":      f"{base_path}2022preEE_UL/hzg_mutrig17_2022_efficiencies.json",
    "2022postEE":    f"{base_path}2022postEE_UL/hzg_mutrig17_2022EE_efficiencies.json",
    "2023preBPix":      f"{base_path}2023preBPix_UL/hzg_mutrig17_2023_efficiencies.json",
    "2023postBPix":  f"{base_path}2023postBPix_UL/hzg_mutrig17_2023BPix_efficiencies.json",
}
DoubleMuon_LowLeg_HLT_FILE = {
    "2016preVFP":   f"{base_path}2016preVFP_UL/hzg_mutrig8_2016APV_efficiencies.json",
    "2016postVFP":      f"{base_path}2016postVFP_UL/hzg_mutrig8_2016_efficiencies.json",
    "2017":      f"{base_path}2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2018":      f"{base_path}2018_UL/hzg_mutrig8_2018_efficiencies.json",
    "2022preEE":      f"{base_path}2022preEE_UL/hzg_mutrig8_2022_efficiencies.json",
    "2022postEE":    f"{base_path}2022postEE_UL/hzg_mutrig8_2022EE_efficiencies.json",
    "2023preBPix":      f"{base_path}2023preBPix_UL/hzg_mutrig8_2023_efficiencies.json",
    "2023postBPix":  f"{base_path}2023postBPix_UL/hzg_mutrig8_2023BPix_efficiencies.json",
}


def get_lepton_probabilities(leptons, single_eval, high_eval, low_eval, is_electron, is_data, central_only, year, single_eval_hole=None, high_eval_hole=None, low_eval_hole=None):
    """
    Computes the probabilities and systematics for each lepton to pass exclusive trigger categories.
    """
    n_leptons = awkward.num(leptons)
    if awkward.all(n_leptons == 0):
        # Return empty arrays with correct shape if there are no leptons
        p_shape = (len(leptons), 0, 4)
        s_shape = (len(leptons), 0, 4)
        return awkward.Array(numpy.empty(p_shape)), awkward.Array(numpy.empty(s_shape))

    leptons_flat = awkward.flatten(leptons)
    lep_pt = awkward.to_numpy(leptons_flat.pt)
    lep_eta = awkward.to_numpy(leptons_flat.eta)
    if not is_electron:
        lep_eta = numpy.abs(lep_eta)

    eff_type = "effdata" if is_data else "effmc"
    syst_type = "systdata" if is_data else "systmc"

    # Default efficiencies and systematics
    eff_s = single_eval[eff_type].evalv(lep_pt, lep_eta)
    eff_h = high_eval[eff_type].evalv(lep_pt, lep_eta)
    eff_l = low_eval[eff_type].evalv(lep_pt, lep_eta)
    
    syst_s = syst_h = syst_l = numpy.zeros_like(lep_pt)
    if not central_only:
        syst_s = single_eval[syst_type].evalv(lep_pt, lep_eta)
        syst_h = high_eval[syst_type].evalv(lep_pt, lep_eta)
        syst_l = low_eval[syst_type].evalv(lep_pt, lep_eta)

    # Handle hole region for 2023postBPix electrons
    if is_electron and year == "2023postBPix" and single_eval_hole is not None:
        lep_phi = awkward.to_numpy(leptons_flat.phi)
        hole_mask = (lep_eta > -1.566) & (lep_eta < 0) & (lep_phi > -1.2) & (lep_phi < -0.8)
        
        eff_s_hole = single_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_h_hole = high_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_l_hole = low_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_s = numpy.where(hole_mask, eff_s_hole, eff_s)
        eff_h = numpy.where(hole_mask, eff_h_hole, eff_h)
        eff_l = numpy.where(hole_mask, eff_l_hole, eff_l)

        if not central_only:
            syst_s_hole = single_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_h_hole = high_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_l_hole = low_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_s = numpy.where(hole_mask, syst_s_hole, syst_s)
            syst_h = numpy.where(hole_mask, syst_h_hole, syst_h)
            syst_l = numpy.where(hole_mask, syst_l_hole, syst_l)

    # Probabilities for exclusive categories for each lepton
    # 0: fail all, 1: pass low only, 2: pass high only, 3: pass single only
    p0 = 1.0 - eff_l
    p1 = eff_l - eff_h
    p2 = eff_h - eff_s
    p3 = eff_s
    
    s0 = syst_l # uncertainty of 1-eff_l is syst_l
    s1 = numpy.sqrt(syst_l**2 + syst_h**2)
    s2 = numpy.sqrt(syst_h**2 + syst_s**2)
    s3 = syst_s

    probs = awkward.unflatten(numpy.stack([p0, p1, p2, p3], axis=1), n_leptons)
    systs = awkward.unflatten(numpy.stack([s0, s1, s2, s3], axis=1), n_leptons)
    
    return probs, systs


def calculate_single_lepton_trigger_prob(probs, systs):
    """
    Calculates trigger probability for single-lepton events.
    """
    # For 1 lepton, only single trigger is possible. The probability is just eff_s.
    prob = probs[:, 0, 3] # p3 = eff_s
    unc_sq = systs[:, 0, 3]**2 # s3 = syst_s
    return prob, unc_sq


def calculate_multi_lepton_trigger_prob(probs, systs):
    """
    Calculates trigger probability for events with >= 2 leptons.
    """
    # Efficiencies for each lepton to pass different trigger legs
    eff_s = probs[:, :, 3]
    eff_h = probs[:, :, 2] + probs[:, :, 3]
    eff_l = probs[:, :, 1] + probs[:, :, 2] + probs[:, :, 3]

    # Systematics for each lepton
    syst_s = systs[:, :, 3]
    syst_h = numpy.sqrt(systs[:, :, 2]**2 + systs[:, :, 3]**2)
    syst_l = numpy.sqrt(systs[:, :, 1]**2 + systs[:, :, 2]**2 + systs[:, :, 3]**2)
    print(probs[:, :, 0])
    print(probs[:, :, 1])
    print(probs[:, :, 2])
    print(probs[:, :, 3])

    # --- Single Lepton Trigger Probability ---
    # P(single) = 1 - product(1 - P_i(single))
    prob_single = 1.0 - awkward.prod(1.0 - eff_s, axis=1)
    
    # Uncertainty propagation for single lepton trigger
    # unc_sq_single = sum( (product_{j!=i}(1-p_j))^2 * unc_i^2 )
    # Avoid division by zero if (1-eff_s) is zero
    one_minus_eff_s = 1.0 - eff_s
    prod_one_minus_eff_s = awkward.prod(one_minus_eff_s, axis=1)
    # Add a small epsilon to avoid division by zero
    safe_one_minus_eff_s = one_minus_eff_s + 1e-12
    unc_sq_single = awkward.sum((prod_one_minus_eff_s / safe_one_minus_eff_s)**2 * (syst_s**2), axis=1)

    # --- Dilepton Trigger Probability ---
    # Create all 2-lepton combinations (pairs)
    eff_h_pairs = awkward.combinations(eff_h, 2, fields=["l1", "l2"])
    eff_l_pairs = awkward.combinations(eff_l, 2, fields=["l1", "l2"])
    syst_h_pairs = awkward.combinations(syst_h, 2, fields=["l1", "l2"])
    syst_l_pairs = awkward.combinations(syst_l, 2, fields=["l1", "l2"])

    # Probability for a given pair to pass the dilepton trigger
    # P(pair) = P(l1_h)*P(l2_l) + P(l1_l)*P(l2_h) - P(l1_h)*P(l2_h)
    # This accounts for either lepton being the high or low leg, and subtracts the overlap (both passing high leg)
    prob_pair_pass = (eff_h_pairs.l1 * eff_l_pairs.l2) + (eff_l_pairs.l1 * eff_h_pairs.l2) - (eff_h_pairs.l1 * eff_h_pairs.l2)
    
    # P(dilep) = 1 - product(1 - P(pair_i))
    prob_dilep = 1.0 - awkward.prod(1.0 - prob_pair_pass, axis=1)

    # Dilepton uncertainty (approximation)
    # Uncertainty for P(pair)
    # unc_sq_pair_pass = (eff_l_pairs.l2 * syst_h_pairs.l1)**2 + (eff_h_pairs.l1 * syst_l_pairs.l2)**2
    unc_sq_pair_pass = (eff_l_pairs.l2 - eff_h_pairs.l2)**2 * syst_h_pairs.l1**2 + \
                       (eff_h_pairs.l1)**2 * syst_l_pairs.l2**2 + \
                       (eff_l_pairs.l1 - eff_h_pairs.l1)**2 * syst_h_pairs.l2**2 + \
                       (eff_h_pairs.l2)**2 * syst_l_pairs.l1**2
    
    one_minus_prob_pair = 1.0 - prob_pair_pass
    prod_one_minus_prob_pair = awkward.prod(one_minus_prob_pair, axis=1)
    safe_one_minus_prob_pair = one_minus_prob_pair + 1e-12
    unc_sq_dilep = awkward.sum((prod_one_minus_prob_pair / safe_one_minus_prob_pair)**2 * unc_sq_pair_pass, axis=1)

    # Total probability: P(total) = P(single) + P(dilep) - P(single) * P(dilep)
    total_prob = prob_single + prob_dilep - prob_single * prob_dilep
    
    # Total uncertainty
    total_unc_sq = (1 - prob_dilep)**2 * unc_sq_single + (1 - prob_single)**2 * unc_sq_dilep

    return total_prob, total_unc_sq


def get_flavor_probability(leptons, single_eval, high_eval, low_eval,
                           single_thresh, high_thresh, low_thresh,
                           is_electron, is_data, central_only, year,
                           single_eval_hole=None, high_eval_hole=None, low_eval_hole=None):
    """
    Calculate the total trigger probability for a single lepton flavor (electron or muon) for each event.
    """
    n_leptons = awkward.num(leptons)
    total_prob = numpy.zeros(len(leptons))
    total_unc_sq = numpy.zeros(len(leptons))

    # Get per-lepton probabilities and systematics for all events
    probs, systs = get_lepton_probabilities(
        leptons, single_eval, high_eval, low_eval,
        is_electron, is_data, central_only, year,
        single_eval_hole, high_eval_hole, low_eval_hole
    )

    # Case 1: One lepton
    mask_1lep = (n_leptons == 1)
    if numpy.any(mask_1lep):
        prob_1lep, unc_sq_1lep = calculate_single_lepton_trigger_prob(probs[mask_1lep], systs[mask_1lep])
        total_prob[mask_1lep] = prob_1lep
        total_unc_sq[mask_1lep] = unc_sq_1lep

    # Case 2 or more: Multiple leptons
    mask_multi_lep = (n_leptons >= 2)
    if numpy.any(mask_multi_lep):
        prob_multi, unc_sq_multi = calculate_multi_lepton_trigger_prob(
            probs[mask_multi_lep], systs[mask_multi_lep]
        )
        total_prob[mask_multi_lep] = awkward.to_numpy(prob_multi)
        total_unc_sq[mask_multi_lep] = awkward.to_numpy(unc_sq_multi)
    
    logger.debug(f"n_lep = 0: {awkward.sum(n_leptons == 0)}, n_lep = 1: {awkward.sum(n_leptons == 1)}, n_lep >= 2: {awkward.sum(n_leptons >= 2)}")
    print("n_lep==1:", awkward.mean(total_prob[n_leptons == 1]))
    print("n_lep>=2:", awkward.mean(total_prob[n_leptons >= 2]))

    return total_prob, numpy.sqrt(total_unc_sq)


def HLT_sf(events, year, central_only):
    """
    Calculate event-level HLT scale factors.
    The logic is adapted from the C++ implementation, calculating a total event probability
    for passing the trigger for data and MC, then taking the ratio.
    """
    # Load electron evaluators
    ele_single_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleElectron_HLT_FILE[year]))
    ele_high_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_HighLeg_HLT_FILE[year]))
    ele_low_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_LowLeg_HLT_FILE[year]))
    
    ele_single_eval_hole = ele_high_eval_hole = ele_low_eval_hole = None
    if year == "2023postBPix":
        ele_single_eval_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleElectron_HLT_FILE_HOLE[year]))
        ele_high_eval_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_HighLeg_HLT_FILE_HOLE[year]))
        ele_low_eval_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleElectron_LowLeg_HLT_FILE_HOLE[year]))

    # Load muon evaluators
    muon_single_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleMuon_HLT_FILE[year]))
    muon_high_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_HighLeg_HLT_FILE[year]))
    muon_low_eval = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_LowLeg_HLT_FILE[year]))

    electrons = events.Electron
    muons = events.Muon
    electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
    muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]

    # Define pt thresholds
    if "2016" in year:
        ele_single_thresh, ele_high_thresh, ele_low_thresh = 27.0, 23.0, 12.0
        muon_single_thresh, muon_high_thresh, muon_low_thresh = 24.0, 17.0, 8.0
    elif year == "2017":
        ele_single_thresh, ele_high_thresh, ele_low_thresh = 32.0, 23.0, 12.0
        muon_single_thresh, muon_high_thresh, muon_low_thresh = 27.0, 17.0, 8.0
    elif year == "2018":
        ele_single_thresh, ele_high_thresh, ele_low_thresh = 32.0, 23.0, 12.0
        muon_single_thresh, muon_high_thresh, muon_low_thresh = 24.0, 17.0, 8.0
    else: # 2022/2023
        ele_single_thresh, ele_high_thresh, ele_low_thresh = 30.0, 23.0, 12.0
        muon_single_thresh, muon_high_thresh, muon_low_thresh = 24.0, 17.0, 8.0

    # Calculate probabilities for data
    prob_e_data, unc_e_data = get_flavor_probability(
            electrons, ele_single_eval, ele_high_eval, ele_low_eval,
            ele_single_thresh, ele_high_thresh, ele_low_thresh,
            True, True, central_only, year,
            ele_single_eval_hole, ele_high_eval_hole, ele_low_eval_hole)
    prob_m_data, unc_m_data = get_flavor_probability(
            muons, muon_single_eval, muon_high_eval, muon_low_eval,
            muon_single_thresh, muon_high_thresh, muon_low_thresh,
            False, True, central_only, year)
    
    # Calculate probabilities for MC
    prob_e_mc, unc_e_mc = get_flavor_probability(
            electrons, ele_single_eval, ele_high_eval, ele_low_eval,
            ele_single_thresh, ele_high_thresh, ele_low_thresh,
            True, False, central_only, year,
            ele_single_eval_hole, ele_high_eval_hole, ele_low_eval_hole)
    prob_m_mc, unc_m_mc = get_flavor_probability(
            muons, muon_single_eval, muon_high_eval, muon_low_eval,
            muon_single_thresh, muon_high_thresh, muon_low_thresh,
            False, False, central_only, year)

    # Combine probabilities: P = P_e + P_m - P_e * P_m
    total_prob_data = prob_e_data + prob_m_data - prob_e_data * prob_m_data
    total_prob_mc = prob_e_mc + prob_m_mc - prob_e_mc * prob_m_mc

    # Calculate SF
    sf = numpy.divide(total_prob_data, total_prob_mc, out=numpy.ones_like(total_prob_data), where=total_prob_mc!=0)
    
    variations = {"central": sf}

    if not central_only:
        # Propagate uncertainties
        unc_data_sq = (1 - prob_m_data)**2 * unc_e_data**2 + (1 - prob_e_data)**2 * unc_m_data**2
        unc_mc_sq = (1 - prob_m_mc)**2 * unc_e_mc**2 + (1 - prob_e_mc)**2 * unc_m_mc**2
        
        rel_unc_sq = numpy.divide(unc_data_sq, total_prob_data**2, out=numpy.zeros_like(unc_data_sq), where=total_prob_data!=0) + \
                     numpy.divide(unc_mc_sq, total_prob_mc**2, out=numpy.zeros_like(unc_mc_sq), where=total_prob_mc!=0)
        
        unc_sf = sf * numpy.sqrt(rel_unc_sq)
        
        variations["up"] = sf + unc_sf
        variations["down"] = sf - unc_sf

    # Clip values to be within a reasonable range, e.g., [0, 5] and convert to awkward array
    for var in variations:
        variations[var] = awkward.Array(numpy.clip(variations[var], 0, 5))

    return variations