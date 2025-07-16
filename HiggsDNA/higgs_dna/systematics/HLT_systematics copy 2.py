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


def get_flavor_probability(leptons, single_eval, high_eval, low_eval,
                           single_thresh, high_thresh, low_thresh,
                           is_electron, is_data, central_only, year,
                           single_eval_hole=None, high_eval_hole=None, low_eval_hole=None):
    """
    Calculate the total trigger probability for a single lepton flavor (electron or muon) for each event.
    This function reproduces the logic from the C++ TriggerWeighter by iterating through
    all possible trigger status combinations for all leptons in an event.
    """
    n_leptons = awkward.num(leptons)
    if awkward.all(n_leptons == 0):
        return numpy.zeros(len(leptons)), numpy.zeros(len(leptons))

    leptons_flat = awkward.flatten(leptons)
    lep_pt = awkward.to_numpy(leptons_flat.pt)
    lep_eta = awkward.to_numpy(leptons_flat.eta)
    if not is_electron:
        lep_eta = numpy.abs(lep_eta)

    # Get efficiencies for each lepton for each trigger leg
    eff_type = "effdata" if is_data else "effmc"
    syst_type = "systdata" if is_data else "systmc"

    hole_mask = None
    if is_electron and year == "2023postBPix":
        lep_phi = awkward.to_numpy(leptons_flat.phi)
        hole_mask = (lep_eta > -1.566) & (lep_eta < 0) & (lep_phi > -1.2) & (lep_phi < -0.8)

    if is_electron and year == "2023postBPix":
        eff_s_normal = single_eval[eff_type].evalv(lep_pt, lep_eta)
        eff_h_normal = high_eval[eff_type].evalv(lep_pt, lep_eta)
        eff_l_normal = low_eval[eff_type].evalv(lep_pt, lep_eta)
        eff_s_hole = single_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_h_hole = high_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_l_hole = low_eval_hole[eff_type].evalv(lep_pt, lep_eta)
        eff_s = numpy.where(hole_mask, eff_s_hole, eff_s_normal)
        eff_h = numpy.where(hole_mask, eff_h_hole, eff_h_normal)
        eff_l = numpy.where(hole_mask, eff_l_hole, eff_l_normal)

        if not central_only:
            syst_s_normal = single_eval[syst_type].evalv(lep_pt, lep_eta)
            syst_h_normal = high_eval[syst_type].evalv(lep_pt, lep_eta)
            syst_l_normal = low_eval[syst_type].evalv(lep_pt, lep_eta)
            syst_s_hole = single_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_h_hole = high_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_l_hole = low_eval_hole[syst_type].evalv(lep_pt, lep_eta)
            syst_s = numpy.where(hole_mask, syst_s_hole, syst_s_normal)
            syst_h = numpy.where(hole_mask, syst_h_hole, syst_h_normal)
            syst_l = numpy.where(hole_mask, syst_l_hole, syst_l_normal)
        else:
            syst_s = syst_h = syst_l = 0

    else:
        eff_s = single_eval[eff_type].evalv(lep_pt, lep_eta)
        eff_h = high_eval[eff_type].evalv(lep_pt, lep_eta)
        eff_l = low_eval[eff_type].evalv(lep_pt, lep_eta)

        syst_s = syst_h = syst_l = 0
        if not central_only:
            syst_s = single_eval[syst_type].evalv(lep_pt, lep_eta)
            syst_h = high_eval[syst_type].evalv(lep_pt, lep_eta)
            syst_l = low_eval[syst_type].evalv(lep_pt, lep_eta)

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

    # Event-level trigger decisions based on pt thresholds
    pass_singlelep = awkward.any(leptons.pt > single_thresh, axis=1)
    
    # For dilepton trigger, need one lepton above high_thresh and another above low_thresh
    # Ensure they are not the same lepton
    pairs = awkward.combinations(leptons[awkward.argsort(leptons.pt, ascending=False, axis=1)], 2, fields = ["slot0", "slot1"])
    pass_dilep = awkward.any(
        ((pairs.slot0.pt > high_thresh) & (pairs.slot1.pt > low_thresh)),
        axis=1
    )

    # Iterate over all combinations of lepton states
    max_leptons = awkward.max(n_leptons)
    if max_leptons > 6: # Cap to prevent excessive computation time
        logger.warning(f"Capping max_leptons from {max_leptons} to 6 to limit HLT SF computation time.")
        max_leptons = 6

    total_prob = numpy.zeros(len(leptons))
    total_unc_sq = numpy.zeros(len(leptons))

    for i in range(4**max_leptons):
        temp_i = i
        
        cat_prob = numpy.ones(len(leptons))
        cat_unc_sq = numpy.zeros(len(leptons))
        
        n_pass_lower = numpy.zeros(len(leptons), dtype=int)
        n_pass_upper = numpy.zeros(len(leptons), dtype=int)
        n_pass_single = numpy.zeros(len(leptons), dtype=int)
        
        relevant_event = (n_leptons > 0)

        for j in range(max_leptons):
            state = temp_i % 4
            temp_i //= 4

            # Mask for events with at least j+1 leptons
            mask = (n_leptons > j)
            relevant_event = relevant_event & mask

            # Get prob and syst for j-th lepton in its current state
            p_j = awkward.to_numpy(probs[mask, j, state])
            s_j = awkward.to_numpy(systs[mask, j, state])
            s_j_sq = s_j**2

            # Propagate uncertainty: (s_new/p_new)^2 = (s_old/p_old)^2 + (s_j/p_j)^2
            # p_new = p_old * p_j
            # s_new^2 = s_old^2 * p_j^2 + p_old^2 * s_j^2
            
            p_old = cat_prob[mask]
            s_old_sq = cat_unc_sq[mask]

            cat_prob[mask] = p_old * p_j
            cat_unc_sq[mask] = s_old_sq * (p_j**2) + (p_old**2) * s_j_sq

            if state == 1: # pass low
                n_pass_lower[mask] += 1
            elif state == 2: # pass high
                n_pass_lower[mask] += 1
                n_pass_upper[mask] += 1
            elif state == 3: # pass single
                n_pass_lower[mask] += 1
                n_pass_upper[mask] += 1
                n_pass_single[mask] += 1
        
        # Check if this category is consistent with the event's trigger decision
        # This logic is from the C++ file
        is_relevant_cat = ((n_pass_single > 0) == pass_singlelep) & \
                          (((n_pass_upper > 0) & (n_pass_lower > 1)) == pass_dilep)
        
        is_relevant_cat_np = awkward.to_numpy(is_relevant_cat)

        # Add probability of this category if it's relevant
        total_prob[relevant_event] += cat_prob[relevant_event] * is_relevant_cat_np[relevant_event]
        total_unc_sq[relevant_event] += cat_unc_sq[relevant_event] * is_relevant_cat_np[relevant_event]

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
        
        # Handle cases where event did not pass trigger
        pass_single_e = awkward.any(electrons.pt > ele_single_thresh, axis=1)
        pairs_e = awkward.combinations(electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)], 2, fields = ["slot0", "slot1"])
        pass_double_e = awkward.any(
            ((pairs_e.slot0.pt > ele_high_thresh) & (pairs_e.slot1.pt > ele_low_thresh)),
            axis=1
        )
        pass_single_m = awkward.any(muons.pt > muon_single_thresh, axis=1)
        pairs_m = awkward.combinations(muons[awkward.argsort(muons.pt, ascending=False, axis=1)], 2, fields = ["slot0", "slot1"])
        pass_double_m = awkward.any(
            ((pairs_m.slot0.pt > muon_high_thresh) & (pairs_m.slot1.pt > muon_low_thresh)),
            axis=1
        )

        pass_trig = awkward.to_numpy(pass_single_e | pass_double_e | pass_single_m | pass_double_m)
        
        variations["up"] = sf + unc_sf
        variations["down"] = sf - unc_sf
        
        variations["up"][~pass_trig] = sf[~pass_trig] - unc_sf[~pass_trig]
        variations["down"][~pass_trig] = sf[~pass_trig] + unc_sf[~pass_trig]

    # Clip values to be within a reasonable range, e.g., [0, 5] and convert to awkward array
    for var in variations:
        variations[var] = awkward.Array(numpy.clip(variations[var], 0, 5))

    return variations