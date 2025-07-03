import awkward
import numpy

from correctionlib import _core
import json

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.utils import systematic_from_bins

SingleEle_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_eltrig27_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_eltrig27_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_eltrig32_2018_efficiencies.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_eltrig30_2022_efficiencies.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_eltrig30_2022EE_efficiencies.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_eltrig30_2023_efficiencies.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig30_2023BPix_efficiencies.json"
}
DoubleEle_LowLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_eltrig12_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_eltrig12_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_eltrig12_2018_efficiencies.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_eltrig12_2022_efficiencies.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_eltrig12_2022EE_efficiencies.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_eltrig12_2023_efficiencies.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig12_2023BPix_efficiencies.json"
}
DoubleEle_HighLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_eltrig23_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_eltrig23_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig23_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_eltrig23_2018_efficiencies.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_eltrig23_2022_efficiencies.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_eltrig23_2022EE_efficiencies.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_eltrig23_2023_efficiencies.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig23_2023BPix_efficiencies.json"
}

# Hole files for 2023postBPix
SingleEle_HLT_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig30_2023BPixHole_efficiencies.json"
}
DoubleEle_LowLeg_HLT_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig12_2023BPixHole_efficiencies.json"
}
DoubleEle_HighLeg_HLT_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_eltrig23_2023BPixHole_efficiencies.json"
}

SingleMuon_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_mutrig24_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_mutrig24_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig27_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_mutrig24_2018_efficiencies.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_mutrig24_2022_efficiencies.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_mutrig24_2022EE_efficiencies.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_mutrig24_2023_efficiencies.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_mutrig24_2023BPix_efficiencies.json"

}
DoubleMuon_LowLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_mutrig8_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_mutrig8_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_mutrig8_2018_efficiencies.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_mutrig8_2022_efficiencies.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_mutrig8_2022EE_efficiencies.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_mutrig8_2023_efficiencies.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_mutrig8_2023BPix_efficiencies.json"
}
DoubleMuon_HighLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_mutrig17_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_mutrig17_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_mutrig17_2018_efficiencies.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_mutrig17_2022_efficiencies.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_mutrig17_2022EE_efficiencies.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_mutrig17_2023_efficiencies.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_mutrig17_2023BPix_efficiencies.json"
}

def Ele_HLT_sf(events, year, central_only):
    """
    Calculate electron HLT scale factors based on pt thresholds.
    For electrons in different pt ranges, use different efficiency calculations:
    - pt < doubleele_lowleg_threshold: 1 - eff(doubleele_lowleg)
    - doubleele_lowleg_threshold <= pt < doubleele_highleg_threshold: eff(doubleele_highleg) - eff(doubleele_lowleg)
    - doubleele_highleg_threshold <= pt < singleele_threshold: eff(singleele) - eff(doubleele_highleg)
    - pt >= singleele_threshold: eff(singleele)
    
    SF = EFF(data) / EFF(mc)
    """
    required_fields = [
        ("Electron", "pt"), ("Electron", "eta")
    ]
    
    # For 2023postBPix, we also need phi to determine which SF file to use
    if year == "2023postBPix":
        required_fields.append(("Electron", "phi"))

    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning(f"Missing fields for electron HLT SF: {missing_fields}")
        return None

    # Load evaluators for all three trigger types
    single_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleEle_HLT_FILE[year]))
    double_high_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_HighLeg_HLT_FILE[year]))
    double_low_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_LowLeg_HLT_FILE[year]))
    
    # For 2023postBPix, also load the Hole files
    single_evaluator_hole = None
    double_high_evaluator_hole = None
    double_low_evaluator_hole = None
    if year == "2023postBPix":
        single_evaluator_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleEle_HLT_FILE_HOLE[year]))
        double_high_evaluator_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_HighLeg_HLT_FILE_HOLE[year]))
        double_low_evaluator_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_LowLeg_HLT_FILE_HOLE[year]))

    electrons = events["Electron"]

    # Flatten electrons then convert to numpy for compatibility with correctionlib
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    ele_eta = numpy.clip(
        awkward.to_numpy(electrons_flattened.eta),
        -2.49999,
        2.49999  # SFs only valid up to eta 2.5
    )

    ele_pt = numpy.clip(
        awkward.to_numpy(electrons_flattened.pt),
        7.0,     # SFs only valid for pT >= 7.0
        499.999  # and pT < 500.
    )
    
    # Get phi for 2023postBPix
    ele_phi = None
    hole_mask = None
    if year == "2023postBPix":
        ele_phi = awkward.to_numpy(electrons_flattened.phi)
        # Hole region: -1.566 < eta < 0 and -1.2 < phi < -0.8
        hole_mask = (ele_eta > -1.566) & (ele_eta < 0) & (ele_phi > -1.2) & (ele_phi < -0.8)

    # Get efficiencies from all three evaluators for both MC and data
    if year == "2023postBPix":
        # For 2023postBPix, combine normal and hole regions
        # Single electron efficiencies
        eff_single_mc_normal = single_evaluator["effmc"].evalv(ele_pt, ele_eta)
        eff_single_mc_hole = single_evaluator_hole["effmc"].evalv(ele_pt, ele_eta)
        eff_single_mc = numpy.where(hole_mask, eff_single_mc_hole, eff_single_mc_normal)
        
        eff_single_data_normal = single_evaluator["effdata"].evalv(ele_pt, ele_eta)
        eff_single_data_hole = single_evaluator_hole["effdata"].evalv(ele_pt, ele_eta)
        eff_single_data = numpy.where(hole_mask, eff_single_data_hole, eff_single_data_normal)
        
        # Double electron high leg efficiencies
        eff_double_high_mc_normal = double_high_evaluator["effmc"].evalv(ele_pt, ele_eta)
        eff_double_high_mc_hole = double_high_evaluator_hole["effmc"].evalv(ele_pt, ele_eta)
        eff_double_high_mc = numpy.where(hole_mask, eff_double_high_mc_hole, eff_double_high_mc_normal)
        
        eff_double_high_data_normal = double_high_evaluator["effdata"].evalv(ele_pt, ele_eta)
        eff_double_high_data_hole = double_high_evaluator_hole["effdata"].evalv(ele_pt, ele_eta)
        eff_double_high_data = numpy.where(hole_mask, eff_double_high_data_hole, eff_double_high_data_normal)
        
        # Double electron low leg efficiencies
        eff_double_low_mc_normal = double_low_evaluator["effmc"].evalv(ele_pt, ele_eta)
        eff_double_low_mc_hole = double_low_evaluator_hole["effmc"].evalv(ele_pt, ele_eta)
        eff_double_low_mc = numpy.where(hole_mask, eff_double_low_mc_hole, eff_double_low_mc_normal)
        
        eff_double_low_data_normal = double_low_evaluator["effdata"].evalv(ele_pt, ele_eta)
        eff_double_low_data_hole = double_low_evaluator_hole["effdata"].evalv(ele_pt, ele_eta)
        eff_double_low_data = numpy.where(hole_mask, eff_double_low_data_hole, eff_double_low_data_normal)
    else:
        # For other years, use the standard approach
        eff_single_mc = single_evaluator["effmc"].evalv(ele_pt, ele_eta)
        eff_double_high_mc = double_high_evaluator["effmc"].evalv(ele_pt, ele_eta)
        eff_double_low_mc = double_low_evaluator["effmc"].evalv(ele_pt, ele_eta)
        
        eff_single_data = single_evaluator["effdata"].evalv(ele_pt, ele_eta)
        eff_double_high_data = double_high_evaluator["effdata"].evalv(ele_pt, ele_eta)
        eff_double_low_data = double_low_evaluator["effdata"].evalv(ele_pt, ele_eta)

    # Define pt thresholds (these should be adjusted based on your trigger configuration)
    # For now using typical CMS values
    if "2016" in year:
        single_threshold = 27.0
        double_high_threshold = 23.0
        double_low_threshold = 12.0
    elif year == "2017":
        single_threshold = 32.0
        double_high_threshold = 23.0
        double_low_threshold = 12.0
    elif year == "2018":
        single_threshold = 32.0
        double_high_threshold = 23.0
        double_low_threshold = 12.0
    elif "2022" in year or "2023" in year:
        single_threshold = 30.0
        double_high_threshold = 23.0
        double_low_threshold = 12.0
    else:
        # Default values
        single_threshold = 30.0
        double_high_threshold = 23.0
        double_low_threshold = 12.0

    # Calculate effective efficiency based on pt ranges for MC and data
    eff_combined_mc = numpy.where(
        ele_pt < double_low_threshold,
        1.0 - eff_double_low_mc,
        numpy.where(
            ele_pt < double_high_threshold,
            eff_double_high_mc - eff_double_low_mc,
            numpy.where(
                ele_pt < single_threshold,
                eff_single_mc - eff_double_high_mc,
                eff_single_mc
            )
        )
    )
    
    eff_combined_data = numpy.where(
        ele_pt < double_low_threshold,
        1.0 - eff_double_low_data,
        numpy.where(
            ele_pt < double_high_threshold,
            eff_double_high_data - eff_double_low_data,
            numpy.where(
                ele_pt < single_threshold,
                eff_single_data - eff_double_high_data,
                eff_single_data
            )
        )
    )
    
    # Calculate scale factor: SF = EFF(data) / EFF(mc)
    sf_combined = numpy.where(
        eff_combined_mc > 0,
        eff_combined_data / eff_combined_mc,
        1.0  # Set SF = 1 when MC efficiency is 0
    )

    variations = {}
    variations["central"] = awkward.unflatten(sf_combined, n_electrons)

    if not central_only:
        # Get systematic uncertainties for both MC and data
        if year == "2023postBPix":
            # For 2023postBPix, combine normal and hole regions
            # Single electron uncertainties
            syst_single_mc_normal = single_evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst_single_mc_hole = single_evaluator_hole["systmc"].evalv(ele_pt, ele_eta)
            syst_single_mc = numpy.where(hole_mask, syst_single_mc_hole, syst_single_mc_normal)
            
            syst_single_data_normal = single_evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_single_data_hole = single_evaluator_hole["systdata"].evalv(ele_pt, ele_eta)
            syst_single_data = numpy.where(hole_mask, syst_single_data_hole, syst_single_data_normal)
            
            # Double electron high leg uncertainties
            syst_double_high_mc_normal = double_high_evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst_double_high_mc_hole = double_high_evaluator_hole["systmc"].evalv(ele_pt, ele_eta)
            syst_double_high_mc = numpy.where(hole_mask, syst_double_high_mc_hole, syst_double_high_mc_normal)
            
            syst_double_high_data_normal = double_high_evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_double_high_data_hole = double_high_evaluator_hole["systdata"].evalv(ele_pt, ele_eta)
            syst_double_high_data = numpy.where(hole_mask, syst_double_high_data_hole, syst_double_high_data_normal)
            
            # Double electron low leg uncertainties
            syst_double_low_mc_normal = double_low_evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst_double_low_mc_hole = double_low_evaluator_hole["systmc"].evalv(ele_pt, ele_eta)
            syst_double_low_mc = numpy.where(hole_mask, syst_double_low_mc_hole, syst_double_low_mc_normal)
            
            syst_double_low_data_normal = double_low_evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_double_low_data_hole = double_low_evaluator_hole["systdata"].evalv(ele_pt, ele_eta)
            syst_double_low_data = numpy.where(hole_mask, syst_double_low_data_hole, syst_double_low_data_normal)
        else:
            # For other years, use the standard approach
            syst_single_mc = single_evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst_double_high_mc = double_high_evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst_double_low_mc = double_low_evaluator["systmc"].evalv(ele_pt, ele_eta)
            
            syst_single_data = single_evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_double_high_data = double_high_evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_double_low_data = double_low_evaluator["systdata"].evalv(ele_pt, ele_eta)

        # Propagate uncertainties for combined efficiency - MC
        syst_combined_mc = numpy.where(
            ele_pt < double_low_threshold,
            syst_double_low_mc,  # For (1 - eff_double_low)
            numpy.where(
                ele_pt < double_high_threshold,
                numpy.sqrt(syst_double_high_mc**2 + syst_double_low_mc**2),  # For (eff_double_high - eff_double_low)
                numpy.where(
                    ele_pt < single_threshold,
                    numpy.sqrt(syst_single_mc**2 + syst_double_high_mc**2),  # For (eff_single - eff_double_high)
                    syst_single_mc  # For eff_single
                )
            )
        )
        
        # Propagate uncertainties for combined efficiency - data
        syst_combined_data = numpy.where(
            ele_pt < double_low_threshold,
            syst_double_low_data,  # For (1 - eff_double_low)
            numpy.where(
                ele_pt < double_high_threshold,
                numpy.sqrt(syst_double_high_data**2 + syst_double_low_data**2),  # For (eff_double_high - eff_double_low)
                numpy.where(
                    ele_pt < single_threshold,
                    numpy.sqrt(syst_single_data**2 + syst_double_high_data**2),  # For (eff_single - eff_double_high)
                    syst_single_data  # For eff_single
                )
            )
        )
        
        # Calculate systematic uncertainty for scale factor: d(SF)/SF = sqrt((d_data/eff_data)^2 + (d_mc/eff_mc)^2)
        # where d_data and d_mc are the absolute uncertainties
        rel_syst_data = numpy.where(eff_combined_data > 0, syst_combined_data / eff_combined_data, 0)
        rel_syst_mc = numpy.where(eff_combined_mc > 0, syst_combined_mc / eff_combined_mc, 0)
        rel_syst_sf = numpy.sqrt(rel_syst_data**2 + rel_syst_mc**2)
        syst_sf = sf_combined * rel_syst_sf

        variations["up"] = awkward.unflatten(sf_combined + syst_sf, n_electrons)
        variations["down"] = awkward.unflatten(sf_combined - syst_sf, n_electrons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
            (electrons.pt < 7.0) | (electrons.pt >= 500.0) | (abs(electrons.eta) >= 2.5),
            awkward.ones_like(variations[var], dtype=float),
            variations[var]
        )

    return variations

def Muon_HLT_sf(events, year, central_only):
    """
    Calculate muon HLT scale factors based on pt thresholds.
    For muons in different pt ranges, use different efficiency calculations:
    - pt < doublemuon_lowleg_threshold: 1 - eff(doublemuon_lowleg)
    - doublemuon_lowleg_threshold <= pt < doublemuon_highleg_threshold: eff(doublemuon_highleg) - eff(doublemuon_lowleg)
    - doublemuon_highleg_threshold <= pt < singlemuon_threshold: eff(singlemuon) - eff(doublemuon_highleg)
    - pt >= singlemuon_threshold: eff(singlemuon)
    
    SF = EFF(data) / EFF(mc)
    """
    required_fields = [
        ("Muon", "pt"), ("Muon", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning(f"Missing fields for muon HLT SF: {missing_fields}")
        return None

    # Load evaluators for all three trigger types
    single_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleMuon_HLT_FILE[year]))
    double_high_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_HighLeg_HLT_FILE[year]))
    double_low_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_LowLeg_HLT_FILE[year]))

    muons = events["Muon"]

    # Flatten muons then convert to numpy for compatibility with correctionlib
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    mu_eta = numpy.clip(
        awkward.to_numpy(muons_flattened.eta),
        -2.39999,
        2.39999  # SFs only valid up to eta 2.4
    )

    mu_pt = numpy.clip(
        awkward.to_numpy(muons_flattened.pt),
        5.0,     # SFs only valid for pT >= 5.0
        499.999  # and pT < 500.
    )

    # Get efficiencies from all three evaluators for both MC and data
    eff_single_mc = single_evaluator["effmc"].evalv(mu_pt, mu_eta)
    eff_double_high_mc = double_high_evaluator["effmc"].evalv(mu_pt, mu_eta)
    eff_double_low_mc = double_low_evaluator["effmc"].evalv(mu_pt, mu_eta)
    
    eff_single_data = single_evaluator["effdata"].evalv(mu_pt, mu_eta)
    eff_double_high_data = double_high_evaluator["effdata"].evalv(mu_pt, mu_eta)
    eff_double_low_data = double_low_evaluator["effdata"].evalv(mu_pt, mu_eta)

    # Define pt thresholds based on trigger configuration
    if "2016" in year:
        single_threshold = 24.0
        double_high_threshold = 17.0
        double_low_threshold = 8.0
    elif year == "2017":
        single_threshold = 27.0
        double_high_threshold = 17.0
        double_low_threshold = 8.0
    elif year == "2018":
        single_threshold = 24.0
        double_high_threshold = 17.0
        double_low_threshold = 8.0
    elif "2022" in year or "2023" in year:
        single_threshold = 24.0
        double_high_threshold = 17.0
        double_low_threshold = 8.0
    else:
        # Default values
        single_threshold = 24.0
        double_high_threshold = 17.0
        double_low_threshold = 8.0

    # Calculate effective efficiency based on pt ranges for MC and data
    eff_combined_mc = numpy.where(
        mu_pt < double_low_threshold,
        1.0 - eff_double_low_mc,
        numpy.where(
            mu_pt < double_high_threshold,
            eff_double_high_mc - eff_double_low_mc,
            numpy.where(
                mu_pt < single_threshold,
                eff_single_mc - eff_double_high_mc,
                eff_single_mc
            )
        )
    )
    
    eff_combined_data = numpy.where(
        mu_pt < double_low_threshold,
        1.0 - eff_double_low_data,
        numpy.where(
            mu_pt < double_high_threshold,
            eff_double_high_data - eff_double_low_data,
            numpy.where(
                mu_pt < single_threshold,
                eff_single_data - eff_double_high_data,
                eff_single_data
            )
        )
    )
    
    # Calculate scale factor: SF = EFF(data) / EFF(mc)
    sf_combined = numpy.where(
        eff_combined_mc > 0,
        eff_combined_data / eff_combined_mc,
        1.0  # Set SF = 1 when MC efficiency is 0
    )

    variations = {}
    variations["central"] = awkward.unflatten(sf_combined, n_muons)

    if not central_only:
        # Get systematic uncertainties for both MC and data
        syst_single_mc = single_evaluator["systmc"].evalv(mu_pt, mu_eta)
        syst_double_high_mc = double_high_evaluator["systmc"].evalv(mu_pt, mu_eta)
        syst_double_low_mc = double_low_evaluator["systmc"].evalv(mu_pt, mu_eta)
        
        syst_single_data = single_evaluator["systdata"].evalv(mu_pt, mu_eta)
        syst_double_high_data = double_high_evaluator["systdata"].evalv(mu_pt, mu_eta)
        syst_double_low_data = double_low_evaluator["systdata"].evalv(mu_pt, mu_eta)

        # Propagate uncertainties for combined efficiency - MC
        syst_combined_mc = numpy.where(
            mu_pt < double_low_threshold,
            syst_double_low_mc,  # For (1 - eff_double_low)
            numpy.where(
                mu_pt < double_high_threshold,
                numpy.sqrt(syst_double_high_mc**2 + syst_double_low_mc**2),  # For (eff_double_high - eff_double_low)
                numpy.where(
                    mu_pt < single_threshold,
                    numpy.sqrt(syst_single_mc**2 + syst_double_high_mc**2),  # For (eff_single - eff_double_high)
                    syst_single_mc  # For eff_single
                )
            )
        )
        
        # Propagate uncertainties for combined efficiency - data
        syst_combined_data = numpy.where(
            mu_pt < double_low_threshold,
            syst_double_low_data,  # For (1 - eff_double_low)
            numpy.where(
                mu_pt < double_high_threshold,
                numpy.sqrt(syst_double_high_data**2 + syst_double_low_data**2),  # For (eff_double_high - eff_double_low)
                numpy.where(
                    mu_pt < single_threshold,
                    numpy.sqrt(syst_single_data**2 + syst_double_high_data**2),  # For (eff_single - eff_double_high)
                    syst_single_data  # For eff_single
                )
            )
        )
        
        # Calculate systematic uncertainty for scale factor: d(SF)/SF = sqrt((d_data/eff_data)^2 + (d_mc/eff_mc)^2)
        # where d_data and d_mc are the absolute uncertainties
        rel_syst_data = numpy.where(eff_combined_data > 0, syst_combined_data / eff_combined_data, 0)
        rel_syst_mc = numpy.where(eff_combined_mc > 0, syst_combined_mc / eff_combined_mc, 0)
        rel_syst_sf = numpy.sqrt(rel_syst_data**2 + rel_syst_mc**2)
        syst_sf = sf_combined * rel_syst_sf

        variations["up"] = awkward.unflatten(sf_combined + syst_sf, n_muons)
        variations["down"] = awkward.unflatten(sf_combined - syst_sf, n_muons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
            (muons.pt < 5.0) | (muons.pt >= 500.0) | (abs(muons.eta) >= 2.4),
            awkward.ones_like(variations[var], dtype=float),
            variations[var]
        )

    return variations