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

MUON_ID_SF_FILE = {
    "2016" : "higgs_dna/systematics/data/2016postVFP_UL/muid_2016_2016APV.json",
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/muid_2016_2016APV.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/muid_2016_2016APV.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/muid_2017.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/muid_2018.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/muid_2022.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/muid_2022EE.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_muid_2023_scalefactors.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_muid_2023BPix_scalefactors.json"
}

def muon_LooseID_sf(events, year, central_only, input_collection):

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(MUON_ID_SF_FILE[year]))

    muons = events[input_collection]

    # Flatten muons then convert to numpy for compatibility with correctionlib
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    muon_eta = numpy.clip(
        awkward.to_numpy(muons_flattened.eta),
        -2.39999,
        2.39999 # SFs only valid up to eta 2.5
    )

    muon_pt = numpy.clip(
        awkward.to_numpy(muons_flattened.pt),
        5.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    # sf = evaluator["Muon_LooseID_MCeff"].evalv(
    #         "effmc",
    #         abs(muon_eta),
    #         muon_pt
    # )
    sf = evaluator["sf_pass"].evalv(
            muon_pt,
            muon_eta
    )
    variations["central"] = awkward.unflatten(sf, n_muons)

    if not central_only:
        # syst_var = "systmc"
        # syst = evaluator["Muon_LooseID_MCeff"].evalv(
        #     syst_var,
        #     abs(muon_eta),
        #     muon_pt
        #     )
        syst = evaluator["unc_pass"].evalv(
            muon_pt,
            muon_eta
            )
        variations["up"] = awkward.unflatten(syst+sf, n_muons)
        variations["down"] = awkward.unflatten(sf-syst, n_muons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (muons.pt < 5.0) | (muons.pt >= 500.0) | (abs(muons.eta) >= 2.4),
                awkward.ones_like(variations[var], dtype=float),
                variations[var]
        )

    return variations

MUON_ISO_SF_FILE = {
    "2016" : ["higgs_dna/systematics/data/2016postVFP_UL/hzg_muiso0p1_2016_efficiencies.json", "higgs_dna/systematics/data/2016postVFP_UL/hzg_muiso0p15_2016_efficiencies.json"],
    "2016preVFP" : ["higgs_dna/systematics/data/2016preVFP_UL/hzg_muiso0p1_2016APV_efficiencies.json", "higgs_dna/systematics/data/2016preVFP_UL/hzg_muiso0p15_2016APV_efficiencies.json"],
    "2016postVFP" : ["higgs_dna/systematics/data/2016postVFP_UL/hzg_muiso0p1_2016_efficiencies.json", "higgs_dna/systematics/data/2016postVFP_UL/hzg_muiso0p15_2016_efficiencies.json"],
    "2017" : ["higgs_dna/systematics/data/2017_UL/hzg_muiso0p1_2017_efficiencies.json", "higgs_dna/systematics/data/2017_UL/hzg_muiso0p15_2017_efficiencies.json"],
    "2018" : ["higgs_dna/systematics/data/2018_UL/hzg_muiso0p1_2018_efficiencies.json", "higgs_dna/systematics/data/2018_UL/hzg_muiso0p15_2018_efficiencies.json"],
    "2022preEE" : ["higgs_dna/systematics/data/2022preEE_UL/hzg_muiso0p1_2022_efficiencies.json", "higgs_dna/systematics/data/2022preEE_UL/hzg_muiso0p15_2022_efficiencies.json"],
    "2022postEE" : ["higgs_dna/systematics/data/2022postEE_UL/hzg_muiso0p1_2022EE_efficiencies.json", "higgs_dna/systematics/data/2022postEE_UL/hzg_muiso0p15_2022EE_efficiencies.json"],
    "2023preBPix" : ["higgs_dna/systematics/data/2023preBPix_UL/hzg_muiso0p1_2023_efficiencies.json", "higgs_dna/systematics/data/2023preBPix_UL/hzg_muiso0p15_2023_efficiencies.json"],
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/hzg_muiso0p1_2023BPix_efficiencies.json", "higgs_dna/systematics/data/2023postBPix_UL/hzg_muiso0p15_2023BPix_efficiencies.json"]
}

def muon_ISO_sf(events, year, central_only, input_collection):
    boundaries = [0.1, 0.15]
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "miniPFRelIso_all")
    ]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluators = []
    for file in MUON_ISO_SF_FILE[year]:
        evaluators.append(_core.CorrectionSet.from_file(misc_utils.expand_path(file)))
    
    muons = events[input_collection]
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    variations_list = []

    muon_eta = numpy.clip(
        numpy.abs(awkward.to_numpy(muons_flattened.eta)),
        0,
        2.39999 # SFs only valid up to eta 2.4
    )
    muon_pt = numpy.clip(
        awkward.to_numpy(muons_flattened.pt),
        5.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )

    for evaluator in evaluators:
        # Calculate SF and syst
        sub_variations = {}
        sf_data = evaluator["effdata"].evalv(
                muon_pt,
                muon_eta
        )
        sf_mc = evaluator["effmc"].evalv(
                muon_pt,
                muon_eta
        )
        sf = sf_data / sf_mc
        sub_variations["central"] = awkward.unflatten(sf, n_muons)

        if not central_only:
            syst_data = evaluator["systdata"].evalv(
                muon_pt,
                muon_eta
            )
            syst_mc = evaluator["systmc"].evalv(
                muon_pt,
                muon_eta
            )
            syst = sf * numpy.sqrt(
                (syst_data / sf_data)**2 + (syst_mc / sf_mc)**2
            )
            sub_variations["up"] = awkward.unflatten(syst+sf, n_muons)
            sub_variations["down"] = awkward.unflatten(sf-syst, n_muons)

        variations_list.append(sub_variations)

    evaluator = evaluators[1]
    sub_variations = {}
    sf_data = evaluator["effdata"].evalv(
            muon_pt,
            muon_eta
    )
    sf_mc = evaluator["effmc"].evalv(
            muon_pt,
            muon_eta
    )
    sf = (1-sf_data)/(1-sf_mc)
    sub_variations["central"] = awkward.unflatten(sf, n_muons)
    if not central_only:
        syst_data = evaluator["systdata"].evalv(
            muon_pt,
            muon_eta
        )
        syst_mc = evaluator["systmc"].evalv(
            muon_pt,
            muon_eta
        )
        syst = numpy.sqrt((syst_data/(1-sf_mc))**2 + (syst_mc/(1-sf_mc)*sf)**2)
        sub_variations["up"] = awkward.unflatten(syst+sf, n_muons)
        sub_variations["down"] = awkward.unflatten(sf-syst, n_muons)
    variations_list.append(sub_variations)

    variations = {}
    for var in variations_list[0].keys():
        variations[var] = awkward.where(
            (muons.pt < 5.0) | (muons.pt >= 500.0) | (abs(muons.eta) >= 2.4) | (muons.miniPFRelIso_all >= boundaries[-1]),
            awkward.ones_like(variations_list[0][var], dtype=float),
            awkward.where(muons.miniPFRelIso_all < boundaries[0], variations_list[0][var], awkward.where(muons.miniPFRelIso_all > boundaries[1], variations_list[2][var], variations_list[1][var]))
        )
    return variations


ELECTRON_ID_SF_FILE = {
    "2016" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_elid_2016_scalefactors.json",
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_elid_2016APV_scalefactors.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_elid_2016_scalefactors.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_elid_2017_scalefactors.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_elid_2018_scalefactors.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/hzg_elid_2022_scalefactors.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/hzg_elid_2022EE_scalefactors.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/hzg_elid_2023_scalefactors.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_elid_2023BPix_scalefactors.json"
}

ELECTRON_ID_SF_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hzg_elid_2023BPixHole_scalefactors.json"
}

def electron_WPL_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]
    
    # For 2023postBPix, we also need phi to determine which SF file to use
    if year == "2023postBPix":
        required_fields.append((input_collection, "phi"))

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(ELECTRON_ID_SF_FILE[year]))
    
    # For 2023postBPix, also load the Hole file
    evaluator_hole = None
    if year == "2023postBPix":
        evaluator_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(ELECTRON_ID_SF_FILE_HOLE[year]))

    electrons = events[input_collection]

    # Flatten electrons then convert to numpy for compatibility with correctionlib
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    ele_eta = numpy.clip(
        awkward.to_numpy(electrons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )

    ele_pt = numpy.clip(
        awkward.to_numpy(electrons_flattened.pt),
        7.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )
    
    # Get phi for 2023postBPix
    ele_phi = None
    if year == "2023postBPix":
        ele_phi = awkward.to_numpy(electrons_flattened.phi)

    # Calculate SF and syst
    variations = {}
    
    if year == "2023postBPix":
        # For 2023postBPix, use Hole file for electrons in the hole region
        # Hole region: -1.566 < eta < 0 and -1.2 < phi < -0.8
        hole_mask = (ele_eta > -1.566) & (ele_eta < 0) & (ele_phi > -1.2) & (ele_phi < -0.8)
        
        # Calculate SF for normal regions
        sf_normal = evaluator["sf_pass"].evalv(ele_pt, ele_eta)
        
        # Calculate SF for hole regions
        sf_hole = evaluator_hole["sf_pass"].evalv(ele_pt, ele_eta)
        
        # Combine SFs based on the hole mask
        sf = numpy.where(hole_mask, sf_hole, sf_normal)
        
        variations["central"] = awkward.unflatten(sf, n_electrons)

        if not central_only:
            # Calculate uncertainties for normal regions
            syst_normal = evaluator["unc_pass"].evalv(ele_pt, ele_eta)
            
            # Calculate uncertainties for hole regions
            syst_hole = evaluator_hole["unc_pass"].evalv(ele_pt, ele_eta)
            
            # Combine uncertainties based on the hole mask
            syst = numpy.where(hole_mask, syst_hole, syst_normal)
            
            variations["up"] = awkward.unflatten(syst+sf, n_electrons)
            variations["down"] = awkward.unflatten(sf-syst, n_electrons)
    else:
        # For other years, use the standard approach
        sf = evaluator["sf_pass"].evalv(ele_pt, ele_eta)
        variations["central"] = awkward.unflatten(sf, n_electrons)

        if not central_only:
            syst = evaluator["unc_pass"].evalv(ele_pt, ele_eta)
            variations["up"] = awkward.unflatten(syst+sf, n_electrons)
            variations["down"] = awkward.unflatten(sf-syst, n_electrons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (electrons.pt < 7.0) | (electrons.pt >= 500.0) | (abs(electrons.eta) >= 2.5),
                awkward.ones_like(variations[var], dtype=float),
                variations[var]
        )

    return variations

ELECTRON_ISO_SF_FILE = {
    "2016" : ["higgs_dna/systematics/data/2016postVFP_UL/hzg_eliso0p1_2016_efficiencies.json", "higgs_dna/systematics/data/2016postVFP_UL/hzg_eliso0p15_2016_efficiencies.json"],
    "2016preVFP" : ["higgs_dna/systematics/data/2016preVFP_UL/hzg_eliso0p1_2016APV_efficiencies.json", "higgs_dna/systematics/data/2016preVFP_UL/hzg_eliso0p15_2016APV_efficiencies.json"],
    "2016postVFP" : ["higgs_dna/systematics/data/2016postVFP_UL/hzg_eliso0p1_2016_efficiencies.json", "higgs_dna/systematics/data/2016postVFP_UL/hzg_eliso0p15_2016_efficiencies.json"],
    "2017" : ["higgs_dna/systematics/data/2017_UL/hzg_eliso0p1_2017_efficiencies.json", "higgs_dna/systematics/data/2017_UL/hzg_eliso0p15_2017_efficiencies.json"],
    "2018" : ["higgs_dna/systematics/data/2018_UL/hzg_eliso0p1_2018_efficiencies.json", "higgs_dna/systematics/data/2018_UL/hzg_eliso0p15_2018_efficiencies.json"],
    "2022preEE" : ["higgs_dna/systematics/data/2022preEE_UL/hzg_eliso0p1_2022_efficiencies.json", "higgs_dna/systematics/data/2022preEE_UL/hzg_eliso0p15_2022_efficiencies.json"],
    "2022postEE" : ["higgs_dna/systematics/data/2022postEE_UL/hzg_eliso0p1_2022EE_efficiencies.json", "higgs_dna/systematics/data/2022postEE_UL/hzg_eliso0p15_2022EE_efficiencies.json"],
    "2023preBPix" : ["higgs_dna/systematics/data/2023preBPix_UL/hzg_eliso0p1_2023_efficiencies.json", "higgs_dna/systematics/data/2023preBPix_UL/hzg_eliso0p15_2023_efficiencies.json"],
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/hzg_eliso0p1_2023BPix_efficiencies.json", "higgs_dna/systematics/data/2023postBPix_UL/hzg_eliso0p15_2023BPix_efficiencies.json"]
}

ELECTRON_ISO_SF_FILE_HOLE = {
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/hzg_eliso0p1_2023BPixHole_efficiencies.json", "higgs_dna/systematics/data/2023postBPix_UL/hzg_eliso0p15_2023BPixHole_efficiencies.json"]
}

def electron_ISO_sf(events, year, central_only, input_collection, is_data=False):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    boundaries = [0.1, 0.15]
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "miniPFRelIso_all")
    ]
    
    # For 2023postBPix, we also need phi to determine which SF file to use
    if year == "2023postBPix":
        required_fields.append((input_collection, "phi"))
        
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluators = []
    evaluators_hole = []
    
    for file in ELECTRON_ISO_SF_FILE[year]:
        evaluators.append(_core.CorrectionSet.from_file(misc_utils.expand_path(file)))
    
    # For 2023postBPix, also load the Hole files
    if year == "2023postBPix":
        for file in ELECTRON_ISO_SF_FILE_HOLE[year]:
            evaluators_hole.append(_core.CorrectionSet.from_file(misc_utils.expand_path(file)))
    
    electrons = events[input_collection]
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    variations_list = []

    ele_eta = numpy.clip(
        awkward.to_numpy(electrons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )
    ele_pt = numpy.clip(
        awkward.to_numpy(electrons_flattened.pt),
        7.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )
    
    # Get phi for 2023postBPix
    ele_phi = None
    hole_mask = None
    if year == "2023postBPix":
        ele_phi = awkward.to_numpy(electrons_flattened.phi)
        # Hole region: -1.566 < eta < 0 and -1.2 < phi < -0.8
        hole_mask = (ele_eta > -1.566) & (ele_eta < 0) & (ele_phi > -1.2) & (ele_phi < -0.8)

    for i, evaluator in enumerate(evaluators):
        # Calculate SF and syst
        sub_variations = {}
        
        if year == "2023postBPix":
            # For 2023postBPix, combine normal and hole regions
            evaluator_hole = evaluators_hole[i]
            
            # Calculate SF for normal regions
            sf_data_normal = evaluator["effdata"].evalv(ele_pt, ele_eta)
            sf_mc_normal = evaluator["effmc"].evalv(ele_pt, ele_eta)
            sf_normal = sf_data_normal / sf_mc_normal
            
            # Calculate SF for hole regions
            sf_data_hole = evaluator_hole["effdata"].evalv(ele_pt, ele_eta)
            sf_mc_hole = evaluator_hole["effmc"].evalv(ele_pt, ele_eta)
            sf_hole = sf_data_hole / sf_mc_hole
            
            # Combine SFs based on the hole mask
            sf = numpy.where(hole_mask, sf_hole, sf_normal)
            sub_variations["central"] = awkward.unflatten(sf, n_electrons)

            if not central_only:
                # Calculate uncertainties for normal regions
                syst_data_normal = evaluator["systdata"].evalv(ele_pt, ele_eta)
                syst_mc_normal = evaluator["systmc"].evalv(ele_pt, ele_eta)
                syst_normal = sf_normal * numpy.sqrt(
                    (syst_data_normal / sf_data_normal)**2 + (syst_mc_normal / sf_mc_normal)**2
                )
                
                # Calculate uncertainties for hole regions
                syst_data_hole = evaluator_hole["systdata"].evalv(ele_pt, ele_eta)
                syst_mc_hole = evaluator_hole["systmc"].evalv(ele_pt, ele_eta)
                syst_hole = sf_hole * numpy.sqrt(
                    (syst_data_hole / sf_data_hole)**2 + (syst_mc_hole / sf_mc_hole)**2
                )
                
                # Combine uncertainties based on the hole mask
                syst = numpy.where(hole_mask, syst_hole, syst_normal)
                
                sub_variations["up"] = awkward.unflatten(syst+sf, n_electrons)
                sub_variations["down"] = awkward.unflatten(sf-syst, n_electrons)
        else:
            # For other years, use the standard approach
            sf_data = evaluator["effdata"].evalv(ele_pt, ele_eta)
            sf_mc = evaluator["effmc"].evalv(ele_pt, ele_eta)
            sf = sf_data / sf_mc
            sub_variations["central"] = awkward.unflatten(sf, n_electrons)

            if not central_only:
                syst_data = evaluator["systdata"].evalv(ele_pt, ele_eta)
                syst_mc = evaluator["systmc"].evalv(ele_pt, ele_eta)
                syst = sf * numpy.sqrt(
                    (syst_data / sf_data)**2 + (syst_mc / sf_mc)**2
                )
                sub_variations["up"] = awkward.unflatten(syst+sf, n_electrons)
                sub_variations["down"] = awkward.unflatten(sf-syst, n_electrons)

        variations_list.append(sub_variations)
    
    # Handle the third evaluator (for isolation > 0.15)
    if year == "2023postBPix":
        evaluator = evaluators[1]
        evaluator_hole = evaluators_hole[1]
        
        sub_variations = {}
        
        # Calculate SF for normal regions
        sf_data_normal = evaluator["effdata"].evalv(ele_pt, ele_eta)
        sf_mc_normal = evaluator["effmc"].evalv(ele_pt, ele_eta)
        sf_normal = (1-sf_data_normal)/(1-sf_mc_normal)
        
        # Calculate SF for hole regions
        sf_data_hole = evaluator_hole["effdata"].evalv(ele_pt, ele_eta)
        sf_mc_hole = evaluator_hole["effmc"].evalv(ele_pt, ele_eta)
        sf_hole = (1-sf_data_hole)/(1-sf_mc_hole)
        
        # Combine SFs based on the hole mask
        sf = numpy.where(hole_mask, sf_hole, sf_normal)
        sub_variations["central"] = awkward.unflatten(sf, n_electrons)
        
        if not central_only:
            # Calculate uncertainties for normal regions
            syst_data_normal = evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_mc_normal = evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst_normal = numpy.sqrt((syst_data_normal/(1-sf_mc_normal))**2 + (syst_mc_normal/(1-sf_mc_normal)*sf_normal)**2)
            
            # Calculate uncertainties for hole regions
            syst_data_hole = evaluator_hole["systdata"].evalv(ele_pt, ele_eta)
            syst_mc_hole = evaluator_hole["systmc"].evalv(ele_pt, ele_eta)
            syst_hole = numpy.sqrt((syst_data_hole/(1-sf_mc_hole))**2 + (syst_mc_hole/(1-sf_mc_hole)*sf_hole)**2)
            
            # Combine uncertainties based on the hole mask
            syst = numpy.where(hole_mask, syst_hole, syst_normal)
            
            sub_variations["up"] = awkward.unflatten(syst+sf, n_electrons)
            sub_variations["down"] = awkward.unflatten(sf-syst, n_electrons)
    else:
        # For other years, use the standard approach
        evaluator = evaluators[1]

        sub_variations = {}
        sf_data = evaluator["effdata"].evalv(ele_pt, ele_eta)
        sf_mc = evaluator["effmc"].evalv(ele_pt, ele_eta)
        sf = (1-sf_data)/(1-sf_mc)
        sub_variations["central"] = awkward.unflatten(sf, n_electrons)

        if not central_only:
            syst_data = evaluator["systdata"].evalv(ele_pt, ele_eta)
            syst_mc = evaluator["systmc"].evalv(ele_pt, ele_eta)
            syst = numpy.sqrt((syst_data/(1-sf_mc))**2 + (syst_mc/(1-sf_mc)*sf)**2)
            sub_variations["up"] = awkward.unflatten(syst+sf, n_electrons)
            sub_variations["down"] = awkward.unflatten(sf-syst, n_electrons)
            
    variations_list.append(sub_variations)

    variations = {}
    for var in variations_list[0].keys():
        variations[var] = awkward.where(
            (electrons.pt < 7) | (electrons.pt >= 500.0) | (abs(electrons.eta) >= 2.5) | (electrons.miniPFRelIso_all >= boundaries[-1]),
            awkward.ones_like(variations_list[0][var], dtype=float),
            awkward.where(electrons.miniPFRelIso_all < boundaries[0], variations_list[0][var], awkward.where(electrons.miniPFRelIso_all > boundaries[1], variations_list[2][var], variations_list[1][var]))
        )

    return variations

ELECTRON_RECO_SF_FILE = {
    "2016" : "higgs_dna/systematics/data/2016preVFP_UL/electron_recoSF2016postVFP.json",
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/electron_recoSF2016preVFP.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/electron_recoSF2016postVFP.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/electron_recoSF2017.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/electron_recoSF2018.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/electron_recoSF2022.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/electron_recoSF2022EE.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/electron_recoSF2023.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/electron_recoSF2023BPix.json"
}

def electron_reco_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(ELECTRON_RECO_SF_FILE[year]))

    electrons = events[input_collection]

    # Flatten electrons then convert to numpy for compatibility with correctionlib
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    ele_eta = numpy.clip(
        awkward.to_numpy(electrons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )

    ele_pt = numpy.clip(
        awkward.to_numpy(electrons_flattened.pt),
        10.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    
    sf = evaluator["sf_pass"].evalv(
        ele_pt,
        ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator["unc_pass"].evalv(
            ele_pt,
            ele_eta
        )
        variations["up"] = awkward.unflatten(syst+sf, n_electrons)
        variations["down"] = awkward.unflatten(sf-syst, n_electrons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
            (electrons.pt < 10.0) | (electrons.pt >= 500.0) | (abs(electrons.eta) >= 2.5),
            awkward.ones_like(variations[var], dtype=float),
            variations[var]
        )

    return variations

#################################
### Electron Scale and Smear ###
#################################

electron_scale_FILE = {
    "2022preEE" : "jsonpog-integration/POG/EGM/2022_Summer22/electronSS_EtDependent.json",
    "2022postEE" : "jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS_EtDependent.json",
    "2023preBPix" : "jsonpog-integration/POG/EGM/2023_Summer23/electronSS_EtDependent.json",
    "2023postBPix" : "jsonpog-integration/POG/EGM/2023_Summer23BPix/electronSS_EtDependent.json"
}

electron_year_names = {
    "2022preEE" : "EGMSmearAndSyst_ElePTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_ElePTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_ElePTsplit_2023preBPIX", 
    "2023postBPix" : "EGMSmearAndSyst_ElePTsplit_2023postBPIX"
}

def electron_scale_smear_run3(events, year):
    """
    This function applies electron scale and smear corrections on the electron pt.
    It returns events with the corrected electron pt values.
    Based on the photon_scale_smear_run3 function but adapted for electrons.
    """
    logger.info("[Lepton Systematics] Applying electron scale and smear corrections for year %s", year)

    required_fields = [
        ("Electron", "eta"), ("Electron", "pt"), ("Electron", "r9"), ("Electron", "deltaEtaSC")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[lepton_systematics : electron_scale_smear_run3] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(electron_scale_FILE[year]))

    electrons = events["Electron"]
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    electrons_AbsEta = numpy.abs(awkward.to_numpy(electrons_flattened.eta+electrons_flattened.deltaEtaSC))
    electrons_pt = awkward.to_numpy(electrons_flattened.pt)
    electrons_r9 = awkward.to_numpy(electrons_flattened.r9)

    # Apply smear corrections first
    smear = evaluator[electron_year_names[year]].evalv("smear", electrons_pt, electrons_r9, electrons_AbsEta)
    rng = numpy.random.default_rng(seed=123)
    smear_val = awkward.where(
        (electrons_AbsEta > 3.0) | (electrons_pt < 20.0),
        awkward.ones_like(electrons_pt, dtype=float),
        rng.normal(loc=1., scale=numpy.abs(smear))
    )
    
    # Apply central smear correction to electron pt
    corrected_pt = awkward.to_numpy(electrons_pt * smear_val)

    # Calculate smear systematics
    for syst in ["smear_up", "smear_down"]:
        smear_syst = evaluator[electron_year_names[year]].evalv(syst, electrons_pt, electrons_r9, electrons_AbsEta)

        if "up" in syst:
            smear_up = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear+smear_syst)), n_electrons)
            events["Electron", "dEsigmaUp"] = electrons.pt * awkward.where(
                (abs(electrons.eta+electrons.deltaEtaSC) > 3.0) | (electrons.pt < 20.0),
                awkward.ones_like(smear_up, dtype=float),
                smear_up
            )
        elif "down" in syst:
            smear_down = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear-smear_syst)), n_electrons)
            events["Electron", "dEsigmaDown"] = electrons.pt * awkward.where(
                (abs(electrons.eta+electrons.deltaEtaSC) > 3.0) | (electrons.pt < 20.0),
                awkward.ones_like(smear_down, dtype=float),
                smear_down
            )

    print("Electron pt before scale correction:", electrons.pt)
    events["Electron", "corrected_pt"] = awkward.unflatten(corrected_pt, n_electrons)
    electrons = events["Electron"]

    # Calculate scale systematics
    for syst in ["scale_up", "scale_down"]:
        scale = evaluator[electron_year_names[year]].evalv(syst, corrected_pt, electrons_r9, electrons_AbsEta)
        
        if "up" in syst:
            scale_up = awkward.unflatten(scale, n_electrons)
            events["Electron", "dEscaleUp"] = electrons.corrected_pt * awkward.where(
                (abs(electrons.eta+electrons.deltaEtaSC) > 3.0) | (electrons.pt < 20.0),
                awkward.ones_like(scale_up, dtype=float),
                scale_up
            )
        elif "down" in syst:
            scale_down = awkward.unflatten(scale, n_electrons)
            events["Electron", "dEscaleDown"] = electrons.corrected_pt * awkward.where(
                (abs(electrons.eta+electrons.deltaEtaSC) > 3.0) | (electrons.pt < 20.0),
                awkward.ones_like(scale_down, dtype=float),
                scale_down
            )

    logger.info("[Lepton Systematics] Electron scale and smear corrections applied successfully")
    print("Electron pt after scale and smear corrections:", events["Electron", "corrected_pt"])
    print("Electron pt Scale Up:", events["Electron", "dEscaleUp"])
    print("Electron pt Scale Down:", events["Electron", "dEscaleDown"])
    print("Electron pt Smear Up:", events["Electron", "dEsigmaUp"])
    print("Electron pt Smear Down:", events["Electron", "dEsigmaDown"])

    return events
