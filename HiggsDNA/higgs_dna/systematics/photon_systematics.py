import awkward
import numpy

from correctionlib import _core

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.utils import systematic_from_bins, ic_systematic_from_bins

########################
##### Photon ID SF #####
########################

PHOTON_ID_SF_FILE = {
    "2016" : ["higgs_dna/systematics/data/2016postVFP_UL/photon.json", "higgs_dna/systematics/data/2016preVFP_UL/hzg_phidvalidate_2016_scalefactors.json"],
    "2016preVFP" : ["higgs_dna/systematics/data/2016preVFP_UL/photon.json", "higgs_dna/systematics/data/2016preVFP_UL/hzg_phidvalidate_2016APV_scalefactors.json"],
    "2016postVFP" : ["higgs_dna/systematics/data/2016postVFP_UL/photon.json", "higgs_dna/systematics/data/2016postVFP_UL/hzg_phidvalidate_2016_scalefactors.json"],
    "2017" : ["higgs_dna/systematics/data/2017_UL/photon.json", "higgs_dna/systematics/data/2017_UL/hzg_phidvalidate_2017_scalefactors.json"],
    "2018" : ["higgs_dna/systematics/data/2018_UL/photon.json", "higgs_dna/systematics/data/2018_UL/hzg_phidvalidate_2018_scalefactors.json"],
    "2022preEE" : ["higgs_dna/systematics/data/2022preEE_UL/photon.json", "higgs_dna/systematics/data/2022preEE_UL/hzg_phidvalidate_2022_scalefactors.json"],
    "2022postEE" : ["higgs_dna/systematics/data/2022postEE_UL/photon.json", "higgs_dna/systematics/data/2022postEE_UL/hzg_phidvalidate_2022EE_scalefactors.json"],
    "2023preBPix" : ["higgs_dna/systematics/data/2023preBPix_UL/photon.json"],
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/photon.json"],
}

PHOTON_ID_SF = {
    "2016" : "2016postVFP",
    "2016preVFP" : "2016preVFP",
    "2016postVFP" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018",
    "2022preEE" : "2022Re-recoBCD",
    "2022postEE" : "2022Re-recoE+PromptFG",
    "2023preBPix" : "2023PromptC",
    "2023postBPix" : "2023PromptD"
}

PHOTON_ID_EVAL = {
    "2016" : "UL-Photon-ID-SF",
    "2016preVFP" : "UL-Photon-ID-SF",
    "2016postVFP" : "UL-Photon-ID-SF",
    "2017" : "UL-Photon-ID-SF",
    "2018" : "UL-Photon-ID-SF",
    "2022preEE" : "Photon-ID-SF",
    "2022postEE" : "Photon-ID-SF",
    "2023preBPix" : "Photon-ID-SF",
    "2023postBPix" : "Photon-ID-SF"
}

def photon_id_sf(events, year, central_only, input_collection, working_point = "wp80"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "phi")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluators = []
    for file in PHOTON_ID_SF_FILE[year]:
        evaluators.append(_core.CorrectionSet.from_file(misc_utils.expand_path(file)))

    photons = events[input_collection]

    # Flatten photons then convert to numpy for compatibility with correctionlib
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    variations_list = []

    evaluator = evaluators[0]

    pho_eta = numpy.clip(
        awkward.to_numpy(photons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )

    pho_pt = numpy.clip(
        awkward.to_numpy(photons_flattened.pt),
        20.0, # SFs only valid for pT >= 20.0
        499.999 # and pT < 500.
    )

    pho_phi = numpy.clip(
        awkward.to_numpy(photons_flattened.phi),
        -999.0,
        999.0 # SFs only valid for phi in [-pi, pi]
    )

    variations = {}
    if int(year[:4]) < 2023:
        sf = evaluator[PHOTON_ID_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                working_point,
                pho_eta,
                pho_pt
        )
    else:
        # For 2023, need to handle pt < 20 photons separately
        pho_pt_orig = awkward.to_numpy(photons_flattened.pt)
        pt_below_20_mask = pho_pt_orig < 20.0
        
        # Create separate pt arrays for below/above 20 GeV
        pho_pt_below_20 = numpy.clip(pho_pt_orig, 15.0, 19.999)
        pho_pt_above_20 = numpy.clip(pho_pt_orig, 20.0, 499.999)
        
        # Calculate SF for pt >= 20
        sf_above_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                working_point,
                pho_eta,
                pho_pt_above_20,
                pho_phi
        )
        
        # Calculate SF for pt < 20 with modified working point
        sf_below_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                f"{working_point}Below20",
                pho_eta,
                pho_pt_below_20,
                pho_phi
        )
        
        # Combine results based on pt threshold
        sf = numpy.where(pt_below_20_mask, sf_below_20, sf_above_20)
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        syst_vars = ["sfup", "sfdown"] 
        for syst_var in syst_vars:
            if int(year[:4]) < 2023:
                syst = evaluator[PHOTON_ID_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        working_point,
                        pho_eta,
                        pho_pt
                )
            else:
                # For 2023, need to handle pt < 20 photons separately
                pho_pt_orig = awkward.to_numpy(photons_flattened.pt)
                pt_below_20_mask = pho_pt_orig < 20.0
                
                # Create separate pt arrays for below/above 20 GeV
                pho_pt_below_20 = numpy.clip(pho_pt_orig, 15.0, 19.999)
                pho_pt_above_20 = numpy.clip(pho_pt_orig, 20.0, 499.999)
                
                # Calculate syst for pt >= 20
                syst_above_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        working_point,
                        pho_eta,
                        pho_pt_above_20,
                        pho_phi
                )
                
                # Calculate syst for pt < 20 with modified working point
                syst_below_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        f"{working_point}Below20",
                        pho_eta,
                        pho_pt_below_20,
                        pho_phi
                )
                
                # Combine results based on pt threshold
                syst = numpy.where(pt_below_20_mask, syst_below_20, syst_above_20)
                
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_photons)
    variations_list.append(variations)

    if len(evaluators) == 2:
        evaluator = evaluators[1]

        pho_eta = numpy.clip(
            awkward.to_numpy(photons_flattened.eta),
            -2.49999,
            2.49999 # SFs only valid up to eta 2.5
        )

        pho_pt = numpy.clip(
            awkward.to_numpy(photons_flattened.pt),
            15.0, # SFs only valid for pT >= 15.0
            20.0,
        )

        variations = {}
        sf = evaluator["sf_pass"].evalv(
                pho_pt,
                pho_eta
        )
        variations["central"] = awkward.unflatten(sf, n_photons)

        if not central_only:
            syst = evaluator["unc_pass"].evalv(
                    pho_pt,
                    pho_eta
            )
            variations["up"] = awkward.unflatten(sf + syst, n_photons)
            variations["down"] = awkward.unflatten(sf - syst, n_photons)
        variations_list.append(variations)

    if len(variations_list) > 1:
        for var in variations_list[0].keys():
            variations[var] = awkward.where(
                (photons.pt < 15.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations_list[0][var], dtype=float),
                awkward.where(photons.pt < 20.0, variations_list[1][var], variations_list[0][var])
            )
    else:
        for var in variations_list[0].keys():
            variations[var] = awkward.where(
                (photons.pt < 15.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations_list[0][var], dtype=float),
                variations_list[0][var]
            )

    return variations

########################
#### Photon CSEV SF ####
########################

PHOTON_CSEV_EVAL = {
    "2016preVFP" : "UL-Photon-CSEV-SF",
    "2016postVFP" : "UL-Photon-CSEV-SF",
    "2017" : "UL-Photon-CSEV-SF",
    "2018" : "UL-Photon-CSEV-SF",
    "2022preEE" : "Photon-CSEV-SF",
    "2022postEE" : "Photon-CSEV-SF",
    "2023preBPix" : "Photon-CSEV-SF",
    "2023postBPix" : "Photon-CSEV-SF"
}

def photon_CSEV_sf(events, year, central_only, input_collection, working_point = "MVA"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """
    
    # working_point run3: MVA80 
    # working_point run2: MVA
    is_run2 = not (year.startswith("2022") or year.startswith("2023"))
    working_point = "MVA" if is_run2 else "MVA80"
    
    required_fields = [
        (input_collection, "eta"), (input_collection, "r9")
    ]
    
    # For Run2, also need isScEtaEB and isScEtaEE fields
    if int(year[:4]) < 2022:
        required_fields.extend([
            (input_collection, "isScEtaEB"), 
            (input_collection, "isScEtaEE")
        ])

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_ID_SF_FILE[year][0]))

    photons = events[input_collection]

    # Flatten photons then convert to numpy for compatibility with correctionlib
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    pho_eta = numpy.clip(
        awkward.to_numpy(photons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )

    pho_r9 = numpy.clip(
        awkward.to_numpy(photons_flattened.r9),
        0.0, # SFs only valid for R9 >= 0.0
        99999.0
    )

    # Calculate SF and syst
    variations = {}
    
    # Different treatment for Run2 vs Run3
    if int(year[:4]) < 2022:  # Run2: use CSEV bins with predefined values
        # Get CSEV SF values from the evaluator
        csev_correction = evaluator[PHOTON_CSEV_EVAL[year]]
        
        # Get SF values for each bin and working point
        sf_eb_high_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], "sf", working_point, "EBHighR9")
        sf_eb_low_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], "sf", working_point, "EBLowR9")
        sf_ee_high_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], "sf", working_point, "EEHighR9")
        sf_ee_low_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], "sf", working_point, "EELowR9")
        
        # Determine which bin each photon belongs to
        pho_isScEtaEB = photons_flattened.isScEtaEB
        pho_isScEtaEE = photons_flattened.isScEtaEE
        
        # Create masks for each bin (R9 boundary = 0.96)
        eb_high_r9_mask = pho_isScEtaEB & (photons_flattened.r9 > 0.96)
        eb_low_r9_mask = pho_isScEtaEB & (photons_flattened.r9 <= 0.96)
        ee_high_r9_mask = pho_isScEtaEE & (photons_flattened.r9 > 0.96)
        ee_low_r9_mask = pho_isScEtaEE & (photons_flattened.r9 <= 0.96)
        
        # Assign SF values using awkward.where
        sf = awkward.where(
            eb_high_r9_mask,
            sf_eb_high_r9,
            awkward.where(
                eb_low_r9_mask,
                sf_eb_low_r9,
                awkward.where(
                    ee_high_r9_mask,
                    sf_ee_high_r9,
                    awkward.where(
                        ee_low_r9_mask,
                        sf_ee_low_r9,
                        1.0  # Default value for photons not in any bin
                    )
                )
            )
        )
    else:  # Run3: use eta and R9 directly
        sf = evaluator[PHOTON_CSEV_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                working_point,
                pho_eta,
                pho_r9
        )
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        syst_vars = ["sfup", "sfdown"] 
        for syst_var in syst_vars:
            if int(year[:4]) < 2022:  # Run2: use CSEV bins with predefined values
                # Get systematic SF values for each bin
                syst_eb_high_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], syst_var, working_point, "EBHighR9")
                syst_eb_low_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], syst_var, working_point, "EBLowR9")
                syst_ee_high_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], syst_var, working_point, "EEHighR9")
                syst_ee_low_r9 = csev_correction.evaluate(PHOTON_ID_SF[year], syst_var, working_point, "EELowR9")
                
                # Assign systematic SF values using awkward.where
                syst = awkward.where(
                    eb_high_r9_mask,
                    syst_eb_high_r9,
                    awkward.where(
                        eb_low_r9_mask,
                        syst_eb_low_r9,
                        awkward.where(
                            ee_high_r9_mask,
                            syst_ee_high_r9,
                            awkward.where(
                                ee_low_r9_mask,
                                syst_ee_low_r9,
                                1.0  # Default value
                            )
                        )
                    )
                )
            else:  # Run3: use eta and R9 directly
                syst = evaluator[PHOTON_CSEV_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        working_point,
                        pho_eta,
                        pho_r9
                )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_photons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations[var]),
                variations[var]
        )

    return variations

########################
##### Fake Photon SF ###
########################

PHOTON_Fake_Photon_SF_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2016postVFP" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2017" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2018" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2022preEE" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2022postEE" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2023preBPix" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2023postBPix" : "higgs_dna/systematics/data/fakephoton/fakephoton.json"
}

def photon_fake_photon_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "genPartFlav")
    ]
    
    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning(f"Missing fields for photon_fake_photon_sf: {missing_fields}")

    is_run2 = not (year.startswith("2022") or year.startswith("2023"))
    run_str = "run2" if is_run2 else "run3"
    
    photons = events[input_collection]
    n_photons = awkward.num(photons)
    
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_Fake_Photon_SF_FILE[year]))
    photons_flattened = awkward.flatten(photons)

    photons_abseta = numpy.clip(
            numpy.abs(awkward.to_numpy(photons_flattened.eta)),
            0.0,
            2.4999
    )
    photons_pt = numpy.clip(
            awkward.to_numpy(photons_flattened.pt),
            15.0,
            499.999
    )

    isjet = awkward.where(photons_flattened.genPartFlav == 1, 0, 1)    
    photons_isjet = awkward.to_numpy(isjet).astype(int)
    
    sf = evaluator["fakephoton_corrections"].evalv(
        run_str,
        photons_isjet,
        photons_pt,
        photons_abseta
    )

    variations = {}
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        # No systematics for fake photon SF, so up/down are the same as central
        variations["up"] = variations["central"]
        variations["down"] = variations["central"]

    for var in variations.keys():
        variations[var] = awkward.where(
            (photons.pt < 15.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
            awkward.ones_like(variations[var], dtype=float),
            variations[var]
        )

    return variations

########################
### Electron veto SF ###
########################

from higgs_dna.systematics.data.electron_veto_sf import PHOTON_ELECTRON_VETO_SF_2016, PHOTON_ELECTRON_VETO_SF_2017, PHOTON_ELECTRON_VETO_SF_2018
photon_electron_veto_sf_bins = {
    "2016" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016preVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016postVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2017" : PHOTON_ELECTRON_VETO_SF_2017,
    "2018" : PHOTON_ELECTRON_VETO_SF_2018,
    "2022preEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2022postEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023preBPix" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023postBPix" : PHOTON_ELECTRON_VETO_SF_2018 #FIXME
}

def photon_electron_veto_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins = photon_electron_veto_sf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        central_only = central_only
    )
    print(variations)
    return variations

##################
### Trigger SF ###
##################
# Note: since the trigger sf applies separate SF for the lead/sublead photons,
# it's easiest to just cast this as an EventWeightSystematic (rather than ObjectWeightSystematic as we would typically do)
# and just multiply the lead/sublead variations manually by hand

from higgs_dna.systematics.data.trigger_sf import LEAD_TRIGGER_SF_2016, SUBLEAD_TRIGGER_SF_2016, LEAD_TRIGGER_SF_2017, SUBLEAD_TRIGGER_SF_2017, LEAD_TRIGGER_SF_2018, SUBLEAD_TRIGGER_SF_2018 
lead_trigger_sf_bins = {
    "2016" : LEAD_TRIGGER_SF_2016,
    "2016preVFP" : LEAD_TRIGGER_SF_2016,
    "2016postVFP" : LEAD_TRIGGER_SF_2016,
    "2017" : LEAD_TRIGGER_SF_2017,
    "2018" : LEAD_TRIGGER_SF_2018,
    "2022preEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : LEAD_TRIGGER_SF_2018 #FIXME
}
sublead_trigger_sf_bins = {
    "2016" : SUBLEAD_TRIGGER_SF_2016,
    "2016preVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2016postVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2017" : SUBLEAD_TRIGGER_SF_2017,
    "2018" : SUBLEAD_TRIGGER_SF_2018, 
    "2022preEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : SUBLEAD_TRIGGER_SF_2018 #FIXME
}

def trigger_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9"), ("Photon", "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations_lead = systematic_from_bins(
        bins = lead_trigger_sf_bins[year], 
        variables = {
            "photon_r9" : events.LeadPhoton.r9,
            "photon_eta" : abs(events.LeadPhoton.eta),
            "photon_pt" : events.LeadPhoton.pt
        },
        central_only = central_only
    )

    variations_sublead = systematic_from_bins(
        bins = sublead_trigger_sf_bins[year],
        variables = {
            "photon_r9" : events.SubleadPhoton.r9,
            "photon_eta" : abs(events.SubleadPhoton.eta),
            "photon_pt" : events.SubleadPhoton.pt
        },
        central_only = central_only
    )
 
    # Multiply up/down/central variations together, following this treatment in flashgg:
    # https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Systematics/interface/DiPhotonFromSinglePhotonViewBase.h#L85-L87
    variations = {}
    for key in variations_lead.keys():
        variations[key] = variations_lead[key] * variations_sublead[key]

    return variations

############
### FNUF ###
############

from higgs_dna.systematics.data.fnuf import FNUF_2016, FNUF_2017, FNUF_2018
fnuf_bins = {
    "2016" : FNUF_2016,
    "2016preVFP" : FNUF_2016,
    "2016postVFP" : FNUF_2016,
    "2017" : FNUF_2017,
    "2018" : FNUF_2018,
    "2022preEE" : FNUF_2018,
    "2022postEE" : FNUF_2018,
    "2023preBPix" : FNUF_2018,
    "2023postBPix" : FNUF_2018,
}

def fnuf_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    """

    """
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message) 

    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "eb":
        mask = photons.isScEtaEB == True
    elif loc == "ee":
        mask = photons.isScEtaEE == True

    variations = ic_systematic_from_bins(
        bins = fnuf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations

################
### Material ###
################

from higgs_dna.systematics.data.material import MATERIAL_2016, MATERIAL_2017, MATERIAL_2018
material_bins = {
    "2016" : MATERIAL_2016,
    "2016preVFP" : MATERIAL_2016,
    "2016postVFP" : MATERIAL_2016,
    "2017" : MATERIAL_2017,
    "2018" : MATERIAL_2018,
    "2022preEE" : MATERIAL_2018,
    "2022postEE" : MATERIAL_2018,
    "2023preBPix" : MATERIAL_2018,
    "2023postBPix" : MATERIAL_2018
}

def material_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "central_barrel":
        mask = abs(photons.eta) <= 1.0
    elif loc == "outer_barrel":
        mask = (abs(photons.eta) > 1.0) & (abs(photons.eta) <= 1.5)
    elif loc == "forward":
        mask = abs(photons.eta) > 1.5
    
    variations = ic_systematic_from_bins(
        bins = material_bins[year],
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations  


def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations


#######################
### Photon MC Smear ###
#######################

def photon_mc_smear(events, r9, loc):
    photons = events.Photon

    # Split into high/low r9 and EE/EB
    if r9 == "high":
        r9_cut = photons.r9 > 0.94
    elif r9 == "low":
        r9_cut = photons.r9 <= 0.94

    if loc == "eb":
        loc_cut = photons.isScEtaEB == True
    elif loc == "ee":
        loc_cut = photons.isScEtaEE == True

    mask = r9_cut & loc_cut

    variations = {}
    variations["up"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaUp,
            photons.pt
    )
    variations["down"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaDown,
            photons.pt
    )

    return variations

photon_scale_FILE = {"2022preEE" : "jsonpog-integration/POG/EGM/2022_Summer22/photonSS_EtDependent.json",
    "2022postEE" : "jsonpog-integration/POG/EGM/2022_Summer22EE/photonSS_EtDependent.json",
    "2023preBPix" : "jsonpog-integration/POG/EGM/2023_Summer23/photonSS_EtDependent.json",
    "2023postBPix" : "jsonpog-integration/POG/EGM/2023_Summer23BPix/photonSS_EtDependent.json"
}

year_names = {
    "2022preEE" : "EGMSmearAndSyst_PhoPTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_PhoPTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_PhoPTsplit_2023preBPIX",
    "2023postBPix" : "EGMSmearAndSyst_PhoPTsplit_2023postBPIX"
}

def photon_scale_smear_run3(events, year):
    """
    This function applies photon scale and smear corrections on the photon pt.
    It returns events with the corrected photon pt values.
    """
    logger.info("[Photon Systematics] Applying photon scale and smear corrections for year %s", year)

    required_fields = [
        ("Photon", "x_calo"), ("Photon", "y_calo"), ("Photon", "z_calo"), ("Photon", "eta"), ("Photon", "pt"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))

    photons = events["Photon"]
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    # Calculate eta from calo coordinates
    # rho = numpy.sqrt(awkward.to_numpy(photons_flattened.x_calo)**2 + awkward.to_numpy(photons_flattened.y_calo)**2)
    # theta = numpy.arctan(rho / awkward.to_numpy(photons_flattened.z_calo))
    photons_AbsEta = numpy.abs(awkward.to_numpy(photons_flattened.eta))

    photons_pt = awkward.to_numpy(photons_flattened.pt)
    photons_r9 = awkward.to_numpy(photons_flattened.r9)

    # Apply smear corrections first
    smear = evaluator[year_names[year]].evalv("smear", photons_pt, photons_r9, photons_AbsEta)
    rng = numpy.random.default_rng(seed=123)
    smear_val = awkward.where(
        (photons_AbsEta > 3.0) | (photons_pt < 20.0),
        awkward.ones_like(photons_pt, dtype=float),
        rng.normal(loc=1., scale=numpy.abs(smear))
    )
    
    # Apply central smear correction to photon pt
    corrected_pt = awkward.to_numpy(photons_pt * smear_val)

    # Calculate smear systematics
    for syst in ["smear_up", "smear_down"]:
        smear_syst = evaluator[year_names[year]].evalv(syst, photons_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            smear_up = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear+smear_syst)), n_photons)
            events["Photon", "dEsigmaUp"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_up, dtype=float),
                smear_up
            )
        elif "down" in syst:
            smear_down = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear-smear_syst)), n_photons)
            events["Photon", "dEsigmaDown"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_down, dtype=float),
                smear_down
            )
    
    print("Photon pt before scale and smear corrections:", photons.pt)
    events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
    photons = events["Photon"]

    # Calculate scale systematics
    for syst in ["scale_up", "scale_down"]:
        scale = evaluator[year_names[year]].evalv(syst, corrected_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            scale_up = awkward.unflatten(scale, n_photons)
            events["Photon", "dEscaleUp"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_up, dtype=float),
                scale_up
            )
        elif "down" in syst:
            scale_down = awkward.unflatten(scale, n_photons)
            events["Photon", "dEscaleDown"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_down, dtype=float),
                scale_down
            )

    print("Photon pt after scale and smear corrections:", events["Photon", "corrected_pt"])
    print("Photon pt Scale Up:", events["Photon", "dEscaleUp"])
    print("Photon pt Scale Down:", events["Photon", "dEscaleDown"])
    print("Photon pt Smear Up:", events["Photon", "dEsigmaUp"])
    print("Photon pt Smear Down:", events["Photon", "dEsigmaDown"])

    return events

############################
### Electron veto SF ###
############################

from higgs_dna.systematics.data.electron_veto_sf import PHOTON_ELECTRON_VETO_SF_2016, PHOTON_ELECTRON_VETO_SF_2017, PHOTON_ELECTRON_VETO_SF_2018
photon_electron_veto_sf_bins = {
    "2016" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016preVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016postVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2017" : PHOTON_ELECTRON_VETO_SF_2017,
    "2018" : PHOTON_ELECTRON_VETO_SF_2018,
    "2022preEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2022postEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023preBPix" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023postBPix" : PHOTON_ELECTRON_VETO_SF_2018 #FIXME
}

def photon_electron_veto_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins = photon_electron_veto_sf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        central_only = central_only
    )
    print(variations)
    return variations

##################
### Trigger SF ###
##################
# Note: since the trigger sf applies separate SF for the lead/sublead photons,
# it's easiest to just cast this as an EventWeightSystematic (rather than ObjectWeightSystematic as we would typically do)
# and just multiply the lead/sublead variations manually by hand

from higgs_dna.systematics.data.trigger_sf import LEAD_TRIGGER_SF_2016, SUBLEAD_TRIGGER_SF_2016, LEAD_TRIGGER_SF_2017, SUBLEAD_TRIGGER_SF_2017, LEAD_TRIGGER_SF_2018, SUBLEAD_TRIGGER_SF_2018 
lead_trigger_sf_bins = {
    "2016" : LEAD_TRIGGER_SF_2016,
    "2016preVFP" : LEAD_TRIGGER_SF_2016,
    "2016postVFP" : LEAD_TRIGGER_SF_2016,
    "2017" : LEAD_TRIGGER_SF_2017,
    "2018" : LEAD_TRIGGER_SF_2018,
    "2022preEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : LEAD_TRIGGER_SF_2018 #FIXME
}
sublead_trigger_sf_bins = {
    "2016" : SUBLEAD_TRIGGER_SF_2016,
    "2016preVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2016postVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2017" : SUBLEAD_TRIGGER_SF_2017,
    "2018" : SUBLEAD_TRIGGER_SF_2018, 
    "2022preEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : SUBLEAD_TRIGGER_SF_2018 #FIXME
}

def trigger_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9"), ("Photon", "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations_lead = systematic_from_bins(
        bins = lead_trigger_sf_bins[year], 
        variables = {
            "photon_r9" : events.LeadPhoton.r9,
            "photon_eta" : abs(events.LeadPhoton.eta),
            "photon_pt" : events.LeadPhoton.pt
        },
        central_only = central_only
    )

    variations_sublead = systematic_from_bins(
        bins = sublead_trigger_sf_bins[year],
        variables = {
            "photon_r9" : events.SubleadPhoton.r9,
            "photon_eta" : abs(events.SubleadPhoton.eta),
            "photon_pt" : events.SubleadPhoton.pt
        },
        central_only = central_only
    )
 
    # Multiply up/down/central variations together, following this treatment in flashgg:
    # https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Systematics/interface/DiPhotonFromSinglePhotonViewBase.h#L85-L87
    variations = {}
    for key in variations_lead.keys():
        variations[key] = variations_lead[key] * variations_sublead[key]

    return variations

############
### FNUF ###
############

from higgs_dna.systematics.data.fnuf import FNUF_2016, FNUF_2017, FNUF_2018
fnuf_bins = {
    "2016" : FNUF_2016,
    "2016preVFP" : FNUF_2016,
    "2016postVFP" : FNUF_2016,
    "2017" : FNUF_2017,
    "2018" : FNUF_2018,
    "2022preEE" : FNUF_2018,
    "2022postEE" : FNUF_2018,
    "2023preBPix" : FNUF_2018,
    "2023postBPix" : FNUF_2018,
}

def fnuf_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    """

    """
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message) 

    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "eb":
        mask = photons.isScEtaEB == True
    elif loc == "ee":
        mask = photons.isScEtaEE == True

    variations = ic_systematic_from_bins(
        bins = fnuf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations

################
### Material ###
################

from higgs_dna.systematics.data.material import MATERIAL_2016, MATERIAL_2017, MATERIAL_2018
material_bins = {
    "2016" : MATERIAL_2016,
    "2016preVFP" : MATERIAL_2016,
    "2016postVFP" : MATERIAL_2016,
    "2017" : MATERIAL_2017,
    "2018" : MATERIAL_2018,
    "2022preEE" : MATERIAL_2018,
    "2022postEE" : MATERIAL_2018,
    "2023preBPix" : MATERIAL_2018,
    "2023postBPix" : MATERIAL_2018
}

def material_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "central_barrel":
        mask = abs(photons.eta) <= 1.0
    elif loc == "outer_barrel":
        mask = (abs(photons.eta) > 1.0) & (abs(photons.eta) <= 1.5)
    elif loc == "forward":
        mask = abs(photons.eta) > 1.5
    
    variations = ic_systematic_from_bins(
        bins = material_bins[year],
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations  


def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations


#######################
### Photon MC Smear ###
#######################

def photon_mc_smear(events, r9, loc):
    photons = events.Photon

    # Split into high/low r9 and EE/EB
    if r9 == "high":
        r9_cut = photons.r9 > 0.94
    elif r9 == "low":
        r9_cut = photons.r9 <= 0.94

    if loc == "eb":
        loc_cut = photons.isScEtaEB == True
    elif loc == "ee":
        loc_cut = photons.isScEtaEE == True

    mask = r9_cut & loc_cut

    variations = {}
    variations["up"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaUp,
            photons.pt
    )
    variations["down"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaDown,
            photons.pt
    )

    return variations

photon_scale_FILE = {"2022preEE" : "jsonpog-integration/POG/EGM/2022_Summer22/photonSS_EtDependent.json",
    "2022postEE" : "jsonpog-integration/POG/EGM/2022_Summer22EE/photonSS_EtDependent.json",
    "2023preBPix" : "jsonpog-integration/POG/EGM/2023_Summer23/photonSS_EtDependent.json",
    "2023postBPix" : "jsonpog-integration/POG/EGM/2023_Summer23BPix/photonSS_EtDependent.json"
}

year_names = {
    "2022preEE" : "EGMSmearAndSyst_PhoPTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_PhoPTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_PhoPTsplit_2023preBPIX",
    "2023postBPix" : "EGMSmearAndSyst_PhoPTsplit_2023postBPIX"
}

def photon_scale_smear_run3(events, year):
    """
    This function applies photon scale and smear corrections on the photon pt.
    It returns events with the corrected photon pt values.
    """
    logger.info("[Photon Systematics] Applying photon scale and smear corrections for year %s", year)

    required_fields = [
        ("Photon", "x_calo"), ("Photon", "y_calo"), ("Photon", "z_calo"), ("Photon", "eta"), ("Photon", "pt"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))

    photons = events["Photon"]
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    # Calculate eta from calo coordinates
    # rho = numpy.sqrt(awkward.to_numpy(photons_flattened.x_calo)**2 + awkward.to_numpy(photons_flattened.y_calo)**2)
    # theta = numpy.arctan(rho / awkward.to_numpy(photons_flattened.z_calo))
    photons_AbsEta = numpy.abs(awkward.to_numpy(photons_flattened.eta))

    photons_pt = awkward.to_numpy(photons_flattened.pt)
    photons_r9 = awkward.to_numpy(photons_flattened.r9)

    # Apply smear corrections first
    smear = evaluator[year_names[year]].evalv("smear", photons_pt, photons_r9, photons_AbsEta)
    rng = numpy.random.default_rng(seed=123)
    smear_val = awkward.where(
        (photons_AbsEta > 3.0) | (photons_pt < 20.0),
        awkward.ones_like(photons_pt, dtype=float),
        rng.normal(loc=1., scale=numpy.abs(smear))
    )
    
    # Apply central smear correction to photon pt
    corrected_pt = awkward.to_numpy(photons_pt * smear_val)

    # Calculate smear systematics
    for syst in ["smear_up", "smear_down"]:
        smear_syst = evaluator[year_names[year]].evalv(syst, photons_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            smear_up = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear+smear_syst)), n_photons)
            events["Photon", "dEsigmaUp"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_up, dtype=float),
                smear_up
            )
        elif "down" in syst:
            smear_down = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear-smear_syst)), n_photons)
            events["Photon", "dEsigmaDown"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_down, dtype=float),
                smear_down
            )
    
    print("Photon pt before scale and smear corrections:", photons.pt)
    events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
    photons = events["Photon"]

    # Calculate scale systematics
    for syst in ["scale_up", "scale_down"]:
        scale = evaluator[year_names[year]].evalv(syst, corrected_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            scale_up = awkward.unflatten(scale, n_photons)
            events["Photon", "dEscaleUp"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_up, dtype=float),
                scale_up
            )
        elif "down" in syst:
            scale_down = awkward.unflatten(scale, n_photons)
            events["Photon", "dEscaleDown"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_down, dtype=float),
                scale_down
            )

    print("Photon pt after scale and smear corrections:", events["Photon", "corrected_pt"])
    print("Photon pt Scale Up:", events["Photon", "dEscaleUp"])
    print("Photon pt Scale Down:", events["Photon", "dEscaleDown"])
    print("Photon pt Smear Up:", events["Photon", "dEsigmaUp"])
    print("Photon pt Smear Down:", events["Photon", "dEsigmaDown"])

    return events

###############################
### Photon Electron Veto SF ### 
###############################

from higgs_dna.systematics.data.electron_veto_sf import PHOTON_ELECTRON_VETO_SF_2016, PHOTON_ELECTRON_VETO_SF_2017, PHOTON_ELECTRON_VETO_SF_2018
photon_electron_veto_sf_bins = {
    "2016" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016preVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016postVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2017" : PHOTON_ELECTRON_VETO_SF_2017,
    "2018" : PHOTON_ELECTRON_VETO_SF_2018,
    "2022preEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2022postEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023preBPix" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023postBPix" : PHOTON_ELECTRON_VETO_SF_2018 #FIXME
}

def photon_electron_veto_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins = photon_electron_veto_sf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        central_only = central_only
    )
    print(variations)
    return variations

##################
### Trigger SF ###
##################
# Note: since the trigger sf applies separate SF for the lead/sublead photons,
# it's easiest to just cast this as an EventWeightSystematic (rather than ObjectWeightSystematic as we would typically do)
# and just multiply the lead/sublead variations manually by hand

from higgs_dna.systematics.data.trigger_sf import LEAD_TRIGGER_SF_2016, SUBLEAD_TRIGGER_SF_2016, LEAD_TRIGGER_SF_2017, SUBLEAD_TRIGGER_SF_2017, LEAD_TRIGGER_SF_2018, SUBLEAD_TRIGGER_SF_2018 
lead_trigger_sf_bins = {
    "2016" : LEAD_TRIGGER_SF_2016,
    "2016preVFP" : LEAD_TRIGGER_SF_2016,
    "2016postVFP" : LEAD_TRIGGER_SF_2016,
    "2017" : LEAD_TRIGGER_SF_2017,
    "2018" : LEAD_TRIGGER_SF_2018,
    "2022preEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : LEAD_TRIGGER_SF_2018 #FIXME
}
sublead_trigger_sf_bins = {
    "2016" : SUBLEAD_TRIGGER_SF_2016,
    "2016preVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2016postVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2017" : SUBLEAD_TRIGGER_SF_2017,
    "2018" : SUBLEAD_TRIGGER_SF_2018, 
    "2022preEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : SUBLEAD_TRIGGER_SF_2018 #FIXME
}

def trigger_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9"), ("Photon", "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations_lead = systematic_from_bins(
        bins = lead_trigger_sf_bins[year], 
        variables = {
            "photon_r9" : events.LeadPhoton.r9,
            "photon_eta" : abs(events.LeadPhoton.eta),
            "photon_pt" : events.LeadPhoton.pt
        },
        central_only = central_only
    )

    variations_sublead = systematic_from_bins(
        bins = sublead_trigger_sf_bins[year],
        variables = {
            "photon_r9" : events.SubleadPhoton.r9,
            "photon_eta" : abs(events.SubleadPhoton.eta),
            "photon_pt" : events.SubleadPhoton.pt
        },
        central_only = central_only
    )
 
    # Multiply up/down/central variations together, following this treatment in flashgg:
    # https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Systematics/interface/DiPhotonFromSinglePhotonViewBase.h#L85-L87
    variations = {}
    for key in variations_lead.keys():
        variations[key] = variations_lead[key] * variations_sublead[key]

    return variations

############
### FNUF ###
############

from higgs_dna.systematics.data.fnuf import FNUF_2016, FNUF_2017, FNUF_2018
fnuf_bins = {
    "2016" : FNUF_2016,
    "2016preVFP" : FNUF_2016,
    "2016postVFP" : FNUF_2016,
    "2017" : FNUF_2017,
    "2018" : FNUF_2018,
    "2022preEE" : FNUF_2018,
    "2022postEE" : FNUF_2018,
    "2023preBPix" : FNUF_2018,
    "2023postBPix" : FNUF_2018,
}

def fnuf_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    """

    """
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message) 

    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "eb":
        mask = photons.isScEtaEB == True
    elif loc == "ee":
        mask = photons.isScEtaEE == True

    variations = ic_systematic_from_bins(
        bins = fnuf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations

################
### Material ###
################

from higgs_dna.systematics.data.material import MATERIAL_2016, MATERIAL_2017, MATERIAL_2018
material_bins = {
    "2016" : MATERIAL_2016,
    "2016preVFP" : MATERIAL_2016,
    "2016postVFP" : MATERIAL_2016,
    "2017" : MATERIAL_2017,
    "2018" : MATERIAL_2018,
    "2022preEE" : MATERIAL_2018,
    "2022postEE" : MATERIAL_2018,
    "2023preBPix" : MATERIAL_2018,
    "2023postBPix" : MATERIAL_2018
}

def material_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "central_barrel":
        mask = abs(photons.eta) <= 1.0
    elif loc == "outer_barrel":
        mask = (abs(photons.eta) > 1.0) & (abs(photons.eta) <= 1.5)
    elif loc == "forward":
        mask = abs(photons.eta) > 1.5
    
    variations = ic_systematic_from_bins(
        bins = material_bins[year],
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations  


def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations


#######################
### Photon MC Smear ###
#######################

def photon_mc_smear(events, r9, loc):
    photons = events.Photon

    # Split into high/low r9 and EE/EB
    if r9 == "high":
        r9_cut = photons.r9 > 0.94
    elif r9 == "low":
        r9_cut = photons.r9 <= 0.94

    if loc == "eb":
        loc_cut = photons.isScEtaEB == True
    elif loc == "ee":
        loc_cut = photons.isScEtaEE == True

    mask = r9_cut & loc_cut

    variations = {}
    variations["up"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaUp,
            photons.pt
    )
    variations["down"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaDown,
            photons.pt
    )

    return variations

photon_scale_FILE = {"2022preEE" : "jsonpog-integration/POG/EGM/2022_Summer22/photonSS_EtDependent.json",
    "2022postEE" : "jsonpog-integration/POG/EGM/2022_Summer22EE/photonSS_EtDependent.json",
    "2023preBPix" : "jsonpog-integration/POG/EGM/2023_Summer23/photonSS_EtDependent.json",
    "2023postBPix" : "jsonpog-integration/POG/EGM/2023_Summer23BPix/photonSS_EtDependent.json"
}

year_names = {
    "2022preEE" : "EGMSmearAndSyst_PhoPTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_PhoPTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_PhoPTsplit_2023preBPIX",
    "2023postBPix" : "EGMSmearAndSyst_PhoPTsplit_2023postBPIX"
}

def photon_scale_smear_run3(events, year):
    """
    This function applies photon scale and smear corrections on the photon pt.
    It returns events with the corrected photon pt values.
    """
    logger.info("[Photon Systematics] Applying photon scale and smear corrections for year %s", year)

    required_fields = [
        ("Photon", "x_calo"), ("Photon", "y_calo"), ("Photon", "z_calo"), ("Photon", "eta"), ("Photon", "pt"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))

    photons = events["Photon"]
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    # Calculate eta from calo coordinates
    # rho = numpy.sqrt(awkward.to_numpy(photons_flattened.x_calo)**2 + awkward.to_numpy(photons_flattened.y_calo)**2)
    # theta = numpy.arctan(rho / awkward.to_numpy(photons_flattened.z_calo))
    photons_AbsEta = numpy.abs(awkward.to_numpy(photons_flattened.eta))

    photons_pt = awkward.to_numpy(photons_flattened.pt)
    photons_r9 = awkward.to_numpy(photons_flattened.r9)

    # Apply smear corrections first
    smear = evaluator[year_names[year]].evalv("smear", photons_pt, photons_r9, photons_AbsEta)
    rng = numpy.random.default_rng(seed=123)
    smear_val = awkward.where(
        (photons_AbsEta > 3.0) | (photons_pt < 20.0),
        awkward.ones_like(photons_pt, dtype=float),
        rng.normal(loc=1., scale=numpy.abs(smear))
    )
    
    # Apply central smear correction to photon pt
    corrected_pt = awkward.to_numpy(photons_pt * smear_val)

    # Calculate smear systematics
    for syst in ["smear_up", "smear_down"]:
        smear_syst = evaluator[year_names[year]].evalv(syst, photons_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            smear_up = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear+smear_syst)), n_photons)
            events["Photon", "dEsigmaUp"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_up, dtype=float),
                smear_up
            )
        elif "down" in syst:
            smear_down = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear-smear_syst)), n_photons)
            events["Photon", "dEsigmaDown"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_down, dtype=float),
                smear_down
            )
    
    print("Photon pt before scale and smear corrections:", photons.pt)
    events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
    photons = events["Photon"]

    # Calculate scale systematics
    for syst in ["scale_up", "scale_down"]:
        scale = evaluator[year_names[year]].evalv(syst, corrected_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            scale_up = awkward.unflatten(scale, n_photons)
            events["Photon", "dEscaleUp"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_up, dtype=float),
                scale_up
            )
        elif "down" in syst:
            scale_down = awkward.unflatten(scale, n_photons)
            events["Photon", "dEscaleDown"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_down, dtype=float),
                scale_down
            )

    print("Photon pt after scale and smear corrections:", events["Photon", "corrected_pt"])
    print("Photon pt Scale Up:", events["Photon", "dEscaleUp"])
    print("Photon pt Scale Down:", events["Photon", "dEscaleDown"])
    print("Photon pt Smear Up:", events["Photon", "dEsigmaUp"])
    print("Photon pt Smear Down:", events["Photon", "dEsigmaDown"])

    return events

#######################
###Electron veto SF ###
########################

from higgs_dna.systematics.data.electron_veto_sf import PHOTON_ELECTRON_VETO_SF_2016, PHOTON_ELECTRON_VETO_SF_2017, PHOTON_ELECTRON_VETO_SF_2018
photon_electron_veto_sf_bins = {
    "2016" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016preVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2016postVFP" : PHOTON_ELECTRON_VETO_SF_2016,
    "2017" : PHOTON_ELECTRON_VETO_SF_2017,
    "2018" : PHOTON_ELECTRON_VETO_SF_2018,
    "2022preEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2022postEE" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023preBPix" : PHOTON_ELECTRON_VETO_SF_2018, #FIXME
    "2023postBPix" : PHOTON_ELECTRON_VETO_SF_2018 #FIXME
}

def photon_electron_veto_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations = systematic_from_bins(
        bins = photon_electron_veto_sf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        central_only = central_only
    )
    print(variations)
    return variations

##################
### Trigger SF ###
##################
# Note: since the trigger sf applies separate SF for the lead/sublead photons,
# it's easiest to just cast this as an EventWeightSystematic (rather than ObjectWeightSystematic as we would typically do)
# and just multiply the lead/sublead variations manually by hand

from higgs_dna.systematics.data.trigger_sf import LEAD_TRIGGER_SF_2016, SUBLEAD_TRIGGER_SF_2016, LEAD_TRIGGER_SF_2017, SUBLEAD_TRIGGER_SF_2017, LEAD_TRIGGER_SF_2018, SUBLEAD_TRIGGER_SF_2018 
lead_trigger_sf_bins = {
    "2016" : LEAD_TRIGGER_SF_2016,
    "2016preVFP" : LEAD_TRIGGER_SF_2016,
    "2016postVFP" : LEAD_TRIGGER_SF_2016,
    "2017" : LEAD_TRIGGER_SF_2017,
    "2018" : LEAD_TRIGGER_SF_2018,
    "2022preEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : LEAD_TRIGGER_SF_2018 #FIXME
}
sublead_trigger_sf_bins = {
    "2016" : SUBLEAD_TRIGGER_SF_2016,
    "2016preVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2016postVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2017" : SUBLEAD_TRIGGER_SF_2017,
    "2018" : SUBLEAD_TRIGGER_SF_2018, 
    "2022preEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : SUBLEAD_TRIGGER_SF_2018 #FIXME
}

def trigger_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9"), ("Photon", "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations_lead = systematic_from_bins(
        bins = lead_trigger_sf_bins[year], 
        variables = {
            "photon_r9" : events.LeadPhoton.r9,
            "photon_eta" : abs(events.LeadPhoton.eta),
            "photon_pt" : events.LeadPhoton.pt
        },
        central_only = central_only
    )

    variations_sublead = systematic_from_bins(
        bins = sublead_trigger_sf_bins[year],
        variables = {
            "photon_r9" : events.SubleadPhoton.r9,
            "photon_eta" : abs(events.SubleadPhoton.eta),
            "photon_pt" : events.SubleadPhoton.pt
        },
        central_only = central_only
    )
 
    # Multiply up/down/central variations together, following this treatment in flashgg:
    # https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Systematics/interface/DiPhotonFromSinglePhotonViewBase.h#L85-L87
    variations = {}
    for key in variations_lead.keys():
        variations[key] = variations_lead[key] * variations_sublead[key]

    return variations

############
### FNUF ###
############

from higgs_dna.systematics.data.fnuf import FNUF_2016, FNUF_2017, FNUF_2018
fnuf_bins = {
    "2016" : FNUF_2016,
    "2016preVFP" : FNUF_2016,
    "2016postVFP" : FNUF_2016,
    "2017" : FNUF_2017,
    "2018" : FNUF_2018,
    "2022preEE" : FNUF_2018,
    "2022postEE" : FNUF_2018,
    "2023preBPix" : FNUF_2018,
    "2023postBPix" : FNUF_2018,
}

def fnuf_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    """

    """
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message) 

    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "eb":
        mask = photons.isScEtaEB == True
    elif loc == "ee":
        mask = photons.isScEtaEE == True

    variations = ic_systematic_from_bins(
        bins = fnuf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations

################
### Material ###
################

from higgs_dna.systematics.data.material import MATERIAL_2016, MATERIAL_2017, MATERIAL_2018
material_bins = {
    "2016" : MATERIAL_2016,
    "2016preVFP" : MATERIAL_2016,
    "2016postVFP" : MATERIAL_2016,
    "2017" : MATERIAL_2017,
    "2018" : MATERIAL_2018,
    "2022preEE" : MATERIAL_2018,
    "2022postEE" : MATERIAL_2018,
    "2023preBPix" : MATERIAL_2018,
    "2023postBPix" : MATERIAL_2018
}

def material_unc(events, year, nominal_only, modify_nominal, loc = "all"):
    photons = events.Photon

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "central_barrel":
        mask = abs(photons.eta) <= 1.0
    elif loc == "outer_barrel":
        mask = (abs(photons.eta) > 1.0) & (abs(photons.eta) <= 1.5)
    elif loc == "forward":
        mask = abs(photons.eta) > 1.5
    
    variations = ic_systematic_from_bins(
        bins = material_bins[year],
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations  


def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations


#######################
### Photon MC Smear ###
#######################

def photon_mc_smear(events, r9, loc):
    photons = events.Photon

    # Split into high/low r9 and EE/EB
    if r9 == "high":
        r9_cut = photons.r9 > 0.94
    elif r9 == "low":
        r9_cut = photons.r9 <= 0.94

    if loc == "eb":
        loc_cut = photons.isScEtaEB == True
    elif loc == "ee":
        loc_cut = photons.isScEtaEE == True

    mask = r9_cut & loc_cut

    variations = {}
    variations["up"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaUp,
            photons.pt
    )
    variations["down"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaDown,
            photons.pt
    )

    return variations

photon_scale_FILE = {"2022preEE" : "jsonpog-integration/POG/EGM/2022_Summer22/photonSS_EtDependent.json",
    "2022postEE" : "jsonpog-integration/POG/EGM/2022_Summer22EE/photonSS_EtDependent.json",
    "2023preBPix" : "jsonpog-integration/POG/EGM/2023_Summer23/photonSS_EtDependent.json",
    "2023postBPix" : "jsonpog-integration/POG/EGM/2023_Summer23BPix/photonSS_EtDependent.json"
}

year_names = {
    "2022preEE" : "EGMSmearAndSyst_PhoPTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_PhoPTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_PhoPTsplit_2023preBPIX",
    "2023postBPix" : "EGMSmearAndSyst_PhoPTsplit_2023postBPIX"
}

def photon_scale_smear_run3(events, year):
    """
    This function applies photon scale and smear corrections on the photon pt.
    It returns events with the corrected photon pt values.
    """
    logger.info("[Photon Systematics] Applying photon scale and smear corrections for year %s", year)

    required_fields = [
        ("Photon", "x_calo"), ("Photon", "y_calo"), ("Photon", "z_calo"), ("Photon", "eta"), ("Photon", "pt"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))

    photons = events["Photon"]
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    # Calculate eta from calo coordinates
    # rho = numpy.sqrt(awkward.to_numpy(photons_flattened.x_calo)**2 + awkward.to_numpy(photons_flattened.y_calo)**2)
    # theta = numpy.arctan(rho / awkward.to_numpy(photons_flattened.z_calo))
    photons_AbsEta = numpy.abs(awkward.to_numpy(photons_flattened.eta))

    photons_pt = awkward.to_numpy(photons_flattened.pt)
    photons_r9 = awkward.to_numpy(photons_flattened.r9)

    # Apply smear corrections first
    smear = evaluator[year_names[year]].evalv("smear", photons_pt, photons_r9, photons_AbsEta)
    rng = numpy.random.default_rng(seed=123)
    smear_val = awkward.where(
        (photons_AbsEta > 3.0) | (photons_pt < 20.0),
        awkward.ones_like(photons_pt, dtype=float),
        rng.normal(loc=1., scale=numpy.abs(smear))
    )
    
    # Apply central smear correction to photon pt
    corrected_pt = awkward.to_numpy(photons_pt * smear_val)

    # Calculate smear systematics
    for syst in ["smear_up", "smear_down"]:
        smear_syst = evaluator[year_names[year]].evalv(syst, photons_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            smear_up = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear+smear_syst)), n_photons)
            events["Photon", "dEsigmaUp"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_up, dtype=float),
                smear_up
            )
        elif "down" in syst:
            smear_down = awkward.unflatten(rng.normal(loc=1., scale=numpy.abs(smear-smear_syst)), n_photons)
            events["Photon", "dEsigmaDown"] = photons.pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(smear_down, dtype=float),
                smear_down
            )
    
    print("Photon pt before scale and smear corrections:", photons.pt)
    events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
    photons = events["Photon"]

    # Calculate scale systematics
    for syst in ["scale_up", "scale_down"]:
        scale = evaluator[year_names[year]].evalv(syst, corrected_pt, photons_r9, photons_AbsEta)
        
        if "up" in syst:
            scale_up = awkward.unflatten(scale, n_photons)
            events["Photon", "pt_ScaleUp"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_up, dtype=float),
                scale_up
            )
        elif "down" in syst:
            scale_down = awkward.unflatten(scale, n_photons)
            events["Photon", "pt_ScaleDown"] = photons.corrected_pt * awkward.where(
                (abs(photons.eta) > 3.0) | (photons.pt < 20.0),
                awkward.ones_like(scale_down, dtype=float),
                scale_down
            )

    print("Photon pt after scale and smear corrections:", events["Photon", "corrected_pt"])
    print("Photon pt Scale Up:", events["Photon", "pt_ScaleUp"])
    print("Photon pt Scale Down:", events["Photon", "pt_ScaleDown"])
    print("Photon pt Smear Up:", events["Photon", "dEsigmaUp"])
    print("Photon pt Smear Down:", events["Photon", "dEsigmaDown"])

    return events
