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
    "2016" : "higgs_dna/systematics/data/2016postVFP_UL/photon_wp80mceff_2016.json",
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/photon_wp80mceff_2016.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/photon_wp80mceff_2016APV.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/photon_wp80mceff_2017.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/photon_wp80mceff_2018.json",
    "2022preEE" : "higgs_dna/systematics/data/2018_UL/photon_wp80mceff_2018.json",
    "2022postEE" : "higgs_dna/systematics/data/2018_UL/photon_wp80mceff_2018.json",
    "2023preBPix" : "higgs_dna/systematics/data/2018_UL/photon_wp80mceff_2018.json",
    "2023postBPix" : "higgs_dna/systematics/data/2018_UL/photon_wp80mceff_2018.json",
}

PHOTON_ID_SF = {
    "2016" : "2016postVFP",
    "2016preVFP" : "2016preVFP",
    "2016postVFP" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018",
    "2022preEE" : "2018",
    "2022postEE" : "2018",
    "2023preBPix" : "2018",
    "2023postBPix" : "2018",
}

def photon_id_sf(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_ID_SF_FILE[year]))

    photons = events[input_collection]

    # Flatten photons then convert to numpy for compatibility with correctionlib
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    pho_eta = numpy.clip(
        awkward.to_numpy(photons_flattened.eta),
        -4.99999,
        4.99999 # SFs only valid up to eta 2.5
    )

    pho_pt = numpy.clip(
        awkward.to_numpy(photons_flattened.pt),
        20.0, # SFs only valid for pT >= 20.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    sf = evaluator["UL-Photon-ID-SF"].evalv(
            PHOTON_ID_SF[year],
            "sf",
            working_point,
            pho_eta,
            pho_pt
    )
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        syst_vars = ["sfup", "sfdown"] 
        for syst_var in syst_vars:
            syst = evaluator["UL-Photon-ID-SF"].evalv(
                    PHOTON_ID_SF[year],
                    syst_var,
                    working_point,
                    pho_eta,
                    pho_pt
            )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_photons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (photons.pt < 20.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations[var], dtype=float),
                variations[var]
        )

    return variations

########################
#### Photon CSEV SF ####
########################

def photon_CSEV_sf(events, year, central_only, input_collection, working_point = "none", Bin = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_ID_SF_FILE[year]))

    photons = events[input_collection]

    # # Flatten photons then convert to numpy for compatibility with correctionlib
    # n_photons = awkward.num(photons)
    # photons_flattened = awkward.flatten(photons)

    # Calculate SF and syst
    variations = {}
    sf = evaluator["UL-Photon-CSEV-SF"].evalv(
            PHOTON_ID_SF[year],
            "sf",
            working_point,
            Bin
    )
    # variations["central"] = awkward.unflatten(sf, n_photons)
    variations["central"] = awkward.unflatten(sf, 1)

    if not central_only:
        syst_vars = ["sfup", "sfdown"] 
        for syst_var in syst_vars:
            syst = evaluator["UL-Photon-CSEV-SF"].evalv(
                    PHOTON_ID_SF[year],
                    syst_var,
                    working_point,
                    Bin
            )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            # variations[syst_var_name] = awkward.unflatten(syst, n_photons)
            variations[syst_var_name] = awkward.unflatten(syst, 1)

    # for var in variations.keys():
    #     # Set SFs = 1 for leptons which are not applicable
    #     variations[var] = awkward.where(
    #             (photons.pt < 20.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
    #             awkward.ones_like(variations[var]),
    #             variations[var]
    #     )

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
def photon_scale_run3(events, year, central_only, input_collection, working_point = "none"):

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "r9"),(input_collection, "run"),(input_collection,"seedGain")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))

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
