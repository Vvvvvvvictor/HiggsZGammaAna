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

muon_ID_SF_FILE = {
    "2016" : "higgs_dna/systematics/data/2017_UL/muon_mceff.json",
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/muon_mceff.json",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/muon_mceff.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/muon_mceff.json",
    "2018" : "higgs_dna/systematics/data/2017_UL/muon_mceff.json"
}

muon_ID_SF = {
    "2016" : "2016postVFP",
    "2016preVFP" : "2016preVFP",
    "2016postVFP" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018"
}

def muon_LooseID_sf(events, year, central_only, input_collection, working_point = "none"):

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(muon_ID_SF_FILE[year]))

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
    sf = evaluator["Muon_LooseID_MCeff"].evalv(
            "effmc",
            abs(muon_eta),
            muon_pt
    )
    variations["central"] = awkward.unflatten(sf, n_muons)

    if not central_only:
        syst_var = "systmc"
        syst = evaluator["Muon_LooseID_MCeff"].evalv(
            syst_var,
            abs(muon_eta),
            muon_pt
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


ELECTRON_ID_SF_FILE = {
    "2016" : "higgs_dna/systematics/data/2017_UL/electron_WPL.json",
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/electron_WPL.json",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/electron_WPL.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/electron_WPL.json",
    "2018" : "higgs_dna/systematics/data/2017_UL/electron_WPL.json"
}

ELECTRON_ID_SF = {
    "2016" : "2016postVFP",
    "2016preVFP" : "2016preVFP",
    "2016postVFP" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018"
}

def electron_WPL_sf(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(ELECTRON_ID_SF_FILE[year]))

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

    # Calculate SF and syst
    variations = {}
    sf = evaluator["ElectronWPL"].evalv(
            "effmc",
            ele_eta,
            ele_pt
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst_var = "systmc"
        syst = evaluator["ElectronWPL"].evalv(
            syst_var,
            ele_eta,
            ele_pt
            )
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



MUON_ID_SF_FILE = {
    "2016" : "jsonpog-integration/POG/MUO/2016postVFP_UL/muon.json",
    "2016preVFP" : "jsonpog-integration/POG/MUO/2016preVFP_UL/muon.json",
    "2016postVFP" : "jsonpog-integration/POG/MUO/2016postVFP_UL/muon.json",
    "2017" : "jsonpog-integration/POG/MUO/2017_UL/muon.json",
    "2018" : "jsonpog-integration/POG/MUO/2018_UL/muon.json"
}

MUON_ID_SF = {
    "2016" : "2016postVFP_UL",
    "2016preVFP" : "2016preVFP_UL",
    "2016postVFP" : "2016postVFP_UL",
    "2017" : "2017_UL",
    "2018" : "2018_UL"
}

DUMMY_LEPTON_SFs = {
    "variables" : ["lepton_pt", "lepton_eta"],
    "bins" : [
        {
            "lepton_pt" : [25.0, 50.0],
            "lepton_eta" : [0.0, 1.5],
            "value" : 0.5,
            "uncertainty" : 0.25
        },
        {
            "lepton_pt" : [25.0, 50.0],
            "lepton_eta" : [1.5, 999],
            "value" : 1.5,
            "uncertainty" : 0.25
        },
        {
            "lepton_pt" : [50.0, 999999.],
            "lepton_eta" : [0.0, 999999.],
            "value" : 0.5,
            "uncertainty" : 0.25
        },
 
    ]
}

def dummy_lepton_sf(events, central_only, input_collection):
    """
    Dummy function illustrating a lepton scale factor.
    Should be deleted once real examples are implemented.
    """
    required_fields = [
            (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    leptons = events[input_collection]
    variations = systematic_from_bins(
            bins = DUMMY_LEPTON_SFs,
            variables = {
                "lepton_pt" : leptons.pt,
                "lepton_eta" : leptons.eta
            },
            central_only = central_only
    )

    return variations 

def muon_id_sfSTAT(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ID_SF

    if not central_only:
        variations["up"] = leptons.ID_SF+leptons.ID_SFstat
        variations["down"] = leptons.ID_SF-leptons.ID_SFstat

    return variations

def muon_id_sfSYS(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ID_SF

    if not central_only:
        variations["up"] = leptons.ID_SF+leptons.ID_SFsyst
        variations["down"] = leptons.ID_SF-leptons.ID_SFsyst

    return variations

def muon_iso_sfSTAT(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ISO_SF

    if not central_only:
        variations["up"] = leptons.ISO_SF+leptons.ISO_SFstat
        variations["down"] = leptons.ISO_SF-leptons.ISO_SFstat

    return variations

def muon_iso_sfSYS(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.ISO_SF

    if not central_only:
        variations["up"] = leptons.ISO_SF+leptons.ISO_SFsyst
        variations["down"] = leptons.ISO_SF-leptons.ISO_SFsyst

    return variations

def tauJ_sf(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.sfDeepTau2017v2p1VSjet_Loose

    if not central_only:
        variations["up"] = leptons.sfDeepTau2017v2p1VSjet_LooseUp
        variations["down"] = leptons.sfDeepTau2017v2p1VSjet_LooseDown

    return variations

def tauM_sf(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.sfDeepTau2017v2p1VSmu_VLoose

    if not central_only:
        variations["up"] = leptons.sfDeepTau2017v2p1VSmu_VLooseUp
        variations["down"] = leptons.sfDeepTau2017v2p1VSmu_VLooseDown

    return variations

def tauE_sf(events, year, central_only, input_collection, working_point = "none"):
    leptons = events[input_collection]
    variations = {}
    variations["central"] = leptons.sfDeepTau2017v2p1VSe_VVLoose

    if not central_only:
        variations["up"] = leptons.sfDeepTau2017v2p1VSe_VVLooseUp
        variations["down"] = leptons.sfDeepTau2017v2p1VSe_VVLooseDown #NB according to Tau POG we should be using Tight for the vsJ scale factor to be valid!! check twiki 

    return variations

SingleEle_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies",
    "2018" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies"
}
DoubleEle_LowLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig12_2017_efficiencies.json"
}
DoubleEle_HighLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies.json"
}

SingleMuon_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig27_2017_efficiencies.jsons",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig27_2017_efficiencies.jsons",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig27_2017_efficiencies.jsons",
    "2018" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig27_2017_efficiencies.jsons"
}
DoubleMuon_LowLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig8_2017_efficiencies.json"
}
DoubleMuon_HighLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig17_2017_efficiencies.json"
}


def singleEle_MCcorrection(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleEle_HLT_FILE[year]))

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

    # Calculate SF and syst
    variations = {}
    sf = evaluator["effmc"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator["systmc"].evalv(
            ele_pt,
            ele_eta
            )
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
def singleEle_Datacorrection(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleEle_HLT_FILE[year]))

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

    # Calculate SF and syst
    variations = {}
    sf = evaluator["effdata"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator["systdata"].evalv(
            ele_pt,
            ele_eta
            )
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

def doubleEleLowLeg_MCcorrection(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_LowLeg_HLT_FILE[year]))

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

    # Calculate SF and syst
    variations = {}
    sf = evaluator["effmc"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator["systmc"].evalv(
            ele_pt,
            ele_eta
            )
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
def doubleEleLowLeg_Datacorrection(events, year, central_only, input_collection, working_point = "none"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "pt"), (input_collection, "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_LowLeg_HLT_FILE[year]))

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

    # Calculate SF and syst
    variations = {}
    sf = evaluator["effdata"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator["systdata"].evalv(
            ele_pt,
            ele_eta
            )
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
