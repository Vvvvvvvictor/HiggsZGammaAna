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
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_eltrig23_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig32_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_eltrig32_2018_efficiencies.json"
}
DoubleEle_LowLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_eltrig12_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_eltrig12_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig12_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_eltrig12_2018_efficiencies.json"
}
DoubleEle_HighLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_eltrig23_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_eltrig23_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_eltrig23_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_eltrig23_2018_efficiencies.json"
}

SingleMuon_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_mutrig24_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_mutrig24_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig27_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_mutrig24_2018_efficiencies.json"
}
DoubleMuon_LowLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_mutrig8_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_mutrig8_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig8_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_mutrig8_2018_efficiencies.json"
}
DoubleMuon_HighLeg_HLT_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/hzg_mutrig17_2016APV_efficiencies.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/hzg_mutrig17_2016_efficiencies.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/hzg_mutrig17_2017_efficiencies.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/hzg_mutrig17_2018_efficiencies.json"
}


def singleEle_correction(events, year, central_only, type="mc"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        ("Electron", "pt"), ("Electron", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleEle_HLT_FILE[year]))

    electrons = events["Electron"]
    electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
    electrons = awkward.firsts(electrons)

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
    sf = evaluator[f"eff{type}"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator[f"syst{type}"].evalv(
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

def doubleEleLowLeg_correction(events, year, central_only, type="mc"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        ("Electron", "pt"), ("Electron", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_LowLeg_HLT_FILE[year]))

    electrons = events["Electron"]
    electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
    electrons = electrons[:,1:2]

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
    sf = evaluator[f"eff{type}"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator[f"syst{type}"].evalv(
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


def doubleEleHighLeg_correction(events, year, central_only, type="mc"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        ("Electron", "pt"), ("Electron", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleEle_HighLeg_HLT_FILE[year]))

    electrons = events["Electron"]
    electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
    electrons = awkward.firsts(electrons)

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
    sf = evaluator[f"eff{type}"].evalv(
            ele_pt,
            ele_eta
    )
    variations["central"] = awkward.unflatten(sf, n_electrons)

    if not central_only:
        syst = evaluator[f"syst{type}"].evalv(
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

def singleMuon_correction(events, year, central_only, type="mc"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        ("Muon", "pt"), ("Muon", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(SingleMuon_HLT_FILE[year]))

    muons = events["Muon"]
    muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]
    muons = awkward.firsts(muons)

    # Flatten muons then convert to numpy for compatibility with correctionlib
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    mu_abseta = numpy.clip(
        awkward.to_numpy(abs(muons_flattened.eta)),
        0,
        2.39999 # SFs only valid up to eta 2.4
    )

    mu_pt = numpy.clip(
        awkward.to_numpy(muons_flattened.pt),
        5.0, # SFs only valid for pT >= 3.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    sf = evaluator[f"eff{type}"].evalv(
            mu_pt,
            mu_eta
    )
    variations["central"] = awkward.unflatten(sf, n_muons)

    if not central_only:
        syst = evaluator[f"syst{type}"].evalv(
            mu_pt,
            mu_eta
            )
        variations["up"] = awkward.unflatten(syst+sf, n_muons)
        variations["down"] = awkward.unflatten(sf-syst, n_muons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (muons.pt < 5.0) | (muons.pt >= 500.0) | (abs(muons.eta) >= 2.4),
                awkward.ones_like(variations[var], dtype=float, highlevel=False),
                variations[var]
        )

    return variations

def doubleMuonLowLeg_correction(events, year, central_only, type="mc"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        ("Muon", "pt"), ("Muon", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_LowLeg_HLT_FILE[year]))

    muons = events["Muon"]
    muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]
    muons = muons[:,1:2]

    # Flatten muons then convert to numpy for compatibility with correctionlib
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    mu_abseta = numpy.clip(
        awkward.to_numpy(abs(muons_flattened.eta)),
        0,
        2.39999 # SFs only valid up to eta 2.4
    )

    mu_pt = numpy.clip(
        awkward.to_numpy(muons_flattened.pt),
        5.0, # SFs only valid for pT >= 3.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    sf = evaluator[f"eff{type}"].evalv(
            mu_pt,
            mu_eta
    )
    variations["central"] = awkward.unflatten(sf, n_muons)

    if not central_only:
        syst = evaluator[f"syst{type}"].evalv(
            mu_pt,
            mu_eta
            )
        variations["up"] = awkward.unflatten(syst+sf, n_muons)
        variations["down"] = awkward.unflatten(sf-syst, n_muons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (muons.pt < 5.0) | (muons.pt >= 500.0) | (abs(muons.eta) >= 2.4),
                awkward.ones_like(variations[var], dtype=float, highlevel=False),
                variations[var]
        )

    return variations
    
def doubleMuonHighLeg_correction(events, year, central_only, type="mc"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        ("Muon", "pt"), ("Muon", "eta")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(DoubleMuon_HighLeg_HLT_FILE[year]))

    muons = events["Muon"]
    muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]
    muons = awkward.firsts(muons)

    # Flatten muons then convert to numpy for compatibility with correctionlib
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    mu_abseta = numpy.clip(
        awkward.to_numpy(abs(muons_flattened.eta)),
        0,
        2.39999 # SFs only valid up to eta 2.4
    )

    mu_pt = numpy.clip(
        awkward.to_numpy(muons_flattened.pt),
        5.0, # SFs only valid for pT >= 3.0
        499.999 # and pT < 500.
    )

    # Calculate SF and syst
    variations = {}
    sf = evaluator[f"eff{type}"].evalv(
            mu_pt,
            mu_eta
    )
    variations["central"] = awkward.unflatten(sf, n_muons)

    if not central_only:
        syst = evaluator[f"syst{type}"].evalv(
            mu_pt,
            mu_eta
            )
        variations["up"] = awkward.unflatten(syst+sf, n_muons)
        variations["down"] = awkward.unflatten(sf-syst, n_muons)

    for var in variations.keys():
        # Set SFs = 1 for leptons which are not applicable
        variations[var] = awkward.where(
                (muons.pt < 5.0) | (muons.pt >= 500.0) | (abs(muons.eta) >= 2.4),
                awkward.ones_like(variations[var],
                dtype=float, highlevel=False),
                variations[var]
        )

    return variations