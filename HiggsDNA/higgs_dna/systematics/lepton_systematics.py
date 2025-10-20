import awkward
import numpy
import json

from correctionlib import _core
import correctionlib
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

    # 1. Handle reconstructed muons
    reco_muons = events[input_collection]
    n_reco_muons = awkward.num(reco_muons)
    reco_muons_flattened = awkward.flatten(reco_muons)

    reco_muon_eta = numpy.clip(
        awkward.to_numpy(reco_muons_flattened.eta),
        -2.39999,
        2.39999 # SFs only valid up to eta 2.5
    )
    reco_muon_pt = numpy.clip(
        awkward.to_numpy(reco_muons_flattened.pt),
        5.0, # SFs only valid for pT >= 10.0
        499.999 # and pT < 500.
    )

    reco_variations = {}
    sf = evaluator["sf_pass"].evalv(reco_muon_pt, reco_muon_eta)
    reco_variations["central"] = awkward.unflatten(sf, n_reco_muons)

    if not central_only:
        syst = evaluator["unc_pass"].evalv(reco_muon_pt, reco_muon_eta)
        reco_variations["up"] = awkward.unflatten(syst + sf, n_reco_muons)
        reco_variations["down"] = awkward.unflatten(sf - syst, n_reco_muons)

    for var in reco_variations.keys():
        reco_variations[var] = awkward.where(
            (reco_muons.pt < 5.0) | (reco_muons.pt >= 500.0) | (abs(reco_muons.eta) >= 2.4),
            awkward.ones_like(reco_variations[var], dtype=float),
            reco_variations[var]
        )

    return reco_variations


def muon_LooseID_sf_nomatch(events, year, central_only, input_collection):
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "pdgId")
    ]
    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning(f"Missing fields for muon_LooseID_sf_nomatch: {missing_fields}")
        # return None

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(MUON_ID_SF_FILE[year]))
    unmatched_gen_muons = events[input_collection]
    n_unmatched_gen = awkward.num(unmatched_gen_muons)
    unmatched_gen_flattened = awkward.flatten(unmatched_gen_muons)

    unmatched_variations = {}
    if len(unmatched_gen_flattened) > 0:
        gen_muon_eta = numpy.clip(awkward.to_numpy(unmatched_gen_flattened.eta), -2.39999, 2.39999)
        gen_muon_pt = numpy.clip(awkward.to_numpy(unmatched_gen_flattened.pt), 5.0, 499.999)
        sf = evaluator["sf_fail"].evalv(gen_muon_pt, gen_muon_eta)
        unmatched_variations["central"] = awkward.unflatten(sf, n_unmatched_gen)

        if not central_only:
            syst = evaluator["unc_fail"].evalv(gen_muon_pt, gen_muon_eta)
            unmatched_variations["up"] = awkward.unflatten(syst + sf, n_unmatched_gen)
            unmatched_variations["down"] = awkward.unflatten(sf - syst, n_unmatched_gen)

        for var in unmatched_variations.keys():
            unmatched_variations[var] = awkward.where(
                (unmatched_gen_muons.pt < 5.0) | (unmatched_gen_muons.pt >= 500.0) | (abs(unmatched_gen_muons.eta) >= 2.4),
                awkward.ones_like(unmatched_variations[var], dtype=float),
                unmatched_variations[var]
            )
    else:
        event_counts = awkward.num(events, axis=0)
        dummy_event_structure = awkward.ArrayBuilder()
        for i in range(event_counts):
            dummy_event_structure.begin_list()
            dummy_event_structure.end_list()
        
        dummy_event_structure = dummy_event_structure.snapshot()

        unmatched_variations = { "central": dummy_event_structure }
        if not central_only:
            unmatched_variations["up"] = dummy_event_structure
            unmatched_variations["down"] = dummy_event_structure

    return unmatched_variations


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

# run2: low: pT 8-40 med: pT 40-60 high: P > 50
# run3: high P > 50
MUON_RECO_SF_FILE = {
    "2016preVFP" : ["higgs_dna/systematics/data/2016preVFP_UL/muon_lowpT_recoSF2016preVFP.json", "higgs_dna/systematics/data/2016preVFP_UL/muon_medpT_recoSF2016preVFP.json", "higgs_dna/systematics/data/2016preVFP_UL/muon_highpT_recoSF2016preVFP.json"],
    "2016postVFP" : ["higgs_dna/systematics/data/2016postVFP_UL/muon_lowpT_recoSF2016postVFP.json","higgs_dna/systematics/data/2016postVFP_UL/muon_medpT_recoSF2016postVFP.json","higgs_dna/systematics/data/2016postVFP_UL/muon_highpT_recoSF2016postVFP.json"],
    "2017" : ["higgs_dna/systematics/data/2017_UL/muon_lowpT_recoSF2017.json","higgs_dna/systematics/data/2017_UL/muon_medpT_recoSF2017.json","higgs_dna/systematics/data/2017_UL/muon_highpT_recoSF2017.json"],
    "2018" : ["higgs_dna/systematics/data/2018_UL/muon_lowpT_recoSF2018.json","higgs_dna/systematics/data/2018_UL/muon_medpT_recoSF2018.json","higgs_dna/systematics/data/2018_UL/muon_highpT_recoSF2018.json"],
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/muon_highpT_recoSF2022preEE.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/muon_highpT_recoSF2022postEE.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/muon_highpT_recoSF2023preBPix.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/muon_highpT_recoSF2023postBPix.json"
}

def get_run2_reco_sf_from_json(sf_dict, eta, pt, syst):
    # Find eta bin
    eta_bin = None
    for key in sf_dict["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"]:
        if "abseta" not in key:
            continue
        eta_range = key.split(":")[1][1:-1].split(",")
        low_eta, high_eta = float(eta_range[0]), float(eta_range[1])
        if low_eta <= abs(eta) < high_eta:
            eta_bin = sf_dict["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"][key]
            break
    if eta_bin is None:
        return 1.0

    # Find pt bin
    pt_bin = None
    for key in eta_bin:
        if "pt" not in key:
            continue
        pt_range = key.split(":")[1][1:-1].split(",")
        low_pt, high_pt = float(pt_range[0]), float(pt_range[1])
        if low_pt <= pt < high_pt:
            pt_bin = eta_bin[key]
            break
    if pt_bin is None:
        return 1.0

    sf = pt_bin.get("value", 1.0)
    if syst == "central":
        return sf
    
    # For low pT, error is 'error'. For med pT, it's 'stat' and 'syst'.
    err = 0.0
    if "error" in pt_bin: # low pT
        err = pt_bin["error"]
    elif "stat" in pt_bin and "syst" in pt_bin: # med pT
        err = numpy.sqrt(pt_bin["stat"]**2 + pt_bin["syst"]**2)

    if syst == "up":
        return sf + err
    elif syst == "down":
        return sf - err
    
    return sf


def muon_reco_sf(events, year, central_only, input_collection):
    """
    See:
        - https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016
        - https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017
        - https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018
    """
    
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning("[muon_reco_sf] The events array is missing the following fields: %s." % str(missing_fields))
        if not central_only:
            return {"central" : awkward.ones_like(events.event), "up" : awkward.ones_like(events.event), "down" : awkward.ones_like(events.event)}
        else:
            return {"central" : awkward.ones_like(events.event)}

    muons = events[input_collection]
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    is_run2 = not (year.startswith("2022") or year.startswith("2023"))

    if is_run2:
        with open(misc_utils.expand_path(MUON_RECO_SF_FILE[year][0])) as f:
            sf_dict_low = json.load(f)
        with open(misc_utils.expand_path(MUON_RECO_SF_FILE[year][1])) as f:
            sf_dict_med = json.load(f)
        evaluator_high = _core.CorrectionSet.from_file(misc_utils.expand_path(MUON_RECO_SF_FILE[year][2]))
    else:
        evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(MUON_RECO_SF_FILE[year]))

    muon_eta = numpy.clip(
            awkward.to_numpy(muons_flattened.eta),
            -2.3999,
            2.3999
    )

    variations = {}
    var_names = ["central", "up", "down"] if not central_only else ["central"]

    for var in var_names:
        syst = "sf"
        if var == "up":
            syst = "sfup"
        elif var == "down":
            syst = "sfdown"

        if is_run2:
            muon_pt_flattened = awkward.to_numpy(muons_flattened.pt)
            muon_eta_flattened = awkward.to_numpy(muons_flattened.eta)
            
            syst_arg = "central"
            if var == "up":
                syst_arg = "up"
            elif var == "down":
                syst_arg = "down"

            sf_low = numpy.array([get_run2_reco_sf_from_json(sf_dict_low, eta, pt, syst_arg) for eta, pt in zip(muon_eta_flattened, muon_pt_flattened)])
            sf_med = numpy.array([get_run2_reco_sf_from_json(sf_dict_med, eta, pt, syst_arg) for eta, pt in zip(muon_eta_flattened, muon_pt_flattened)])

            muon_p_high = muon_pt_flattened * numpy.cosh(muon_eta_flattened)
            muon_p_high_clipped = numpy.clip(muon_p_high, 50.0, 4999.999)
            muon_eta_abs_high = numpy.clip(numpy.abs(muon_eta_flattened), 0.0, 2.3999)
            
            syst_arg_high = "nominal"
            if var == "up":
                syst_arg_high = "systup"
            elif var == "down":
                syst_arg_high = "systdown"

            sf_high = evaluator_high["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evalv(muon_eta_abs_high, muon_p_high_clipped, syst_arg_high)

            sf = numpy.where(
                    muon_pt_flattened < 40.0,
                    sf_low,
                    numpy.where(
                        muon_pt_flattened < 60.0,
                        sf_med,
                        sf_high
                    )
            )
        else:
            muon_pt_np = awkward.to_numpy(muons_flattened.pt)
            muon_eta_np = awkward.to_numpy(muons_flattened.eta)
            muon_p_np = muon_pt_np * numpy.cosh(muon_eta_np)

            muon_eta_abs = numpy.clip(numpy.abs(muon_eta_np), 0, 2.3999)
            muon_p_clipped = numpy.clip(muon_p_np, 50.0, 4999.999) # SFs are in p, not pt

            syst_arg = syst
            if syst == "sf":
                syst_arg = "nominal"
            elif syst == "sfup":
                syst_arg = "systup"
            elif syst == "sfdown":
                syst_arg = "systdown"
            
            valid_mask = (muon_p_clipped >= 50.0) & (muon_p_clipped < 5000.0) & (muon_eta_abs >= 0.0) & (muon_eta_abs < 2.4)

            # 默认值设为 1.0
            sf = numpy.ones_like(muon_eta_abs)

            # 对合法区域进行计算
            sf[valid_mask] = evaluator["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evalv(
                muon_eta_abs[valid_mask], muon_p_clipped[valid_mask], syst_arg
            )

        variations[var] = awkward.unflatten(sf, n_muons)

    is_run2_pt_low = awkward.ones_like(muons.pt, dtype=bool)
    if is_run2:
        is_run2_pt_low = (muons.pt >= 5.0)
    else:
        is_run2_pt_low = (muons.pt >= 50.0)

    for var in variations.keys():
        variations[var] = awkward.where(
            (is_run2_pt_low) & (muons.pt < 500.0) & (abs(muons.eta) < 2.4),
            variations[var],
            awkward.ones_like(variations[var], dtype=float)
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

    # 1. Handle reconstructed electrons
    reco_electrons = events[input_collection]
    n_reco_electrons = awkward.num(reco_electrons)
    reco_electrons_flattened = awkward.flatten(reco_electrons)

    reco_ele_eta = numpy.clip(
        awkward.to_numpy(reco_electrons_flattened.eta),
        -2.49999,
        2.49999    )

    reco_ele_pt = numpy.clip(
        awkward.to_numpy(reco_electrons_flattened.pt),
        7.0,        499.999    )
    
    reco_ele_phi = None
    if year == "2023postBPix":
        reco_ele_phi = awkward.to_numpy(reco_electrons_flattened.phi)

    reco_variations = {}
    
    if year == "2023postBPix":
        hole_mask = (reco_ele_eta > -1.566) & (reco_ele_eta < 0) & (reco_ele_phi > -1.2) & (reco_ele_phi < -0.8)
        sf_normal = evaluator["sf_pass"].evalv(reco_ele_pt, reco_ele_eta)
        sf_hole = evaluator_hole["sf_pass"].evalv(reco_ele_pt, reco_ele_eta)
        sf = numpy.where(hole_mask, sf_hole, sf_normal)
        reco_variations["central"] = awkward.unflatten(sf, n_reco_electrons)

        if not central_only:
            syst_normal = evaluator["unc_pass"].evalv(reco_ele_pt, reco_ele_eta)
            syst_hole = evaluator_hole["unc_pass"].evalv(reco_ele_pt, reco_ele_eta)
            syst = numpy.where(hole_mask, syst_hole, syst_normal)
            reco_variations["up"] = awkward.unflatten(syst + sf, n_reco_electrons)
            reco_variations["down"] = awkward.unflatten(sf - syst, n_reco_electrons)
    else:
        sf = evaluator["sf_pass"].evalv(reco_ele_pt, reco_ele_eta)
        reco_variations["central"] = awkward.unflatten(sf, n_reco_electrons)

        if not central_only:
            syst = evaluator["unc_pass"].evalv(reco_ele_pt, reco_ele_eta)
            reco_variations["up"] = awkward.unflatten(syst + sf, n_reco_electrons)
            reco_variations["down"] = awkward.unflatten(sf - syst, n_reco_electrons)

    for var in reco_variations.keys():
        reco_variations[var] = awkward.where(
            (reco_electrons.pt < 7.0) | (reco_electrons.pt >= 500.0) | (abs(reco_electrons.eta) >= 2.5),
            awkward.ones_like(reco_variations[var], dtype=float),
            reco_variations[var]
        )

    return reco_variations


def electron_WPL_sf_nomatch(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "pdgId")
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

    unmatched_gen_electrons = events[input_collection]
    n_unmatched_gen = awkward.num(unmatched_gen_electrons)
    unmatched_gen_flattened = awkward.flatten(unmatched_gen_electrons)

    unmatched_variations = {}
    if len(unmatched_gen_flattened) > 0:
        gen_ele_eta = numpy.clip(awkward.to_numpy(unmatched_gen_flattened.eta), -2.49999, 2.49999)
        gen_ele_pt = numpy.clip(awkward.to_numpy(unmatched_gen_flattened.pt), 7.0, 499.999)
        
        gen_ele_phi = None
        if year == "2023postBPix":
            gen_ele_phi = awkward.to_numpy(unmatched_gen_flattened.phi)

        if year == "2023postBPix":
            hole_mask = (gen_ele_eta > -1.566) & (gen_ele_eta < 0) & (gen_ele_phi > -1.2) & (gen_ele_phi < -0.8)
            sf_normal = evaluator["sf_fail"].evalv(gen_ele_pt, gen_ele_eta)
            sf_hole = evaluator_hole["sf_fail"].evalv(gen_ele_pt, gen_ele_eta)
            sf = numpy.where(hole_mask, sf_hole, sf_normal)
            unmatched_variations["central"] = awkward.unflatten(sf, n_unmatched_gen)

            if not central_only:
                syst_normal = evaluator["unc_fail"].evalv(gen_ele_pt, gen_ele_eta)
                syst_hole = evaluator_hole["unc_fail"].evalv(gen_ele_pt, gen_ele_eta)
                syst = numpy.where(hole_mask, syst_hole, syst_normal)
                unmatched_variations["up"] = awkward.unflatten(syst + sf, n_unmatched_gen)
                unmatched_variations["down"] = awkward.unflatten(sf - syst, n_unmatched_gen)
        else:
            sf = evaluator["sf_fail"].evalv(gen_ele_pt, gen_ele_eta)
            unmatched_variations["central"] = awkward.unflatten(sf, n_unmatched_gen)

            if not central_only:
                syst = evaluator["unc_fail"].evalv(gen_ele_pt, gen_ele_eta)
                unmatched_variations["up"] = awkward.unflatten(syst + sf, n_unmatched_gen)
                unmatched_variations["down"] = awkward.unflatten(sf - syst, n_unmatched_gen)

        for var in unmatched_variations.keys():
            unmatched_variations[var] = awkward.where(
                (unmatched_gen_electrons.pt < 7.0) | (unmatched_gen_electrons.pt >= 500.0) | (abs(unmatched_gen_electrons.eta) >= 2.5),
                awkward.ones_like(unmatched_variations[var], dtype=float),
                unmatched_variations[var]
            )
    else:
        event_counts = awkward.num(events, axis=0)
        dummy_event_structure = awkward.ArrayBuilder()
        for i in range(event_counts):
            dummy_event_structure.begin_list()
            dummy_event_structure.end_list()
        
        dummy_event_structure = dummy_event_structure.snapshot()

        unmatched_variations = { "central": dummy_event_structure }
        if not central_only:
            unmatched_variations["up"] = dummy_event_structure
            unmatched_variations["down"] = dummy_event_structure

    return unmatched_variations


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

def electron_ISO_sf(events, year, central_only, input_collection):
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
            sf_data = evaluator["effdata"].evalv(
                    ele_pt,
                    ele_eta
            )
            sf_mc = evaluator["effmc"].evalv(
                    ele_pt,
                    ele_eta
            )
            sf = sf_data / sf_mc
            sub_variations["central"] = awkward.unflatten(sf, n_electrons)

            if not central_only:
                syst_data = evaluator["systdata"].evalv(
                    ele_pt,
                    ele_eta
                )
                syst_mc = evaluator["systmc"].evalv(
                    ele_pt,
                    ele_eta
                )
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
        evaluator = evaluators[1]
        sub_variations = {}
        sf_data = evaluator["effdata"].evalv(
                ele_pt,
                ele_eta
        )
        sf_mc = evaluator["effmc"].evalv(
                ele_pt,
                ele_eta
        )
        sf = (1-sf_data)/(1-sf_mc)
        sub_variations["central"] = awkward.unflatten(sf, n_electrons)
        if not central_only:
            syst_data = evaluator["systdata"].evalv(
                ele_pt,
                ele_eta
            )
            syst_mc = evaluator["systmc"].evalv(
                ele_pt,
                ele_eta
            )
            syst = numpy.sqrt((syst_data/(1-sf_mc))**2 + (syst_mc/(1-sf_mc)*sf)**2)
            sub_variations["up"] = awkward.unflatten(syst+sf, n_electrons)
            sub_variations["down"] = awkward.unflatten(sf-syst, n_electrons)
            
    variations_list.append(sub_variations)

    variations = {}
    for var in variations_list[0].keys():
        variations[var] = awkward.where(
            (electrons.pt < 7.0) | (electrons.pt >= 500.0) | (abs(electrons.eta) >= 2.5) | (electrons.miniPFRelIso_all >= boundaries[-1]),
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

electron_smear_names = {
    "2022preEE" : "EGMSmearAndSyst_ElePTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_ElePTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_ElePTsplit_2023preBPIX", 
    "2023postBPix" : "EGMSmearAndSyst_ElePTsplit_2023postBPIX"
}

electron_scale_names = {
    "2022preEE": "EGMScale_Compound_Ele_2022preEE",
    "2022postEE": "EGMScale_Compound_Ele_2022postEE",
    "2023preBPix": "EGMScale_Compound_Ele_2023preBPIX",
    "2023postBPix": "EGMScale_Compound_Ele_2023postBPIX"
}

def electron_scale_smear_run3(events, year, is_data):
    """
    See:
        - https://gitlab.cern.ch/cms-egamma/EgammaPostRecoTools/-/tree/master/EgammaAnalysis/ElectronTools/data/Run3
    """
    logger.info("[Lepton Systematics] Applying electron scale and smear corrections for year %s", year)

    required_fields = [
        ("Electron", "eta"), ("Electron", "pt"), ("Electron", "r9"), ("Electron", "deltaEtaSC")
    ]
    if is_data:
        required_fields.append(("Electron", "seedGain"))
        required_fields.append("run")

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        logger.warning(
            "[Lepton Systematics] Electron scale and smear corrections will not be applied, because the following fields are missing: %s" % str(missing_fields)
        )
        return events

    if is_data:
        evaluator = correctionlib.CorrectionSet.from_file(misc_utils.expand_path(electron_scale_FILE[year]))
    else:
        evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(electron_scale_FILE[year]))

    electrons = events["Electron"]
    n_electrons = awkward.num(electrons)
    electrons_flattened = awkward.flatten(electrons)

    electrons_scEta = awkward.to_numpy(electrons_flattened.eta + electrons_flattened.deltaEtaSC)
    electrons_AbsScEta = numpy.abs(electrons_scEta)
    electrons_pt = awkward.to_numpy(electrons_flattened.pt)
    electrons_r9 = awkward.to_numpy(electrons_flattened.r9)
    if is_data:
        electrons_seedGain = awkward.to_numpy(electrons_flattened.seedGain)
        run_arr_flattened = numpy.repeat(awkward.to_numpy(events["run"]), n_electrons)

    if is_data:
        scale = evaluator.compound[electron_scale_names[year]].evaluate("scale", run_arr_flattened, electrons_scEta, electrons_r9, electrons_AbsScEta, electrons_pt, electrons_seedGain)
        scale = awkward.where(
            (electrons_AbsScEta > 3.0) | (electrons_pt < 20.0),
            awkward.ones_like(electrons_pt, dtype=float),
            scale
        )
        corrected_pt = awkward.to_numpy(electrons_pt * scale)
        events["Electron", "corrected_pt"] = awkward.unflatten(corrected_pt, n_electrons)
        logger.info("[Lepton Systematics] Electron pt before scale correction (data): %s", electrons.pt)
        logger.info("[Lepton Systematics] Electron pt after scale correction (data): %s", corrected_pt)
        return events

    # Apply smear corrections first
    smear = evaluator[electron_smear_names[year]].evalv("smear", electrons_pt, electrons_r9, electrons_AbsScEta)
    rng = numpy.random.default_rng(seed=8011)
    smear_val = awkward.where(
        (electrons_AbsScEta > 3.0) | (electrons_pt < 20.0),
        awkward.ones_like(electrons_pt, dtype=float),
        rng.normal(loc=1., scale=numpy.abs(smear))
    )
    
    # Apply central smear correction to electron pt
    corrected_pt = awkward.to_numpy(electrons_pt * smear_val)

    # Calculate smear systematics
    for syst in ["smear_up", "smear_down"]:
        smear_syst = evaluator[electron_smear_names[year]].evalv(syst, electrons_pt, electrons_r9, electrons_AbsScEta)
        smear_val_syst = awkward.where(
            (electrons_AbsScEta > 3.0) | (electrons_pt < 20.0),
            awkward.ones_like(electrons_pt, dtype=float),
            rng.normal(loc=1., scale=numpy.abs(smear_syst))
        )
        events["Electron", "dEsigma" + syst.replace("smear_", "").capitalize()] = awkward.unflatten(corrected_pt * (smear_val_syst - 1), n_electrons)

    print("Electron pt before smear correction:", electrons.pt)
    events["Electron", "corrected_pt"] = awkward.unflatten(corrected_pt, n_electrons)
    electrons = events["Electron"]

    # Calculate scale systematics
    for syst in ["scale_up", "scale_down"]:
        scale_syst = evaluator[electron_smear_names[year]].evalv(syst, electrons_pt, electrons_r9, electrons_AbsScEta)
        events["Electron", "dEscale" + syst.replace("scale_", "").capitalize()] = awkward.unflatten(corrected_pt * scale_syst, n_electrons)

    logger.info("[Lepton Systematics] Electron scale and smear corrections applied successfully")
    print("Electron pt after scale and smear corrections:", events["Electron", "corrected_pt"])
    print("Electron pt Scale Up:", events["Electron", "dEscaleUp"])
    print("Electron pt Scale Down:", events["Electron", "dEscaleDown"])
    print("Electron pt Smear Up:", events["Electron", "dEsigmaUp"])
    print("Electron pt Smear Down:", events["Electron", "dEsigmaDown"])

    return events

##########################
### Muon Scale Run3 ###
##########################

MUON_SCALE_FILE = {
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/muon_scale_2022preEE.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/muon_scale_2022postEE.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/muon_scale_2023preBPix.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/muon_scale_2023postBPix.json"
}

def muon_scale_run3(events, year, is_data):
    """
    Apply muon scale corrections for Run3 data only.
    Based on MuonScaRe package.
    
    See:
        - https://github.com/CMS-MUO-POG/MuonScaRe
    """
    logger.info("[Lepton Systematics] Applying muon scale corrections for year %s", year)

    required_fields = [
        ("Muon", "eta"), ("Muon", "pt"), ("Muon", "phi"), ("Muon", "charge")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        logger.warning(
            "[Lepton Systematics] Muon scale corrections will not be applied, because the following fields are missing: %s" % str(missing_fields)
        )
        return events

    # Only apply to data
    if not is_data:
        logger.info("[Lepton Systematics] Muon scale corrections are only applied to data, skipping MC")
        return events

    evaluator = correctionlib.CorrectionSet.from_file(misc_utils.expand_path(MUON_SCALE_FILE[year]))

    muons = events["Muon"]
    n_muons = awkward.num(muons)
    muons_flattened = awkward.flatten(muons)

    muon_pt = awkward.to_numpy(muons_flattened.pt)
    muon_eta = awkward.to_numpy(muons_flattened.eta)
    muon_phi = awkward.to_numpy(muons_flattened.phi)
    muon_charge = awkward.to_numpy(muons_flattened.charge)

    # Get scale correction parameters from correctionlib
    a_data = evaluator["a_data"].evaluate(muon_eta, muon_phi, "nom")
    m_data = evaluator["m_data"].evaluate(muon_eta, muon_phi, "nom")

    # Apply scale correction: pt_corr = 1 / (m/pt + charge * a)
    pt_corr = 1.0 / (m_data / muon_pt + muon_charge * a_data)

    # Filter boundaries: only apply corrections for pt in [26, 200] GeV
    pt_corr = numpy.where((muon_pt < 26) | (muon_pt > 200), muon_pt, pt_corr)
    
    # Filter out NaN entries
    pt_corr = numpy.where(numpy.isnan(pt_corr), muon_pt, pt_corr)
    
    # # Additional safety check: reject unreasonable corrections
    # pt_ratio = pt_corr / muon_pt
    # pt_corr = numpy.where((pt_ratio > 2) | (pt_ratio < 0.1) | (pt_corr < 0), muon_pt, pt_corr)

    # Store corrected pt
    events["Muon", "corrected_pt"] = awkward.unflatten(pt_corr, n_muons)
    
    logger.info("[Lepton Systematics] Muon pt before scale correction (data): mean = %.2f", awkward.mean(muons.pt))
    logger.info("[Lepton Systematics] Muon pt after scale correction (data): mean = %.2f", awkward.mean(events["Muon", "corrected_pt"]))

    return events
