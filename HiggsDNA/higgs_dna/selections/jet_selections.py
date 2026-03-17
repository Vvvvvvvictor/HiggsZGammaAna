import awkward
import numpy

import logging
logger = logging.getLogger(__name__)

from correctionlib import _core
from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

JET_VETO_MAP_FILE = {
    "2016" : "higgs_dna/systematics/data/2016postVFP_UL/jetvetomaps.json",
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/jetvetomaps.json",
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/jetvetomaps.json",
    "2017" : "higgs_dna/systematics/data/2017_UL/jetvetomaps.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/jetvetomaps.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/jetvetomaps.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/jetvetomaps.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/jetvetomaps.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/jetvetomaps.json"
}

DEFAULT_JETS = {
    "pt" : 30.0,
    "eta" : 4.7,
    "looseID" : True
}

def select_jets(jets, options, clean, year, name = "none", tagger = None, event_runs = None, event_numbers = None):
    """

    """
    # Replace DEFAULT_JETS to the jet selection in zgamma_tagger_run2
    options = misc_utils.update_dict(
        original = DEFAULT_JETS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name 

    standard_cuts, photon_removal, lepton_removal = object_selections.select_objects(jets, options, clean, name, tagger)
    clean_eta_cuts, _, _ = object_selections.select_objects(
        jets,
        {"eta": options["eta"]},
        clean,
        name,
        None
    )
    eta_cut = abs(jets.eta) < options["eta"]
    standard_cuts = standard_cuts & eta_cut
    clean_eta_cuts = clean_eta_cuts & eta_cut

    jet_pt_for_veto = jets.raw_pt if "raw_pt" in jets.fields else jets.pt

    # Pei-Zhu, Jet Horn 
    jets_horn_year = {"2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"}
    horn_cut = jets.pt > 0 # default: keep all
    if year in jets_horn_year and "jets_horn" in options:
        horn_eta_min, horn_eta_max = options["jets_horn"]["eta"]
        horn_pt_threshold = options["jets_horn"]["pt"]
        # Identify jets to exclude: 2.5 < |eta| < 3.0 and pt < 40.0
        horn_cut = ~((abs(jets.eta) > horn_eta_min) & (abs(jets.eta) < horn_eta_max) & (jets.pt < horn_pt_threshold))
        logger.info(f"Excluding jets with {horn_eta_min} < |eta| < {horn_eta_max} and pt < {horn_pt_threshold} "
                    f"for era {year}")
    
    # TODO: jet ID
    if options["looseID"]:
        if int(year[:4]) >= 2022:
            id_cut = (
                ((jets.jetId & 0b010) > 0)
                & (
                    (abs(jets.eta) <= 2.7)
                    | (
                        (abs(jets.eta) > 2.7)
                        & (abs(jets.eta) <= 3.0)
                        & (jets.neHEF < 0.99)
                    )
                    | (
                        (abs(jets.eta) > 3.0)
                        & (jets.neEmEF < 0.4)
                    )
                )
            )
        else:
            id_cut = jets.jetId >= 1 & ((jets.puId >= 1) | (jets.pt > 50)) # jetID stored bitwise for loose/tight/tightLepVeto
    else:
        id_cut = jets.pt > 0

    jet_veto_map_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(JET_VETO_MAP_FILE[year]))
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)
    jet_eta = numpy.clip(
        awkward.to_numpy(jets_flattened.eta),
        -5.19099,
        5.19099
    )
    jet_pt_flattened_for_veto = awkward.to_numpy(awkward.flatten(jet_pt_for_veto))
    jet_phi = numpy.clip(
        awkward.to_numpy(jets_flattened.phi),
        -3.1415925,
        3.1415925
    )

    jet_veto_cut = jets.pt > 0  # default: keep all
    jet_event_veto_cut = jets.pt > 0
    if int(year[:4]) >= 2022:
        jet_veto_values = jet_veto_map_evaluator["jetvetomap"].evalv(
            "jetvetomap",
            jet_eta,
            jet_phi
        )
        jet_veto_mask = (
            (jet_pt_flattened_for_veto > 15.0)
            & (abs(jet_eta) < 5.191)
            & (jet_veto_values > 0)
        )
        jet_veto_sf = numpy.where(
            jet_veto_mask,
            False,
            True
        )
        jet_veto_cut = awkward.unflatten(jet_veto_sf, n_jets)
        isgood_min_flat = awkward.to_numpy(awkward.flatten(clean_eta_cuts & id_cut))
        jet_event_veto_cut = awkward.unflatten(
            numpy.where(jet_veto_mask & isgood_min_flat, False, True),
            n_jets
        )

    # ---- 2018 HEM jet-level cleaning (moved from event-level) ----
    hem_mask = jets.pt > 0  # default: keep all
    if year == "2018":
        if tagger is not None and event_runs is not None:
            if len(jets) == len(event_runs):
                run_broadcast = awkward.broadcast_arrays(event_runs, jets.pt)[0]
                region = (
                    (jet_pt_for_veto > 15.0)
                    & (jets.phi > -1.57)
                    & (jets.phi < -0.87)
                    & (jets.eta > -3.2)
                    & (jets.eta < -1.3)
                )
                if tagger.is_data:
                    run_region = run_broadcast >= 319077
                    hem_mask = ~(region & run_region)
                else:
                    if event_numbers is not None:
                        event_broadcast = awkward.broadcast_arrays(event_numbers, jets.pt)[0]
                        hem_region = (event_broadcast % 10000) > 3564
                    else:
                        fraction = 0.6515623538907509
                        rand = numpy.random.random(len(event_runs))
                        hem_region = awkward.broadcast_arrays(rand < fraction, jets.pt)[0]
                    hem_mask = ~(region & hem_region)
                if awkward.any(~hem_mask):
                    logger.debug(f"[HEM] Removed jets: {awkward.sum(~hem_mask, axis=1)[:10]}")
            else:
                logger.warning("[HEM] event_runs length mismatch; skip HEM cleaning.")

    all_cuts = standard_cuts & horn_cut & id_cut & jet_veto_cut & hem_mask

    jet_veto = jet_event_veto_cut & hem_mask

    standard_cuts = awkward.flatten(standard_cuts)
    id_cut = standard_cuts & awkward.flatten(id_cut)
    horn_cut = id_cut & awkward.flatten(horn_cut)
    jet_veto_cut = horn_cut & awkward.flatten(jet_veto_cut)
    hem_event_cut = jet_veto_cut & awkward.flatten(hem_mask)
    all_record_cut = hem_event_cut

    if tagger is not None:
        tagger.register_cuts(
            names = ["std cuts", "id cut", "horn cut", "jet veto cut", "hem cut", "all cuts"],
            results = [standard_cuts, id_cut, horn_cut, jet_veto_cut, hem_event_cut, all_record_cut],
            cut_type = name
        )

    return all_cuts, jet_veto, photon_removal, lepton_removal
