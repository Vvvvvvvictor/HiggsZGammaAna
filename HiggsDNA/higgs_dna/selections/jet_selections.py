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

def select_jets(jets, options, clean, year, name = "none", tagger = None, event_runs = None):
    """

    """
    # Replace DEFAULT_JETS to the jet selection in zgamma_tagger_run2
    options = misc_utils.update_dict(
        original = DEFAULT_JETS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name 

    standard_cuts = object_selections.select_objects(jets, options, clean, name, tagger)

    # Pei-Zhu, Jet Horn 
    jets_horn_year = {"2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"}
    if year in jets_horn_year and "jets_horn" in options:
        horn_eta_min, horn_eta_max = options["jets_horn"]["eta"]
        horn_pt_threshold = options["jets_horn"]["pt"]
        # Identify jets to exclude: 2.5 < |eta| < 3.0 and pt < 40.0
        horn_exclusion = (abs(jets.eta) > horn_eta_min) & (abs(jets.eta) < horn_eta_max) & (jets.pt < horn_pt_threshold)
        # Apply the exclusion by inverting it
        standard_cuts = standard_cuts & ~horn_exclusion
        logger.info(f"Excluding jets with {horn_eta_min} < |eta| < {horn_eta_max} and pt < {horn_pt_threshold} "
                    f"for era {year}")

    # TODO: jet ID
    if options["looseID"]:
        id_cut = jets.jetId >= 1 # jetID stored bitwise for loose/tight/tightLepVeto
    else:
        id_cut = jets.pt > 0
        
    if "photons" in clean:
        photons = clean["photons"]["objects"]
        jet_idx = awkward.local_index(jets.pt, axis=1)
        new_jet = awkward.unflatten(awkward.unflatten(awkward.flatten(jet_idx), [1]*awkward.sum(awkward.num(jet_idx))), awkward.num(jet_idx, axis=1))
        # Replace empty lists in photons.jetIdx with [-1]
        photons_jetIdx = awkward.where(awkward.num(photons.jetIdx, axis=1) == 0, awkward.ones_like(awkward.num(photons.jetIdx, axis=1))*-1, photons.jetIdx)
        new_pho = awkward.broadcast_arrays(photons_jetIdx[:, None], new_jet, depth_limit=2)[0]
        photon_veto_cut = ~awkward.flatten(awkward.any(new_jet[:, :, None] == new_pho, axis=2), axis=-1)
    else:
        photon_veto_cut = jets.pt > 0

    jet_veto_map_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(JET_VETO_MAP_FILE[year]))
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)
    jet_eta = numpy.clip(
        awkward.to_numpy(jets_flattened.eta),
        -5.19099,
        5.19099
    )
    jet_phi = numpy.clip(
        awkward.to_numpy(jets_flattened.phi),
        -3.1415925,
        3.1415925
    )
    jet_veto_sf = numpy.where(
        jet_veto_map_evaluator["jetvetomap"].evalv(
            "jetvetomap",
            jet_eta,
            jet_phi
        ) > 0,
        False,
        True
    )
    jet_veto_sf = numpy.where(
        (abs(jet_eta) >= 5.191) | (abs(jet_phi) >= 3.1415926),
        True,
        jet_veto_sf
    )
    jet_veto_cut = awkward.unflatten(jet_veto_sf, n_jets)

    # ---- 2018 HEM jet-level cleaning (moved from event-level) ----
    hem_mask = jets.pt > 0  # default: keep all
    if year == "2018":
        if tagger is not None and event_runs is not None:
            # broadcast event-level info to jet dimension
            if len(jets) == len(event_runs):
                run_broadcast = awkward.broadcast_arrays(event_runs, jets.pt)[0]
                region = ((jets.phi > -1.57) & (jets.phi < -0.87) &
                          (jets.eta > -3.0) & (jets.eta < -1.3))
                if tagger.is_data:
                    run_region = (run_broadcast > 319077)
                    remove = region & run_region
                else:
                    fraction = 0.6515623538907509
                    # 產生一次 per-event 隨機數 (可重複執行時保持非決定性；若需可加種子)
                    rand = numpy.random.random(len(event_runs))
                    hem_run = rand < fraction
                    hem_run_broadcast = awkward.broadcast_arrays(hem_run, jets.pt)[0]
                    remove = region & hem_run_broadcast
                hem_mask = ~remove
                if awkward.any(~hem_mask):
                    logger.debug(f"[HEM] Removed jets: {awkward.sum(~hem_mask, axis=1)[:10]}")
            else:
                logger.warning("[HEM] event_runs length mismatch; skip HEM cleaning.")

    all_cuts = standard_cuts & id_cut & photon_veto_cut & jet_veto_cut & hem_mask

    standard_cuts = awkward.sum(standard_cuts, axis=1) > 0
    id_cut = awkward.sum(standard_cuts & id_cut, axis=1) > 0
    photon_veto_cut = awkward.sum(standard_cuts & id_cut & photon_veto_cut, axis=1) > 0
    jet_veto_cut = awkward.sum(standard_cuts & id_cut & photon_veto_cut & jet_veto_cut, axis=1) > 0
    hem_event_cut = awkward.sum(standard_cuts & id_cut & photon_veto_cut & jet_veto_cut & hem_mask, axis=1) > 0

    if tagger is not None:
        tagger.register_cuts(
            names = ["std cuts", "id cut", "photon veto cut", "jet veto cut", "hem cut", "all cuts"],
            results = [standard_cuts, id_cut, photon_veto_cut, jet_veto_cut, hem_event_cut, all_cuts],
            cut_type = name
        )

    return all_cuts
