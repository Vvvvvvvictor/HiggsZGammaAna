import awkward

import logging
logger = logging.getLogger(__name__)

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils

DEFAULT_JETS = {
    "pt" : 30.0,
    "eta" : 4.7,
    "looseID" : True
}

def select_jets(jets, options, clean, name = "none", tagger = None, year=None):
    """

    """
    # Replace DEFAULT_JETS to the jet selection in zgamma_tagger_run2
    options = misc_utils.update_dict(
        original = DEFAULT_JETS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name         

    standard_cuts = object_selections.select_objects(jets, options, clean, name, tagger)

    jets_horn_year = {"2022preEE", "2022postEE", "2023preBPix", "2023postBPix"}
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
        # Replace empty lists in photons.jetIdx with [-1], safely handling depth
        # print("Photon fields:", photons.fields)
        if "jetIdx" in photons.fields:
            photons_jetIdx = awkward.where(awkward.num(photons.jetIdx, axis=1) == 0, awkward.ones_like(awkward.num(photons.jetIdx, axis=1))*-1, photons.jetIdx)
        else:
            # Default to -1 (no jet association) if jetIdx is missing
            photons_jetIdx = awkward.ones_like(awkward.num(photons.pt, axis=1)) * -1
        new_pho = awkward.broadcast_arrays(photons_jetIdx[:, None], new_jet, depth_limit=2)[0]
        photon_veto_cut = ~awkward.flatten(awkward.any(new_jet[:, :, None] == new_pho, axis=2), axis=-1)
    else:
        photon_veto_cut = jets.pt > 0

    all_cuts = standard_cuts & id_cut & photon_veto_cut
    standard_cuts = awkward.sum(standard_cuts, axis=1) > 0
    id_cut = awkward.sum(standard_cuts & id_cut, axis=1) > 0
    photon_veto_cut = awkward.sum(standard_cuts & id_cut & photon_veto_cut, axis=1) > 0

    if tagger is not None:
        tagger.register_cuts(
            names = ["std cuts", "id cut", "photon veto cut", "all cuts"],
            results = [standard_cuts, id_cut, photon_veto_cut, all_cuts],
            cut_type = name
        )

    return all_cuts
