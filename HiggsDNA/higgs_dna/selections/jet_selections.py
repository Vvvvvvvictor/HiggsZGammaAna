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

def select_jets(jets, options, clean, name = "none", tagger = None):
    """

    """
    options = misc_utils.update_dict(
        original = DEFAULT_JETS,
        new = options
    )

    tagger_name = "none" if tagger is None else tagger.name 

    standard_cuts = object_selections.select_objects(jets, options, clean, name, tagger)

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
