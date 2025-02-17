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
    "2018" : "higgs_dna/systematics/data/2018_UL/jetvetomaps.json"
}

DEFAULT_JETS = {
    "pt" : 30.0,
    "eta" : 4.7,
    "looseID" : True
}

def select_jets(jets, options, clean, year, name = "none", tagger = None):
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

    all_cuts = standard_cuts & id_cut & photon_veto_cut & jet_veto_cut
    standard_cuts = awkward.sum(standard_cuts, axis=1) > 0
    id_cut = awkward.sum(standard_cuts & id_cut, axis=1) > 0
    photon_veto_cut = awkward.sum(standard_cuts & id_cut & photon_veto_cut, axis=1) > 0
    jet_veto_cut = awkward.sum(standard_cuts & id_cut & photon_veto_cut & jet_veto_cut, axis=1) > 0

    if tagger is not None:
        tagger.register_cuts(
            names = ["std cuts", "id cut", "photon veto cut", "jet veto cut", "all cuts"],
            results = [standard_cuts, id_cut, photon_veto_cut, jet_veto_cut, all_cuts],
            cut_type = name
        )

    return all_cuts
