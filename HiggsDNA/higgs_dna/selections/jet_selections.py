import awkward
import numpy
import vector

vector.register_awkward()

import logging
logger = logging.getLogger(__name__)

from correctionlib import _core
from higgs_dna.utils import misc_utils

JET_VETO_MAP_FILE = {
    "2016": "higgs_dna/systematics/data/2016postVFP_UL/jetvetomaps.json",
    "2016preVFP": "higgs_dna/systematics/data/2016preVFP_UL/jetvetomaps.json",
    "2016postVFP": "higgs_dna/systematics/data/2016postVFP_UL/jetvetomaps.json",
    "2017": "higgs_dna/systematics/data/2017_UL/jetvetomaps.json",
    "2018": "higgs_dna/systematics/data/2018_UL/jetvetomaps.json",
    "2022preEE": "higgs_dna/systematics/data/2022preEE_UL/jetvetomaps.json",
    "2022postEE": "higgs_dna/systematics/data/2022postEE_UL/jetvetomaps.json",
    "2023preBPix": "higgs_dna/systematics/data/2023preBPix_UL/jetvetomaps.json",
    "2023postBPix": "higgs_dna/systematics/data/2023postBPix_UL/jetvetomaps.json",
}

JET_ID_FILE = {
    "2022preEE": "jsonpog-integration/POG/JME/2022_Summer22/jetid.json",
    "2022postEE": "jsonpog-integration/POG/JME/2022_Summer22EE/jetid.json",
    "2023preBPix": "jsonpog-integration/POG/JME/2023_Summer23/jetid.json",
    "2023postBPix": "jsonpog-integration/POG/JME/2023_Summer23BPix/jetid.json",
}

DEFAULT_JETS = {
    "pt": 30.0,
    "eta": 4.7,
    "looseID": True,
}

JET_ID_EVALUATORS = {}


def _empty_jet_mask(jets):
    return jets.pt < 0


def _as_p4(objects):
    return awkward.Array(objects, with_name="Momentum4D")


def _get_photon_jet_idx(photons):
    if photons is None:
        return None
    if "jetIdx_11p9" in photons.fields:
        return awkward.fill_none(photons.jetIdx_11p9, -1)
    if "jetIdx" in photons.fields:
        return awkward.fill_none(photons.jetIdx, -1)
    return None


def _jet_lepton_overlap_mask(jets, leptons, min_dr):
    if leptons is None or awkward.count(leptons.pt) == 0:
        return _empty_jet_mask(jets)

    jet_p4 = awkward.unflatten(_as_p4(jets), counts=1, axis=-1)
    lepton_p4 = awkward.unflatten(_as_p4(leptons), counts=1, axis=0)
    dR = jet_p4.deltaR(lepton_p4)
    jet_pt = awkward.unflatten(jets.pt, counts=1, axis=-1)
    lepton_pt = awkward.unflatten(leptons.pt, counts=1, axis=0)
    pt_match = (numpy.abs(jet_pt - lepton_pt) / lepton_pt) < 1.0
    return awkward.any((dR < min_dr) & pt_match, axis=-1)


def _jet_photon_overlap_mask(jets, photons, min_dr):
    if photons is None or awkward.count(photons.pt) == 0:
        return _empty_jet_mask(jets)

    jet_p4 = awkward.unflatten(_as_p4(jets), counts=1, axis=-1)
    photon_p4 = awkward.unflatten(_as_p4(photons), counts=1, axis=0)
    dr_overlap = awkward.any(jet_p4.deltaR(photon_p4) < min_dr, axis=-1)

    photon_jet_idx = _get_photon_jet_idx(photons)
    if photon_jet_idx is None:
        return dr_overlap

    jet_idx = awkward.local_index(jets.pt, axis=1)
    jetidx_overlap = awkward.any(jet_idx[:, :, None] == photon_jet_idx[:, None, :], axis=-1)
    return dr_overlap | jetidx_overlap


def _get_run3_jet_id_bits(jets):
    if "jetId_11p9" in jets.fields:
        return jets.jetId_11p9
    if "jetId" in jets.fields:
        return jets.jetId
    return None


def _evaluate_run3_tight_id(year, jets):
    required_fields = [
        "eta",
        "chHEF",
        "neHEF",
        "chEmEF",
        "neEmEF",
        "muEF",
        "chMultiplicity",
        "neMultiplicity",
    ]
    if not all(field in jets.fields for field in required_fields):
        return None

    if year not in JET_ID_EVALUATORS:
        JET_ID_EVALUATORS[year] = _core.CorrectionSet.from_file(
            misc_utils.expand_path(JET_ID_FILE[year])
        )

    evaluator = JET_ID_EVALUATORS[year]["AK4PUPPI_Tight"]
    n_jets = awkward.num(jets)
    jets_flat = awkward.flatten(jets)
    tight_id = evaluator.evalv(
        awkward.to_numpy(jets_flat.eta),
        awkward.to_numpy(jets_flat.chHEF),
        awkward.to_numpy(jets_flat.neHEF),
        awkward.to_numpy(jets_flat.chEmEF),
        awkward.to_numpy(jets_flat.neEmEF),
        awkward.to_numpy(jets_flat.muEF),
        awkward.to_numpy(jets_flat.chMultiplicity).astype(numpy.int32),
        awkward.to_numpy(jets_flat.neMultiplicity).astype(numpy.int32),
        awkward.to_numpy(jets_flat.chMultiplicity + jets_flat.neMultiplicity).astype(numpy.int32),
    )
    return awkward.unflatten(tight_id > 0, n_jets)


def _get_run3_tight_id(year, jets):
    jet_id_bits = _get_run3_jet_id_bits(jets)
    if jet_id_bits is not None and "jetId_11p9" in jets.fields:
        return (
            ((jet_id_bits & 0b010) > 0)
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

    tight_id = _evaluate_run3_tight_id(year, jets)
    if tight_id is not None:
        return tight_id

    if jet_id_bits is None:
        return jets.pt > 0

    return (
        ((jet_id_bits & 0b010) > 0)
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


def _get_run2_tight_id(jets):
    return (jets.jetId >= 1) & ((jets.puId >= 1) | (jets.pt > 50))


def select_jets(jets, options, clean, year, name="none", tagger=None, event_runs=None, event_numbers=None):
    """
    Apply n2p-like jet keep/remove logic:
    - signal jet overlap removal with photons/leptons
    - Run3 tight jet ID using Jet_jetId_11p9 when available
    - eta horn veto, jet veto map, and HEM veto
    - final pT requirement only after isgood_min / veto decisions
    """

    options = misc_utils.update_dict(original=DEFAULT_JETS, new=options)

    eta_cut = abs(jets.eta) < options["eta"]
    pt_cut = jets.pt > options["pt"]

    photon_overlap = _empty_jet_mask(jets)
    electron_overlap = _empty_jet_mask(jets)
    muon_overlap = _empty_jet_mask(jets)
    if "photons" in clean:
        photon_overlap = _jet_photon_overlap_mask(
            jets, clean["photons"]["objects"], clean["photons"]["min_dr"]
        )
    if "electrons" in clean:
        electron_overlap = _jet_lepton_overlap_mask(
            jets, clean["electrons"]["objects"], clean["electrons"]["min_dr"]
        )
    if "muons" in clean:
        muon_overlap = _jet_lepton_overlap_mask(
            jets, clean["muons"]["objects"], clean["muons"]["min_dr"]
        )

    photon_keep = ~photon_overlap
    lepton_keep = ~(electron_overlap | muon_overlap)
    standard_cuts = pt_cut & eta_cut & photon_keep & lepton_keep
    photon_removal = pt_cut & eta_cut & photon_keep
    lepton_removal = standard_cuts
    clean_eta_cuts = eta_cut & photon_keep & lepton_keep

    jet_pt_for_veto = jets.raw_pt if "raw_pt" in jets.fields else jets.pt

    jets_horn_year = {"2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"}
    horn_cut = jets.pt > 0
    if year in jets_horn_year and "jets_horn" in options:
        horn_eta_min, horn_eta_max = options["jets_horn"]["eta"]
        horn_pt_threshold = options["jets_horn"]["pt"]
        horn_cut = ~(
            (abs(jets.eta) > horn_eta_min)
            & (abs(jets.eta) < horn_eta_max)
            & (jets.pt < horn_pt_threshold)
        )
        logger.info(
            "Excluding jets with %s < |eta| < %s and pt < %s for era %s",
            horn_eta_min,
            horn_eta_max,
            horn_pt_threshold,
            year,
        )

    if options["looseID"]:
        if int(year[:4]) >= 2022:
            id_cut = _get_run3_tight_id(year, jets)
        else:
            id_cut = _get_run2_tight_id(jets)
    else:
        id_cut = jets.pt > 0

    isgood_min = clean_eta_cuts & id_cut

    jet_veto_map_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(JET_VETO_MAP_FILE[year]))
    n_jets = awkward.num(jets)
    jets_flattened = awkward.flatten(jets)
    jet_eta = numpy.clip(awkward.to_numpy(jets_flattened.eta), -5.19099, 5.19099)
    jet_phi = numpy.clip(awkward.to_numpy(jets_flattened.phi), -3.1415925, 3.1415925)
    jet_pt_flattened_for_veto = awkward.to_numpy(awkward.flatten(jet_pt_for_veto))

    jet_veto_mask = _empty_jet_mask(jets)
    if int(year[:4]) >= 2022:
        jet_veto_values = jet_veto_map_evaluator["jetvetomap"].evalv("jetvetomap", jet_eta, jet_phi)
        jet_veto_mask_flat = (
            (jet_pt_flattened_for_veto > 15.0)
            & (abs(jet_eta) < 5.191)
            & (jet_veto_values > 0)
        )
        isgood_min_flat = awkward.to_numpy(awkward.flatten(isgood_min))
        jet_veto_mask = awkward.unflatten(jet_veto_mask_flat & isgood_min_flat, n_jets)

    hem_veto_mask = _empty_jet_mask(jets)
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
                    & isgood_min
                )
                if tagger.is_data:
                    hem_veto_mask = region & (run_broadcast >= 319077)
                else:
                    if event_numbers is not None:
                        event_broadcast = awkward.broadcast_arrays(event_numbers, jets.pt)[0]
                        hem_region = (event_broadcast % 10000) > 3564
                    else:
                        fraction = 0.6515623538907509
                        rand = numpy.random.random(len(event_runs))
                        hem_region = awkward.broadcast_arrays(rand < fraction, jets.pt)[0]
                    hem_veto_mask = region & hem_region
                if awkward.any(hem_veto_mask):
                    logger.debug(f"[HEM] Removed jets: {awkward.sum(hem_veto_mask, axis=1)[:10]}")
            else:
                logger.warning("[HEM] event_runs length mismatch; skip HEM cleaning.")

    jet_isgood_min = isgood_min & horn_cut
    jet_veto_cut = ~jet_veto_mask
    hem_mask = ~hem_veto_mask
    all_cuts = pt_cut & jet_isgood_min & jet_veto_cut & hem_mask
    jet_veto = jet_veto_cut & hem_mask

    standard_cuts_flat = awkward.flatten(standard_cuts)
    id_record_cut = standard_cuts_flat & awkward.flatten(id_cut)
    horn_record_cut = id_record_cut & awkward.flatten(horn_cut)
    jet_veto_record_cut = horn_record_cut & awkward.flatten(jet_veto_cut)
    hem_event_cut = jet_veto_record_cut & awkward.flatten(hem_mask)
    all_record_cut = hem_event_cut

    if tagger is not None:
        tagger.register_cuts(
            names=["std cuts", "id cut", "horn cut", "jet veto cut", "hem cut", "all cuts"],
            results=[standard_cuts_flat, id_record_cut, horn_record_cut, jet_veto_record_cut, hem_event_cut, all_record_cut],
            cut_type=name,
        )

    return all_cuts, jet_veto, photon_removal, lepton_removal
