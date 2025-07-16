import awkward
import numpy

from correctionlib import _core

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils

PU_REWEIGHTING_FILE = {
    "2016preVFP" : "jsonpog-integration/POG/LUM/2016preVFP_UL/puWeights.json",
    "2016postVFP" : "jsonpog-integration/POG/LUM/2016postVFP_UL/puWeights.json",
    "2017" : "jsonpog-integration/POG/LUM/2017_UL/puWeights.json",
    "2018" : "jsonpog-integration/POG/LUM/2018_UL/puWeights.json",
    "2022preEE" : "jsonpog-integration/POG/LUM/2022_Summer22/puWeights.json",
    "2022postEE" : "jsonpog-integration/POG/LUM/2022_Summer22EE/puWeights.json",
    "2023preBPix" : "jsonpog-integration/POG/LUM/2023_Summer23/puWeights.json",
    "2023postBPix" : "jsonpog-integration/POG/LUM/2023_Summer23BPix/puWeights.json"
}

PU_CAMPAIGN = {
    "2016preVFP" : "Collisions16_UltraLegacy_goldenJSON",
    "2016postVFP" : "Collisions16_UltraLegacy_goldenJSON",
    "2017" : "Collisions17_UltraLegacy_goldenJSON",
    "2018" : "Collisions18_UltraLegacy_goldenJSON",
    "2022preEE" : "Collisions2022_355100_357900_eraBCD_GoldenJson",
    "2022postEE" : "Collisions2022_359022_362760_eraEFG_GoldenJson",
    "2023preBPix" : "Collisions2023_366403_369802_eraBC_GoldenJson",
    "2023postBPix" : "Collisions2023_369803_370790_eraD_GoldenJson",
}

def pu_reweight_sf(events, year, central_only):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/LUMI_puWeights_Run2_UL/LUMI_puWeights_2017_UL.html
    """

    required_fields = ["Pileup_nTrueInt"]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PU_REWEIGHTING_FILE[year]))

    num_true_interactions = awkward.to_numpy(events.Pileup_nTrueInt).astype(float)

    # Calculate SF and syst
    variations = {}

    variations["central"] = awkward.from_numpy(evaluator[PU_CAMPAIGN[year]].evalv(
            num_true_interactions,
            "nominal"
    ))

    if not central_only:
        for var in ["up", "down"]:
            variations[var] = awkward.from_numpy(evaluator[PU_CAMPAIGN[year]].evalv(
                    num_true_interactions,
                    var
            ))

    # Cap variations at 10 to avoid extremely large weights
    for key in variations:
        variations[key] = awkward.where(variations[key] > 10, 10, variations[key])

    return variations
