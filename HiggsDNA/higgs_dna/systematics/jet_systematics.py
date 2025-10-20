import awkward

import numpy

from correctionlib import _core

from higgs_dna.utils import misc_utils, awkward_utils

import logging
logger = logging.getLogger(__name__)

###################################
### b-tag continuous reshape SF ###
###################################


BATG_MED = {
    "2016preVFP": 0.2598,
    "2016postVFP": 0.2489,
    "2017": 0.3040,
    "2018": 0.2783,
    "2022preEE": 0.3086,
    "2022postEE": 0.3196,
    "2023preBPix": 0.2431,
    "2023postBPix": 0.2435
}

BTAG_MCEFF_FILE = {
    "2016" : "higgs_dna/systematics/data/2016postVFP_UL/btag_mceff.json",
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/btag_mceff.json", 
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/btag_mceff.json", 
    "2017" : "higgs_dna/systematics/data/2017_UL/btag_mceff.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/btag_mceff.json",
    "2022preEE" : "higgs_dna/systematics/data/2022preEE_UL/btag_mceff.json",
    "2022postEE" : "higgs_dna/systematics/data/2022postEE_UL/btag_mceff.json",
    "2023preBPix" : "higgs_dna/systematics/data/2023preBPix_UL/btag_mceff.json",
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/btag_mceff.json"
}

BTAG_RESHAPE_SF_FILE = {
    "2016" : "jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json", 
    "2016preVFP" : "jsonpog-integration/POG/BTV/2016preVFP_UL/btagging.json", 
    "2016postVFP" : "jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json", 
    "2017" : "jsonpog-integration/POG/BTV/2017_UL/btagging.json",
    "2018" : "jsonpog-integration/POG/BTV/2018_UL/btagging.json",
    "2022preEE" : "jsonpog-integration/POG/BTV/2022_Summer22/btagging.json",
    "2022postEE" : "jsonpog-integration/POG/BTV/2022_Summer22EE/btagging.json",
    "2023preBPix" : "jsonpog-integration/POG/BTV/2023_Summer23/btagging.json",
    "2023postBPix" : "jsonpog-integration/POG/BTV/2023_Summer23BPix/btagging.json"
}

DEEPJET_RESHAPE_SF = {
    "2016" : "deepJet_shape", 
    "2016preVFP" : "deepJet_shape",
    "2016postVFP" : "deepJet_shape",
    "2017" : "deepJet_shape",
    "2018" : "deepJet_shape",
    "2022preEE" : "deepJet_shape", #FIXME
    "2022postEE" : "deepJet_shape", #FIXME
    "2023preBPix" : "deepJet_shape", #FIXME
    "2023postBPix" : "deepJet_shape", #FIXME
}

DEEPJET_VARIATIONS = { # b, c, light
    "up_correlated" : [5, 4, 0], 
    "down_correlated" : [5, 4, 0],
    "up_uncorrelated" : [5, 4, 0],
    "down_uncorrelated" : [5, 4, 0],
    # "up_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm (4)
    # "up_lf" : [5],
    # "up_hfstats1" : [5],
    # "up_hfstats2" : [5],
    # "up_cferr1" : [4],
    # "up_cferr2" : [4],
    # "up_hf" : [0],
    # "up_lfstats1" : [0],
    # "up_lfstats2" : [0],
    # "down_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm(4)
    # "down_lf" : [5],
    # "down_hfstats1" : [5],
    # "down_hfstats2" : [5],
    # "down_cferr1" : [4],
    # "down_cferr2" : [4],
    # "down_hf" : [0],
    # "down_lfstats1" : [0],
    # "down_lfstats2" : [0],
}


def btag_deepjet_wp_sf_heavy(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_bjets_Run2_UL/
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/btvExample.py

    Note: application of SFs should not change the overall normalization of a sample (before any b-tagging selection) and each sample should be adjusted by an overall weight derived in a phase space with no requirements on b-jets such that the normalization is unchanged. TODO: link BTV TWiki that describes this.
    """
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "phi"), (input_collection, "hadronFlavour"), (input_collection, "btagDeepFlavB") 
    ]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
    jet_mc_eff_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_MCEFF_FILE[year]))
   
    jets = events[input_collection]
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end
    logger.debug(f"Number of jets(syst): {n_jets[:10]}")
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    jet_btag_deepjet = awkward.to_numpy(jets_flattened.btagDeepFlavB)
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.39999 # SFs only valid up to eta 2.5
    )
    jet_eta = numpy.clip(
        awkward.to_numpy(jets_flattened.eta),
        -2.39999,
        2.39999
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.99999
    )

    variations_list = ["central"]
    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}

    central_sf = numpy.ones_like(jet_flavor)
    jet_mc_eff, jet_mc_eff_syst = numpy.ones_like(jet_flavor), numpy.ones_like(jet_flavor)
    for f in [4, 5]:
        evaluator_key = "deepJet_mujets" if f > 0 else "deepJet_light" if int(year[:4]) > 2020 else "deepJet_incl"
        central_sf = numpy.where(
            jet_flavor == f,
            evaluator[evaluator_key].evalv(
                "central",
                "M",
                numpy.ones_like(jet_flavor) * f,
                jet_abs_eta,
                jet_pt
            ),
            central_sf
        )
        flavour_name = "b" if f == 5 else "c" if f == 4 else "uds"
        jet_mc_eff = numpy.where(
            jet_flavor == f,
            jet_mc_eff_evaluator["Btag_{}_WPmedium_MCeff".format(flavour_name)].evalv(
                "effmc",
                jet_eta,
                jet_pt
            ),
            jet_mc_eff
        )
        jet_mc_eff_syst = numpy.where(
            jet_flavor == f,
            jet_mc_eff_evaluator["Btag_{}_WPmedium_MCeff".format(flavour_name
            )].evalv(
                "systmc",
                jet_eta,
                jet_pt
            ),
            jet_mc_eff
        )

    jet_mc_eff_up = jet_mc_eff + jet_mc_eff_syst
    jet_mc_eff_down = jet_mc_eff - jet_mc_eff_syst
    
    central_sf = numpy.where(
        jet_btag_deepjet > BATG_MED[year],
        central_sf,
        (1 - central_sf*jet_mc_eff) / (1 - jet_mc_eff)
    )

    variations["central"] = awkward.unflatten(central_sf, n_jets)

    for var in variations_list:
        if var == "central":
            continue
        applicable_flavors = DEEPJET_VARIATIONS[var] # the up/down variations are only applicable to specific flavors of jet
        var_sf = central_sf
        for f in applicable_flavors:
            if f not in [4, 5]:
                continue
            evaluator_key = "deepJet_mujets" if f > 0 else "deepJet_light" if int(year[:4]) > 2020 else "deepJet_incl"
            var_sf = numpy.where(
                jet_flavor == f,
                evaluator[evaluator_key].evalv(
                    var,
                    "M",
                    numpy.ones_like(jet_flavor) * f,
                    jet_abs_eta,
                    jet_pt
                ),
                var_sf
            )

        var_sf = numpy.where(
            jet_btag_deepjet > BATG_MED[year],
            var_sf,
            (1 - var_sf * jet_mc_eff) / (1 - (jet_mc_eff_up if "up_" in var else jet_mc_eff_down if "down_" in var else jet_mc_eff))
        )

        variations[var] = awkward.unflatten(var_sf, n_jets) # make jagged again

    for var in variations.keys():
        # Set SFs = 1 for jets which are not applicable (pt <= 20 or |eta| >= 2.5)
        variations[var] = awkward.where(
                (jets.pt <= 20.0) | (abs(jets.eta) >= 2.4),
                awkward.ones_like(variations[var]),
                variations[var]
        )
    # Swap 'up' or 'down' with the rest of the variation name
    swapped_variations = {
        (var.replace("up_", "") + "_up") if "up" in var else (var.replace("down_", "") + "_down") if "down" in var else var: val
        for var, val in variations.items()
    }
    return swapped_variations


def btag_deepjet_wp_sf_light(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_bjets_Run2_UL/
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/btvExample.py

    Note: application of SFs should not change the overall normalization of a sample (before any b-tagging selection) and each sample should be adjusted by an overall weight derived in a phase space with no requirements on b-jets such that the normalization is unchanged. TODO: link BTV TWiki that describes this.
    """
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "phi"), (input_collection, "hadronFlavour"), (input_collection, "btagDeepFlavB")
    ]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
    jet_mc_eff_evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_MCEFF_FILE[year]))

    jets = events[input_collection]
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end
    logger.debug(f"Number of jets(syst): {n_jets[:10]}")
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    jet_btag_deepjet = awkward.to_numpy(jets_flattened.btagDeepFlavB)
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.39999 # SFs only valid up to eta 2.5
    )
    jet_eta = numpy.clip(
        awkward.to_numpy(jets_flattened.eta),
        -2.39999,
        2.39999
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.99999
    )

    variations_list = ["central"]
    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}

    central_sf = numpy.ones_like(jet_flavor)
    jet_mc_eff, jet_mc_eff_syst = numpy.ones_like(jet_flavor), numpy.ones_like(jet_flavor)
    for f in [0]:
        evaluator_key = "deepJet_mujets" if f > 0 else "deepJet_light" if int(year[:4]) > 2020 else "deepJet_incl"
        central_sf = numpy.where(
            jet_flavor == f,
            evaluator[evaluator_key].evalv(
                "central",
                "M",
                numpy.ones_like(jet_flavor) * f,
                jet_abs_eta,
                jet_pt
            ),
            central_sf
        )
        flavour_name = "b" if f == 5 else "c" if f == 4 else "uds"
        jet_mc_eff = numpy.where(
            jet_flavor == f,
            jet_mc_eff_evaluator["Btag_{}_WPmedium_MCeff".format(flavour_name)].evalv(
                "effmc",
                jet_eta,
                jet_pt
            ),
            jet_mc_eff
        )
        jet_mc_eff_syst = numpy.where(
            jet_flavor == f,
            jet_mc_eff_evaluator["Btag_{}_WPmedium_MCeff".format(flavour_name
            )].evalv(
                "systmc",
                jet_eta,
                jet_pt
            ),
            jet_mc_eff
        )

    jet_mc_eff_up = jet_mc_eff + jet_mc_eff_syst
    jet_mc_eff_down = jet_mc_eff - jet_mc_eff_syst

    central_sf = numpy.where(
        jet_btag_deepjet > BATG_MED[year],
        central_sf,
        (1 - central_sf*jet_mc_eff) / (1 - jet_mc_eff)
    )

    variations["central"] = awkward.unflatten(central_sf, n_jets)

    for var in variations_list:
        if var == "central":
            continue
        applicable_flavors = DEEPJET_VARIATIONS[var] # the up/down variations are only applicable to specific flavors of jet
        var_sf = central_sf
        for f in applicable_flavors:
            if f not in [0]:
                continue
            evaluator_key = "deepJet_mujets" if f > 0 else "deepJet_light" if int(year[:4]) > 2020 else "deepJet_incl"
            var_sf = numpy.where(
                jet_flavor == f,
                evaluator[evaluator_key].evalv(
                    var,
                    "M",
                    numpy.ones_like(jet_flavor) * f,
                    jet_abs_eta,
                    jet_pt
                ),
                var_sf
            )

        var_sf = numpy.where(
            jet_btag_deepjet > BATG_MED[year],
            var_sf,
            (1 - var_sf * jet_mc_eff) / (1 - (jet_mc_eff_up if "up_" in var else jet_mc_eff_down if "down_" in var else jet_mc_eff))
        )

        variations[var] = awkward.unflatten(var_sf, n_jets) # make jagged again

    for var in variations.keys():
        # Set SFs = 1 for jets which are not applicable (pt <= 20 or |eta| >= 2.5)
        variations[var] = awkward.where(
                (jets.pt <= 20.0) | (abs(jets.eta) >= 2.4),
                awkward.ones_like(variations[var]),
                variations[var]
        )
    # Swap 'up' or 'down' with the rest of the variation name
    swapped_variations = {
        (var.replace("up_", "") + "_up") if "up" in var else (var.replace("down_", "") + "_down") if "down" in var else var: val
        for var, val in variations.items()
    }
    return swapped_variations

def btag_deepjet_mujet_sf(events, year, central_only, input_collection, working_point ="M"):
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "hadronFlavour")
    ]    
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
    
    jets = events[input_collection]   
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    is_light = jet_flavor == 0
    
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.49999 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.99
    )
    variations_list = ["central"]

    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}
    print(type(jet_flavor))
    print(type(numpy.where(is_light,4,jet_flavor)))
    central_sf = evaluator["deepJet_mujets"].evalv(
            "central",
            working_point,
            numpy.where(is_light,4,jet_flavor),
            jet_abs_eta,
            jet_pt
    )    
    central_sf = numpy.where(is_light,1,central_sf)
    variations["central"] = awkward.unflatten(central_sf, n_jets)
    if not central_only:
        syst_vars = ["up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated"]
        for syst_var in syst_vars:
            syst = evaluator["deepJet_mujets"].evalv(
                    syst_var,
                    working_point,
                    numpy.where(is_light,4,jet_flavor),
                    jet_abs_eta,
                    jet_pt,
            )
            if "up_correlated" in syst_var:
                syst_var_name = "up_correlated"
            elif "down_correlated" in syst_var:
                syst_var_name = "down_correlated"
            elif "up_uncorrelated" in syst_var:
                syst_var_name = "up_uncorrelated"
            elif "down_uncorrelated" in syst_var:
                syst_var_name = "down_uncorrelated"
            variations[syst_var_name] = awkward.unflatten(syst, n_jets)
    return variations

def btag_deepjet_incl_sf(events, year, central_only, input_collection, working_point ="M"):
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "hadronFlavour")
    ]    
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
    
    jets = events[input_collection]   
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    is_light = jet_flavor == 0

    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.49999 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.99
    )
    variations_list = ["central"]

    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()
    
    variations = {}
    central_sf = evaluator["deepJet_incl"].evalv(
            "central",
            working_point,
            numpy.where(~is_light,0,jet_flavor),
            jet_abs_eta,
            jet_pt
    )   
    central_sf = numpy.where(~is_light,1,central_sf)
    variations["central"] = awkward.unflatten(central_sf, n_jets)
    if not central_only:
        syst_vars = ["up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated"]
        for syst_var in syst_vars:
            syst = evaluator["deepJet_incl"].evalv(
                    syst_var,
                    working_point,
                    numpy.where(is_light,4,jet_flavor),
                    jet_abs_eta,
                    jet_pt,
            )
            if "up_correlated" in syst_var:
                syst_var_name = "up_correlated"
            elif "down_correlated" in syst_var:
                syst_var_name = "down_correlated"
            elif "up_uncorrelated" in syst_var:
                syst_var_name = "up_uncorrelated"
            elif "down_uncorrelated" in syst_var:
                syst_var_name = "down_uncorrelated"
            variations[syst_var_name] = awkward.unflatten(syst, n_jets)
    return variations
