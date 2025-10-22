import awkward

import numpy

from correctionlib import _core
import correctionlib

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

def pt_correction_mc(
    events,
    year,
    jes_unc_source,
    jer_variation,
    input_collection="Jet",
    jec_algo="AK4PFPuppi",
    reapply_jes=True
):
    """
    Applies jet energy scale (JES) corrections and jet energy resolution (JER) smearing to simulated jets.
    This function is a Python replica of the C++ version in jets.cxx.

    Parameters
    ----------
    events : awkward.Array
        The events array.
    year : str
        The year of the data-taking period (e.g., "2018", "2022preEE").
    jes_unc_source : str
        The JES uncertainty source to vary. If empty, no JES uncertainty is applied.
        If "All", a total uncertainty is applied.
    jer_variation : str
        The JER variation ("central", "up", "down").
    input_collection : str, optional
        The name of the jet collection, by default "Jet".
    jec_algo : str, optional
        The jet algorithm, by default "AK4PFchs".
    reapply_jes : bool, optional
        Whether to reapply JES corrections, by default False.
    is_data : bool, optional
        Whether the input is data or MC, by default False.

    Returns
    -------
    awkward.Array
        An array with the corrected jet pT.
    """
    JEC_TAG_MC = {
        "2016preVFP": "Summer19UL16APV_V7_MC",
        "2016postVFP": "Summer19UL16_V7_MC",
        "2017": "Summer19UL17_V5_MC",
        "2018": "Summer19UL18_V5_MC",
        "2022preEE": "Summer22_22Sep2023_V2_MC",
        "2022postEE": "Summer22EE_22Sep2023_V2_MC",
        "2023preBPix": "Summer23Prompt23_V2_MC",
        "2023postBPix": "Summer23BPixPrompt23_V3_MC",
    }
    JEC_FILE = {
        "2016preVFP": "jsonpog-integration/POG/JME/2016preVFP_UL/jet_jerc.json",
        "2016postVFP": "jsonpog-integration/POG/JME/2016postVFP_UL/jet_jerc.json",
        "2017": "jsonpog-integration/POG/JME/2017_UL/jet_jerc.json",
        "2018": "jsonpog-integration/POG/JME/2018_UL/jet_jerc.json",
        "2022preEE": "jsonpog-integration/POG/JME/2022_Summer22/jet_jerc.json",
        "2022postEE": "jsonpog-integration/POG/JME/2022_Summer22EE/jet_jerc.json",
        "2023preBPix": "jsonpog-integration/POG/JME/2023_Summer23/jet_jerc.json",
        "2023postBPix": "jsonpog-integration/POG/JME/2023_Summer23BPix/jet_jerc.json",
    }

    jets = events[input_collection]
    
    lhc_run = 2
    if "2022" in year or "2023" in year:
        lhc_run = 3

    jes_tag = JEC_TAG_MC[year]
    jer_tag = jes_tag.replace("_V", "_JRV") if lhc_run == 2 else jes_tag.replace("_V", "_JR_V")
    jec_file = misc_utils.expand_path(JEC_FILE[year])
    evaluator = _core.CorrectionSet.from_file(jec_file)

    # Flatten jets for easier processing
    n_jets = awkward.num(jets)
    jets_flat = awkward.flatten(jets)

    corrected_pts = jets_flat.pt

    if reapply_jes:
        raw_pt = jets_flat.pt * (1 - jets_flat.rawFactor)
        jes_evaluator = evaluator.compound[f"{jes_tag}_L1L2L3Res_{jec_algo}"]
        corr = jes_evaluator.evaluate(
            jets_flat.area, jets_flat.eta, raw_pt, events.fixedGridRhoFastjetAll
        )
        corrected_pts = raw_pt * corr

    # JER smearing
    jer_resolution_evaluator = evaluator[f"{jer_tag}_PtResolution_{jec_algo}"]
    reso = jer_resolution_evaluator.evalv(jets_flat.eta, corrected_pts, events.fixedGridRhoFastjetAll)

    jer_sf_evaluator = evaluator[f"{jer_tag}_ScaleFactor_{jec_algo}"]
    jer_var_map = {"central": "nom", "up": "up", "down": "down"}
    jer_shift_str = jer_var_map[jer_variation]

    if lhc_run == 2:
        reso_sf = jer_sf_evaluator.evalv(jets_flat.eta, jer_shift_str)
    else: # lhc_run == 3
        reso_sf = jer_sf_evaluator.evalv(jets_flat.eta, corrected_pts, jer_shift_str)

    # Gen jet matching
    gen_jets = events.GenJet
    jets_flat_cartesian = awkward.zip({"pt": jets_flat.pt, "eta": jets_flat.eta, "phi": jets_flat.phi, "mass": jets_flat.mass}, with_name="Momentum4D")
    gen_jets_cartesian = awkward.zip({"pt": gen_jets.pt, "eta": gen_jets.eta, "phi": gen_jets.phi, "mass": gen_jets.mass}, with_name="Momentum4D")
    
    matched_gen_jets = awkward.firsts(gen_jets_cartesian[awkward.singletons(jets_flat_cartesian.delta_r(gen_jets_cartesian) < 0.2)])
    
    gen_jet_pt = awkward.where(
        (matched_gen_jets.pt > 0) & (abs(corrected_pts - matched_gen_jets.pt) < (3.0 * reso * corrected_pts)),
        matched_gen_jets.pt,
        -1.0
    )

    # Apply smearing
    smear_factor = numpy.ones_like(corrected_pts, dtype=numpy.float32)
    
    # Scaling method for matched jets
    matched_mask = gen_jet_pt > 0
    shift = (reso_sf - 1.0) * (corrected_pts - gen_jet_pt) / corrected_pts
    smear_factor = numpy.where(matched_mask, numpy.maximum(0.0, 1.0 + shift), smear_factor)

    # Stochastic method for non-matched jets
    not_matched_mask = ~matched_mask
    seed = events.event
    randm = numpy.random.RandomState(seed=seed)
    gaus_smear = randm.normal(0, reso)
    shift_stochastic = gaus_smear * numpy.sqrt(numpy.maximum(reso_sf * reso_sf - 1.0, 0.0))
    smear_factor = numpy.where(not_matched_mask, numpy.maximum(0.0, 1.0 + shift_stochastic), smear_factor)

    corrected_pts = corrected_pts * smear_factor

    # JES uncertainties
    pt_scale_sf = numpy.ones_like(corrected_pts, dtype=numpy.float32)
    if jes_unc_source:
        jes_shift = 1.0 if "up" in jes_unc_source.lower() else -1.0
        source_name = jes_unc_source.replace("_up", "").replace("_down", "")

        if source_name == "HEMIssue" and year == "2018":
            if jes_shift == -1.0:
                hem_mask = (corrected_pts > 15.0) & (jets_flat.phi > -1.57) & (jets_flat.phi < -0.87) & (jets_flat.jetId == 2)
                eta_mask1 = (jets_flat.eta > -2.5) & (jets_flat.eta < -1.3)
                eta_mask2 = (jets_flat.eta > -3.0) & (jets_flat.eta <= -2.5)
                pt_scale_sf = numpy.where(hem_mask & eta_mask1, 0.8, pt_scale_sf)
                pt_scale_sf = numpy.where(hem_mask & eta_mask2, 0.65, pt_scale_sf)
        else:
            # Assuming 'Total' for combined sources, or a specific source name
            if source_name == "Total":
                # This part requires listing all sources and summing in quadrature.
                # For simplicity, let's assume a 'Total' uncertainty is available.
                # This might need adjustment based on what's in the JEC file.
                try:
                    jes_source_evaluator = evaluator[f"{jes_tag}_{source_name}_{jec_algo}"]
                    unc = jes_source_evaluator.evalv(jets_flat.eta, corrected_pts)
                    pt_scale_sf = 1.0 + jes_shift * unc
                except:
                    logger.warning(f"Could not find total JES uncertainty '{jes_tag}_{source_name}_{jec_algo}'. Using individual sources.")
                    # Fallback to summing individual sources if 'Total' is not present
                    all_sources = [s.name for s in evaluator if f"{jes_tag}_" in s.name and f"_{jec_algo}" in s.name and "_L1" not in s.name] # simplified
                    quad_sum = numpy.zeros_like(corrected_pts, dtype=numpy.float32)
                    for source in all_sources:
                        if "PtResolution" in source or "ScaleFactor" in source: continue
                        s_eval = evaluator[source]
                        quad_sum += numpy.power(s_eval.evalv(jets_flat.eta, corrected_pts), 2.0)
                    pt_scale_sf = 1.0 + jes_shift * numpy.sqrt(quad_sum)

            else: # single source
                jes_source_evaluator = evaluator[f"{jes_tag}_{source_name}_{jec_algo}"]
                unc = jes_source_evaluator.evalv(jets_flat.eta, corrected_pts)
                pt_scale_sf = 1.0 + jes_shift * unc

    corrected_pts = corrected_pts * pt_scale_sf

    return awkward.unflatten(corrected_pts, n_jets)


def pt_correction_data(
    events,
    year,
    run_period=None,
    input_collection="Jet",
    jec_algo="AK4PFPuppi"
):
    """
    Applies jet energy corrections (JEC) to jets in real data.
    This function is a Python replica of the C++ version in jets.cxx.

    Parameters
    ----------
    events : awkward.Array
        The events array.
    year : str
        The year of the data-taking period (e.g., "2018", "2022preEE").
    run_period : str
        The data-taking run period (e.g., "A", "B", "C").
    input_collection : str, optional
        The name of the jet collection, by default "Jet".
    jec_algo : str, optional
        The jet algorithm, by default "AK4PFchs".

    Returns
    -------
    awkward.Array
        An array with the corrected jet pT.
    """
    JEC_TAG_DATA = {
        "2016preVFP": "Summer19UL16APV_V7_DATA",
        "2016postVFP": "Summer19UL16_V7_DATA",
        "2017": "Summer19UL17_V5_DATA",
        "2018": "Summer19UL18_V5_DATA",
        "2022preEE": "Summer22_22Sep2023_RunCD_V2_DATA",
        "2022postEE": "Summer22EE_22Sep2023_V2_DATA",
        "2023preBPix": "Summer23Prompt23_V2_DATA",
        "2023postBPix": "Summer23BPixPrompt23_V3_DATA",
    }
    JEC_FILE = {
        "2016preVFP": "jsonpog-integration/POG/JME/2016preVFP_UL/jet_jerc.json",
        "2016postVFP": "jsonpog-integration/POG/JME/2016postVFP_UL/jet_jerc.json",
        "2017": "jsonpog-integration/POG/JME/2017_UL/jet_jerc.json",
        "2018": "jsonpog-integration/POG/JME/2018_UL/jet_jerc.json",
        "2022preEE": "jsonpog-integration/POG/JME/2022_Summer22/jet_jerc.json",
        "2022postEE": "jsonpog-integration/POG/JME/2022_Summer22EE/jet_jerc.json",
        "2023preBPix": "jsonpog-integration/POG/JME/2023_Summer23/jet_jerc.json",
        "2023postBPix": "jsonpog-integration/POG/JME/2023_Summer23BPix/jet_jerc.json",
    }

    jets = events[input_collection]
    
    if run_period:
        jes_tag = f"{JEC_TAG_DATA[year][:-8]}_Run{run_period}_{JEC_TAG_DATA[year][-7:]}"
    else:
        jes_tag = JEC_TAG_DATA[year]
    jec_file = misc_utils.expand_path(JEC_FILE[year])
    evaluator = correctionlib.CorrectionSet.from_file(jec_file)

    n_jets = awkward.num(jets)
    jets_flat = awkward.flatten(jets)
    jets_area = awkward.to_numpy(jets_flat.area)
    jets_eta = awkward.to_numpy(jets_flat.eta)
    rho_arr = numpy.repeat(awkward.to_numpy(events["Rho_fixedGridRhoFastjetAll"]), n_jets)

    raw_pt = jets_flat.pt * (1 - jets_flat.rawFactor)
    raw_pt = awkward.to_numpy(raw_pt)
    if "2023" in year:
        run_arr = numpy.repeat(awkward.to_numpy(events["run"]), n_jets)
        jets_phi = awkward.to_numpy(jets_flat.phi)
    
    jec_evaluator = evaluator.compound[f"{jes_tag}_L1L2L3Res_{jec_algo}"]
    if "2023post" in year:
        corr = jec_evaluator.evaluate(
            jets_area, jets_eta, raw_pt, rho_arr, jets_phi, run_arr
        )
    elif "2023pre" in year:
        corr = jec_evaluator.evaluate(
            jets_area, jets_eta, raw_pt, rho_arr, run_arr
        )
    else:
        corr = jec_evaluator.evaluate(
            jets_area, jets_eta, raw_pt, rho_arr
        )
    corrected_pts = raw_pt * corr
    events["Jet", "corrected_pt"] = awkward.unflatten(corrected_pts, n_jets)

    return events
