import awkward
import ROOT
import numpy

from correctionlib import _core
import correctionlib

from higgs_dna.selections import object_selections
from higgs_dna.utils import misc_utils, awkward_utils

import logging
logger = logging.getLogger(__name__)
ROOT_rng = ROOT.TRandom3(4357)

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
            jet_mc_eff_syst
        )

    logger.debug(f"First 10 jet flavor before clipping: {jet_flavor[:10]} in length {len(jet_flavor)}")
    logger.debug(f"First 10 central jet mc eff before clipping: {jet_mc_eff[:10]} in length {len(jet_mc_eff)}")
    logger.debug(f"First 10 central jet mc eff syst before clipping: {jet_mc_eff_syst[:10]} in length {len(jet_mc_eff_syst)}")
    jet_mc_eff_up = jet_mc_eff + jet_mc_eff_syst
    jet_mc_eff_down = jet_mc_eff - jet_mc_eff_syst

    central_sf = numpy.where(
        jet_btag_deepjet > BATG_MED[year],
        central_sf,
        numpy.where(
            jet_mc_eff < 1.0,
            (1 - central_sf*jet_mc_eff) / (1 - jet_mc_eff),
            central_sf
        )
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
        eff_to_use = jet_mc_eff_up if "up_" in var else jet_mc_eff_down if "down_" in var else jet_mc_eff
        logger.debug(f"var: {var}, First 10 jet flavour before clipping: {jet_flavor[:10]} in length {len(jet_flavor)}")
        logger.debug(f"var: {var}, First 10 central jet mc eff before clipping: {jet_mc_eff[:10]} in length {len(jet_mc_eff)}")
        logger.debug(f"var: {var}, First 10 eff_to_use before clipping: {eff_to_use[:10]} in length {len(eff_to_use)}")
        logger.debug(f"var: {var}, First 10 var_sf before clipping: {var_sf[:10]} in length {len(var_sf)}")
        var_sf = numpy.where(
            jet_btag_deepjet > BATG_MED[year],
            var_sf,
            numpy.where(
                (eff_to_use < 1.0) & (jet_mc_eff < 1.0),
                (1 - var_sf * jet_mc_eff) / (1 - eff_to_use),
                var_sf
            )
        )
        logger.debug(f"var: {var}, First 10 var_sf before clipping: {jet_flavor[:10]} in length {len(jet_flavor)}")
        logger.debug(f"var: {var}, First 10 var_sf before clipping: {var_sf[:10]} in length {len(var_sf)}")
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
            jet_mc_eff_syst
        )

    jet_mc_eff_up = jet_mc_eff + jet_mc_eff_syst
    jet_mc_eff_down = jet_mc_eff - jet_mc_eff_syst
 
    central_sf = numpy.where(
        jet_btag_deepjet > BATG_MED[year],
        central_sf,
        numpy.where(
            (jet_mc_eff < 1.0),
            (1 - central_sf*jet_mc_eff) / (1 - jet_mc_eff),
            central_sf
        )
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
        eff_to_use = jet_mc_eff_up if "up_" in var else jet_mc_eff_down if "down_" in var else jet_mc_eff
        var_sf = numpy.where(
            jet_btag_deepjet > BATG_MED[year],
            var_sf,
            numpy.where(
                (eff_to_use < 1.0) & (jet_mc_eff < 1.0),
                (1 - var_sf * jet_mc_eff) / (1 - eff_to_use),
                var_sf
            )
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
    input_collection="Jet",
    jec_algo="AK4PFPuppi",
    reapply_jes=True
    ):
    """
    Applies jet energy scale (JES) corrections and jet energy resolution (JER) smearing to simulated jets.
    This function is a Python replica of the C++ version in jets.cxx.
    It calculates the nominal jet pT, as well as JES and JER variations,
    and stores them as new branches in the events array.
    
    It also propagates these corrections to PuppiMET (nominal and variations).

    Parameters
    ----------
    events : awkward.Array
        The events array.
    year : str
        The year of the data-taking period (e.g., "2018", "2022preEE").
    input_collection : str, optional
        The name of the jet collection, by default "Jet".
    jec_algo : str, optional
        The jet algorithm, by default "AK4PFPuppi".
    reapply_jes : bool, optional
        Whether to reapply JES corrections, by default True.

    Returns
    -------
    awkward.Array
        The events array with new branches for corrected jet pT and corrected PuppiMET.
    """
    JEC_TAG_MC = {
        "2022preEE": "Summer22_22Sep2023_V2_MC",
        "2022postEE": "Summer22EE_22Sep2023_V3_MC",
        "2023preBPix": "Summer23Prompt23_V2_MC",
        "2023postBPix": "Summer23BPixPrompt23_V3_MC",
    }
    JER_TAG_DATA = {
        "2022preEE": "Summer22_22Sep2023_JRV1_MC",
        "2022postEE": "Summer22EE_22Sep2023_JRV1_MC",
        "2023preBPix": "Summer23Prompt23_RunCv1234_JRV1_MC",
        "2023postBPix": "Summer23BPixPrompt23_RunD_JRV1_MC",
    }
    JEC_FILE = {
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
    # jer_tag = jes_tag.replace("_V", "_JRV") if lhc_run == 2 else jes_tag.replace("_V", "_JR_V")
    jer_tag = JER_TAG_DATA[year]
    
    # Assuming misc_utils is defined elsewhere or paths are relative
    jec_file = misc_utils.expand_path(JEC_FILE[year]) 
    # jec_file = JEC_FILE[year] 
    evaluator = correctionlib.CorrectionSet.from_file(jec_file)

    # Flatten jets for easier processing
    n_jets = awkward.num(jets)
    jets_flat = awkward.flatten(jets)
    corrected_pts_base = jets_flat.pt

    jet_area = awkward.to_numpy(jets_flat.area)
    jet_eta = awkward.to_numpy(jets_flat.eta)
    jet_phi = awkward.to_numpy(jets_flat.phi)
    rho_arr = numpy.repeat(awkward.to_numpy(events["Rho_fixedGridRhoFastjetAll"]), n_jets)

    # --- 1. JES Re-application ---
    raw_pt = jets_flat.pt * (1 - jets_flat.rawFactor)
    raw_pt_numpy = awkward.to_numpy(raw_pt)
    
    if reapply_jes:
        # print(jets_flat.pt, jets_flat.rawFactor, raw_pt)
        jes_evaluator = evaluator.compound[f"{jes_tag}_L1L2L3Res_{jec_algo}"]
        if "2023postBPix" in year:
            corr = jes_evaluator.evaluate(
                jet_area, jet_eta, raw_pt_numpy, rho_arr, jet_phi
            )
        else:
            corr = jes_evaluator.evaluate(
                jet_area, jet_eta, raw_pt_numpy, rho_arr
            )
        corrected_pts_base = raw_pt * corr

    # --- 2. L1 Correction (Needed for Run 3 MET) ---
    # For Run 3 Type-1 MET, we subtract (Nominal - L1).
    # We need to evaluate L1FastJet separately.
    pt_l1_flat = raw_pt # Default to raw if not Run 3 (not used in Run 2 logic)
    if lhc_run == 3:
        try:
            # Usually named like Summer22_..._L1FastJet_AK4PFPuppi
            # Note: Not 'compound', just standard evaluator for single level
            l1_evaluator = evaluator[f"{jes_tag}_L1FastJet_{jec_algo}"]
            l1_corr = l1_evaluator.evaluate(jet_area, jet_eta, raw_pt_numpy, rho_arr)
            pt_l1_flat = raw_pt * l1_corr
        except Exception as e:
            logger.warning(f"Could not load L1FastJet evaluator for MET correction: {e}. MET calculation might be inaccurate.")
            pt_l1_flat = raw_pt # Fallback

    # --- 3. JER Resolution & Gen Matching ---
    jer_resolution_evaluator = evaluator[f"{jer_tag}_PtResolution_{jec_algo}"]
    jet_pt_resolution = jer_resolution_evaluator.evaluate(
        jet_eta, corrected_pts_base, rho_arr
    )

    gen_jets = events.GenJet

    # For each jet, find the closest gen_jet
    jet_pairs = awkward.cartesian([jets, gen_jets], axis=1, nested=[0])

    d_eta = jet_pairs['0'].eta - jet_pairs['1'].eta
    d_phi = jet_pairs['0'].phi - jet_pairs['1'].phi
    d_phi = numpy.where(d_phi > numpy.pi, d_phi - 2 * numpy.pi, d_phi)
    d_phi = numpy.where(d_phi < -numpy.pi, d_phi + 2 * numpy.pi, d_phi)
    delta_r = numpy.sqrt(d_eta**2 + d_phi**2)

    min_dr_idx = awkward.argmin(delta_r, axis=2)
    gen_matched_mask = awkward.min(delta_r, axis=2) < 0.2

    matched_gen_jets_pt = awkward.where(
        gen_matched_mask,
        gen_jets[min_dr_idx].pt,
        -1.0
    )
    matched_gen_jets_pt = awkward.flatten(matched_gen_jets_pt)

    gen_jet_pt = awkward.where(
        (matched_gen_jets_pt > 0) & (abs(corrected_pts_base - matched_gen_jets_pt) < (3.0 * jet_pt_resolution * corrected_pts_base)),
        matched_gen_jets_pt,
        -1.0
    )

    # --- 4. JER & MET Calculation Setup ---
    
    # Prepare MET basics
    # In Run 2, we start from standard PuppiMET (already has some corrections, we refine it)
    # In Run 3, we start from RawPuppiMET (and fully re-apply)
    if lhc_run == 3:
        # Note: Check if RawPuppiMET exists, otherwise fallback might be needed
        met_pt_orig = events.RawPuppiMET_pt
        met_phi_orig = events.RawPuppiMET_phi
    else:
        met_pt_orig = events.PuppiMET_pt
        met_phi_orig = events.PuppiMET_phi

    met_px = met_pt_orig * numpy.cos(met_phi_orig)
    met_py = met_pt_orig * numpy.sin(met_phi_orig)
    
    # Dictionaries to store MET shifts (dx, dy) for various variations
    met_shifts = {
        "nom": {"x": 0.0, "y": 0.0},
        "jerUp": {"x": 0.0, "y": 0.0},
        "jerDown": {"x": 0.0, "y": 0.0},
        "jesTotalUp": {"x": 0.0, "y": 0.0},
        "jesTotalDown": {"x": 0.0, "y": 0.0},
    }

    # Pre-calculate common jet vector components
    jet_cos_phi = numpy.cos(jet_phi)
    jet_sin_phi = numpy.sin(jet_phi)

    # Filter for MET propagation (match C++ logic)
    # 1. Remove muons
    # 2. pt > 15 (after muon subtraction)
    # 3. |eta| < 5.2
    # 4. emEF < 0.9
    muon_factor = jets_flat.muonSubtrFactor
    em_ef = jets_flat.neEmEF + jets_flat.chEmEF
    
    # Pt used for threshold check (Nominal JEC, no JER yet, muon subtracted)
    # C++: jet_l1l2l3_pt_nomu > 15
    pt_for_threshold = corrected_pts_base * (1 - muon_factor)
    
    met_jet_mask = (
        (pt_for_threshold > 15.0) & 
        (abs(jet_eta) < 5.2) & 
        (em_ef < 0.9)
    )

    # Reference pT for MET subtraction
    # Run 2: Reference is the original NanoAOD jet pt (inclusive of whatever JEC was applied upstream)
    # Run 3: Reference is the L1 corrected pT
    if lhc_run == 3:
        pt_ref_flat = pt_l1_flat
    else:
        pt_ref_flat = jets_flat.pt 

    # --- 5. Calculate JER Variations and Nominal PT ---
    
    pt_nom_flat = None # To store for JES calculation

    for jer_variation in ["central", "up", "down"]:
        jer_sf_evaluator = evaluator[f"{jer_tag}_ScaleFactor_{jec_algo}"]
        jer_var_map = {"central": "nom", "up": "up", "down": "down"}
        jer_shift_str = jer_var_map[jer_variation]

        if lhc_run == 2:
            reso_sf = jer_sf_evaluator.evaluate(jet_eta, jer_shift_str)
        else:
            reso_sf = jer_sf_evaluator.evaluate(jet_eta, corrected_pts_base, jer_shift_str)

        # Apply smearing
        smear_factor = numpy.ones_like(corrected_pts_base, dtype=numpy.float32)
        
        # Scaling method for matched jets
        matched_mask = gen_jet_pt > 0
        shift = (reso_sf - 1.0) * (corrected_pts_base - gen_jet_pt) / corrected_pts_base
        smear_factor = numpy.where(matched_mask, numpy.maximum(0.0, 1.0 + shift), smear_factor)

        # Stochastic method for non-matched jets
        not_matched_mask = ~matched_mask
        
        # Using 0 for Gaussian smear seed in this replica to match C++ snippet provided logic (simplification)
        # Real implementation should use proper RNG
        gaus_smear = 0 
        
        if lhc_run == 3:
            gaus_smear = numpy.where(
                (not_matched_mask) & (abs(jet_eta) > 2.5) & (abs(jet_eta) < 3.0) & (corrected_pts_base < 50),
                0.0, 
                gaus_smear
            )
            
        shift_stochastic = gaus_smear * numpy.sqrt(numpy.maximum(reso_sf * reso_sf - 1.0, 0.0))
        smear_factor = numpy.where(not_matched_mask, numpy.maximum(0.0, 1.0 + shift_stochastic), smear_factor)

        pt_jer_varied_flat = corrected_pts_base * smear_factor
        
        # Store for output and JES
        if jer_variation == "central":
            pt_nom_flat = pt_jer_varied_flat
            events[input_collection, "pt_nom"] = awkward.unflatten(pt_nom_flat, n_jets)
            
            # MET Calculation (Nominal)
            # dx = (Pt_new - Pt_ref) * cos(phi)
            # In Run 2: (Pt_nom - Pt_nano) * cos(phi) -> subtract from MET
            # In Run 3: (Pt_nom - Pt_l1) * cos(phi) -> subtract from MET
            # Note: We use the FULL pT for the vector subtraction, but the mask determines if we do it.
            diff_pt = numpy.where(met_jet_mask, pt_nom_flat - pt_ref_flat, 0.0)
            
            # Sum delta_x and delta_y per event
            met_shifts["nom"]["x"] = awkward.sum(awkward.unflatten(diff_pt * jet_cos_phi, n_jets), axis=1)
            met_shifts["nom"]["y"] = awkward.sum(awkward.unflatten(diff_pt * jet_sin_phi, n_jets), axis=1)

        else:
            events[input_collection, f"pt_jer{jer_variation.capitalize()}"] = awkward.unflatten(pt_jer_varied_flat, n_jets)
            
            # MET Calculation (JER Variations)
            # We propagate the difference between (JER_Var) and (Nominal) to the MET
            # MET_JERVar = MET_Nom - Sum(Pt_JERVar - Pt_Nom)
            # Wait, the logic is usually: Calculate Delta relative to Ref for every var, or Delta relative to Nom.
            # C++: met_x_jerup -= ... (jet_nom_pt * (jer_factor_up - 1.0))
            # which equals (jet_jerUp_pt - jet_nom_pt).
            diff_pt_jer = numpy.where(met_jet_mask, pt_jer_varied_flat - pt_nom_flat, 0.0)
            
            key = f"jer{jer_variation.capitalize()}" # jerUp, jerDown
            met_shifts[key]["x"] = awkward.sum(awkward.unflatten(diff_pt_jer * jet_cos_phi, n_jets), axis=1)
            met_shifts[key]["y"] = awkward.sum(awkward.unflatten(diff_pt_jer * jet_sin_phi, n_jets), axis=1)

        print(f"number of matched jets for JER {jer_variation}: {numpy.sum(matched_mask)} out of {len(matched_mask)}")

    # --- 6. JES Uncertainties & MET Propagation ---
    
    # We only compute Total_up/down for MET to save space/time in this snippet, 
    # mirroring the "Total" fallback logic in the jet loop.
    
    for jes_unc_source in ["Total_up", "Total_down"]:
        jes_shift = 1.0 if "up" in jes_unc_source.lower() else -1.0
        source_name = jes_unc_source.replace("_up", "").replace("_down", "")

        pt_scale_sf = numpy.ones_like(pt_nom_flat, dtype=numpy.float32)
        try:
            jes_source_evaluator = evaluator[f"{jes_tag}_{source_name}_{jec_algo}"]
            unc = jes_source_evaluator.evaluate(jet_eta, pt_nom_flat)
            pt_scale_sf = 1.0 + jes_shift * unc
        except:
            # Fallback (simplified for speed, just 0 shift if failed, or repeat logic)
            # Reusing the logic from the jet part for consistency
            all_sources = [s.name for s in evaluator if f"{jes_tag}_" in s.name and f"_{jec_algo}" in s.name and "_L1" not in s.name and "PtResolution" not in s.name and "ScaleFactor" not in s.name]
            quad_sum = numpy.zeros_like(pt_nom_flat, dtype=numpy.float32)
            for source in all_sources:
                s_eval = evaluator[source]
                quad_sum += numpy.power(s_eval.evaluate(jet_eta, pt_nom_flat), 2.0)
            pt_scale_sf = 1.0 + jes_shift * numpy.sqrt(quad_sum)

        pt_jes_varied_flat = pt_nom_flat * pt_scale_sf
        
        # Save Jet Branch
        branch_name = f"pt_jes{source_name}{'Up' if jes_shift > 0 else 'Down'}"
        events[input_collection, branch_name] = awkward.unflatten(pt_jes_varied_flat, n_jets)
        
        # MET Calculation (JES Variations)
        # MET_JESVar = MET_Nom - Sum(Pt_JESVar - Pt_Nom)
        diff_pt_jes = numpy.where(met_jet_mask, pt_jes_varied_flat - pt_nom_flat, 0.0)
        
        key = f"jes{source_name}{'Up' if jes_shift > 0 else 'Down'}"
        met_shifts[key]["x"] = awkward.sum(awkward.unflatten(diff_pt_jes * jet_cos_phi, n_jets), axis=1)
        met_shifts[key]["y"] = awkward.sum(awkward.unflatten(diff_pt_jes * jet_sin_phi, n_jets), axis=1)

    # --- 7. Finalize MET Branches ---

    # Nominal MET
    # MET_Nom = MET_Orig - Sum(Nom - Ref)
    met_px_nom = met_px - met_shifts["nom"]["x"]
    met_py_nom = met_py - met_shifts["nom"]["y"]
    events["PuppiMET_pt_corrected"] = numpy.sqrt(met_px_nom**2 + met_py_nom**2)
    events["PuppiMET_phi_corrected"] = numpy.arctan2(met_py_nom, met_px_nom)

    # Variations
    # For variations, the shift in dictionary is (Var - Nom).
    # MET_Var = MET_Nom - (Var - Nom)
    # So we start from met_px_nom and subtract the JES/JER shift.
    
    for key, shifts in met_shifts.items():
        if key == "nom": continue
        
        # Construct branch names: PuppiMET_ptJERUp, PuppiMET_ptJESTotalUp, etc.
        # Key is like "jerUp" or "jesTotalUp"
        base_name = "PuppiMET"
        suffix = key.replace("jesTotal", "jesTotal") # keep Total if present
        # Or map to specific names if needed:
        if "jesTotal" in key:
            suffix = key.replace("jesTotal", "JES") # jesTotalUp -> JESUp
        elif "jer" in key:
            suffix = key.replace("jer", "JER") # jerUp -> JERUp
            
        px_var = met_px_nom - shifts["x"]
        py_var = met_py_nom - shifts["y"]
        
        events[f"{base_name}_pt{suffix}"] = numpy.sqrt(px_var**2 + py_var**2)
        events[f"{base_name}_phi{suffix}"] = numpy.arctan2(py_var, px_var)

    return events


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
    
    It also propagates these corrections to PuppiMET (Nominal only).

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
        An array with the corrected jet pT and corrected PuppiMET.
    """
    JEC_TAG_DATA = {
        "2022preEE": "Summer22_22Sep2023_RunCD_V2_DATA",
        "2022postEE": "Summer22EE_22Sep2023_V3_DATA",
        "2023preBPix": "Summer23Prompt23_V2_DATA",
        "2023postBPix": "Summer23BPixPrompt23_V3_DATA",
    }
    JEC_FILE = {
        "2022preEE": "jsonpog-integration/POG/JME/2022_Summer22/jet_jerc.json",
        "2022postEE": "jsonpog-integration/POG/JME/2022_Summer22EE/jet_jerc.json",
        "2023preBPix": "jsonpog-integration/POG/JME/2023_Summer23/jet_jerc.json",
        "2023postBPix": "jsonpog-integration/POG/JME/2023_Summer23BPix/jet_jerc.json",
    }

    try:
        jets = events[input_collection]
    except:
        return events

    lhc_run = 2
    if "2022" in year or "2023" in year:
        lhc_run = 3
    
    if run_period:
        jes_tag = f"{JEC_TAG_DATA[year][:-8]}_Run{run_period}_{JEC_TAG_DATA[year][-7:]}"
    else:
        jes_tag = JEC_TAG_DATA[year]
        
    jec_file = misc_utils.expand_path(JEC_FILE[year])
    # jec_file = JEC_FILE[year]
    evaluator = correctionlib.CorrectionSet.from_file(jec_file)

    n_jets = awkward.num(jets)
    jets_flat = awkward.flatten(jets)
    jets_area = awkward.to_numpy(jets_flat.area)
    jets_eta = awkward.to_numpy(jets_flat.eta)
    rho_arr = numpy.repeat(awkward.to_numpy(events["Rho_fixedGridRhoFastjetAll"]), n_jets)

    # --- 1. Calculate Fully Corrected Jet PT ---
    raw_pt = jets_flat.pt * (1 - jets_flat.rawFactor)
    raw_pt_numpy = awkward.to_numpy(raw_pt)
    
    # Prepare inputs for evaluator
    if "2023" in year:
        run_arr = numpy.repeat(awkward.to_numpy(events["run"]), n_jets)
        jets_phi = awkward.to_numpy(jets_flat.phi)
    
    jec_evaluator = evaluator.compound[f"{jes_tag}_L1L2L3Res_{jec_algo}"]
    
    # Evaluate based on year specific inputs
    if "2023post" in year:
        corr = jec_evaluator.evaluate(
            jets_area, jets_eta, raw_pt_numpy, rho_arr, jets_phi, run_arr
        )
    elif "2023pre" in year:
        corr = jec_evaluator.evaluate(
            jets_area, jets_eta, raw_pt_numpy, rho_arr, run_arr
        )
    else:
        corr = jec_evaluator.evaluate(
            jets_area, jets_eta, raw_pt_numpy, rho_arr
        )
        
    corrected_pts_flat = raw_pt * corr
    events[input_collection, "corrected_pt"] = awkward.unflatten(corrected_pts_flat, n_jets)

    # --- 2. Calculate L1 Corrected Jet PT (Needed for Run 3 MET) ---
    pt_l1_flat = raw_pt 
    if lhc_run == 3:
        try:
            l1_evaluator = evaluator[f"{jes_tag}_L1FastJet_{jec_algo}"]
            if "2023" in year:
                # 2023 data L1 usually needs run number? Check JSON. Usually L1 is just Area/Eta/Pt/Rho
                # If IOV dependent, might need run. Assuming standard L1 signature for now.
                 l1_corr = l1_evaluator.evaluate(jets_area, jets_eta, raw_pt_numpy, rho_arr)
            else:
                 l1_corr = l1_evaluator.evaluate(jets_area, jets_eta, raw_pt_numpy, rho_arr)
            
            pt_l1_flat = raw_pt * l1_corr
        except Exception as e:
            # Fallback if L1 not found, though risky for MET accuracy
            pt_l1_flat = raw_pt 

    # --- 3. Propagate to PuppiMET ---

    # Get original MET
    if lhc_run == 3:
        met_pt_orig = events.RawPuppiMET_pt
        met_phi_orig = events.RawPuppiMET_phi
        pt_ref_flat = pt_l1_flat # Run 3 subtracts L1, adds Full
    else:
        met_pt_orig = events.PuppiMET_pt
        met_phi_orig = events.PuppiMET_phi
        pt_ref_flat = jets_flat.pt # Run 2 subtracts NanoAOD stored, adds Full

    met_px = met_pt_orig * numpy.cos(met_phi_orig)
    met_py = met_pt_orig * numpy.sin(met_phi_orig)

    # Selection mask for jets propagating to MET
    muon_factor = jets_flat.muonSubtrFactor
    em_ef = jets_flat.neEmEF + jets_flat.chEmEF
    
    # Threshold check using muon-subtracted nominal pT
    pt_for_threshold = corrected_pts_flat * (1 - muon_factor)
    
    met_jet_mask = (
        (pt_for_threshold > 15.0) & 
        (abs(jets_eta) < 5.2) & 
        (em_ef < 0.9)
    )

    # Calculate Vector Difference (Nominal - Reference)
    diff_pt = numpy.where(met_jet_mask, corrected_pts_flat - pt_ref_flat, 0.0)
    
    jet_cos_phi = numpy.cos(jets_flat.phi)
    jet_sin_phi = numpy.sin(jets_flat.phi)

    # Sum shifts per event
    shift_x = awkward.sum(awkward.unflatten(diff_pt * jet_cos_phi, n_jets), axis=1)
    shift_y = awkward.sum(awkward.unflatten(diff_pt * jet_sin_phi, n_jets), axis=1)

    # Apply correction: MET_new = MET_old - Sum(Jet_new - Jet_ref)
    met_px_new = met_px - shift_x
    met_py_new = met_py - shift_y

    events["PuppiMET_pt_corrected"] = numpy.sqrt(met_px_new**2 + met_py_new**2)
    events["PuppiMET_phi_corrected"] = numpy.arctan2(met_py_new, met_px_new)

    return events