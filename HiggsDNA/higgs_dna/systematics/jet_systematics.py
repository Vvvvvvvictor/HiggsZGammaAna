import awkward

import numpy

from correctionlib import _core

from higgs_dna.utils import misc_utils, awkward_utils

###################################
### b-tag continuous reshape SF ###
###################################

BTAG_RESHAPE_SF_FILE = {
    "2016" : "jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json", 
    "2016UL_preVFP" : "jsonpog-integration/POG/BTV/2016preVFP_UL/btagging.json", 
    "2016UL_postVFP" : "jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json", 
    "2017" : "jsonpog-integration/POG/BTV/2017_UL/btagging.json",
    "2018" : "jsonpog-integration/POG/BTV/2018_UL/btagging.json"
}


DEEPJET_VARIATIONS = { 
    "up_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm (4)
    "up_lf" : [5],
    "up_hfstats1" : [5],
    "up_hfstats2" : [5],
    "up_cferr1" : [4],
    "up_cferr2" : [4],
    "up_hf" : [0],
    "up_lfstats1" : [0],
    "up_lfstats2" : [0],
    "down_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm(4)
    "down_lf" : [5],
    "down_hfstats1" : [5],
    "down_hfstats2" : [5],
    "down_cferr1" : [4],
    "down_cferr2" : [4],
    "down_hf" : [0],
    "down_lfstats1" : [0],
    "down_lfstats2" : [0],
}

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
    print(jet_flavor[jet_flavor==0])
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
    central_sf = evaluator["deepJet_mujets"].evalv(
            "central",
            working_point,
            jet_flavor,
            jet_abs_eta,
            jet_pt
    )    
    variations["central"] = awkward.unflatten(central_sf, n_jets)
    if not central_only:
        syst_vars = ["up", "down"] 
        for syst_var in syst_vars:
            syst = evaluator["deepJet_mujets"].evalv(
                    syst_var,
                    working_point,
                    jet_flavor,
                    jet_abs_eta,
                    jet_pt,
            )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
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
    jet_flavor = 0 * jet_flavor # set all jets to light flavor
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.00001,
        2.499 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.9
    )
    variations_list = ["central"]

    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}
    central_sf = evaluator["deepJet_incl"].evalv(
            "central",
            "M",
            jet_flavor,
            jet_abs_eta,
            jet_pt
    )    
    variations["central"] = awkward.unflatten(central_sf, n_jets)
    if not central_only:
        syst_vars = ["up", "down"] 
        for syst_var in syst_vars:
            syst = evaluator["deepJet_incl"].evalv(
                    syst_var,
                    working_point,
                    jet_flavor,
                    jet_abs_eta,
                    jet_pt,
            )
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_jets)
    return variations

def dummy_jes_syst(events, is_data):
    """
    Dummy function illustrating a jet energy scale uncertainty that results in new jet collections with Jet.pt varied.
    Should be deleted once real examples are implemented.
    """
    jets = events.Jet 

    variations = {}
    variations["central"] = jets.pt + (2 * awkward.ones_like(jets.pt))
    if not is_data:
        variations["up"] = jets.pt + (12 * awkward.ones_like(jets.pt))
        variations["down"] = jets.pt - (8 * awkward.ones_like(jets.pt))
    return variations
