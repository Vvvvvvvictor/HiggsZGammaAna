import awkward
import time
import numpy
import numba
import vector

from pdb import set_trace

vector.register_awkward()

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG 
from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.selections import object_selections, lepton_selections, jet_selections, tau_selections, physics_utils
from higgs_dna.selections import gen_selections

DUMMY_VALUE = -999.
DEFAULT_OPTIONS = {
    "photons" : {
        "use_central_nano" : True,
        "pt" : 15.0,
        "eta" : [
            [0.0, 1.4442],
            [1.566, 2.5]
        ],
        "mvaID_barrel" : -0.4,
        "mvaID_endcap" : -0.58,
        "e_veto" : 0.5
    },
    "zgammas" : {
        "relative_pt_gamma" : 15.0/110.,
        "mass_h" : [100., 180.],
        "mass_sum" : 185,
        "select_highest_pt_sum" : True
    },
    "single_muon_trigger":{
        "2016":["HLT_IsoMu24", "HLT_IsoTkMu24"],
        "2017":["HLT_IsoMu27"],
        "2018":["HLT_IsoMu24"],
        "2022":["HLT_IsoMu24"],
        "2023":["HLT_IsoMu24"]
    },
    "single_ele_trigger":{
        "2016":["HLT_Ele27_WPTight_Gsf"],
        "2017":["HLT_Ele32_WPTight_Gsf_L1DoubleEG"],
        "2018":["HLT_Ele32_WPTight_Gsf"],
        "2022":["HLT_Ele30_WPTight_Gsf"],
        "2023":["HLT_Ele30_WPTight_Gsf"]
    },
    "double_muon_trigger":{
        "2016":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"],
        "2017":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2018":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2022":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2023":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"]
    },
    "double_ele_trigger":{
        "2016":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
        "2017":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2018":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2022":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2023":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
    },
    # "trigger" : {
    #     "2016" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele27_WPTight_Gsf"],
    #     "2016UL_preVFP" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele27_WPTight_Gsf"],
    #     "2016UL_postVFP" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele27_WPTight_Gsf"],
    #     "2017" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Ele27_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_IsoMu24_eta2p1"],
    #     "2018" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Mu37_TkMu27", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_IsoMu24_eta2p1", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_DoubleEle25_CaloIdL_MW", "HLT_Ele27_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_Ele20_WPLoose_Gsf"],
    #     "2022" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Ele30_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"]
    # }, 
    "electrons" : {
        "pt" : 7.0
    },
    "muons" : {
        "pt" : 5.0
    },
    "lead_ele_pt":{
        "2016": 30,
        "2017": 35,
        "2018": 35,
        "2022": 35,
        "2023": 35
    },
    "lead_mu_pt":{
        "2016": 25,
        "2017": 28,
        "2018": 25,
        "2022": 25,
        "2023": 25
    },
    "jets" : {
        "pt" : 30.0,
        "eta" : 4.7,
        "eta" : 4.7,
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
        "jets_horn" : {
            "pt" : 50.0,
            "eta" : [2.5, 3.0]
        }
    },
    "btag_med": {
        "2016preVFP": 0.2598,
        "2016postVFP": 0.2489,
        "2017": 0.3040,
        "2018": 0.2783,
        "2022preEE": 0.3086,
        "2022postEE": 0.3196,
        "2023preBPix": 0.2431,
        "2023postBPix": 0.2435
    },
    "FSRphotons" : {
        "iso" : 1.8,
        "eta" : 2.4,
        "pt" : 2.0,
        "dROverEt2" : 0.012
    },
    "gen_info" : {
        "calculate" : False,
        "max_dr" : 0.2,
        "max_pt_diff" : 15.
    }  
}


# Diphoton preselection below synced with flashgg, see details in:
#   - https://indico.cern.ch/event/1071721/contributions/4551056/attachments/2320292/3950844/HiggsDNA_DiphotonPreselectionAndSystematics_30Sep2021.pdf

class ZGammaTaggerRun2(Tagger):
    def __init__(self, name = "default_zgamma_tagger", options = {}, is_data = None, year = None):
        super(ZGammaTaggerRun2, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, events):
        """
        Select photons and create diphoton pairs.
        Add a record "Diphoton" to events array with relevant information about each diphoton pair.
        In principle, there can be more than one Diphoton pair per event.
        """

        # Determine what type of rho variable is available in nanoAOD
        # To be deleted once a standard rho variable is added to central nanoAOD
        if not self.options["photons"]["use_central_nano"]:
            if "fixedGridRhoAll" in events.fields:
                rho = events.fixedGridRhoAll
            elif "fixedGridRhoFastjetAll" in events.fields:
                rho = events.fixedGridRhoFastjetAll
            elif "Rho_fixedGridRhoAll" in events.fields:
                rho = events.Rho_fixedGridRhoAll
        else:
            rho = awkward.ones_like(events.Photon)

        if not self.is_data:
            self.overlap_removal(events=events)

        zgamma_selection, zgammas = self.produce_and_select_zgammas(
                events = events,
                rho = rho,
                options = self.options["zgammas"]
        )

        if not self.is_data and self.options["gen_info"]["calculate"]:
            zgammas = self.calculate_gen_info(zgammas, self.options["gen_info"])

        return zgamma_selection, zgammas 

    def overlap_removal(self, events):
        """
        Select isolation photons in events
        Add number of isolation photons (n_iso_photons) in output .parquet file to indetify thé overlap events
        """
        
        """
        statusFlags usage: (events.GenPart.statusFlags // numpy.power(2, i)) % 2 == 1
        "statusFlags" is a number with 14 bits. 
        Filling "1" on corresponding digit when the particle meets one of the 14 conditions, else, remaining "0".
        Echo particles can meet more than one kind of condition, thus, more than one digit in "statusFlags" is "1".
        """
        iso_photons_cut = (events.GenPart.pdgId == 22) & (events.GenPart.pt > 15) & (abs(events.GenPart.eta) < 2.6) & (( (events.GenPart.statusFlags // numpy.power(2, 0)) % 2 == 1 ) | ( (events.GenPart.statusFlags // numpy.power(2, 8)) % 2 == 1 ))
        iso_photons = events.GenPart[iso_photons_cut]

        truth_objects_cut =  (events.GenPart.pdgId != 22) & (events.GenPart.pt > 5) & ( (events.GenPart.statusFlags // numpy.power(2, 8)) % 2 == 1 ) 
        truth_objects = events.GenPart[truth_objects_cut]

        iso_cut = object_selections.delta_R(iso_photons, truth_objects, 0.05)
        iso_photons = iso_photons[iso_cut]

        n_iso_photons = awkward.num(iso_photons)
        awkward_utils.add_field(events, "n_iso_photons", n_iso_photons, overwrite=True)

    def produce_and_select_zgammas(self, events, rho, options):
        """
        Perform diphoton preselection.
        For events with more than 2 photons, more than 1 diphoton candidate
        per event is possible.

        :param events: events array to calculate diphoton candidates from
        :type events: awkward.highlevel.Array
        :param photons: array of selected photons from events
        :type photons: awkward.highlevel.Array
        :param options: dictionary containing configurable options for the diphoton preselection
        :type options: dict
        :return: boolean array indicating which events pass the diphoton preselection
        :rtype: awkward.highlevel.Array
        """

        start = time.time()

        # events = events[(events.run == 356077) & (events.luminosityBlock == 158) & (events.event == 196518696)]
        # Electrons
        electron_cut = lepton_selections.select_electrons(
            electrons = events.Electron,
            options = self.options["electrons"],
            clean = {
            },
            name = "SelectedElectron",
            tagger = self,
            year = self.year[:4]
        )
        
        electrons = awkward_utils.add_field(
            events = events,
            name = "SelectedElectron",
            data = events.Electron[electron_cut]
        )
        
        # generate the index in the original array and add to electrons
        arr = awkward.local_index(events.Electron["pt"], axis=1)[electron_cut]
        electron_idx = awkward.mask(arr, awkward.num(arr) > 0)
        awkward_utils.add_field(events = electrons, name = "Idx", data = electron_idx)

        # Muons
        muon_cut = lepton_selections.select_muons(
            muons = events.Muon,
            options = self.options["muons"],
            clean = {
            },
            name = "SelectedMuon",
            tagger = self
        )

        muons = awkward_utils.add_field(
            events = events,
            name = "SelectedMuon",
            data = events.Muon[muon_cut]
        )

        # gen mu not reco
        gen_muons = events.GenPart[(abs(events.GenPart.pdgId) == 13)]
        reco_muons = events.Muon

        unmatched_gen_mask = awkward.fill_none(object_selections.delta_R(gen_muons, reco_muons, 0.4), True)
        
        # reco_gen_muon_indices = reco_muons.genPartIdx
        # all_gen_muon_indices = awkward.local_index(gen_muons.pt)

        # # For each event, check which of the gen muon indices are not in the list of reco-matched gen muon indices
        # unmatched_gen_mask = awkward.all(all_gen_muon_indices[:, :, None] != reco_gen_muon_indices[:, None, :], axis=-1)
        
        unmatched_gen_muons = gen_muons[unmatched_gen_mask]
        unmatched_gen_muons = unmatched_gen_muons[unmatched_gen_muons.pt > self.options["muons"]["pt"]]

        unmatched_gen_muons = awkward_utils.add_field(
            events = events,
            name = "SelectedGenNoRecoMuon",
            data = unmatched_gen_muons
        )

        # gen ele not reco
        gen_eles = events.GenPart[(abs(events.GenPart.pdgId) == 11)]
        reco_eles = events.Electron

        unmatched_gen_mask = awkward.fill_none(object_selections.delta_R(gen_eles, reco_eles, 0.4), True)
        
        # reco_gen_ele_indices = reco_eles.genPartIdx
        # all_gen_ele_indices = awkward.local_index(gen_eles.pt)

        # # For each event, check which of the gen electron indices are not in the list of reco-matched gen electron indices
        # unmatched_gen_mask = awkward.all(all_gen_ele_indices[:, :, None] != reco_gen_ele_indices[:, None, :], axis=-1)

        unmatched_gen_eles = gen_eles[unmatched_gen_mask]
        unmatched_gen_eles = unmatched_gen_eles[unmatched_gen_eles.pt > self.options["electrons"]["pt"]]

        unmatched_gen_eles = awkward_utils.add_field(
            events = events,
            name = "SelectedGenNoRecoElectron",
            data = unmatched_gen_eles
        )

        # Photons
        photon_selection = self.select_photons(
                photons = events.Photon,
                electrons = electrons,
                rho = rho,
                options = self.options["photons"]
        )

        photons = events.Photon[photon_selection]
        
        # lepton-photon overlap removal 
        clean_photon_mask = awkward.fill_none(object_selections.delta_R(photons, muons, 0.3), True) & awkward.fill_none(object_selections.delta_R(photons, electrons, 0.3), True) # FIXME: 0.4 -> 0.3(baseline)
        # object_selections.delta_R(photons, muons, 0.3) & object_selections.delta_R(photons, electrons, 0.3)
        photons = photons[clean_photon_mask]
        
        # Jets
        jet_cut = jet_selections.select_jets(
            jets = events.Jet,
            options = self.options["jets"],
            clean = {
                "photons" : {
                    "objects" : photons,
                    "min_dr" : self.options["jets"]["dr_photons"]
                },
                "electrons" : {
                    
                    "objects" : electrons,
                    "min_dr" : self.options["jets"]["dr_electrons"]
                },
                "muons" : {
                    "objects" : muons,
                    "min_dr" : self.options["jets"]["dr_muons"]
                }
            },
            year = self.year,
            name = "SelectedJet",
            tagger = self
        )

        jets = awkward_utils.add_field(
            events = events,
            name = "SelectedJet",
            data = events.Jet[jet_cut]
        )
        
        photon = awkward_utils.add_field(
                events = events,
                name = "Photon",
                data = events.Photon[photon_selection],
        )

        # if "2017" not in self.year and "2018" not in self.year:
        #     b_jet_cut = jets.btagDeepFlavB > self.options["btag_med"][self.year[:4]]
        # else:
        b_jet_cut = jets.btagDeepFlavB > self.options["btag_med"][self.year]
        jets = awkward.with_field(jets, b_jet_cut, "is_med_bjet") 

        # Add object fields to events array
        for objects, name in zip([electrons, muons, jets, unmatched_gen_muons, unmatched_gen_eles], ["electron", "muon", "jet", "gen_no_reco_muon", "gen_no_reco_electron"]):
            for var in objects.fields:
                if var in ["charge", "pt", "eta", "phi", "mass", "id", "ptE_error"]:
                    awkward_utils.add_field(
                        events,
                        f"{name}_{var}",
                        awkward.fill_none(objects[var], DUMMY_VALUE)
                    )
                else:
                    awkward_utils.add_field(
                        events,
                        f"{name}_{var}",
                        awkward.fill_none(objects[var], 0)
                    )

        if not self.is_data:
            dZ = events.GenVtx_z - events.PV_z
            awkward_utils.add_field(events, "dZ", dZ, overwrite=True)

        n_electrons = awkward.fill_none(awkward.num(electrons), 0)
        # N_e_cut = n_electrons>=2
        awkward_utils.add_field(events, "n_electrons", n_electrons, overwrite=True)

        n_muons = awkward.num(muons)
        # N_mu_cut = n_muons>=2
        awkward_utils.add_field(events, "n_muons", n_muons, overwrite=True)

        n_leptons = n_electrons + n_muons
        # N_e_mu_cut = N_e_cut | N_mu_cut
        awkward_utils.add_field(events, "n_leptons", n_leptons, overwrite=True)

        n_jets = awkward.num(jets)
        logger.debug(f"Number of jets(tagger): {n_jets[:10]}")
        awkward_utils.add_field(events, "n_jets", n_jets, overwrite=True)

        n_b_jets = awkward.sum(b_jet_cut, axis=1)
        awkward_utils.add_field(events, "n_b_jets", n_b_jets, overwrite=True)

        # PDG ID
        electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 11, "id")
        # electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 0.00051099895, "mass")
        muons = awkward.with_field(muons, awkward.ones_like(muons.pt) * 13, "id")

        # leptons ptE_error
        electrons = awkward.with_field(electrons, electrons.energyErr, "ptE_error")
        muons = awkward.with_field(muons, muons.ptErr, "ptE_error")

        # Sort objects by pt
        photons = photons[awkward.argsort(photons.pt, ascending=False, axis=1)]
        electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
        muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]

        # self.select_fake_and_medium_photons(events=events, photons=photons)

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        photons = awkward.Array(photons, with_name = "Momentum4D")
        electrons = awkward.Array(electrons, with_name = "Momentum4D")
        muons = awkward.Array(muons, with_name = "Momentum4D")


        ee_pairs = awkward.combinations(electrons, 2, fields = ["LeadLepton", "SubleadLepton"])
        os_cut = (ee_pairs.LeadLepton.charge * ee_pairs.SubleadLepton.charge) == -1
        ee_pairs = ee_pairs[os_cut]
        mm_pairs = awkward.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
        os_cut = (mm_pairs.LeadLepton.charge * mm_pairs.SubleadLepton.charge) == -1
        mm_pairs = mm_pairs[os_cut]
        z_cands = awkward.concatenate([ee_pairs, mm_pairs], axis = 1)
        z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton
        z_cands = z_cands[awkward.argsort(abs(z_cands.ZCand.mass - 91.1876), axis = 1)]
        z_ee_cut = awkward.fill_none(awkward.firsts(z_cands).LeadLepton.id == 11, False)
        z_mumu_cut = awkward.fill_none(awkward.firsts(z_cands).LeadLepton.id == 13, False)

        # Make trigger cuts 
        if self.year is not None:
            year = self.year[:4]
            single_ele_trigger_cut = awkward.num(events.Photon) < 0 # dummy cut, all False
            double_ele_trigger_cut = awkward.num(events.Photon) < 0
            single_mu_trigger_cut = awkward.num(events.Photon) < 0
            double_mu_trigger_cut = awkward.num(events.Photon) < 0 
            for hlt in self.options["single_ele_trigger"][year]:
                if hasattr(events, hlt):
                    single_ele_trigger_cut = (single_ele_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
            for hlt in self.options["double_ele_trigger"][year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    double_ele_trigger_cut = (double_ele_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
            for hlt in self.options["single_muon_trigger"][year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    single_mu_trigger_cut = (single_mu_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
            for hlt in self.options["double_muon_trigger"][year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    double_mu_trigger_cut = (double_mu_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
        else:
            single_ele_trigger_cut = awkward.num(events.Photon) >= 0 # dummy cut, all True
            double_ele_trigger_cut = awkward.num(events.Photon) >= 0
            single_mu_trigger_cut = awkward.num(events.Photon) >= 0
            double_mu_trigger_cut = awkward.num(events.Photon) >= 0
            
        if "2017" in self.year:
            single_ele_trigger_cut = awkward.any((events.TrigObj.id == 11) & ((events.TrigObj.filterBits & 0x400) != 0), axis=1) & single_ele_trigger_cut

        trigger_cut = single_ele_trigger_cut | double_ele_trigger_cut | single_mu_trigger_cut  | double_mu_trigger_cut
        ele_trigger_cut = single_ele_trigger_cut | double_ele_trigger_cut
        mu_trigger_cut = single_mu_trigger_cut  | double_mu_trigger_cut
        # HLT lepton status: [FIXME]
        # 0: failed all triggers ->  double_ele_trigger_cut = False, single_ele_trigger_cut = False
        # 1: passed lower dilepton trigger -> double_ele_trigger_cut = False, single_ele_trigger_cut = True
        # 2: passed upper dilepton trigger -> double_ele_trigger_cut = True, single_ele_trigger_cut = False
        # 3: passed single lepton trigger ->  single_ele_trigger_cut = True
        # HLT_ele_cat0 = (~double_ele_trigger_cut) & (~single_ele_trigger_cut)
        # HLT_mu_cat0 = (~double_mu_trigger_cut) & (~single_mu_trigger_cut)
        # HLT_ele_cat1 = (~double_ele_trigger_cut) & single_ele_trigger_cut
        # HLT_mu_cat1 = (~double_mu_trigger_cut) & single_mu_trigger_cut
        # HLT_ele_cat2 = double_ele_trigger_cut & (~single_ele_trigger_cut)
        # HLT_mu_cat2 = double_mu_trigger_cut & (~single_mu_trigger_cut)
        # HLT_ele_cat3 = single_ele_trigger_cut
        # HLT_mu_cat3 = single_mu_trigger_cut
        # def GetLeptonProbability(lepton_pt, lepton_eta, is_data, is_electron, trigger_leg) :
        #     if (is_electron):
        #         if (trigger_leg == pass_lowerdilep):
        #             if is_data:
        #                 prob_map = diele12_correction_data
        #                 unc_map = diele12_uncertainty_data
        #             else :
        #                 prob_map = diele12_correction_mc
        #                 unc_map = diele12_uncertainty_mc
        #         elif (trigger_leg == pass_upperdilep):
        #             if is_data:
        #                 prob_map = diele23_correction_data
        #                 unc_map = diele23_uncertainty_data
        #             else :
        #                 prob_map = diele23_correction_mc
        #                 unc_map = diele23_uncertainty_mc
        #         else :  
        #             if is_data:
        #                 prob_map = singleele_correction_data
        #                 unc_map = singleele_uncertainty_data
        #             else :
        #                 prob_map = singleele_correction_mc
        #                 unc_map = singleele_uncertainty_mc
        #     else:
        #         if (trigger_leg == pass_lowerdilep):
        #             if is_data:
        #                 prob_map = dimu12_correction_data
        #                 unc_map = dimu12_uncertainty_data
        #             else :
        #                 prob_map = dimu12_correction_mc
        #                 unc_map = dimu12_uncertainty_mc
        #         elif (trigger_leg == pass_upperdilep):
        #             if is_data:
        #                 prob_map = dimu23_correction_data
        #                 unc_map = dimu23_uncertainty_data
        #             else :
        #                 prob_map = dimu23_correction_mc
        #                 unc_map = dimu23_uncertainty_mc
        #         else :  
        #             if is_data:
        #                 prob_map = singlemu_correction_data
        #                 unc_map = singlemu_uncertainty_data
        #             else :
        #                 prob_map = singlemu_correction_mc
        #                 unc_map = singlemu_uncertainty_mc
        #     prob = prob_map
        #     uncr = unc_map
        #     return prob, uncr

        # def GetFlavorProbability(lepton_pt, lepton_eta, pass_singlelep, pass_dilep, is_data, is_electron):
        #     n_pass_lower = awkward.num(1*HLT_ele_cat1[1*HLT_ele_cat1>0])+awkward.num(1*HLT_mu_cat1[1*HLT_mu_cat1>0])
        #     n_pass_upper = awkward.num(1*HLT_ele_cat2[1*HLT_ele_cat2>0])+awkward.num(1*HLT_mu_cat2[1*HLT_mu_cat2>0])
        #     n_pass_single = awkward.num(1*HLT_ele_cat3[1*HLT_ele_cat3>0])+awkward.num(1*HLT_mu_cat3[1*HLT_mu_cat3>0])
        #     relevant_cat = awkward.ones_like(lepton_pt) * True
        #     relevant_cat = awkward.where(((n_pass_single>0)!=pass_singlelep), awkward.ones_like(relevant_cat) * False, relevant_cat)
        #     relevant_cat = awkward.where((((n_pass_upper>0)&(n_pass_lower>1))!=pass_dilep), awkward.ones_like(relevant_cat) * False, relevant_cat)
        #     lep_prob = awkward.zero_like(lepton_pt) 
        #     lep_unc = awkward.zero_like(lepton_pt)

        # def GetTotalProbability(electron_pt,muon_pt, electron_eta,muon_eta, pass_singleel,pass_singlemu, pass_diel, pass_dimu, is_data):
        #     electron_prob = GetFlavorProbability(electron_pt, electron_eta,pass_singleel, pass_diel, is_data, True)
        #     muon_prob = GetFlavorProbability(muon_pt, muon_eta,pass_singlemu, pass_dimu, is_data, False)
        
        # mc_prob = GetTotalProbability(electron_pt, muon_pt, electron_eta, muon_eta, pass_singleel, pass_singlemu, pass_diel, pass_dimu, False)


        # 0: lep_prob = 1 - probability(1), lep_unc = uncertainty(1)
        # 1: lep_prob = prob(1) - prob(2), lep_unc = sqrt(uncertainty(1)**2 + uncertainty(2)**2)
        # 2: lep_prob = prob(2) - prob(3), lep_unc = sqrt(uncertainty(2)**2 + uncertainty(3)**2)
        # 3: lep_prob = prob(3), lep_unc = uncertainty(3)
        
        #if mc_prob < 0.001: mc_prob = 0.
        # sf = 1.0;
        # unc = 0.0;
        # propagate_uncertainty_ratio(data_prob, data_unc, mc_prob, mc_unc, sf, unc);
        # sf,unc=1
        # if mc_prob !=0: sf = data_prob/mc_prob; unc = sqrt((data_unc/mc_prob)**2 + ((mc_unc*data_prob)/(mc_prob*mc*prob))**2)
        # elif mc_prob == 0 sf = 1; unc = 0
        #[FIXME]
        # a trick that give 0 to the empty array    
        if self.year is not None:
            year = self.year[:4]
            e_cut = awkward.fill_none(awkward.pad_none(electrons.pt, 1, axis=1)[:, 0], 0) > self.options["lead_ele_pt"][year]
            m_cut = awkward.fill_none(awkward.pad_none(muons.pt, 1, axis=1)[:, 0], 0) > self.options["lead_mu_pt"][year]
        else:
            e_cut = awkward.fill_none(awkward.pad_none(electrons.pt, 1, axis=1)[:, 0], 0) > 25
            m_cut = awkward.fill_none(awkward.pad_none(muons.pt, 1, axis=1)[:, 0], 0) > 20
        ee_cut = (awkward.fill_none(awkward.pad_none(electrons.pt, 1, axis=1)[:, 0], 0) > 25) & (awkward.fill_none(awkward.pad_none(electrons.pt, 2, axis=1)[:, 1], 0) > 15)
        mm_cut = (awkward.fill_none(awkward.pad_none(muons.pt, 1, axis=1)[:, 0], 0) > 20) & (awkward.fill_none(awkward.pad_none(muons.pt, 2, axis=1)[:, 1], 0) > 10)
        
        ele_trigger_pt_cut = (single_ele_trigger_cut & e_cut) | (double_ele_trigger_cut & ee_cut)
        mu_trigger_pt_cut = (single_mu_trigger_cut & m_cut) | (double_mu_trigger_cut & mm_cut)
        trigger_pt_cut = (single_ele_trigger_cut & e_cut) | (double_ele_trigger_cut & ee_cut) | (single_mu_trigger_cut & m_cut) | (double_mu_trigger_cut & mm_cut)
        
        z_mumu = z_mumu_cut 
        z_ee = z_ee_cut
        awkward_utils.add_field(events, "z_mumu", z_mumu, overwrite=True)
        awkward_utils.add_field(events, "z_ee", z_ee, overwrite=True)

        logger.info(f"Z candidates lead lepton pt: {awkward.flatten(z_cands.LeadLepton.pt[:50])}"
        f" Z candidates sublead lepton pt: {awkward.flatten(z_cands.SubleadLepton.pt[:50])}")
        logger.info(f"Z candidates mass: {awkward.flatten(z_cands.ZCand.mass[:50])}")
        mass_cut = (z_cands.ZCand.mass > 80.) & (z_cands.ZCand.mass < 100.)
        # mass_cut = z_cands.ZCand.mass > 50.
        z_cands = z_cands[mass_cut] # OSSF lepton pairs with m_ll > 50.
        
        # HEM cut
        if self.year=="2018" and self.is_data:
            hem_run=events.run > 319077        
            # checked 65.15623538907509% events in data could pass this run cut
            hem_jet=awkward.num(events.Jet[(events.Jet.phi>-1.57) & (events.Jet.phi<-0.87) & (events.Jet.eta>-3) & (events.Jet.eta<-1.3)])>0
            hem_fatjet=awkward.num(events.FatJet[(events.FatJet.phi>-1.57) & (events.FatJet.phi<-0.87) & (events.FatJet.eta>-3) & (events.FatJet.eta<-1.3)])>0
            hem_cut=~((hem_run & hem_jet) | (hem_run & hem_fatjet))        
        elif self.year=="2018" and not self.is_data:
            #random number generator from 0 to 1
            fraction=0.6515623538907509
            events['random'] = numpy.random.rand(len(events))
            hem_run=events.random < fraction
            hem_jet=awkward.num(events.Jet[(events.Jet.phi>-1.57) & (events.Jet.phi<-0.87) & (events.Jet.eta>-3) & (events.Jet.eta<-1.3)])>0
            hem_fatjet=awkward.num(events.FatJet[(events.FatJet.phi>-1.57) & (events.FatJet.phi<-0.87) & (events.FatJet.eta>-3) & (events.FatJet.eta<-1.3)])>0
            hem_cut=~((hem_run & hem_jet) | (hem_run & hem_fatjet))
        else:
            hem_cut=awkward.num(events.Photon) >= 0 
        events = events[hem_cut]

        has_z_cand = awkward.num(z_cands) >= 1
        z_cand = awkward.firsts(z_cands)

        # Add Z-related fields to array
        for field in ["pt", "eta", "phi", "mass", "charge", "id", "ptE_error"]:
            if not field in ["charge", "id", "ptE_error"]:
                awkward_utils.add_field(
                        events,
                        "Z_%s" % field,
                        awkward.fill_none(getattr(z_cand.ZCand, field), DUMMY_VALUE)
                )
            awkward_utils.add_field(
                    events,
                    "Z_lead_lepton_%s" % field,
                    awkward.fill_none(z_cand.LeadLepton[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                    events,
                    "Z_sublead_lepton_%s" % field,
                    awkward.fill_none(z_cand.SubleadLepton[field], DUMMY_VALUE)
            )

        FSRphoton_selection = self.select_FSRphotons(
                FSRphotons = events.FsrPhoton,
                electrons = electrons,
                photons = photons,
                z_cand = z_cand,
                options = self.options["FSRphotons"]
        )
        FSRphotons = awkward_utils.add_field(
            events = events,
            name = "SelectedFSRPhotons",
            data = events.FsrPhoton[FSRphoton_selection]
        )
        FSRphotons = awkward.with_field(FSRphotons, awkward.ones_like(FSRphotons.pt) * 0.0, "mass")        
        FSRphotons = awkward.Array(FSRphotons, with_name = "Momentum4D")
        awkward_utils.add_object_fields(
                events = events,
                name = "gamma_fsr",
                objects = FSRphotons,
                n_objects = 1,
                dummy_value = DUMMY_VALUE
            )
        awkward_utils.add_field(events, "n_fsr", awkward.num(FSRphotons), overwrite=True)
        logger.debug(f"Number of FSR photons: {awkward.num(FSRphotons)}")
        logger.debug(f"Total number of FSR photons: {sum(awkward.num(FSRphotons)>0)}")

        # Make gamma candidate-level cuts
        has_gamma_cand = (awkward.num(photons) >= 1) #& (events.n_iso_photons == 0) # only for dy samples
        gamma_cand = awkward.firsts(photons)
        gamma_mvaID_WPL = ((gamma_cand.isScEtaEB & (gamma_cand.mvaID > self.options["photons"]["mvaID_barrel"])) | (gamma_cand.isScEtaEE & (gamma_cand.mvaID > self.options["photons"]["mvaID_endcap"])))
        gamma_e_veto = gamma_cand.electronVeto > self.options["photons"]["e_veto"]

        awkward_utils.add_field(gamma_cand, "mass", awkward.ones_like(gamma_cand.pt) * 0) #TODO: run3 BUG

        # Add gamma-related fields to array
        for field in ["pt", "eta", "phi", "mass", "mvaID", "energyErr", "sieie", "hoe", "r9", "mvaID_WP80", "mvaID_WP90"]:
            awkward_utils.add_field(
                events,
                "gamma_%s" % field,
                awkward.fill_none(getattr(gamma_cand, field), DUMMY_VALUE)
            )
        awkward_utils.add_field(events, "gamma_mvaID_WPL",  gamma_mvaID_WPL)
        awkward_utils.add_field(events, "gamma_e_veto",  gamma_e_veto)
        if int(self.year[:4]) < 2020:
            awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg) #run2
            awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all) #run2
        elif int(self.year[:4]) > 2020:
            awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg_quadratic) #run3
            awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all_quadratic) #run3
        #awkward_utils.add_field(events, "gamma_mvaID_17",  gamma_cand.mvaID_Fall17V2) #run3

        # Make Higgs candidate-level cuts
        h_cand = (z_cand.ZCand + gamma_cand)
        sel_h_1 = (gamma_cand.pt / h_cand.mass) > options["relative_pt_gamma"]
        sel_h_2 = (z_cand.ZCand.mass + h_cand.mass) > options["mass_sum"]
        sel_h_3 = (h_cand.mass > options["mass_h"][0]) & (h_cand.mass < options["mass_h"][1])

        sel_h_1 = awkward.fill_none(sel_h_1, value = False)
        sel_h_2 = awkward.fill_none(sel_h_2, value = False)
        sel_h_3 = awkward.fill_none(sel_h_3, value = False)
        
        print(f'!!!!has H: {sum(has_z_cand & has_gamma_cand & z_mumu_cut)} | {sum(has_z_cand & has_gamma_cand & z_ee_cut)}')

        # Add Higgs-related fields to array
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(
                events,
                "H_%s" % field,
                awkward.fill_none(getattr(h_cand, field), DUMMY_VALUE)
            )

        # additional leptons
        leptons = awkward.concatenate([electrons, muons], axis = 1)
        max_I_mini = awkward.fill_none(awkward.max(leptons.miniPFRelIso_all, axis = 1), 9999)
        awkward_utils.add_field(events, "max_I_mini", max_I_mini)
        
        veto_Z_leptons = (leptons.pt != events.Z_lead_lepton_pt) & (leptons.pt != events.Z_sublead_lepton_pt)
        additional_leptons = leptons[veto_Z_leptons]       
        additional_leptons = additional_leptons[awkward.argsort(additional_leptons.pt, ascending=False, axis=1)]

        for objects, name in zip([additional_leptons], ["additional_lepton"]):
            awkward_utils.add_object_fields(
                events = events,
                name = name,
                objects = objects,
                n_objects = 2,
                dummy_value = DUMMY_VALUE
            )

        event_filter = (events.Flag_goodVertices & 
                        events.Flag_globalSuperTightHalo2016Filter & 
                        ((awkward.num(events.Photon) >= 0) if "202" in self.year else events.Flag_HBHENoiseFilter) & 
                        ((awkward.num(events.Photon) >= 0) if "202" in self.year else events.Flag_HBHENoiseIsoFilter) & 
                        events.Flag_EcalDeadCellTriggerPrimitiveFilter & 
                        events.Flag_BadPFMuonFilter & 
                        events.Flag_BadPFMuonDzFilter & 
                        events.Flag_hfNoisyHitsFilter & 
                        events.Flag_eeBadScFilter & 
                        ((awkward.num(events.Photon) >= 0) if "2016" in self.year else events.Flag_ecalBadCalibFilter) # 2016 dummy cut, all True
                        )
        
        all_cuts = trigger_pt_cut & has_z_cand & has_gamma_cand & sel_h_1 & sel_h_2 & event_filter & sel_h_3 #& awkward.fill_none((h_cand.mass>80) & (h_cand.mass < options["mass_h"][1]), False)

        for cut_type in ["zgammas", "zgammas_ele", "zgammas_mu", "zgammas_w", "zgammas_ele_w", "zgammas_mu_w"]:
            if "_w" in cut_type:
                if hasattr(events, 'Generator_weight'):
                    weighted = True
                else:
                    continue
            else:
                weighted = False

            cut0 = awkward.num(events.Photon) >= 0

            cut1 = z_ee_cut | z_mumu_cut
            if "ele" in cut_type:
                cut1 = z_ee_cut
            elif "mu" in cut_type:
                cut1 = z_mumu_cut
            cut2 = cut1 & trigger_cut
            if "ele" in cut_type:
                cut2 = cut1 & ele_trigger_cut
            elif "mu" in cut_type:
                cut2 = cut1 & mu_trigger_cut
            cut3 = cut2 & trigger_pt_cut
            if "ele" in cut_type:
                cut3 = cut2 & ele_trigger_pt_cut
            elif "mu" in cut_type:
                cut3 = cut2 & mu_trigger_pt_cut
            cut4 = cut3 & has_gamma_cand
            cut5 = cut4 & has_z_cand
            cut6 = cut5 & sel_h_1
            cut7 = cut6 & sel_h_2
            cut8 = cut7 & sel_h_3
            cut9 = cut8 & event_filter
            
            if cut_type == "zgammas_ele":
                ee_all_cut = cut9
            if cut_type == "zgammas_mu":
                mm_all_cut = cut9
            
            # if cut_type == "zgammas_ele":
            #     print(f"!!!start check events tag({cut_type})!!!")
            #     for i in events[cut1]:
            #         print(f"{i.run} {i.luminosityBlock} {i.event}")
            #     print(f"!!!end check events tag({cut_type})!!!")

            self.register_event_cuts(
                # names = ["all", "N_lep_sel", "trig_cut", "lead_lep_pt_cut", "sub_lep_pt_cut", "has_g_cand", "has_z_cand", "sel_h_1", "sel_h_2", "sel_h_3"],
                # results = [cut0, cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9],
                names = ["all", "N_lep_sel", "trig_cut", "lep_pt_cut", "has_g_cand", "has_z_cand", "sel_h_1", "sel_h_2", "sel_h_3", "event", "all cuts"],
                results = [cut0, cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9, all_cuts],
                events = events,
                cut_type = cut_type,
                weighted = weighted
            )
            # cut_names = ["N_lep_sel", "trig_cut", "lep_pt_cut", "has_g_cand", "os_cut", "has_z_cand", "sel_h_1", "sel_h_2", "sel_h_3"]
            # for cut, cut_name in zip([cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9], cut_names):
            #     awkward_utils.add_field(events, f"{cut_type}_{cut_name}", cut)

        all_cuts = ee_all_cut | mm_all_cut
        
        # print(f"Sum of all_cuts: {sum(all_cuts)}")
        # all_cuts = ee_all_cut | mm_all_cut
        # print(f"Sum of all_cuts: {sum(all_cuts)}")

        # checked_cut = (z_ee_cut | z_mumu_cut) & pair_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(inclusive)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(inclusive)!!!")

        # checked_cut = z_ee_cut & ele_trigger_pt_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(electron)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(electron)!!!")

        # checked_cut = z_mumu_cut & mm_trigger_pt_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(muon)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(muon)!!!")

        elapsed_time = time.time() - start
        logger.debug("[ZGammaTagger] %s, syst variation : %s, total time to execute select_zgammas: %.6f s" % (self.name, self.current_syst, elapsed_time))

        #dummy_cut =  awkward.num(events.Photon) >= 0
        return all_cuts, events 

    
    def calculate_gen_info(self, zgammas, options):
        """
        Calculate gen info, adding the following fields to the events array:
            GenHggHiggs : [pt, eta, phi, mass, dR]
            GenHggLeadPhoton : [pt, eta, phi, mass, dR, pt_diff]
            GenHggSubleadPhoton : [pt, eta, phi, mass, dR, pt_diff]
            LeadPhoton : [gen_dR, gen_pt_diff]
            SubleadPhoton : [gen_dR, gen_pt_diff]

        Perform both matching of
            - closest gen photons from Higgs to reco lead/sublead photons from diphoton candidate
            - closest reco photons to gen photons from Higgs

        If no match is found for a given reco/gen photon, it will be given values of -999. 
        """
        gen_hzg = gen_selections.select_x_to_yz(zgammas.GenPart, 25, 23, 22)

        print("DEBUG: bing", gen_hzg.fields)
        
        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgHiggs",
                objects = gen_hzg.GenParent,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgLeadGenChild",
                objects = gen_hzg.LeadGenChild,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgSubleadGenChild",
                objects = gen_hzg.SubleadGenChild,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgLeadGenChildChild1",
                objects = gen_hzg.LeadGenChildChild1,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgLeadGenChildChild2",
                objects = gen_hzg.LeadGenChildChild2,
                n_objects = 1
        )

        return zgammas 
        

    def select_photons(self, photons, electrons, rho, options):
        """
        Enforces all photon cuts that are commmon to both
        leading and subleading photons in the diphoton preselection.
        Cuts specific to a diphoton pair are not enforced here.

        :param photons: input collection of photons to use for calculating selected photons
        :type photons: awkward.highlevel.Array
        :param rho: energy density in each event, used for corrections to photon isolation cuts
        :type rho: awkward.highlevel.Array
        :param options: dictionary containing configurable options for the photon selection
        :type options: dict
        :return: boolean array indicating which photons pass the photon selection
        :rtype: awkward.highlevel.Array
        """
        # pt
        pt_cut = photons.pt > options["pt"]

        # eta
        #eta_cut = Tagger.get_range_cut(abs(photons.eta), options["eta"]) | (photons.isScEtaEB | photons.isScEtaEE)
        eta_cut = (photons.isScEtaEB | photons.isScEtaEE) 
        id_cut = photons.mvaID_WP80
        # eta_cut = ((photons.isScEtaEB & (photons.mvaID > options["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > options["mvaID_endcap"])))

        # electron veto
        e_veto_cut = (photons.electronVeto > options["e_veto"])
        
        photon_ele_idx = awkward.where(awkward.num(photons.electronIdx, axis=1) == 0, awkward.ones_like(photons.pt)*-1, photons.electronIdx)
        new_pho = awkward.unflatten(awkward.unflatten(awkward.flatten(photon_ele_idx), [1]*awkward.sum(awkward.num(photon_ele_idx))), awkward.num(photon_ele_idx, axis=1))
        new_ele = awkward.broadcast_arrays(electrons.Idx[:,None], new_pho, depth_limit=2)[0]
        eg_overlap_cut = ~awkward.where(
            awkward.is_none(electrons.Idx),
            awkward.broadcast_arrays(photons.electronIdx, False)[1],
            awkward.flatten(awkward.any(new_pho[:, :, None] == new_ele, axis=-2), axis=-1)
        ) # some events may have no electrons, so we need to replace None with False

        # use_central_nano = options["use_central_nano"] # indicates whether we are using central nanoAOD (with some branches that are necessary for full diphoton preselection missing) or custom nanoAOD (with these branches added)

        all_cuts = pt_cut & eta_cut & id_cut & e_veto_cut & eg_overlap_cut
        # all_cuts = pt_cut & eta_cut & e_veto_cut & eg_overlap_cut # bing for CR selection

        self.register_cuts(
                names = ["pt", "eta", "id", "e_veto", "ele_pho_overlap", "all"], #"pt", "eta", "id", "e_veto", "ele_pho_overlap", "all"
                results = [pt_cut, eta_cut, id_cut, e_veto_cut, eg_overlap_cut, all_cuts], #pt_cut, eta_cut, id_cut, e_veto_cut, eg_overlap_cut, all_cuts
                cut_type = "photon"
        )

        return all_cuts

    def select_FSRphotons(self, FSRphotons, electrons, photons, z_cand, options):
        # Basic kinematic cuts: pt > 2 GeV and |eta| < 2.4
        FSR_pt_cut = FSRphotons.pt > options["pt"]
        FSR_eta_cut = abs(FSRphotons.eta) < options["eta"]
        FSR_iso_cut = FSRphotons.relIso03 < options["iso"]
        FSR_dROverEt2_cut = FSRphotons.dROverEt2 < options["dROverEt2"]
        
        # Clean from electrons and photons - remove FSR photons too close to electrons or photons
        FSRphoton_clean = object_selections.delta_R(FSRphotons, electrons, 0.001) & object_selections.delta_R(FSRphotons, photons, 0.001)
        
        # Require ΔR > 0.2 from highest pt photon
        highest_pt_photon = photons[awkward.argmax(photons.pt, axis=1, keepdims=True)]
        FSR_photon_separation = object_selections.delta_R(FSRphotons, highest_pt_photon, 0.2)
        
        # Apply lepton-type specific selections
        is_muon_channel = awkward.fill_none(z_cand.LeadLepton.id == 13, False)
        # is_electron_channel = awkward.fill_none(abs(z_cand.LeadLepton.id) == 11, False)
        
        FSR_all_cuts = FSR_pt_cut & FSR_eta_cut & FSR_iso_cut & FSR_dROverEt2_cut & FSRphoton_clean & FSR_photon_separation & is_muon_channel
        
        # For multiple FSR photons, select the one with smallest ΔR/ET²
        # This follows the criteria: for events with multiple FSR photons passing criteria 
        # for the same lepton, choose FSR photon with smallest ΔR(lepton,γ)/E²ₜ,γ
        best_fsr_idx = awkward.argmin(
            awkward.where(FSR_all_cuts, FSRphotons.dROverEt2, float('inf')), 
            axis=1, 
            keepdims=True
        )
        
        # Create final selection mask: only one FSR photon per event
        event_range = awkward.local_index(FSRphotons.pt, axis=1)
        best_fsr_mask = (event_range == best_fsr_idx)
        all_cut = best_fsr_mask & FSR_all_cuts

        self.register_cuts(
            names = ["FSR_pt", "FSR_eta", "FSR_iso", "FSR_dROverEt2", "FSR_clean", "FSR_separation", "FSR_all_cuts", "best_fsr_mask", "all"],
            results = [FSR_pt_cut, FSR_eta_cut, FSR_iso_cut, FSR_dROverEt2_cut, FSRphoton_clean, FSR_photon_separation, FSR_all_cuts, best_fsr_mask, all_cut],
            cut_type = "FSRphotons"
        )

        FSR_all_cuts = awkward.fill_none(FSR_all_cuts, False) & awkward.fill_none(best_fsr_mask, False)
        
        return FSR_all_cuts

# Below is an example of how the diphoton preselection could be performed with an explicit loop (C++ style) 
# that is compiled with numba for increased performance.

# NOTE: pre-compiled numba functions should be defined outside the class,
# as numba does not like when a class instance is passed to a function
@numba.njit
def produce_diphotons(photons, n_photons, lead_pt_cut, lead_pt_mgg_cut, sublead_pt_mgg_cut):
    n_events = len(photons)

    diphotons_offsets = numpy.zeros(n_events + 1, numpy.int64)
    diphotons_contents = []
    lead_photons_idx = []
    sublead_photons_idx = []

    # Loop through events and select diphoton pairs
    for i in range(n_events):
        n_diphotons_event = 0
        # Enumerate photon pairs
        if n_photons[i] >= 2: # only try if there are at least 2 photons
            sum_pt = 0
            for j in range(n_photons[i]):
                for k in range(j+1, n_photons[i]):
                    # Choose who is the leading photon
                    lead_idx = j if photons[i][j].pt > photons[i][k].pt else k
                    sublead_idx = k if photons[i][j].pt > photons[i][k].pt else j
                    lead_photon = photons[i][lead_idx]
                    sublead_photon = photons[i][sublead_idx]

                    # Lead photon must satisfy lead pt cut
                    if lead_photon.pt < lead_pt_cut:
                        continue

                    # Construct four vectors
                    lead_photon_vec = vector.obj(
                            pt = lead_photon.pt,
                            eta = lead_photon.eta,
                            phi = lead_photon.phi,
                            mass = lead_photon.mass
                    )
                    sublead_photon_vec = vector.obj(
                            pt = sublead_photon.pt,
                            eta = sublead_photon.eta,
                            phi = sublead_photon.phi,
                            mass = sublead_photon.mass
                    )
                    diphoton = vector.obj(px = 0., py = 0., pz = 0., E = 0.) # IMPORTANT NOTE: you need to initialize this to an empty vector first. Otherwise, you will get ZeroDivisionError exceptions for like 1 out of a million events (seemingly only with numba). 
                    diphoton = diphoton + lead_photon_vec + sublead_photon_vec

                    if (diphoton.mass < 100) | (diphoton.mass > 180):
                        continue

                    if lead_photon.pt / diphoton.mass < lead_pt_mgg_cut:
                        continue

                    if sublead_photon.pt / diphoton.mass < sublead_pt_mgg_cut:
                        continue

                    # This diphoton candidate passes
                    n_diphotons_event += 1

                    diphotons_contents.append([
                        diphoton.pt,
                        diphoton.eta,
                        diphoton.phi,
                        diphoton.mass
                    ])

                    lead_photons_idx.append(lead_idx)
                    sublead_photons_idx.append(sublead_idx)

        diphotons_offsets[i+1] = diphotons_offsets[i] + n_diphotons_event

    return diphotons_offsets, numpy.array(diphotons_contents), numpy.array(lead_photons_idx), numpy.array(sublead_photons_idx)
