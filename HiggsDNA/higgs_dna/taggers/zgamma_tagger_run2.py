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
        "2022": 30,
        "2023": 30
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
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
    },
    "btag_med": {
        "2016preVFP": 0.2598,
        "2016postVFP": 0.2489,
        "2017": 0.3040,
        "2018": 0.2783,
        "2022": 0.3086,
        "2022EE": 0.3196,
        "2023": 0.2431,
        "2023BPix": 0.2435
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
        Filling "1" on corresponding digit when the particle meets one of the 14 conditions, else remaining "0".
        Echo paticles can meet more than one kind of condition, thus, more than one digit in "statusFlags" is "1".
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

        # Electrons
        electron_cut = lepton_selections.select_electrons(
            electrons = events.Electron,
            options = self.options["electrons"],
            clean = {
            },
            name = "SelectedElectron",
            tagger = self
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

        # Photons
        photon_selection, sr_cut = self.select_photons(
                photons = events.Photon,
                electrons = electrons,
                rho = rho,
                options = self.options["photons"]
        )

        photons = events.Photon[photon_selection]
        
        # lepton-photon overlap removal 
        clean_photon_mask = object_selections.delta_R(photons, muons, 0.3) & object_selections.delta_R(photons, electrons, 0.3)
        photons = photons[clean_photon_mask]
<<<<<<< HEAD


=======
        
>>>>>>> main
        # Jets
        jet_cut = jet_selections.select_jets(
            jets = events.Jet,
            options = self.options["jets"],
            clean = {
                "photons" : {
                    "objects" : photons[(photons.electronVeto > self.options["photons"]["e_veto"]) & (photons.mvaID_WP80 == 1)],
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
            name = "SelectedJet",
            tagger = self
        )

        jets = awkward_utils.add_field(
            events = events,
            name = "SelectedJet",
            data = events.Jet[jet_cut]
        )

        FSRphoton_selection = self.select_FSRphotons(
                FSRphotons = events.FsrPhoton,
                electrons = electrons,
                photons = photons,
                options = self.options["FSRphotons"]
        )
        FSRphotons = awkward_utils.add_field(
            events = events,
            name = "SelectedFSRPhotons",
            data = events.FsrPhoton[FSRphoton_selection]
        )
        FSRphotons = awkward.with_field(FSRphotons, awkward.ones_like(FSRphotons.pt) * 0.0, "mass")

        if "2016" not in self.year:
            year = self.year[:4]
            b_jet_cut = jets.btagDeepFlavB > self.options["btag_med"][year]
        else:
            b_jet_cut = jets.btagDeepFlavB > self.options["btag_med"][self.year]
        jets = awkward.with_field(jets, b_jet_cut, "is_med_bjet") 

        # Add object fields to events array
        for objects, name in zip([electrons, muons, jets], ["electron", "muon", "jet"]):
            awkward_utils.add_object_fields(
                events = events,
                name = name,
                objects = objects,
                n_objects = 4,
                dummy_value = DUMMY_VALUE
            )

        n_electrons = awkward.num(electrons)
        # N_e_cut = n_electrons>=2
        awkward_utils.add_field(events, "n_electrons", n_electrons, overwrite=True)

        n_muons = awkward.num(muons)
        # N_mu_cut = n_muons>=2
        awkward_utils.add_field(events, "n_muons", n_muons, overwrite=True)

        n_leptons = n_electrons + n_muons
        # N_e_mu_cut = N_e_cut | N_mu_cut
        awkward_utils.add_field(events, "n_leptons", n_leptons, overwrite=True)

        n_jets = awkward.num(jets)
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
        FSRphotons = awkward.Array(FSRphotons, with_name = "Momentum4D")

        #awkward_utils.add_field(events, "gamma_fsr_pt",  awkward.fill_none(FSRphotons.pt, DUMMY_VALUE))

        # bing with control regions
        #var_CR1 = ((photons.isScEtaEB & (photons.mvaID > self.options["photons"]["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > self.options["photons"]["mvaID_endcap"])))
        var_CR1 = awkward.fill_none(photons.mvaID_WP80, False)
        var_CR2 = awkward.fill_none(photons.mvaID_WP90, False)
        var_CR3 = awkward.fill_none(photons.electronVeto, False)
        SR_cut = (var_CR1) & (var_CR3)
        CR1_cut = (~var_CR1) & (var_CR2) & (var_CR3)
        CR2_cut = (~var_CR2) & (var_CR3)
        CR3_cut = (var_CR1) & (~var_CR3)
        CR4_cut = (~var_CR1) & (var_CR2) & (~var_CR3)
        CR5_cut = (~var_CR2) & (~var_CR3)
        
        SR = awkward.fill_none(awkward.sum(SR_cut, axis=1) > 0, value = False)
        CR1 = awkward.fill_none(awkward.sum(CR1_cut, axis=1) > 0, value=False)
        CR2 = awkward.fill_none(awkward.sum(CR2_cut, axis=1) > 0, value=False)
        CR3 = awkward.fill_none(awkward.sum(CR3_cut, axis=1) > 0, value=False)
        CR4 = awkward.fill_none(awkward.sum(CR4_cut, axis=1) > 0, value=False)
        CR5 = awkward.fill_none(awkward.sum(CR5_cut, axis=1) > 0, value=False)
        regions = awkward.where(SR, 0, awkward.where(CR1, 1, awkward.where(CR2, 2, awkward.where(CR3, 3, awkward.where(CR4, 4, awkward.where(CR5, 5, -1))))))

        masks, photon_masks, non_overlap_masks, non_overlap_photon_masks = [SR, CR1, CR2, CR3, CR4, CR5], [SR_cut, CR1_cut, CR2_cut, CR3_cut, CR4_cut, CR5_cut], [], []
        recent_mask = awkward.zeros_like(SR)
        for mask, photon_mask in zip(masks, photon_masks):
            non_overlap_mask = mask & ~recent_mask
            non_overlap_masks.append(non_overlap_mask)
            non_overlap_photon_mask = non_overlap_mask & photon_mask
            non_overlap_photon_masks.append(non_overlap_photon_mask)
            recent_mask = mask | recent_mask
        photons = awkward.concatenate([photons[non_overlap_photon_mask] for non_overlap_photon_mask in non_overlap_photon_masks], axis = 1)
        
            
        #for i in range(5000):
        #   print("SR:", SR[i], "CR1:", CR1[i], "CR2", CR2[i], "CR3", CR3[i], "CR4", CR4[i], "CR5", CR5[i], "var_CR1", var_CR1[i], "var_CR2", var_CR2[i], "var_CR3", var_CR3[i], "regions:", regions[i])
        awkward_utils.add_field(events, "regions",  regions)

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
        ele_trigger_cut = single_ele_trigger_cut #| double_ele_trigger_cut
        mu_trigger_cut = single_mu_trigger_cut  #| double_mu_trigger_cut
        
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
        
        ele_trigger_pt_cut = (single_ele_trigger_cut & e_cut) #| (double_ele_trigger_cut & ee_cut)
        mu_trigger_pt_cut = (single_mu_trigger_cut & m_cut) #| (double_mu_trigger_cut & mm_cut)
        trigger_pt_cut = (single_ele_trigger_cut & e_cut) | (double_ele_trigger_cut & ee_cut) | (single_mu_trigger_cut & m_cut) | (double_mu_trigger_cut & mm_cut)

        mass_cut = (z_cands.ZCand.mass > 80.) & (z_cands.ZCand.mass < 100.)
        # mass_cut = z_cands.ZCand.mass > 50.
        z_cands = z_cands[mass_cut] # OSSF lepton pairs with m_ll > 50.

        # # Construct di-electron/di-muon pairs
        # ee_pairs = awkward.combinations(electrons, 2, fields = ["LeadLepton", "SubleadLepton"])
        # if self.year is not None:
        #     year = self.year[:4]
        #     e_cut = ee_pairs.LeadLepton.pt > self.options["lead_ele_pt"][year]
        # else:
        #     e_cut = ee_pairs.LeadLepton.pt > 25
        # ee_cut = (ee_pairs.LeadLepton.pt > 25) & (ee_pairs.SubleadLepton.pt > 15)
        # ee_pairs = awkward.concatenate([ee_pairs[double_ele_trigger_cut & ee_cut], ee_pairs[single_ele_trigger_cut & e_cut]], axis = 1)
        # ele_pair_cut = awkward.num(ee_pairs) >= 1

        # mm_pairs = awkward.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
        # if self.year is not None:
        #     year = self.year[:4]
        #     m_cut = mm_pairs.LeadLepton.pt > self.options["lead_mu_pt"][year]
        # else:
        #     m_cut = mm_pairs.LeadLepton.pt > 20
        # mm_cut = (mm_pairs.LeadLepton.pt > 20) & (mm_pairs.SubleadLepton.pt > 10)
        # mm_pairs = awkward.concatenate([mm_pairs[double_mu_trigger_cut & mm_cut], mm_pairs[single_mu_trigger_cut & m_cut]], axis = 1)
        # muon_pair_cut = awkward.num(mm_pairs) >= 1

        # # Concatenate these together
        # z_cands = awkward.concatenate([ee_pairs, mm_pairs], axis = 1)

        # # Make Z candidate-level cuts
        # # print('!!!!charge of lead lepton{}, charge of sublead lepton {}, flavour of lepton {}'.format(awkward.flatten(z_cands[os_cut]).LeadLepton.charge, awkward.flatten(z_cands[os_cut]).SubleadLepton.charge, awkward.flatten(z_cands[os_cut]).LeadLepton.id))
        # os_cut = (z_cands.LeadLepton.charge * z_cands.SubleadLepton.charge) == -1
        # z_cands = z_cands[os_cut]
        # os_cut = awkward.num(z_cands) >= 1

        # z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton # these add as 4-vectors since we registered them as "Momentum4D" objects
        # mass_cut = (z_cands.ZCand.mass > 80.) & (z_cands.ZCand.mass < 100.)
        # # mass_cut = (z_cands.ZCand.mass > 50.)
        # z_cands = z_cands[mass_cut] # OSSF lepton pairs with m_ll > 50.

        # z_cands = z_cands[awkward.argsort(abs(z_cands.ZCand.mass - 91.1876), axis = 1)] # take the one with mass closest to mZ

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

        # Make gamma candidate-level cuts
        has_gamma_cand = (awkward.num(photons) >= 1)#& (events.n_iso_photons == 0) # only for dy samples
        gamma_cand = awkward.firsts(photons)
        gamma_mvaID_WPL = ((gamma_cand.isScEtaEB & (gamma_cand.mvaID > self.options["photons"]["mvaID_barrel"])) | (gamma_cand.isScEtaEE & (gamma_cand.mvaID > self.options["photons"]["mvaID_endcap"])))
        gamma_e_veto = gamma_cand.electronVeto > self.options["photons"]["e_veto"]

        awkward_utils.add_field(gamma_cand, "mass", awkward.ones_like(gamma_cand.pt) * 0) #TODO: run3 BUG

        # Add gamma-related fields to array
        for field in ["pt", "eta", "phi", "mass", "mvaID", "energyErr", "sieie", "hoe", "r9", "mvaID_WP80"]:
            awkward_utils.add_field(
                events,
                "gamma_%s" % field,
                awkward.fill_none(getattr(gamma_cand, field), DUMMY_VALUE)
            )
        awkward_utils.add_field(events, "gamma_mvaID_WPL",  gamma_mvaID_WPL)
        awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg) #run2
        awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all) #run2
        awkward_utils.add_field(events, "gamma_e_veto",  gamma_e_veto)
        awkward_utils.add_field(events, "gamma_mvaID_WP80",  gamma_cand.mvaID_WP80)
        awkward_utils.add_field(events, "gamma_mvaID_WP90",  gamma_cand.mvaID_WP90)
        # awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg_quadratic) #run3
        # awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all_quadratic) #run3
        #awkward_utils.add_field(events, "gamma_mvaID_17",  gamma_cand.mvaID_Fall17V2) #run3

        # Make Higgs candidate-level cuts
        h_cand = (z_cand.ZCand + gamma_cand)
        sel_h_1 = (gamma_cand.pt / h_cand.mass) > options["relative_pt_gamma"]
        sel_h_2 = (z_cand.ZCand.mass + h_cand.mass) > options["mass_sum"]
        sel_h_3 = (h_cand.mass > options["mass_h"][0]) & (h_cand.mass < options["mass_h"][1])

        sel_h_1 = awkward.fill_none(sel_h_1, value = False)
        sel_h_2 = awkward.fill_none(sel_h_2, value = False)
        sel_h_3 = awkward.fill_none(sel_h_3, value = False)

        # Add Higgs-related fields to array
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(
                events,
                "H_%s" % field,
                awkward.fill_none(getattr(h_cand, field), DUMMY_VALUE)
            )

        # additional leptons
        additional_leptons = awkward.concatenate([electrons, muons], axis = 1)
        veto_Z_leptons = (additional_leptons.pt != events.Z_lead_lepton_pt) & (additional_leptons.pt != events.Z_sublead_lepton_pt)
        additional_leptons = additional_leptons[veto_Z_leptons]       
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
                        events.Flag_HBHENoiseFilter & 
                        events.Flag_HBHENoiseIsoFilter & 
                        events.Flag_EcalDeadCellTriggerPrimitiveFilter & 
                        events.Flag_BadPFMuonFilter & 
                        events.Flag_BadPFMuonDzFilter & 
                        events.Flag_hfNoisyHitsFilter & 
                        events.Flag_eeBadScFilter & 
                        ((awkward.num(events.Photon) >= 0) if "2016" in self.year else events.Flag_ecalBadCalibFilter) # 2016 dummy cut, all True
                        )
        
        all_cuts = trigger_pt_cut & has_z_cand & has_gamma_cand & sel_h_1 & sel_h_2 & sel_h_3 & event_filter

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

        # all_cuts = cut0

        # checked_cut = (z_ee_cut | z_mumu_cut) & pair_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(inclusive)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(inclusive)!!!")

        # checked_cut = z_ee_cut & ee_trigger_pt_cut
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

        # self.register_cuts(
        #     names = ["has_z_cand", "has_gamma_cand", "sel_h_1", "sel_h_2", "sel_h_3", "all cuts"],
        #     results = [has_z_cand, has_gamma_cand, sel_h_1, sel_h_2, sel_h_3, all_cuts],
        #     cut_type = "zgammas_unweighted"
        # )
        ################## fsr recovery
        #clean_FSRphoton_mask = object_selections.delta_R_fsrlep(FSRphotons, z_cand.LeadLepton, 0.001) & object_selections.delta_R_fsrlep(FSRphotons, z_cand.SubleadLepton, 0.001) & object_selections.delta_R_fsrGamma(FSRphotons, gamma_cand, 0.001)
        #FSRphotons = FSRphotons[clean_FSRphoton_mask]
        
        FSRphotons = awkward.pad_none(FSRphotons,1)
        for field in ["pt", "eta", "phi", "mass", "relIso03", "dROverEt2"]:
            awkward_utils.add_field(
                events = events,
                name = "gamma_fsr_%s" % field,
                data = awkward.fill_none(FSRphotons[field][:,0], DUMMY_VALUE)
            )


        # # bing with control regions
        # #var_CR1 = ((photons.isScEtaEB & (photons.mvaID > self.options["photons"]["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > self.options["photons"]["mvaID_endcap"])))
        # var_CR1 = gamma_cand.mvaID_WP80
        # var_CR2 = gamma_cand.mvaID_WP90
        # var_CR3 = gamma_e_veto
        # SR = (var_CR1) & (var_CR3)
        # CR1 = (~var_CR1) & (var_CR2) & (var_CR3)
        # CR2 = (~var_CR2) & (var_CR3)
        # CR3 = (var_CR1) & (~var_CR3)
        # CR4 = (~var_CR1) & (var_CR2) & (~var_CR3)
        # CR5 = (~var_CR2) & (~var_CR3)
        
        # SR = awkward.fill_none(SR, value = False)
        # CR1 = awkward.fill_none(CR1, value=False)
        # CR2 = awkward.fill_none(CR2, value=False)
        # CR3 = awkward.fill_none(CR3, value=False)
        # CR4 = awkward.fill_none(CR4, value=False)
        # CR5 = awkward.fill_none(CR5, value=False)
        # regions = awkward.where(SR, 0, awkward.where(CR1, 1, awkward.where(CR2, 2, awkward.where(CR3, 3, awkward.where(CR4, 4, awkward.where(CR5, 5, -1))))))
<<<<<<< HEAD
        # # for i in range(5000):
        # #    print("SR:", SR[i], "CR1:", CR1[i], "CR2", CR2[i], "CR3", CR3[i], "CR4", CR4[i], "CR5", CR5[i], "var_CR1", var_CR1[i], "var_CR2", var_CR2[i], "var_CR3", var_CR3[i], "regions:", regions[i])
=======
        # for i in range(5000):
        #    print("SR:", SR[i], "CR1:", CR1[i], "CR2", CR2[i], "CR3", CR3[i], "CR4", CR4[i], "CR5", CR5[i], "var_CR1", var_CR1[i], "var_CR2", var_CR2[i], "var_CR3", var_CR3[i], "regions:", regions[i])
>>>>>>> main
        # awkward_utils.add_field(events, "regions",  regions)


        elapsed_time = time.time() - start
        logger.debug("[ZGammaTagger] %s, syst variation : %s, total time to execute select_zgammas: %.6f s" % (self.name, self.current_syst, elapsed_time))

        #dummy_cut =  awkward.num(events.Photon) >= 0
        return all_cuts, events 


    # def select_fake_and_medium_photons(self, events, photons):
    #     # | pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho | 
    #     # listed from the right side
    #     mask1 = 0b10101010101010  # full medium ID
    #     mask2 = 0b00101010101010  # remove Isopho
    #     mask3 = 0b10001010101010  # remove IsoNeu
    #     mask4 = 0b10100010101010  # remove Isoch 
    #     mask5 = 0b10101000101010  # remove sigma ieie
    #     mask6 = 0b10100000101010  # remove the Isoch and sigma ieie

    #     # photons = photons[photons.pixelSeed]
    #     bitmap = photons.vidNestedWPBitmap

    #     # select medium and control photons
    #     # after adding the photons that pass the full ID, add the photons that pass the inverted ID
    #     # select control photons that don't pass the full ID but pass ID that one of cut inverted, which means this cut is inverted
    #     # also the fake photon enriched region
    #     # use photon_selection to identified the type of events
    #     medium_and_control_cut = ((bitmap & mask1) == mask1) | ((bitmap & (mask2 + (3<<12))) == mask2) | ((bitmap & (mask3 + (3<<10))) == mask3) | ((bitmap & (mask4 + (3<<8))) == mask4) | ((bitmap & mask6) == mask6)
    #     selected_medium_or_control_photons = photons[medium_and_control_cut] # append the medium and control photons

    #     pass_selection1 = awkward.num(selected_medium_or_control_photons) >= 1  # select medium and control photons without fake photon
    #     awkward_utils.add_field(events, "pass_selection1", pass_selection1)  # has selected medium and control photons

    #     med_cand = awkward.firsts(selected_medium_or_control_photons) #similar to the gamma_cand

    #     bitmap = med_cand.vidNestedWPBitmap 
    #     photon_selection = (
    #         (((bitmap & mask1) == mask1) << 0) + 
    #         (((bitmap & (mask2 + (3<<12))) == mask2) << 1) + 
    #         ((((bitmap & (mask3 + (3<<10)))) == mask3) << 2) + 
    #         (((bitmap & (mask4 + (3<<8))) == mask4) << 3) + 
    #         (((bitmap & (mask5 + (3<<6))) == mask5) << 4) + 
    #         (((bitmap & mask5) == mask5) << 5) + 
    #         ((((bitmap & mask6) == mask6) & ((bitmap & mask5) != mask5)) << 6)
    #     )
    #     awkward_utils.add_field(events, "photon_selection", awkward.fill_none(photon_selection, -1))
    #     #pass_selection1 && photon_selection==1 or 5 or 7 -> build ture template from MC and data template from data
    #     #pass_selection1 && photon_selection==4 or 6 or 8  && chiso side band -> build fake tempalte from data
    #     #pass_selection1 && ((photon_selection!=1 && photon_selection==2) || (!1 && ==3) || (!1 && ==4) || (!=1 && ==5)) -> build non-prompt photon sample from data

    #     awkward_utils.add_field(events, "photon_is_barrel", awkward.fill_none(med_cand.isScEtaEB, -1))
    #     awkward_utils.add_field(events, "photon_is_endcap", awkward.fill_none(med_cand.isScEtaEE, -1))
    #     for field in ["pt", "eta", "phi", "sieie"]:
    #         awkward_utils.add_field(
    #             events,
    #             "photon_%s" % field,
    #             awkward.fill_none(getattr(med_cand, field), DUMMY_VALUE)
    #         )
    #     awkward_utils.add_field(events, "photon_chiso",  awkward.fill_none(med_cand.pfRelIso03_chg, DUMMY_VALUE))
    
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
        id_cut = photons.mvaID_WP80 == 1
        # eta_cut = ((photons.isScEtaEB & (photons.mvaID > options["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > options["mvaID_endcap"])))

        # electron veto
        e_veto_cut = (photons.electronVeto > options["e_veto"])
        
        new_pho = awkward.unflatten(awkward.unflatten(awkward.flatten(photons.electronIdx), [1]*awkward.sum(awkward.num(photons.electronIdx))), awkward.num(photons.electronIdx, axis=1))
        new_ele = awkward.broadcast_arrays(electrons.Idx[:,None], new_pho, depth_limit=2)[0]
        eg_overlap_cut = ~awkward.where(
            awkward.is_none(electrons.Idx),
            awkward.broadcast_arrays(photons.electronIdx, False)[1],
            awkward.flatten(awkward.any(new_pho[:, :, None] == new_ele, axis=-2), axis=-1)
        )

        # use_central_nano = options["use_central_nano"] # indicates whether we are using central nanoAOD (with some branches that are necessary for full diphoton preselection missing) or custom nanoAOD (with these branches added)

        all_cuts = pt_cut & eta_cut & id_cut & e_veto_cut & eg_overlap_cut
        # all_cuts = pt_cut & eta_cut # bing for CR selection

        self.register_cuts(
                names = ["pt", "eta", "id", "e_veto", "ele_pho_overlap", "all"],
                results = [pt_cut, eta_cut, id_cut, e_veto_cut, eg_overlap_cut, all_cuts],
                cut_type = "photon"
        )

        return all_cuts, id_cut & e_veto_cut

    def select_FSRphotons(self, FSRphotons, electrons, photons, options):
        FSR_pt_cut = FSRphotons.pt > options["pt"]
        FSR_eta_cut = abs(FSRphotons.eta) < options["eta"]

        FSR_iso_cut = FSRphotons.relIso03 < options["iso"]
        FSR_dROverEt2_cut = FSRphotons.dROverEt2 < options["dROverEt2"]

        FSR_all_cuts = FSR_pt_cut & FSR_eta_cut & FSR_iso_cut & FSR_dROverEt2_cut

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
