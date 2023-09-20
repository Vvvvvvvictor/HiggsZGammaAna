import awkward
import time
import numpy
import numba
import vector

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
        "mvaID_endcap" : -0.59,
        "e_veto" : 0.5
    },
    "zgammas" : {
        "relative_pt_gamma" : 0.136363636364,
        "mass_h" : [100., 180.],
        "mass_sum" : 185,
        "select_highest_pt_sum" : True
    },
    "trigger" : {
        "2016" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2016UL_preVFP" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2016UL_postVFP" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2017" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2018" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"]
    },
    "electrons" : {
        "pt" : 7.0
    },
    "muons" : {
        "pt" : 5.0
    },
    "jets" : {
        "pt" : 30.0,
        "eta" : 4.7,
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
    },
    "gen_info" : {
        "calculate" : False,
        "max_dr" : 0.2,
        "max_pt_diff" : 15.
    }  
}


# Diphoton preselection below synced with flashgg, see details in:
#   - https://indico.cern.ch/event/1071721/contributions/4551056/attachments/2320292/3950844/HiggsDNA_DiphotonPreselectionAndSystematics_30Sep2021.pdf

class ZGammaTaggerDY(Tagger):
    def __init__(self, name = "default_zgamma_tagger", options = {}, is_data = None, year = None):
        super(ZGammaTaggerDY, self).__init__(name, options, is_data, year)

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

        photon_selection = self.select_photons(
                photons = events.Photon,
                rho = rho,
                options = self.options["photons"]
        )

        self.overlap_removal(events=events)

        zgamma_selection, zgammas = self.produce_and_select_zgammas(
                events = events,
                photons = events.Photon[photon_selection],
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


    def produce_and_select_zgammas(self, events, photons, options):
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
            name = "SelectedJet",
            tagger = self
        )

        jets = awkward_utils.add_field(
            events = events,
            name = "SelectedJet",
            data = events.Jet[jet_cut]
        )

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
        awkward_utils.add_field(events, "n_electrons", n_electrons, overwrite=True)

        n_muons = awkward.num(muons)
        awkward_utils.add_field(events, "n_muons", n_muons, overwrite=True)

        n_leptons = n_electrons + n_muons
        awkward_utils.add_field(events, "n_leptons", n_leptons, overwrite=True)

        n_jets = awkward.num(jets)
        awkward_utils.add_field(events, "n_jets", n_jets, overwrite=True)


        # lepton-photon overlap removal 
        clean_photon_mask = object_selections.delta_R(photons, muons, 0.4) & object_selections.delta_R(photons, electrons, 0.4)
        photons = photons[clean_photon_mask]

        # PDG ID
        electrons = awkward.with_field(electrons, awkward.ones_like(electrons.pt) * 11, "id")
        muons = awkward.with_field(muons, awkward.ones_like(muons.pt) * 13, "id")

        # Sort objects by pt
        photons = photons[awkward.argsort(photons.pt, ascending=False, axis=1)]
        electrons = electrons[awkward.argsort(electrons.pt, ascending=False, axis=1)]
        muons = muons[awkward.argsort(muons.pt, ascending=False, axis=1)]

        # self.select_fake_and_medium_photons(events=events, photons=photons)

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        photons = awkward.Array(photons, with_name = "Momentum4D")
        electrons = awkward.Array(electrons, with_name = "Momentum4D")
        muons = awkward.Array(muons, with_name = "Momentum4D")
        
        # Construct di-electron/di-muon pairs
        ee_pairs = awkward.combinations(electrons, 2, fields = ["LeadLepton", "SubleadLepton"])
        ee_cut = (ee_pairs.LeadLepton.pt > 25) & (ee_pairs.SubleadLepton.pt > 15)
        ee_pairs = ee_pairs[ee_cut]

        mm_pairs = awkward.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
        mm_cut = (mm_pairs.LeadLepton.pt > 20) & (mm_pairs.SubleadLepton.pt > 10)
        mm_pairs = mm_pairs[mm_cut]

        # Concatenate these together
        z_cands = awkward.concatenate([ee_pairs, mm_pairs], axis = 1)
        z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton # these add as 4-vectors since we registered them as "Momentum4D" objects

        # Make Z candidate-level cuts
        os_cut = z_cands.LeadLepton.charge * z_cands.SubleadLepton.charge == -1
        mass_cut = (z_cands.ZCand.mass > 50.)
        z_cut = os_cut & mass_cut
        z_cands = z_cands[z_cut] # OSSF lepton pairs with m_ll > 50.
        
        has_z_cand = awkward.num(z_cands) >= 1
        z_cands = z_cands[awkward.argsort(abs(z_cands.ZCand.mass - 91.), axis = 1)] # take the one with mass closest to mZ
        z_cand = awkward.firsts(z_cands)

        # Add Z-related fields to array
        for field in ["pt", "eta", "phi", "mass", "charge", "id"]:
            if not field in ["charge", "id"]:
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
        has_gamma_cand = (awkward.num(photons) >= 1) & (events.n_iso_photons == 0) # only for dy samples
        gamma_cand = awkward.firsts(photons)
        gamma_mvaID_WPL = ((gamma_cand.isScEtaEB & (gamma_cand.mvaID > self.options["photons"]["mvaID_barrel"])) | (gamma_cand.isScEtaEE & (gamma_cand.mvaID > self.options["photons"]["mvaID_endcap"])))

        # Add gamma-related fields to array
        for field in ["pt", "eta", "phi", "mass", "mvaID", "energyErr", "sieie", "hoe", "r9", "genPartFlav", "mvaID_WP80"]:
            awkward_utils.add_field(
                events,
                "gamma_%s" % field,
                awkward.fill_none(getattr(gamma_cand, field), DUMMY_VALUE)
            )
        awkward_utils.add_field(events, "gamma_mvaID_WPL",  gamma_mvaID_WPL)
        awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg)
        awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all)

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

        # Make trigger cuts 
        if self.year is not None:
            trigger_cut = awkward.num(events.Photon) < 0 # dummy cut, all False
            for hlt in self.options["trigger"][self.year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    trigger_cut = (trigger_cut) | (events[hlt] == True)
        else:
            trigger_cut = awkward.num(events.Photon) >= 0 # dummy cut, all True

        all_cuts = trigger_cut & has_z_cand & has_gamma_cand & sel_h_1 & sel_h_2 & sel_h_3

        self.register_cuts(
            names = ["trigger_cut", "has_z_cand", "has_gamma_cand", "sel_h_1", "sel_h_2", "sel_h_3", "all cuts"],
            results = [trigger_cut, has_z_cand, has_gamma_cand, sel_h_1, sel_h_2, sel_h_3, all_cuts],
            cut_type = "zgammas"
        )

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
    #     medium_and_control_cut = ((bitmap & mask1) == mask1) | ((bitmap & (mask2 + (3<<12))) == mask2) | ((bitmap & (mask3 + (3<<10))) == mask3) | ((bitmap & mask6) == mask6)
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
    #         (((bitmap & mask6 == mask6) & ((bitmap & mask5) != mask5)) << 6)
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

    #     isprompt_mask = (1 << 0)  #isPrompt
    #     isdirectprompttaudecayproduct_mask = (1 << 5)  #isDirectPromptTauDecayProduct
    #     isfromhardprocess_mask = (1 << 8)  #isPrompt

    #     if hasattr(events, 'nGenPart') and hasattr(med_cand,'genPartIdx'):
    #         match_gen = awkward.local_index(events.GenPart, axis=1) == med_cand.genPartIdx

    #         genpart_cut = (events.GenPart.pt > 5) & (abs(events.GenPart.pdgId) == 22) 
    #         photon_isprompt = (med_cand.genPartIdx >=0) &  ((events.GenPart[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (events.GenPart[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) | (events.GenPart[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) & (awkward.sum(~object_selections.delta_R(events.GenPart[genpart_cut], med_cand, 0.3), axis=1) >= 1)

    #         photon_isprompt = awkward.sum(photon_isprompt, axis=1)
    #         awkward_utils.add_field(events, "photon_isprompt",  awkward.fill_none(photon_isprompt, -1))

    #         awkward_utils.add_field(events, "photon_genPartFlav",  awkward.fill_none(med_cand.genPartFlav, -1))


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
        
        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgHiggs",
                objects = gen_hzg.GenParent,
                n_objects = 1
        )

        return zgammas 
        

    def select_photons(self, photons, rho, options):
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
        # eta_cut = ((photons.isScEtaEB & (photons.mvaID > options["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > options["mvaID_endcap"])))

        # electron veto
        e_veto_cut = photons.electronVeto > options["e_veto"]

        use_central_nano = options["use_central_nano"] # indicates whether we are using central nanoAOD (with some branches that are necessary for full diphoton preselection missing) or custom nanoAOD (with these branches added)


        all_cuts = pt_cut & eta_cut & e_veto_cut

        self.register_cuts(
                names = ["pt", "eta", "e_veto", "all"],
                results = [pt_cut, eta_cut, e_veto_cut, all_cuts],
                cut_type = "photon"
        )

        return all_cuts

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
