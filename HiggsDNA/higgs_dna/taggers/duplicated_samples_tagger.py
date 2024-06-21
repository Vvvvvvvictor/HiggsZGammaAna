import awkward
import numpy
import json

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.taggers.tagger import Tagger

INPUT = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",  "HLT_Mu37_TkMu27", "HLT_IsoMu20",  "HLT_IsoMu22",  "HLT_IsoMu24", "HLT_IsoMu27",  "HLT_IsoTkMu20",  "HLT_IsoTkMu22", "HLT_IsoTkMu24",  "HLT_Mu50",  "HLT_Mu55", "HLT_TkMu50",  "HLT_IsoMu22_eta2p1",  "HLT_IsoMu24_eta2p1", "HLT_Mu45_eta2p1",  "HLT_Mu15_IsoVVVL_PFHT350",  "HLT_Mu15_IsoVVVL_PFHT400", "HLT_Mu15_IsoVVVL_PFHT450",  "HLT_Mu15_IsoVVVL_PFHT600", "HLT_Mu50_IsoVVVL_PFHT400", "HLT_Mu50_IsoVVVL_PFHT450", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",  "HLT_DoubleEle25_CaloIdL_MW", "HLT_DoublePhoton70",  "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55", "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", "HLT_Ele25_WPTight_Gsf", "HLT_Ele27_WPTight_Gsf", "HLT_Ele28_WPTight_Gsf", "HLT_Ele30_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_Ele20_WPLoose_Gsf", "HLT_Ele45_WPLoose_Gsf", "HLT_Ele105_CaloIdVT_GsfTrkIdT", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Ele135_CaloIdVT_GsfTrkIdT", "HLT_Ele145_CaloIdVT_GsfTrkIdT","HLT_Ele25_eta2p1_WPTight_Gsf",  "HLT_Ele27_eta2p1_WPTight_Gsf", "HLT_Ele20_eta2p1_WPLoose_Gsf",  "HLT_Ele25_eta2p1_WPLoose_Gsf", "HLT_Ele27_eta2p1_WPLoose_Gsf",  "HLT_Ele15_IsoVVVL_PFHT350", "HLT_Ele15_IsoVVVL_PFHT400",  "HLT_Ele15_IsoVVVL_PFHT450", "HLT_Ele15_IsoVVVL_PFHT600",  "HLT_Ele50_IsoVVVL_PFHT450", "HLT_Mu17_Photon30_IsoCaloId", 
# "HLT_PFMET90_PFMHT90_IDTight", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_PFMET100_PFMHT100_IDTight", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", "HLT_PFMET110_PFMHT110_IDTight", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", "HLT_PFMET120_PFMHT120_IDTight", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMET130_PFMHT130_IDTight", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", "HLT_PFMET140_PFMHT140_IDTight", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", "HLT_PFMET100_PFMHT100_IDTight_PFHT60", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", "HLT_PFMET110_PFMHT110_IDTight_PFHT60", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60", "HLT_PFMET120_PFMHT120_IDTight_PFHT60", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "HLT_PFMET130_PFMHT130_IDTight_PFHT60", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60", "HLT_PFMET140_PFMHT140_IDTight_PFHT60", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60", "HLT_PFMET120_PFMHT120_IDTight_HFCleaned", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned", "HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned", "HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1", "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1", "HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1""HLT_PFJet500", "HLT_PFHT125", "HLT_PFHT200", "HLT_PFHT300", "HLT_PFHT400", "HLT_PFHT475", "HLT_PFHT600",  "HLT_PFHT650", "HLT_PFHT800", "HLT_PFHT900", "HLT_PFHT180", "HLT_PFHT370", "HLT_PFHT430", "HLT_PFHT510", "HLT_PFHT590", "HLT_PFHT680", "HLT_PFHT780", "HLT_PFHT890", "HLT_PFHT1050", "HLT_PFHT250", "HLT_PFHT350"
]

class DuplicatedSamplesTagger(Tagger):
    """
    Tagger to remove those duplicated samples in different NanoAOD files
    """

    def __init__(self, name='DuplicatedSamplesTagger', options = {}, is_data = None):
        super(DuplicatedSamplesTagger, self).__init__(name, options, is_data)

        if not self.is_data:
            logger.exception("[DuplicatedSamplesTagger : __init__] Only need to remove duplicated events in data.")
            # raise RuntimeError

        self.input = INPUT
    
    def calculate_selection(self, file, tree):
        # Use HLT to select the events that is not duplicated
        self.hlt_branch = [x for x in self.input if x in tree.keys()]
        self.hlt_branch.append("run")
        # data = awkward.Array([])
        # for array in tree.iterate(self.hlt_branch, library="ak", how='zip', step_size=100000):
        #     data.concatenate(array)
        data = tree.arrays(self.hlt_branch, library="ak", how='zip')

        logger.debug("Duplicated_samples:  events fileds %s" % data.fields)

        cut = data.run>=0
        
        if not self.is_data:
            return cut

        # Only use specific events in different NanoAOD
        if "/DoubleMuon/" in file or "/SingleMuon/" in file or "/Muon" in file:
            double_mu_trig = self.get_double_muon_trig(data)
            if not "Double" in file:
                single_mu_trig = self.get_single_muon_trig(data)

            if "/DoubleMuon/" in file:
                cut = double_mu_trig
            elif "/Muon/" in file:
                cut = single_mu_trig | double_mu_trig
            elif "/SingleMuon/" in file:
                cut = single_mu_trig & numpy.logical_not(double_mu_trig)


        if "/DoubleEG/" in file or "/SingleElectron/" in file or "/EGamma" in file:
            double_mu_trig = self.get_double_muon_trig(data)
            single_mu_trig = self.get_single_muon_trig(data)
            double_eg_trig = self.get_double_eg_trig(data)
            if not "Double" in file:
                single_eg_trig = self.get_single_eg_trig(data)

            if "/DoubleEG/" in file:
                cut = double_eg_trig & numpy.logical_not(single_mu_trig | double_mu_trig)
            elif "/EGamma/" in file:
                cut = (single_eg_trig | double_eg_trig) & numpy.logical_not(single_mu_trig | double_mu_trig)
            elif "/SingleElectron/" in file:
                cut = single_eg_trig & numpy.logical_not(double_eg_trig) & numpy.logical_not(single_mu_trig | double_mu_trig)

        if "/MuonEG/" in file:
            double_mu_trig = self.get_double_muon_trig(data)
            single_mu_trig = self.get_single_muon_trig(data)
            double_eg_trig = self.get_double_eg_trig(data)
            single_eg_trig = self.get_single_eg_trig(data)
            muon_eg_trig   = self.get_muon_eg_trig(data)
        
            cut = muon_eg_trig & numpy.logical_not(single_eg_trig | double_eg_trig) & numpy.logical_not(single_mu_trig | double_mu_trig)
        
        cut = cut>0
        logger.debug("Duplicated_samples:  cut type %s" % cut.type)
        logger.debug("Duplicated_samples:  selected samples: %s(%s)" % (sum(cut), len(cut)))
        del data
        return cut

    def get_double_muon_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",  "HLT_Mu37_TkMu27"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | trig_cut
        return trig_cut
    

    def get_single_muon_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_IsoMu20",  "HLT_IsoMu22",  "HLT_IsoMu24", "HLT_IsoMu27",  "HLT_IsoTkMu20",  "HLT_IsoTkMu22", "HLT_IsoTkMu24",  "HLT_Mu50",  "HLT_Mu55", "HLT_TkMu50",  "HLT_IsoMu22_eta2p1",  "HLT_IsoMu24_eta2p1", "HLT_Mu45_eta2p1",  "HLT_Mu15_IsoVVVL_PFHT350",  "HLT_Mu15_IsoVVVL_PFHT400", "HLT_Mu15_IsoVVVL_PFHT450",  "HLT_Mu15_IsoVVVL_PFHT600", "HLT_Mu50_IsoVVVL_PFHT400", "HLT_Mu50_IsoVVVL_PFHT450"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | trig_cut
        return trig_cut


    def get_double_eg_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",  "HLT_DoubleEle25_CaloIdL_MW", "HLT_DoublePhoton70",  "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55", "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | trig_cut
        return trig_cut


    def get_single_eg_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_Ele25_WPTight_Gsf", "HLT_Ele27_WPTight_Gsf", "HLT_Ele28_WPTight_Gsf", "HLT_Ele30_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_Ele20_WPLoose_Gsf", "HLT_Ele45_WPLoose_Gsf", "HLT_Ele105_CaloIdVT_GsfTrkIdT", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Ele135_CaloIdVT_GsfTrkIdT", "HLT_Ele145_CaloIdVT_GsfTrkIdT","HLT_Ele25_eta2p1_WPTight_Gsf",  "HLT_Ele27_eta2p1_WPTight_Gsf", "HLT_Ele20_eta2p1_WPLoose_Gsf",  "HLT_Ele25_eta2p1_WPLoose_Gsf", "HLT_Ele27_eta2p1_WPLoose_Gsf",  "HLT_Ele15_IsoVVVL_PFHT350", "HLT_Ele15_IsoVVVL_PFHT400",  "HLT_Ele15_IsoVVVL_PFHT450", "HLT_Ele15_IsoVVVL_PFHT600",  "HLT_Ele50_IsoVVVL_PFHT450"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | trig_cut
        return trig_cut

    
    def get_muon_eg_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_Mu17_Photon30_IsoCaloId"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | trig_cut
        return trig_cut


    def get_met_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_PFMET90_PFMHT90_IDTight", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_PFMET100_PFMHT100_IDTight", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", "HLT_PFMET110_PFMHT110_IDTight", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", "HLT_PFMET120_PFMHT120_IDTight", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMET130_PFMHT130_IDTight", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", "HLT_PFMET140_PFMHT140_IDTight", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", "HLT_PFMET100_PFMHT100_IDTight_PFHT60", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", "HLT_PFMET110_PFMHT110_IDTight_PFHT60", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60", "HLT_PFMET120_PFMHT120_IDTight_PFHT60", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "HLT_PFMET130_PFMHT130_IDTight_PFHT60", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60", "HLT_PFMET140_PFMHT140_IDTight_PFHT60", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60", "HLT_PFMET120_PFMHT120_IDTight_HFCleaned", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned", "HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned", "HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1", "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1", "HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | self.hlt_branch
        return trig_cut

    
    def get_jetht_trig(self, data):
        trig_cut = data.run < 0
        for trig in ["HLT_PFJet500", "HLT_PFHT125", "HLT_PFHT200", "HLT_PFHT300", "HLT_PFHT400", "HLT_PFHT475", "HLT_PFHT600",  "HLT_PFHT650", "HLT_PFHT800", "HLT_PFHT900", "HLT_PFHT180", "HLT_PFHT370", "HLT_PFHT430", "HLT_PFHT510", "HLT_PFHT590", "HLT_PFHT680", "HLT_PFHT780", "HLT_PFHT890", "HLT_PFHT1050", "HLT_PFHT250", "HLT_PFHT350"]:
            if trig in self.hlt_branch:
                trig_cut = data[trig] | self.hlt_branch
        return trig_cut


        
"""

=============add in analysis.load_events==================

from higgs_dna.taggers.duplicated_samples_tagger import DuplicatedSamplesTagger

duplicated_sample_remover = DuplicatedSamplesTagger(is_data=is_data)
duplicated_cut = duplicated_sample_remover.calculate_selection(file, tree, config["sample"]["year"])

events_file = events_file[duplicated_cut]

"""

    
