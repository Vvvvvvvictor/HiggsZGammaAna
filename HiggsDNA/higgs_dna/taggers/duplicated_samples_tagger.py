import awkward
import numpy
import json

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.taggers.tagger import Tagger

INPUT = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Mu37_TkMu27", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_IsoMu24_eta2p1", "HLT_Mu15_IsoVVVL_PFHT450", "HLT_Mu15_IsoVVVL_PFHT600", "HLT_Mu50_IsoVVVL_PFHT450", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_DoubleEle25_CaloIdL_MW",  "HLT_DoublePhoton70" ,"HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55",  "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", "HLT_Ele27_WPTight_Gsf", "HLT_Ele28_WPTight_Gsf", "HLT_Ele30_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_Ele20_WPLoose_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Ele135_CaloIdVT_GsfTrkIdT", "HLT_Ele15_IsoVVVL_PFHT450", "HLT_Ele15_IsoVVVL_PFHT600", "HLT_Ele50_IsoVVVL_PFHT450", "HLT_Mu17_Photon30_IsoCaloId"]

class DuplicatedSamplesTagger(Tagger):
    """
    Tagger to remove those duplicated samples in different NanoAOD files
    """

    def __init__(self, name='DuplicatedSamplesTagger', options = {}, is_data = None, year = None):
        super(DuplicatedSamplesTagger, self).__init__(name, options, is_data, year)

        if not self.is_data:
            logger.exception("[DuplicatedSamplesTagger : __init__] Only need to remove duplicated events in data.")
            # raise RuntimeError

        self.input = INPUT
    
    def calculate_selection(self, file, tree, year):
        # Use HLT to select the events that is not duplicated
        hlt_branch = [x for x in self.input if x in tree.keys()]
        hlt_branch.append("run")
        data = tree.arrays(hlt_branch, library='ak', how='zip')

        logger.debug("Duplicated_samples:  events fileds %s" % data.fields)

        if not self.is_data:
            return data.HLT_IsoMu20>=0

        # Only use specific events in different NanoAOD
        if "/DoubleMuon/" in file or "/SingleMuon/" in file or "/Muon/" in file:
            double_mu_trig = self.get_double_muon_trig(data, year)
            if not "Double" in file:
                single_mu_trig = self.get_single_muon_trig(data, year)

            if "/DoubleMuon/" in file:
                cut = double_mu_trig
            elif "/Muon/" in file:
                cut = single_mu_trig | double_mu_trig
            elif "/SingleMuon/" in file:
                cut = single_mu_trig & (1-double_mu_trig)


        if "/DoubleEG/" in file or "/SingleElectron/" in file or "/EGamma/" in file:
            double_mu_trig = self.get_double_muon_trig(data, year)
            single_mu_trig = self.get_single_muon_trig(data, year)
            double_eg_trig = self.get_double_eg_trig(data, year)
            if not "Double" in file:
                single_eg_trig = self.get_single_eg_trig(data, year)

            if "/DoubleEG/" in file:
                cut = double_eg_trig & (1-(single_mu_trig | double_mu_trig))
            elif "/EGamma/" in file:
                cut = (single_eg_trig | double_eg_trig) & (1-(single_mu_trig | double_mu_trig))
            elif "/SingleElectron/" in file:
                cut = single_eg_trig & (1-double_eg_trig) & (1-(single_mu_trig | double_mu_trig))

        if "/MuonEG/" in file:
            double_mu_trig = self.get_double_muon_trig(data, year)
            single_mu_trig = self.get_single_muon_trig(data, year)
            double_eg_trig = self.get_double_eg_trig(data, year)
            single_eg_trig = self.get_single_eg_trig(data, year)
            muon_eg_trig   = self.get_muon_eg_trig(data, year)
        
            cut = muon_eg_trig & (1-(single_eg_trig | double_eg_trig)) & (1-(single_mu_trig | double_mu_trig))
        
        cut = cut>0
        logger.debug("Duplicated_samples:  cut type %s" % cut.type)
        logger.debug("Duplicated_samples:  selected samples: %s(%s)" % (sum(cut), len(cut)))
        return cut

    def get_double_muon_trig(self, data, year):
        if year == "2016":
            return data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
        elif year == "2017":
            return data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8
        elif year == "2018":
            return data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 | data.HLT_Mu37_TkMu27
        else:
            return data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 | data.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 | data.HLT_Mu37_TkMu27
    

    def get_single_muon_trig(self, data, year):
        if year == "2016":
            return data.HLT_IsoMu20 | data.HLT_IsoMu24 | data.HLT_IsoMu27 | data.HLT_Mu50 | data.HLT_Mu55 | data.HLT_Mu15_IsoVVVL_PFHT600
        elif year == "2017":
            return data.HLT_IsoMu20 | data.HLT_IsoMu24 | data.HLT_IsoMu27 | data.HLT_Mu50 | data.HLT_IsoMu24_eta2p1
        elif year == "2018":
            return data.HLT_IsoMu20 | data.HLT_IsoMu24 | data.HLT_IsoMu27 | data.HLT_Mu50 | data.HLT_Mu55 | data.HLT_IsoMu24_eta2p1 | data.HLT_Mu15_IsoVVVL_PFHT450 | data.HLT_Mu15_IsoVVVL_PFHT600 | data.HLT_Mu50_IsoVVVL_PFHT450
        else:
            return data.HLT_IsoMu20 | data.HLT_IsoMu24 | data.HLT_IsoMu27 | data.HLT_Mu50 | data.HLT_Mu55 | data.HLT_IsoMu24_eta2p1 | data.HLT_Mu15_IsoVVVL_PFHT450 | data.HLT_Mu15_IsoVVVL_PFHT600 | data.HLT_Mu50_IsoVVVL_PFHT450


    def get_double_eg_trig(self, data, year):
        if year == "2016":
            return data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ | data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
        elif year == "2017":
            return data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ | data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
        elif year == "2018":
            return data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ | data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL | data.HLT_DoubleEle25_CaloIdL_MW | data.HLT_DoublePhoton70 | data.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 | data.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95
        else:
            return data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ | data.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL | data.HLT_DoubleEle25_CaloIdL_MW | data.HLT_DoublePhoton70 | data.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId | data.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55 | data.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 | data.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95


    def get_single_eg_trig(self, data, year):
        if year == "2016":
            return data.HLT_Ele27_WPTight_Gsf | data.HLT_Ele115_CaloIdVT_GsfTrkIdT | data.HLT_Ele15_IsoVVVL_PFHT600
        elif year == "2017":
            return data.HLT_Ele27_WPTight_Gsf | data.HLT_Ele32_WPTight_Gsf_L1DoubleEG | data.HLT_Ele35_WPTight_Gsf
        elif year == "2018":
            return data.HLT_Ele27_WPTight_Gsf | data.HLT_Ele28_WPTight_Gsf | data.HLT_Ele32_WPTight_Gsf | data.HLT_Ele32_WPTight_Gsf_L1DoubleEG | data.HLT_Ele35_WPTight_Gsf | data.HLT_Ele20_WPLoose_Gsf | data.HLT_Ele115_CaloIdVT_GsfTrkIdT | data.HLT_Ele135_CaloIdVT_GsfTrkIdT | data.HLT_Ele15_IsoVVVL_PFHT450 | data.HLT_Ele15_IsoVVVL_PFHT600 | data.HLT_Ele50_IsoVVVL_PFHT450 | data.HLT_Ele30_WPTight_Gsf 
        else:
            return data.HLT_Ele27_WPTight_Gsf | data.HLT_Ele28_WPTight_Gsf | data.HLT_Ele32_WPTight_Gsf | data.HLT_Ele32_WPTight_Gsf_L1DoubleEG | data.HLT_Ele35_WPTight_Gsf | data.HLT_Ele20_WPLoose_Gsf | data.HLT_Ele115_CaloIdVT_GsfTrkIdT | data.HLT_Ele135_CaloIdVT_GsfTrkIdT | data.HLT_Ele15_IsoVVVL_PFHT450 | data.HLT_Ele15_IsoVVVL_PFHT600 | data.HLT_Ele50_IsoVVVL_PFHT450 | data.HLT_Ele30_WPTight_Gsf 

    
    def get_muon_eg_trig(self, data, year):
        if year == "2016":
            return data.HLT_Mu17_Photon30_IsoCaloId
        elif year == "2017":
            return data.HLT_Mu17_Photon30_IsoCaloId
        elif year == "2018":
            return data.HLT_Mu17_Photon30_IsoCaloId
        else:
            return data.HLT_Mu17_Photon30_IsoCaloId


    def get_met_trig(self, data, year):
        return data.HLT_PFMET90_PFMHT90_IDTight | data.HLT_PFMETNoMu90_PFMHTNoMu90_IDTight | data.HLT_PFMET100_PFMHT100_IDTight | data.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight | data.HLT_PFMET110_PFMHT110_IDTight | data.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight | data.HLT_PFMET120_PFMHT120_IDTight | data.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight | data.HLT_PFMET130_PFMHT130_IDTight | data.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight | data.HLT_PFMET140_PFMHT140_IDTight | data.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight | data.HLT_PFMET100_PFMHT100_IDTight_PFHT60 | data.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 | data.HLT_PFMET110_PFMHT110_IDTight_PFHT60 | data.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60 | data.HLT_PFMET120_PFMHT120_IDTight_PFHT60 | data.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 | data.HLT_PFMET130_PFMHT130_IDTight_PFHT60 | data.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60 | data.HLT_PFMET140_PFMHT140_IDTight_PFHT60 | data.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60 | data.HLT_PFMET120_PFMHT120_IDTight_HFCleaned | data.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned | data.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned | data.HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1 | data.HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1 | data.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1

    
    def get_jetht_trig(self, data, year):
        return data.HLT_PFJet500 | data.HLT_PFHT125 | data.HLT_PFHT200 | data.HLT_PFHT300 | data.HLT_PFHT400 | data.HLT_PFHT475 | data.HLT_PFHT600 | data.HLT_PFHT650 | data.HLT_PFHT800 | data.HLT_PFHT900 | data.HLT_PFHT180 | data.HLT_PFHT370 | data.HLT_PFHT430 | data.HLT_PFHT510 | data.HLT_PFHT590 | data.HLT_PFHT680 | data.HLT_PFHT780 | data.HLT_PFHT890 | data.HLT_PFHT1050 | data.HLT_PFHT250 | data.HLT_PFHT350


        
"""

=============add in analysis.load_events==================

from higgs_dna.taggers.duplicated_samples_tagger import DuplicatedSamplesTagger

duplicated_sample_remover = DuplicatedSamplesTagger()
duplicated_cut = duplicated_sample_remover.calculate_selection(file, tree)

events_file = events_file[duplicated_cut]

"""

    
