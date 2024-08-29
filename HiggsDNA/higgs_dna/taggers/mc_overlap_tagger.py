import awkward
import numpy
import json

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.taggers.tagger import Tagger
from higgs_dna.selections import object_selections

INPUT = ["GenPart_eta", "GenPart_genPartIdxMother", "GenPart_mass", "GenPart_pdgId", "GenPart_phi", "GenPart_pt", "GenPart_status","GenPart_statusFlags"]

class MCOverlapTagger(Tagger):
    """
    Tagger to remove those duplicated samples in different NanoAOD files
    """

    def __init__(self, name='MCOverlapTagger', options = {}, is_data = None, year = None):
        super(MCOverlapTagger, self).__init__(name, options, is_data, year)

        if self.is_data:
            logger.exception("[MCOverlapTagger : __init__] Only need to remove overlap events in MC samples, using the events only with fake photon or true photon in one MC dataset.")
            # raise RuntimeError

        self.input = INPUT
    
    def overlap_selection(self, file, tree):
        # Use HLT to select the events that is not duplicated
        self.branch = [x for x in self.input if x in tree.keys()]
        self.branch.append("run")
        data = tree.arrays(self.branch, library='ak', how='zip')
        # logger.debug("Duplicated_samples:  events fileds %s" % data.fields)

        cut = data.run>=0
        if self.is_data:
            return cut

        # Only use specific events in different MC NanoAOD
        if "DYJetsToLL" in file or "EWKZ2J" in file or "TTTo2L2Nu" in file:
            cut = self.get_n_iso_photon(data) == 0
        elif "ZGToLLG" in file or "ZG2JToG2L2J" in file or "TTGJets" in file:
            cut = self.get_n_iso_photon(data) > 0
        
        # cut = cut>0
        logger.debug("MC_overlap_samples:  cut type: %s" % cut.type)
        logger.debug("MC_overlap_samples: selected samples: %s(%s)" % (sum(cut), len(cut)))
        return cut


    def get_n_iso_photon(self, data):
        """
        Select isolation photons in data
        Return the number of isolation photons(int)
        """
        
        """
        statusFlags usage: (data.GenPart.statusFlags // numpy.power(2, i)) % 2 == 1
        "statusFlags" is a number with 14 bits. 
        Filling "1" on corresponding digit when the particle meets one of the 14 conditions, else remaining "0".
        Echo paticles can meet more than one kind of condition, thus, more than one digit in "statusFlags" is "1".
        """
        iso_photons_cut = (data.GenPart.pdgId == 22) & (data.GenPart.pt > 9) & (( (data.GenPart.statusFlags & 0x1) != 0 ) | ( (data.GenPart.statusFlags & 0x100) != 0 )) #& (abs(data.GenPart.eta) < 2.6)
        iso_photons = data.GenPart[iso_photons_cut]

        truth_objects_cut =  (data.GenPart.pdgId != 22) & (data.GenPart.pt > 5) & ( (data.GenPart.statusFlags & 0x100) != 0 ) 
        truth_objects = data.GenPart[truth_objects_cut]

        iso_cut = object_selections.delta_R(iso_photons, truth_objects, 0.05)
        iso_photons = iso_photons[iso_cut]

        n_iso_photons = awkward.num(iso_photons, axis=1)
        logger.debug("MC_overlap_samples:  variable type of n_iso_photons: %s" % n_iso_photons.type)
        return n_iso_photons

        
"""

=============add in analysis.load_events==================

from higgs_dna.taggers.mc_overlap_tagger import MCOverlapTagger

mc_overlap_remover = MCOverlapTagger(is_data=is_data)
overlap_cut = mc_overlap_remover.overlap_selection(file, tree)

events_file = events_file[overlap_cut]

"""

    
