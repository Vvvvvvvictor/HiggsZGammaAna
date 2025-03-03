import awkward
import numpy

from higgs_dna.samples.file import File
from higgs_dna.utils import awkward_utils
from higgs_dna.constants import CENTRAL_WEIGHT, LUMI

import logging
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

# FIXME need to add option for specifying campaign (e.g. Prompt, UL, etc)

class Sample():
    """

    """
    def __init__(self, process, year, files, is_data = None, xs = None, bf = None, systematics = None, process_id = None, fpo = None, **kwargs):
        self.process = process
        self.year = year
        self.name = process + "_" + year
        self.files = files
        self.xs = xs
        self.bf = bf
        self.fpo = fpo
        self.systematics = systematics
        self.process_id = process_id

        if is_data is None:
            self.is_data = self.xs is None and self.bf is None
        else:
            self.is_data = is_data

        # Normalization factor including cross section and BF (if provided)
        self.norm_factor = 1
        if not self.is_data: 
            self.norm_factor = self.xs
        if self.bf is not None:
            self.norm_factor *= self.bf
        if self.year in LUMI.keys():
            self.lumi = LUMI[self.year]

        self.is_prepped = False


    def prep(self, events):
        if self.is_prepped:
            return events

        awkward_utils.add_field(events, CENTRAL_WEIGHT, awkward.ones_like(events.run))

        if not self.is_data:

            awkward_utils.add_field(
                    events = events,
                    name = CENTRAL_WEIGHT,
                    data = events[CENTRAL_WEIGHT] * events.Generator_weight,
                    overwrite = True
            )
            awkward_utils.add_field(
                    events = events,
                    name = CENTRAL_WEIGHT + "_initial",
                    data = events[CENTRAL_WEIGHT], # central value of weight that will stay unchanged after applying corrections
            )

        self.is_prepped = True

        return events
