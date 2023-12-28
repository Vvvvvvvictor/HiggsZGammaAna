import awkward
import numpy
import json

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.taggers.tagger import Tagger
from higgs_dna.constants import GOLDEN_JSON
from higgs_dna.utils import misc_utils

class GoldenJsonTagger(Tagger):
    """
    Tagger to select events in the golden json.
    Should only be applied to data.
    NOTE: TagSequence will automatically add this for data if a corresponding json file for the relevant year is listed in higgs_dna.constants.GOLDEN_JSON. You should not have to add it yourself.
    """
    def __init__(self, name = "golden_json_tagger", options = {}, is_data = True, year = None):
        super(GoldenJsonTagger, self).__init__(name, options, is_data, year)

        if not self.is_data:
            logger.exception("[GoldenJsonTagger : __init__] Golden json should only be applied on data.")
            raise RuntimeError()

        if self.year not in GOLDEN_JSON.keys():
            logger.warning("[A corresponding json file was not found in higgs_dna.constants.GOLDEN_JSON for year '%s'. This tagger will apply a dummy selection." % (self.year))
            self.json_file = None
            self.golden_json = None

        else:
            self.json_file = misc_utils.expand_path(GOLDEN_JSON[self.year])
            with open(self.json_file, "r") as f_in:
                self.golden_json = json.load(f_in)

    
    def calculate_selection(self, events):
        # If no valid golden json file, return all True
        if self.json_file is None:
            cut = awkward.ones_like(events.run, dtype=bool)#abs(events.run) >= 0
        else:
            # Start with dummy all False cut and mark passing run x lumis as True
            cut = awkward.zeros_like(events.run, dtype=bool) 
 
            # Find all runs present in events
            runs = numpy.unique(awkward.to_numpy(events.run))

            # Loop through each run x lumi pair
            for run in runs:
                # If run is not in json file, all lumis are bad
                if str(run) not in self.golden_json.keys():
                    continue # cut is initialized to all False, so no need to update
                
                for lumi_range in self.golden_json[str(run)]:
                    cut = cut | ((events.run == run) & (events.luminosityBlock >= lumi_range[0]) & (events.luminosityBlock <= lumi_range[1]))

        self.register_cuts(
                names = ["golden json"],
                results = [cut]
        )

        self.Get_ProcessedLumi(events, cut)

        return cut, events


    def Get_ProcessedLumi(self, events, cut):

        events_golden = events[cut]
        lumiBlocks = events_golden.luminosityBlock
        runs = events_golden.run
        events = events_golden.event

        # dic to store the final results
        result_dict = {}

        # find runs
        unique_values_runs = set(awkward.to_list(runs))

        for unique_run in unique_values_runs:
            mask = (runs == unique_run)  
            filtered_lumiBlocks = lumiBlocks[mask] 
            filtered_events = events[mask] 

            filtered_lumiBlocks_list = awkward.to_list(filtered_lumiBlocks)

            #print("unique_run", unique_run, "mask", mask, "filtered_lumiBlocks_list:", filtered_lumiBlocks_list , "event:", awkward.to_list(filtered_events))
            filtered_lumiBlocks_list.sort()
            #print(filtered_lumiBlocks_list)

            # find the lumiBlock range
            ranges = []
            start = filtered_lumiBlocks_list[0]
            for i in range(1, len(filtered_lumiBlocks_list)):
                if filtered_lumiBlocks_list[i] - filtered_lumiBlocks_list[i - 1] > 1:
                    if start == filtered_lumiBlocks_list[i - 1]:
                        ranges.append([start, start])
                    else:
                        ranges.append([start, filtered_lumiBlocks_list[i - 1]])
                    start = filtered_lumiBlocks_list[i]
            #print("ranges,", ranges)
            if start == filtered_lumiBlocks_list[-1]:
                ranges.append([start, start])
            else:
                ranges.append([start, filtered_lumiBlocks_list[-1]])

            if str(unique_run) not in result_dict:
                result_dict[str(unique_run)] = [ranges]
            else:
                result_dict[str(unique_run)].append(ranges)

        for key in result_dict.keys():
            if len(result_dict[key]) == 1:
                result_dict[key] = result_dict[key][0]

        print("[[INFO]] Processed Lumi:",  result_dict)




