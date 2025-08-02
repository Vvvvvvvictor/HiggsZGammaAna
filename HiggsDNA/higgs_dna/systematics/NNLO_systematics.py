import awkward
import numpy

from correctionlib import _core

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.utils import systematic_from_bins, ic_systematic_from_bins

NNLO_SF_FILE = {
    "2016preVFP": "higgs_dna/systematics/data/2016preVFP_UL/nnlo_run2.json",
    "2016postVFP": "higgs_dna/systematics/data/2016preVFP_UL/nnlo_run2.json",
    "2017": "higgs_dna/systematics/data/2016preVFP_UL/nnlo_run2.json",
    "2018": "higgs_dna/systematics/data/2016preVFP_UL/nnlo_run2.json",
    "2022preEE": "higgs_dna/systematics/data/2022preEE_UL/nnlo_run3.json",
    "2022postEE": "higgs_dna/systematics/data/2022preEE_UL/nnlo_run3.json",
    "2023preBPix": "higgs_dna/systematics/data/2022preEE_UL/nnlo_run3.json",
    "2023postBPix": "higgs_dna/systematics/data/2022preEE_UL/nnlo_run3.json",
}

def nnlo_sf(events, year, central_only, input_collection):
    """
    Get the NNLO scale factor for a given year.
    """
    if year not in NNLO_SF_FILE:
        raise ValueError(f"Year {year} is not supported for NNLO scale factors.")
    
    required_fields = [(input_collection, "pt")]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = f"[NNLO_systematics : get_nnlo_sf] The events array is missing the following fields: {str(missing_fields)} which are needed as inputs."
        logger.exception(message)
        raise ValueError(message)

    # Load correction file
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(NNLO_SF_FILE[year]))
    nnlo_sf = evaluator["w_nnlo"]

    # Get input collection
    input_objects = events[input_collection]

    # Flatten objects then convert to numpy for compatibility with correctionlib
    n_objects = awkward.num(input_objects)
    objects_flattened = awkward.flatten(input_objects)

    # Clip pt values to valid range
    pt_values = numpy.clip(
        awkward.to_numpy(objects_flattened.pt),
        0.0,  # Minimum pt
        999999.0  # Maximum pt
    )

    # Calculate scale factor
    variations = {}
    
    # Central value
    sf = nnlo_sf.evalv(pt_values)
    variations["central"] = awkward.unflatten(sf, n_objects)

    # Systematic variations - set uncertainty to 0 as requested
    if not central_only:
        # Create up/down variations with zero uncertainty
        variations["up"] = awkward.unflatten(sf, n_objects)  # Same as central
        variations["down"] = awkward.unflatten(sf, n_objects)  # Same as central

    return variations

    