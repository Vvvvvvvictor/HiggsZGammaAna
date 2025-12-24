import awkward
import numpy

from higgs_dna.utils import misc_utils, awkward_utils

import logging
logger = logging.getLogger(__name__)

def ps_isr_sf(events, year, central_only):
    """
    Apply PS ISR scale factors based on PSWeight variations.
    """

    required_fields = ["PSWeight"]
    missing_fields = awkward_utils.missing_fields(events, required_fields)
    logger.debug(f"Missing fields for PS ISR SF: {missing_fields}")

    variations = {}
    variations["central"] = awkward.ones_like(events.PSWeight[:, 0])

    if not central_only:
        for i, var in enumerate(["up", "down"]):
            variations[var] = events.PSWeight[:, i]

    return variations

def ps_fsr_sf(events, year, central_only):
    """
    Apply PS FSR scale factors based on PSWeight variations.
    """

    required_fields = ["PSWeight"]
    missing_fields = awkward_utils.missing_fields(events, required_fields)
    logger.debug(f"Missing fields for PS FSR SF: {missing_fields}")

    variations = {}
    variations["central"] = awkward.ones_like(events.PSWeight[:, 0])

    if not central_only:
        for i, var in enumerate(["up", "down"]):
            variations[var] = events.PSWeight[:, i * 2 + 1]

    return variations