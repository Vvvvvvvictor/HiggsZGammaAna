import pandas as pd
import uproot
import numpy as np
import os

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)


def get_nominal_data(year, sample):
    data = uproot.open(os.path.join(path, f"{sample}/{year}.root"))["inclusive"].arrays(library="pd")
    return data