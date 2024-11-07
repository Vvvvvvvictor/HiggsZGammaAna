import uproot
import awkward as ak
import numpy
import numba
import vector
import os

from pdb import set_trace

branches = [
    "run",
    "luminosityBlock",
    "event",
    "Electron_pt",
    "Electron_eta",
    "Electron_phi",
    "Electron_mass",
    "Electron_charge"
]
# Load the NanoAOD file
print("Loading NanoAOD file...")
file = "root://xrootd-cms.infn.it//store/data/Run2018A/EGamma/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/50D18A5D-64D2-5A45-A2F5-C5B8DB90514F.root"
os.system(f"xrdcp '{file}' '/tmp/{os.getpid()}/{os.path.basename(file)}'")
f = uproot.open(f'/tmp/{os.getpid()}/{os.path.basename(file)}')
# file = uproot.open("/eos/home-j/jiehan/02466EF4-1D3E-BF4C-B43F-0DCC7C61BCC2.root")
data = f["Events"].arrays(branches, library="ak", how="zip")
data = data[(data["run"]==316470) & (data["luminosityBlock"]==370) & (data["event"]==486186232)]
electrons = data.Electron
electrons = ak.Array(electrons, with_name = "Momentum4D")
ee_pairs = ak.combinations(electrons, 2, fields=["i0", "i1"])
os_cut = (ee_pairs.i0.charge * ee_pairs.i1.charge < 0)
ee_pairs["z"] = ee_pairs.i0 + ee_pairs.i1
print(os_cut)
set_trace()
print(data)
