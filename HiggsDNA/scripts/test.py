import uproot
import time
import os

file = "root://cmsxrootd.fnal.gov//store/data/Run2017B/DoubleEG/NANOAOD/UL2017_MiniAODv2_NanoAODv9-v1/70000/0E5E3581-6C20-9D44-B44C-576C11B72848.root"
# file = "root://cmsxrootd.fnal.gov//store/data/Run2017F/DoubleMuon/NANOAOD/UL2017_MiniAODv2_NanoAODv9-v1/120000/F369A309-E4D4-FA45-B7BD-02C46C1178F0.root"
start = time.time_ns()
try:
    with uproot.open(file, timeout = 15) as f:
        print(f["Events"]["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"])
except Exception:
    raise RuntimeError("xrootd failed")
print((time.time_ns()-start)/1000000000.0)