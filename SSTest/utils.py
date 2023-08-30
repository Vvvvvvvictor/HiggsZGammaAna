import ROOT
from tqdm import trange
import uproot, uproot3
import awkward
import os
import re

def ReadFile(file, tree = "test", selections = []):
    '''
    read file and apply selections
    selections: a list of selection(<str>). If there is '&' or '|' in one item, MUST use "()" to enclose each subitem.(exp. '(H_eta>1.5) & (H_pt>10)')
    '''
    print("Reading {}:{}...".format(file, tree))
    variabels = ["weight", "H_mass"]
    decro_sel = []
    for i in selections:
        var_can = list(set([j for j in re.split('\W+', i) if "_" in j]))
        variabels += var_can
        for j in var_can:
            i = i.replace(j, "arrays."+j)
        decro_sel.append(i)
    arrays = uproot.open(file+':'+tree).arrays(variabels,library='ak')
    for i in decro_sel:
        arrays = arrays[eval(i)]
    return arrays

def AddHist(arrays, variable, hist):
    for i in trange(0, len(arrays[variable])):
        hist.Fill(float(arrays[variable][i]), float(arrays['weight'][i]))
    hist.SetLineWidth(3)
    hist.SetMarkerStyle(0)
    yield_hist = hist.Integral()
    return hist, yield_hist

def GetBDTBinary(file):
    with open(file) as f:
        binaries = f.readlines() #//TODO: still need to check the store pattern in log file from ML.
        binaries = binaries[0].split(" ")
    return binaries