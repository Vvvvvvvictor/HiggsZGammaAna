import ROOT
from tqdm import trange
import uproot, uproot3
import awkward
import os
import re

def read_file(file, tree = "inclusive", selections = []):
    '''
    read file and apply selections
    selections: a list of selection(<str>). If there is '&' or '|' in one item, MUST use "()" to enclose each subitem.(exp. '(H_eta>1.5) & (H_pt>10)')
    '''
    print("Reading {}:{}...".format(file, tree))
    decro_sel = []
    for i in selections:
        var_can = list(set([j for j in re.split('\W+', i) if "_" in j]))
        for j in var_can:
            i = i.replace(j, "arrays."+j)
        decro_sel.append(i)
    arrays = uproot.open(file+':'+tree).arrays(library='ak')
    for i in decro_sel:
        arrays = arrays[eval(i)]
    return arrays

