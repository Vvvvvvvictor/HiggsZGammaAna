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

def get_hist(arrays, variable, name, bins, range, hist=None):
    if hist is None:
        hist = ROOT.TH1D(name,name,bins,range[0],range[1])
        hist.Sumw2()
    for i in trange(0, len(arrays[variable])):
        hist.Fill(float(arrays[variable][i]), float(arrays['weight'][i]))
    hist.SetLineWidth(3)
    hist.SetMarkerStyle(0)
    yield_hist = hist.Integral()
    return hist, yield_hist

def get_ratio_hist(numerator_h, denominator_h):
    ratio_h = denominator_h.Clone("ratio_h")
    ratio_h.SetLineColor(ROOT.kBlack)
    ratio_h.SetMarkerStyle(0)
    ratio_h.SetMarkerSize(1)
    ratio_h.SetTitle("")
    ratio_h.SetMinimum(0.35)
    ratio_h.SetMaximum(1.65)
    
    # Set up plot for markers and errors
    # ratio_h.Sumw2()
    ratio_h.SetStats(0)
    ratio_h.Divide(numerator_h)

    # Adjust y-axis settings
    y = ratio_h.GetYaxis()
    y.SetTitle("")
    y.SetNdivisions(105)
    y.SetTitleSize(20)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(20)

    # Adjust x-axis settings
    x = ratio_h.GetXaxis()
    x.SetTitle("")
    x.SetTitleSize(30)
    x.SetTitleFont(43)
    x.SetTitleOffset(1.00)
    x.SetLabelFont(43)
    x.SetLabelSize(20)

    return ratio_h
