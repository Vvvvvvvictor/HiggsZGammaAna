import ROOT
from tqdm import trange
import uproot
import os
import pandas as pd
import numpy as np
import re
from pdb import set_trace

def read_root_file(file, var=None, tree="inclusive", selections=[]):
    '''
    read file and apply selections
    selections: a list of selection(<str>). If there is '&' or '|' in one item, MUST use "()" to enclose each subitem.(exp. '(H_eta>1.5) & (H_pt>10)')
    '''
    print("Reading {}:{}...".format(file, tree))
    variables = ["weight", "H_mass"]
    if isinstance(var, str):
        variables.append(var)
    else:
        for i in var:
            variables.append(i)
    decro_sel, var_not_exist = [], []
    for i in selections:
        var_can = list(set([j for j in re.split('\W+', i) if not j.isdigit() and j]))
        print(var_can)
        variables += var_can
        for j in var_can:
            i = i.replace(j, "arrays."+j)
        decro_sel.append(i)
    data= uproot.open(file+':'+tree)
    variables = list(set(variables))
    for var in variables:
        if var not in data:
            var_not_exist.append(var)
    variables = [x for x in variables if x not in var_not_exist]
    arrays = data.arrays(variables,library='pd')
    # for i in decro_sel:
    #     print(i)
    #     arrays = arrays[eval(i)]
    for i in selections:
        arrays = arrays.query(i)
    return arrays, var_not_exist

def read_parquet_file(file, var=None, selections=[]):
    '''
    read file and apply selections
    selections: a list of selection(<str>). If there is '&' or '|' in one item, MUST use "()" to enclose each subitem.(exp. '(H_eta>1.5) & (H_pt>10)')
    '''
    print("Reading {}...".format(file))
    variables = ["weight_central"]
    if isinstance(var, str):
        variables.append(var)
    else:
        for i in var:
            variables.append(i)
    decro_sel, var_not_exist = [], []
    for i in selections:
        var_can = list(set([j for j in re.split('\W+', i) if not j.isdigit() and j]))
        print(var_can)
        variables += var_can
    data = pd.read_parquet(file)
    variables = list(set(variables))
    for var in variables:
        if var not in data:
            var_not_exist.append(var)
    variables = [x for x in variables if x not in var_not_exist]
    arrays = data[variables]
    for i in selections:
        arrays = arrays.query(i)
    return arrays, var_not_exist

def get_hist(arrays, variable, ratio, name, bins, range, hist=None, selections=[]):
    if hist is None:
        hist = ROOT.TH1D(name,name,bins,range[0],range[1])
        hist.Sumw2()
    yields = hist.Integral()
    for i in selections:
        var_can = list(set([j for j in re.split('\W+', i) if not j.isdigit()]))
        for j in var_can:
            i = i.replace(j, "arrays."+j)
        print(i)
        arrays = arrays[eval(i)]
    arrays = arrays.reset_index(drop=True)
    for i in trange(0, len(arrays[variable])):
        hist.Fill(float(arrays[variable][i]), float(arrays['weight'][i])*ratio)
    hist.SetLineWidth(3)
    hist.SetMarkerStyle(0)
    yield_inte = hist.Integral()
    yield_sing = (yield_inte - yields)# / ratio
    print(yield_sing)
    return hist, yield_sing, yield_inte

def get_hist_sb(arrays, variable, ratio, name, bins, range, blind_range, hist=None):
    if hist is None:
        hist = ROOT.TH1D(name,name,bins,range[0],range[1])
        hist.Sumw2()
    yields = hist.Integral()
    for i in trange(0, len(arrays[variable])):
        if (arrays["H_mass"][i]>blind_range[1] and arrays["H_mass"][i]<blind_range[0]):
            continue
        hist.Fill(float(arrays[variable][i]), float(arrays['weight'][i])*ratio)
    hist.SetLineWidth(3)
    hist.SetMarkerStyle(0)
    yield_inte = hist.Integral()
    yield_sing = (yield_inte - yields)# / ratio
    print(yield_sing)
    return hist, yield_sing, yield_inte

def get_ratio_hist(numerator_h, denominator_h, range=None):
    ratio_h = numerator_h.Clone("ratio_h")
    ratio_h.SetLineColor(ROOT.kBlack)
    ratio_h.SetMarkerStyle(0)
    ratio_h.SetMarkerSize(1)
    ratio_h.SetTitle("")
    if range:
        ratio_h.SetMinimum(range[0])
        ratio_h.SetMaximum(range[1])
    
    # Set up plot for markers and errors
    # ratio_h.Sumw2()
    ratio_h.SetStats(0)
    ratio_h.Divide(denominator_h)

    # Adjust y-axis settings
    y = ratio_h.GetYaxis()
    y.SetTitle("")
    y.SetNdivisions(505)
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

def get_S_over_sqrtB(hist_s, hist_b, ratio, yrange=None):
    hist_sosb = hist_s.Clone("hist_sosb")
    hist_sosb.Reset()
    hist_sosb.SetLineColor(ROOT.kBlack)
    hist_sosb.SetMarkerStyle(0)
    hist_sosb.SetMarkerSize(1)
    hist_sosb.SetTitle("")
    if yrange:
        hist_sosb.SetMinimum(yrange[0])
        hist_sosb.SetMaximum(yrange[1])

    for i in range(hist_s.GetNbinsX()):
        s_bin_content = hist_s.GetBinContent(i) / ratio
        b_bin_content = hist_b.GetBinContent(i)
        s_bin_content_err = hist_s.GetBinError(i)
        b_bin_content_err = hist_b.GetBinError(i)
        # s_bin_content_err = s_bin_content ** 0.5
        # b_bin_content_err = b_bin_content ** 0.5
        # print(s_bin_content, b_bin_content, sep=", ")
        
        maximum = 0
        # Avoid division by zero
        if b_bin_content > 0 and b_bin_content_err / b_bin_content < 5:
            ssqrtoverb = s_bin_content / (b_bin_content ** 0.5)
            if ssqrtoverb > maximum:
                maximum = ssqrtoverb
            ssqrtoverb_err = ((s_bin_content_err / b_bin_content ** 0.5) ** 2 + (s_bin_content * b_bin_content_err / b_bin_content ** 1.5 / 2) ** 2) ** 0.5
            hist_sosb.SetBinContent(i, ssqrtoverb)
            hist_sosb.SetBinError(i, ssqrtoverb_err)
    
    # hist_sosb.SetMaximum(1.1*maximum)

    # Adjust y-axis settings
    y = hist_sosb.GetYaxis()
    y.SetTitle("")
    y.SetNdivisions(505)
    y.SetTitleSize(20)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(20)

    # Adjust x-axis settings
    x = hist_sosb.GetXaxis()
    x.SetTitle("")
    x.SetTitleSize(30)
    x.SetTitleFont(43)
    x.SetTitleOffset(1.00)
    x.SetLabelFont(43)
    x.SetLabelSize(20)
    
    return hist_sosb
