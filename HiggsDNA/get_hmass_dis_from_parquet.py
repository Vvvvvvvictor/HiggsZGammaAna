import pandas as pd
from tqdm import trange
import re
import ROOT

def ReadFile(file, selections = []):
    '''
    read file and apply selections
    selections: a list of selection(<str>). If there is '&' or '|' in one item, MUST use "()" to enclose each subitem.(exp. '(H_eta>1.5) & (H_pt>10)')
    '''
    print("Reading {}...".format(file))
    variabels = ["weight_central"]
    decro_sel = []
    for i in selections:
        var_can = list(set([j for j in re.split('\W+', i) if "_" in j]))
        variabels += var_can
        for j in var_can:
            i = i.replace(j, "arrays."+j)
        decro_sel.append(i)
    arrays = pd.read_parquet(file, columns=variabels)
    for i in decro_sel:
        arrays = arrays[eval(i)]
    return arrays.reset_index(drop=True)

def AddHist(arrays, variable, hist):
    for i in trange(0, len(arrays[variable])):
        hist.Fill(float(arrays[variable][i]), float(arrays['weight_central'][i]))
    hist.SetLineWidth(3)
    hist.SetMarkerStyle(0)
    yield_hist = hist.Integral()
    return hist, yield_hist

if __name__ == "__main__":
    target_path = "/afs/cern.ch/user/j/jiehan/public/data_H_mass_dis.root"
    input_path = "/eos/home-j/jiehan/parquet/2022/data/Data_E_2022/merged_nominal.parquet"
    
    f = ROOT.TFile(target_path, "RECREATE")
    selections = ["(H_mass<120) | (H_mass>130)"]
    individual_sel = {"ele":["Z_lead_lepton_id==11"], "mu":["Z_lead_lepton_id==13"]}

    for cat in ("ele", "mu"):
        name = "data_{}".format(cat)
        sel = selections + individual_sel[cat]
        data_hist = ROOT.TH1D(name,name, 80, 100, 180)
        data_hist.Sumw2()
        arrays = ReadFile(input_path, selections = sel)
        data_hist, yield_h = AddHist(arrays, "H_mass", data_hist)
        print("Yields for {} channel is: {:.0f}".format(cat, yield_h))
        data_hist.Write(name)