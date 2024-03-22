#############################
from ROOT import *

#############################

import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: background Fit")
parser.add_argument("-d", "--data", dest="data", type=str, default="./samples/DY_deathvalley_v3_untagged.dat, ./samples/SMZg_deathvalley_v3_untagged.dat", help="path to the input dataset")
parser.add_argument("-j", "--json", dest="json", type=str, default="./bin_binaries_two_jet.txt", help="path to the input json file")
parser.add_argument("-o", "--out", dest="out", type=str, default="./output_file", help="path to the output dir")
parser.add_argument("-t", "--tree", dest="tree", type=str, default="test", help="name of the tree")
args = parser.parse_args()

# Creat output dir
############################
import os
if not os.path.exists(args.out):
    os.makedirs(args.out) 


# Get the number of categories and the corresponding boundaries
############################
file_json = open(args.json, "r")
json = file_json.readline().split(' ')
nCat = int(json[0])
boundaries = list(map(float, json[1:6]))

print boundaries



# Read the input dataset
############################

var = ["mllg", "photon_pT", "BDT", "weight", "year", "lep_id", "pho_eta", "Nlep", "Njet"]

dataset = {}
for v in var:
    dataset[v] = []

def read_dat_file(files_path):
    for file_path in files_path:
        file_path = file_path.lstrip()
        with open(file_path, 'r') as file:
            for line in file:
                data = line.strip().split()
                for i in range(len(data)):
                    dataset[var[i]].append(float(data[i]))

    return dataset

datasets = args.data.split(',')

data = read_dat_file(datasets)

print len(data["mllg"])


for c in range(nCat):
    w = RooWorkspace("CMS_hzg_workspace")

    Sqrts = RooRealVar("Sqrts","Sqrts",13)
    CMS_hzg_mass = RooRealVar("CMS_hzg_mass","CMS_hzg_mass",125.,105.,170.)
    CMS_hzg_weight = RooRealVar("CMS_hzg_weight","CMS_hzg_weight",-100000,1000000)

    getattr(w,'import')(Sqrts)
    getattr(w,'import')(CMS_hzg_mass)
    getattr(w,'import')(CMS_hzg_weight)

    ArgSet = RooArgSet("args")
    ArgSet.add(CMS_hzg_mass)
    ArgSet.add(CMS_hzg_weight)

    data_mass_cats = RooDataSet("data_mass_cat0","data_mass_cat0",ArgSet, "CMS_hzg_weight")
    

    for jentry in range(len(data[var[0]])):
        if data["mllg"][jentry]<105. or data["mllg"][jentry]>170.: continue
        if data["photon_pT"][jentry]<10. or data["photon_pT"][jentry]>80.: continue
        if data["Nlep"][jentry]>2: continue
        if data["Njet"][jentry]>1: continue

        CMS_hzg_mass.setVal(data["mllg"][jentry])
        CMS_hzg_weight.setVal(data["weight"][jentry])

        if data["BDT"][jentry] > boundaries[c] and data["BDT"][jentry] < boundaries[c+1]:
            data_mass_cats.add(ArgSet, data["weight"][jentry])


    data_mass_cats.Print("v")
    c1 = TCanvas("c1","With Weight")
    massDist = CMS_hzg_mass.frame(RooFit.Title("CMS_hzg_mss"))
    data_mass_cats.plotOn(massDist)
    massDist.Draw()

    getattr(w,'import')(data_mass_cats)

    w.writeToFile(args.out+"/HZGamma_data_bkg_workspace_cat{0}.root".format(c))
    del w
