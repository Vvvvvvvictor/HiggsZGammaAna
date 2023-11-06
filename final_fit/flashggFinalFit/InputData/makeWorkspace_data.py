#############################
from ROOT import *

#############################

import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: background Fit")
parser.add_argument("-d", "--data", dest="data", type=str, default="./two_jet/data.root", help="path to the input dataset")
parser.add_argument("-j", "--json", dest="json", type=str, default="./significances/bin_binaries_two_jet.txt", help="path to the input json file")
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



# Read the input dataset
############################
myfile = TFile(args.data)
mychain = myfile.Get(args.tree)
entries = mychain.GetEntriesFast()
print entries


w = RooWorkspace("CMS_hzg_workspace")

Sqrts = RooRealVar("Sqrts","Sqrts",13)
CMS_hzg_mass = RooRealVar("CMS_hzg_mass","CMS_hzg_mass",125.,105.,170.)

getattr(w,'import')(Sqrts)
getattr(w,'import')(CMS_hzg_mass)

data_mass_cats = {}
for c in range(nCat):
    data_mass_cats[c] = RooDataSet("data_mass_cat"+str(c),"data_mass_cat"+str(c),RooArgSet(CMS_hzg_mass))

for jentry in range(entries):
    nb = mychain.GetEntry(jentry)
    if mychain.H_mass<105. or mychain.H_mass>170.: continue

    CMS_hzg_mass.setVal(mychain.H_mass)

    if mychain.bdt_score_t > boundaries[0] and mychain.bdt_score_t < boundaries[1]:
        data_mass_cats[3].add(RooArgSet(CMS_hzg_mass))
    elif mychain.bdt_score_t > boundaries[1] and mychain.bdt_score_t < boundaries[2]:
        data_mass_cats[2].add(RooArgSet(CMS_hzg_mass))
    elif mychain.bdt_score_t > boundaries[2] and mychain.bdt_score_t < boundaries[3]:
        data_mass_cats[1].add(RooArgSet(CMS_hzg_mass))
    elif mychain.bdt_score_t > boundaries[3] and mychain.bdt_score_t < boundaries[4]:
        data_mass_cats[0].add(RooArgSet(CMS_hzg_mass))


#getattr(w,'import')(data_mass_cat0)
for c in range(nCat):
    getattr(w,'import')(data_mass_cats[c])


w.writeToFile(args.out+"/HZGamma_data_bkg_workspace.root")
del w

