#############################
from ROOT import *

#############################

from pdb import set_trace

import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: background Fit")
parser.add_argument("-d", "--data", dest="data", type=str, default="./two_jet/data.root", help="path to the input dataset")
parser.add_argument("-j", "--json", dest="json", type=str, default="./significances/bin_binaries_two_jet.txt", help="path to the input json file")
parser.add_argument("-o", "--out", dest="out", type=str, default="./output_file", help="path to the output dir")
parser.add_argument("-t", "--tree", dest="tree", type=str, default="test", help="name of the tree")
parser.add_argument("-mb", "--multibdt", dest="multibdt", type=bool, default=False, help="whether use two BDT for categorization")
args = parser.parse_args()

# Creat output dir
############################
import os
if not os.path.exists(args.out):
    os.makedirs(args.out) 


# Get the number of categories and the corresponding boundaries
############################
file_json = open(args.json, "r")
json = file_json.readlines()
num_bdt = int(json[0])
if num_bdt == 2:
    bdt_div = float(json[1])
nCat, boundaries = [], []
for i in range(num_bdt):
    sub_json = json[i+2].split(" ")
    nCat.append(int(sub_json[0]))
    boundaries.append(list(map(float, sub_json[1:nCat[0]+2])))

print nCat, "\n", boundaries

# Read the input dataset
############################
myfile = TFile(args.data)
mychain = myfile.Get(args.tree)
entries = mychain.GetEntriesFast()
print entries

for c in range(nCat[0]):
    w = RooWorkspace("CMS_hzg_workspace")

    Sqrts = RooRealVar("Sqrts","Sqrts",13)
    CMS_hzg_mass = RooRealVar("CMS_hzg_mass","CMS_hzg_mass",125.,100.,180.)

    getattr(w,'import')(Sqrts)
    getattr(w,'import')(CMS_hzg_mass)

    data_mass_cats = RooDataSet("data_mass_cat0","data_mass_cat0",RooArgSet(CMS_hzg_mass))

    for jentry in range(entries):
        nb = mychain.GetEntry(jentry)
        if mychain.H_mass<100. or mychain.H_mass>180.: continue
        if args.multibdt: 
            if mychain.bdt_score_VBF > bdt_div: continue

        CMS_hzg_mass.setVal(mychain.H_mass)

        if mychain.bdt_score_t > boundaries[0][c] and mychain.bdt_score_t < boundaries[0][c+1]:
            data_mass_cats.add(RooArgSet(CMS_hzg_mass))

    #getattr(w,'import')(data_mass_cat0)
    getattr(w,'import')(data_mass_cats)

    w.writeToFile(args.out+"/HZGamma_data_bkg_workspace_cat{0}.root".format(c))
    del w

if args.multibdt:
    for c in range(nCat[1]):
        w = RooWorkspace("CMS_hzg_workspace")

        Sqrts = RooRealVar("Sqrts","Sqrts",13)
        CMS_hzg_mass = RooRealVar("CMS_hzg_mass","CMS_hzg_mass",125.,100.,180.)

        getattr(w,'import')(Sqrts)
        getattr(w,'import')(CMS_hzg_mass)

        data_mass_cats = RooDataSet("data_mass_cat0","data_mass_cat0",RooArgSet(CMS_hzg_mass))

        for jentry in range(entries):
            nb = mychain.GetEntry(jentry)
            if mychain.H_mass<100. or mychain.H_mass>180.: continue
            if mychain.bdt_score_VBF < bdt_div: continue

            CMS_hzg_mass.setVal(mychain.H_mass)

            if mychain.bdt_score_t > boundaries[1][c] and mychain.bdt_score_t < boundaries[1][c+1]:
                data_mass_cats.add(RooArgSet(CMS_hzg_mass))

        #getattr(w,'import')(data_mass_cat0)
        getattr(w,'import')(data_mass_cats)

        w.writeToFile(args.out+"/HZGamma_data_bkg_workspace_cat{0}.root".format(c+nCat[0]))
        del w
