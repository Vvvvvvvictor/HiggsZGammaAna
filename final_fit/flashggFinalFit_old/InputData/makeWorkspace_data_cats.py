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
boundaries = list(map(float, json[1:nCat+2]))



# Read the input dataset
############################
myfile = TFile(args.data)
mychain = myfile.Get(args.tree)
entries = mychain.GetEntries()
print entries

for c in range(nCat):
    print 'Boundaries:',boundaries[c],boundaries[c+1]
    w = RooWorkspace("CMS_hzg_workspace")

    Sqrts = RooRealVar("Sqrts","Sqrts",13)
    CMS_hzg_mass = RooRealVar("CMS_hzg_mass","CMS_hzg_mass",125.,105.,170.)

    getattr(w,'import')(Sqrts)
    getattr(w,'import')(CMS_hzg_mass)

    data_mass_cats = RooDataSet("data_mass_cat0","data_mass_cat0",RooArgSet(CMS_hzg_mass))

    for jentry in range(entries):
        nb = mychain.GetEntry(jentry)
        if mychain.H_mass<105. or mychain.H_mass>170.: continue

        CMS_hzg_mass.setVal(mychain.H_mass)

        #if mychain.BDT_score > boundaries[c] and mychain.BDT_score <= boundaries[c+1]:
        #    data_mass_cats.add(RooArgSet(CMS_hzg_mass))
        if mychain.bdt_score_t > boundaries[c] and mychain.bdt_score_t <= boundaries[c+1]:
            data_mass_cats.add(RooArgSet(CMS_hzg_mass))

    #getattr(w,'import')(data_mass_cat0)
    getattr(w,'import')(data_mass_cats)

    w.writeToFile(args.out+"/HZGamma_data_bkg_workspace_cat{0}.root".format(c))
    del w

