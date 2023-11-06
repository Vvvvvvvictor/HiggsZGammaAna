
#############################
from ROOT import *

#############################
import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: signal Fit")
parser.add_argument("-d", "--data", dest="data", type=str, default="./two_jet/data.root", help="path to the input dataset")
parser.add_argument("-j", "--json", dest="json", type=str, default="./significances/bin_binaries_two_jet.txt", help="path to the input json file")
parser.add_argument("-o", "--out", dest="out", type=str, default="./output_file", help="path to the output dir")
parser.add_argument("-t", "--tree", dest="tree", type=str, default="test", help="name of the tree")
parser.add_argument("-dm", "--shift", dest="shift", type=float, default=5, help="shift ALP mass")
args = parser.parse_args()

shift = int(args.shift)

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
lumi = {'16':16.81, '16APV':19.52, '17':41.48, '18':59.83}
myfile = TFile(args.data)
if myfile:
    print "open " + args.data + " success!"
else:
    print args.data + " does not exist!"
    sys.exit(0)
mychain = myfile.Get(args.tree)
entries = mychain.GetEntriesFast()
print "Entries: ",entries

for c in range(nCat):
    for mass_H in [125 - shift, 125, 125 + shift]:
        print "[[ INFO ]] prepare Higgs mass: {0}".format(mass_H)
        w = RooWorkspace("CMS_hzg_workspace")

        Sqrts = RooRealVar("Sqrts","Sqrts",13)
        IntLumi = RooRealVar("IntLumi","IntLumi", lumi['17'])
        CMS_hzg_mass = RooRealVar("CMS_hzg_mass","CMS_hzg_mass",mass_H,105.,170.)
        CMS_hzg_weight = RooRealVar("CMS_hzg_weight","CMS_hzg_weight",-100000,1000000)


        getattr(w,'import')(Sqrts)
        getattr(w,'import')(IntLumi)
        getattr(w,'import')(CMS_hzg_mass)
        getattr(w,'import')(CMS_hzg_weight)

        ArgSet = RooArgSet("args")
        ArgSet.add(CMS_hzg_mass)
        ArgSet.add(CMS_hzg_weight)
        print "#"*51
        ArgSet.Print("v")
        print "#"*51

        data_mass_cats = RooDataSet("ggh_{0}_13TeV_cat0".format(mass_H),"ggh_{0}_13TeV_cat0".format(mass_H), ArgSet, "CMS_hzg_weight")

        for jentry in range(entries):
            nb = mychain.GetEntry(jentry)
            if mychain.H_mass>170. or mychain.H_mass<105.: continue


            CMS_hzg_mass.setVal(mychain.H_mass + mass_H - 125.0)
            CMS_hzg_weight.setVal(mychain.weight)

            if mychain.bdt_score_t > boundaries[c] and mychain.bdt_score_t < boundaries[c+1]:
                data_mass_cats.add(ArgSet,mychain.weight)

        data_mass_cats.Print("v")
        #dataset_WithoutWeight.Print("v")
        print type(ArgSet)
        print "#"*51

        
        c1 = TCanvas("c1","With Weight")
        massDist = CMS_hzg_mass.frame(RooFit.Title("CMS_hzg_mss"))
        data_mass_cats.plotOn(massDist)
        massDist.Draw()

        getattr(w,'import')(data_mass_cats)

        w.writeToFile(args.out+"/HZGamma_data_sig_Hm{0}_workspace_cat{1}.root".format(mass_H,c))
        del w

