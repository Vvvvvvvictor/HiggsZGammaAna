#############################
from ROOT import *

from xgboost import XGBClassifier
import pickle

#############################

import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: signal Fit")
parser.add_argument("-m", "--mass", dest="mass", type=float, default=1.0, help="ALP mass")
parser.add_argument('--interp', dest='interp', action='store_true', default=False, help='interpolation?')
args = parser.parse_args()

mass = int(args.mass)

BDT_filename="/publicfs/cms/user/laipeizhu/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
if args.interp:
    mvaCuts = {11:0.99, 12:0.99, 13:0.99, 14:0.99, 16:0.99, 17:0.99, 18:0.99, 19:0.99, 21:0.99, 22:0.99, 23:0.985, 24:0.985, 26:0.985, 27:0.985, 28:0.98, 29:0.98}
else:
    mvaCuts = {1:0.945, 2:0.975, 3:0.985, 4:0.985, 5:0.99, 6:0.99, 7:0.99, 8:0.995, 9:0.995, 10:0.995, 15:0.99, 20:0.99, 25:0.99, 30:0.99}
    # mvaCuts = {1:0.955, 2:0.98, 3:0.985, 4:0.98, 5:0.985, 6:0.99, 7:0.985, 8:0.99, 9:0.99, 10:0.99, 15:0.99, 20:0.99, 25:0.985, 30:0.98}

mvaCut = mvaCuts[mass]
model = pickle.load(open(BDT_filename, 'rb'))

############################
myfile = TFile('/publicfs/cms/user/laipeizhu/ALP/Analysis_output/UL/run2/ALP_data.root')
mychain = myfile.Get('passedEvents')
entries = mychain.GetEntriesFast()

print (entries)


w = RooWorkspace("CMS_hza_workspace")

Sqrts = RooRealVar("Sqrts","Sqrts",13)
CMS_hza_mass = RooRealVar("CMS_hza_mass","CMS_hza_mass",125.,95.,180.)

getattr(w,'import')(Sqrts)
getattr(w,'import')(CMS_hza_mass)


data_mass_cat0 = RooDataSet("data_mass_cat0","data_mass_cat0",RooArgSet(CMS_hza_mass))

# ROOT.gROOT.SetOwnership(data_mass_cat0, False)

for jentry in range(entries):
    nb = mychain.GetEntry(jentry)
    if not mychain.passChaHadIso: continue
    if not mychain.passNeuHadIso: continue
    if not mychain.passdR_gl: continue
    if not mychain.passHOverE: continue
    if mychain.H_m<95. or mychain.H_m>180.: continue


    param = (mychain.ALP_m - float(mass))/mychain.H_m
    MVA_list = [mychain.pho1Pt, mychain.pho1R9, mychain.pho1IetaIeta55, mychain.pho1PIso_noCorr ,mychain.pho2Pt, mychain.pho2R9, mychain.pho2IetaIeta55,mychain.pho2PIso_noCorr,mychain.ALP_calculatedPhotonIso, mychain.var_dR_Za, mychain.var_dR_g1g2, mychain.var_dR_g1Z, mychain.var_PtaOverMh, mychain.H_pt, param]
    MVA_value = model.predict_proba([MVA_list])[:, 1]
    if MVA_value < mvaCut:continue

    CMS_hza_mass.setVal(mychain.H_m)
    data_mass_cat0.add(RooArgSet(CMS_hza_mass))

# c1 = TCanvas("c1","With Weight")
# c1.cd()
# massDist = CMS_hza_mass.frame(RooFit.Title("CMS_hza_mass"))
# data_mass_cat0.plotOn(massDist)
# massDist.Draw()
# c1.SaveAs("ALP_data_bkg_Am{0}_workspace.png".format(mass))

getattr(w,'import')(data_mass_cat0)
# getattr(w, 'import')(data_mass_cat0, ROOT.RooFit.RecycleConflictNodes())

w.writeToFile("./output/data/ALP_data_bkg_Am{0}_workspace.root".format(mass))

del w
