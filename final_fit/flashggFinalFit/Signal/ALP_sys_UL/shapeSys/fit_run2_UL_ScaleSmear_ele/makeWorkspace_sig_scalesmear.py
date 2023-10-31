
#############################
from ROOT import *

from xgboost import XGBClassifier
import pickle


#############################
import argparse
parser = argparse.ArgumentParser(description="Prepare flashggFinalFit workspace: signal Fit")
parser.add_argument("-m", "--mass", dest="mass", type=float, default=1.0, help="ALP mass")
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
args = parser.parse_args()

mass = int(args.mass)

if args.ele:
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param_ele.pkl"
    mvaCuts = {1:0.965, 2:0.98, 3:0.98, 4:0.975, 5:0.975, 6:0.955, 7:0.97, 8:0.975, 9:0.98, 10:0.98, 15:0.98, 20:0.985, 25:0.985, 30:0.98}
elif args.mu:
    BDT_filename="/publicfs/cms/user/wangzebing/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param_mu.pkl"
    mvaCuts = {1:0.95, 2:0.975, 3:0.97, 4:0.98, 5:0.985, 6:0.985, 7:0.985, 8:0.985, 9:0.985, 10:0.99, 15:0.99, 20:0.99, 25:0.985, 30:0.98}
else:
    print "channel error"
    exit(0)

mvaCut = mvaCuts[mass]
model = pickle.load(open(BDT_filename, 'rb'))

####################
lumi = {'16':16.81, '16APV':19.52, '17':41.48, '18':59.83}
for year in ['16', '16APV','17','18']:

    print "prepare year: "+year

    filename = '/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/'+year+'/ALP_M{0}.root'.format(mass)
    myfile = TFile(filename)
    if myfile:
        print "open " + filename + " success!"
    else:
        print filename + " does not exist!"
        sys.exit(0)
    mychain = myfile.Get('passedEvents')
    entries = mychain.GetEntriesFast()
    print "Entries: ",entries


    mass_H = 125

    print "[[ INFO ]] prepare ALP mass: {0}".format(mass_H)
    w = RooWorkspace("CMS_hza_workspace")

    Sqrts = RooRealVar("Sqrts","Sqrts",13)
    IntLumi = RooRealVar("IntLumi","IntLumi", lumi[year])
    CMS_hza_mass = RooRealVar("CMS_hza_mass","CMS_hza_mass",mass_H,110.,180.)
    CMS_hza_weight = RooRealVar("CMS_hza_weight","CMS_hza_weight",-100000,1000000)


    getattr(w,'import')(Sqrts)
    getattr(w,'import')(IntLumi)
    getattr(w,'import')(CMS_hza_mass)
    getattr(w,'import')(CMS_hza_weight)

    ArgSet = RooArgSet("args")
    ArgSet.add(CMS_hza_mass)
    ArgSet.add(CMS_hza_weight)
    print "#"*51
    ArgSet.Print("v")
    print "#"*51

    sysList = ['normal', 'phoScale_up', 'phoScale_dn', 'phoSmear_up', 'phoSmear_dn', 'lepScale_up', 'lepScale_dn', 'lepSmear_up', 'lepSmear_dn']
    for sys in sysList:

        dataset = RooDataSet("ggh_{0}_13TeV_cat0_{1}".format(mass_H,sys),"ggh_{0}_13TeV_cat0_{1}".format(mass_H,sys), ArgSet, "CMS_hza_weight")

        for jentry in range(entries):
            nb = mychain.GetEntry(jentry)
            if not mychain.passChaHadIso: continue
            if not mychain.passNeuHadIso: continue
            if not mychain.passdR_gl: continue
            if not mychain.passHOverE: continue
            if mychain.H_m>180. or mychain.H_m<110.: continue

            if args.ele:
                if abs(mychain.l1_id) == 13: 
                    continue
            if args.mu:
                if abs(mychain.l1_id) == 11: 
                    continue

            param = (mychain.ALP_m - mass)/mychain.H_m
            MVA_list = [mychain.pho1Pt, mychain.pho1R9, mychain.pho1IetaIeta55, mychain.pho1PIso_noCorr ,mychain.pho2Pt, mychain.pho2R9, mychain.pho2IetaIeta55,mychain.pho2PIso_noCorr,mychain.ALP_calculatedPhotonIso, mychain.var_dR_Za, mychain.var_dR_g1g2, mychain.var_dR_g1Z, mychain.var_PtaOverMh, mychain.H_pt, param]
            MVA_value = model.predict_proba(MVA_list)[:, 1]
            if MVA_value < mvaCut:continue

            pho1 = TLorentzVector()
            pho2 = TLorentzVector()
            lep1 = TLorentzVector()
            lep2 = TLorentzVector()
            pho1.SetPtEtaPhiM(mychain.pho1Pt, mychain.pho1eta, mychain.pho1phi, 0.0)
            pho2.SetPtEtaPhiM(mychain.pho2Pt, mychain.pho2eta, mychain.pho2phi, 0.0)
            lep1.SetPtEtaPhiM(mychain.l1_pt, mychain.l1_eta, mychain.l1_phi, mychain.l1_mass)
            lep2.SetPtEtaPhiM(mychain.l2_pt, mychain.l2_eta, mychain.l2_phi, mychain.l2_mass)

            
            pho1_corr = TLorentzVector()
            pho2_corr = TLorentzVector()
            lep1_corr = TLorentzVector()
            lep2_corr = TLorentzVector()
            
            pho1_corrE = {'phoScale_up':mychain.pho1scaleup, 'phoScale_dn':mychain.pho1scaledn, 'phoSmear_up':mychain.pho1smearup, 'phoSmear_dn':mychain.pho1smeardn}
            pho2_corrE = {'phoScale_up':mychain.pho2scaleup, 'phoScale_dn':mychain.pho2scaledn, 'phoSmear_up':mychain.pho2smearup, 'phoSmear_dn':mychain.pho2smeardn}
            
            if abs(mychain.l1_id) == 11:
                lep1_corrE = {'lepScale_up':mychain.l1_scaleup, 'lepScale_dn':mychain.l1_scaledn, 'lepSmear_up':mychain.l1_smearup, 'lepSmear_dn':mychain.l1_smeardn}
                lep2_corrE = {'lepScale_up':mychain.l2_scaleup, 'lepScale_dn':mychain.l2_scaledn, 'lepSmear_up':mychain.l2_smearup, 'lepSmear_dn':mychain.l2_smeardn}
            else:
                lep1_corrE = {'lepScale_up':mychain.l1_scaleup*lep1.E(), 'lepScale_dn':mychain.l1_scaledn*lep1.E(), 'lepSmear_up':mychain.l1_smearup*lep1.E(), 'lepSmear_dn':mychain.l1_smeardn*lep1.E()}
                lep2_corrE = {'lepScale_up':mychain.l2_scaleup*lep2.E(), 'lepScale_dn':mychain.l2_scaledn*lep2.E(), 'lepSmear_up':mychain.l2_smearup*lep2.E(), 'lepSmear_dn':mychain.l2_smeardn*lep2.E()}

            if sys in pho1_corrE.keys():
                pho1_f = pho1_corrE[sys]/mychain.pho1EPostCorr
                pho2_f = pho2_corrE[sys]/mychain.pho2EPostCorr
                pho1_corr.SetPxPyPzE(pho1.Px()*pho1_f, pho1.Py()*pho1_f, pho1.Pz()*pho1_f, pho1_corrE[sys])
                pho2_corr.SetPxPyPzE(pho2.Px()*pho2_f, pho2.Py()*pho2_f, pho2.Pz()*pho2_f, pho2_corrE[sys])

                H_m = (pho1_corr+pho2_corr+lep1+lep2).M()
            elif sys in lep1_corrE.keys():
                lep1_f = lep1_corrE[sys]/lep1.E()
                lep2_f = lep2_corrE[sys]/lep2.E()
                lep1_corr.SetPxPyPzE(lep1.Px()*lep1_f, lep1.Py()*lep1_f, lep1.Pz()*lep1_f, lep1_corrE[sys])
                lep2_corr.SetPxPyPzE(lep2.Px()*lep2_f, lep2.Py()*lep2_f, lep2.Pz()*lep2_f, lep2_corrE[sys])
                
                H_m = (pho1+pho2+lep1_corr+lep2_corr).M()
            else:
                H_m = mychain.H_m
                


            CMS_hza_mass.setVal(H_m)
            CMS_hza_weight.setVal(mychain.factor*mychain.pho1SFs*mychain.pho2SFs)
            dataset.add(ArgSet,mychain.factor*mychain.pho1SFs*mychain.pho2SFs)
            #dataset_WithoutWeight.add(ArgSet)

        dataset.Print("v")
        #dataset_WithoutWeight.Print("v")
        print type(ArgSet)
        print "#"*51

        c1 = TCanvas("c1","With Weight")
        massDist = CMS_hza_mass.frame(RooFit.Title("CMS_hza_mss"))
        dataset.plotOn(massDist)
        massDist.Draw()
        #c1.SaveAs("test_"+sys+"_"+year+".png")
        '''
        c2 = TCanvas("c2","Without Weight")
        massDist2 = CMS_hza_mass.frame(RooFit.Title("CMS_hza_mss"))
        dataset_WithoutWeight.plotOn(massDist2)
        massDist2.Draw()
        '''
        getattr(w,'import')(dataset)
        #getattr(w,'import')(dataset_WithoutWeight)

    w.writeToFile("ALP_sig_Am{0}_Hm{1}_{2}_workspace.root".format(mass,mass_H,year))
    del w

#raw_input()
