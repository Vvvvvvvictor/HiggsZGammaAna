import os
import sys
sys.path.insert(0, '%s/lib' % os.getcwd())
import ZConstrain_Helper as Zrefit
import Analyzer_Configs as AC
from Plot_Helper import LoadNtuples
import math

from ROOT import *
import ROOT
import tdrstyle
gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

import argparse
parser = argparse.ArgumentParser(description="Plotter for HZGamma analysis")
parser.add_argument("-i", "--input", dest="input", default="./Data/2017/ZReFit/data_genZ/ggH.root", help="Path to the input Gen datasets")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
parser.add_argument("-t", "--tree", dest="tree", default="passedEvents", help="Tree name")
parser.add_argument("-c", "--channel", dest="channel", default="", help="ele/mu channel?")
parser.add_argument('--single_check', dest='single_check', action='store_true', default=False, help='single_check refit')
parser.add_argument('--ievent', type=int, dest='ievent', default=1, help="Event id for single_check")
parser.add_argument('--trueLineshapeFit', dest='trueLineshapeFit', action='store_true', default=False, help='Fit Gen Z true lineshape')
parser.add_argument("--genFile", dest="genFile", default="./Data/2017/ZReFit/data_genZ/ggH.root", help="Path of the input Gen datasets")
parser.add_argument('--binningFit', dest='binningFit', action='store_true', default=False, help='Fit with nBins')
parser.add_argument("--nBins", dest="nBins", default=100, help="number of bins")
parser.add_argument('--doReFit', dest='doReFit', action='store_true', default=False, help='make the Z constrain refit')
parser.add_argument("--sampleName", dest="sampleName", default="ggH", help="Name of the refitted sample")
parser.add_argument('--splitJobs', dest='splitJobs', action='store_true', default=False, help='Split jobs')
parser.add_argument('--nJobs', type=int, dest='nJobs', default=10, help="number of total jobs")
parser.add_argument('--iJob', type=int, dest='iJob', default=1, help="number of the current job")
parser.add_argument('--verbose', dest='verbose', action='store_true', default=False, help='Output the fit informatim')
parser.add_argument("--truelineshape", dest="truelineshape", default="./Data/2017/ZReFit", help="Path to the true lineshape")
parser.add_argument('--plotReFit', dest='plotReFit', action='store_true', default=False, help='plot the Z constrain refit')
parser.add_argument("-o", "--output", dest="output", default="./Data/2017/ZReFit", help="Path of the output")
args = parser.parse_args()

if not args.verbose:
    os.close(1)
    os.open(os.devnull, os.O_WRONLY)

def Model_GenZ_lineshape(mZ):
    global meanCB, sigmaCB, alphaCB, nCB
    global meanGauss1, sigmaGauss1, f1
    global meanGauss2, sigmaGauss2, f2
    global meanGauss3, sigmaGauss3, f3
    global singleCB, gaussShape1, gaussShape2, gaussShape3, CBplusGauss, CBplusGaussplusGauss, CBplusGaussplusGaussplusGauss

    global parameters
    parameters = ["meanCB","sigmaCB","alphaCB","nCB","meanGauss1","meanGauss2","meanGauss3","sigmaGauss1","sigmaGauss2","sigmaGauss3","f1","f2","f3"]

    meanCB = RooRealVar("meanCB","meanCB",90., 60., 120.)
    sigmaCB = RooRealVar("sigmaCB","sigmaCB",4., 0.1, 10.)
    alphaCB = RooRealVar("alphaCB","alphaCB",0.8, 0.01, 5.)
    nCB = RooRealVar("nCB","nCB",9., 0.1, 20.)
    meanGauss1 = RooRealVar("meanGauss1","meanGauss1",91., 60., 120.)
    sigmaGauss1 = RooRealVar("sigmaGauss1","sigmaGauss1",6, 0.5, 20.0)
    f1 = RooRealVar("f1","f1",0.6, 0., 1.)
    meanGauss2 = RooRealVar("meanGauss2","meanGauss2",91., 60., 120.)
    sigmaGauss2 = RooRealVar("sigmaGauss2","sigmaGauss2",2., 0.5, 20.)
    f2 = RooRealVar("f2","f2",0.6, 0., 1.)
    meanGauss3 = RooRealVar("meanGauss3","meanGauss3",91., 60., 120.)
    sigmaGauss3 = RooRealVar("sigmaGauss3","sigmaGauss3",1., 0.5, 20.)
    f3 = RooRealVar("f3","f3",0.8, 0., 1.)

    singleCB = RooCBShape("singleCB", "singleCB", mZ, meanCB, sigmaCB, alphaCB, nCB)
    gaussShape1 = RooGaussian("gaussShape1", "gaussShape1", mZ, meanGauss1, sigmaGauss1)
    CBplusGauss = RooAddPdf("CBplusGauss", "CBplusGauss", RooArgList(singleCB, gaussShape1), f1)
    gaussShape2 = RooGaussian("gaussShape2", "gaussShape2", mZ, meanGauss2, sigmaGauss2)
    CBplusGaussplusGauss = RooAddPdf("CBplusGaussplusGauss", "CBplusGaussplusGauss", RooArgList(CBplusGauss, gaussShape2), f2)
    gaussShape3 = RooGaussian("gaussShape3", "gaussShape3", mZ, meanGauss3, sigmaGauss3)
    CBplusGaussplusGaussplusGauss = RooAddPdf("CBplusGaussplusGaussplusGauss", "CBplusGaussplusGaussplusGauss", RooArgList(CBplusGaussplusGauss, gaussShape3), f3)

    return CBplusGaussplusGaussplusGauss


def get_GenZ_lineshape(filename, treeName, outpath):
    file = TFile(filename)
    tree = file.Get(treeName)

    massZ = RooRealVar("GenHzgLeadGenChild_mass", "GenHzgLeadGenChild_mass", 91., 60., 120.)
    lep_id = RooRealVar("GenHzgLeadGenChildChild1_pdgId", "GenHzgLeadGenChildChild1_pdgId", 11, -1000, 1000)
    massZ.SetTitle("m_{ll}")
    massZ.setUnit("GeV")

    if args.channel == "ele":
        dataset_GenZ = RooDataSet("dataset_GenZ", "dataset_GenZ", {massZ,lep_id}, Import=tree)
        dataset_GenZ = dataset_GenZ.reduce(massZ, "abs(GenHzgLeadGenChildChild1_pdgId)==11")
    elif args.channel == "mu":
        dataset_GenZ = RooDataSet("dataset_GenZ", "dataset_GenZ", {massZ,lep_id}, Import=tree)
        dataset_GenZ = dataset_GenZ.reduce(massZ, "abs(GenHzgLeadGenChildChild1_pdgId)==13")
    else:
        dataset_GenZ = RooDataSet("dataset_GenZ", "dataset_GenZ", massZ, Import=tree)
    
    model = Model_GenZ_lineshape(massZ)

    if args.binningFit:
        nbins = int(args.nBins)
        massZ.setBins(nbins)
        dataset_GenZ_binned = dataset_GenZ.binnedClone()

        fitResult = model.fitTo(dataset_GenZ_binned,RooFit.SumW2Error(kFALSE),RooFit.Save(kTRUE))

        canv = TCanvas("c1","c1",1500,1000)
        frame = massZ.frame(nbins)

        dataset_GenZ_binned.plotOn(frame,RooFit.Binning(nbins))
    else:
        fitResult = model.fitTo(dataset_GenZ,RooFit.SumW2Error(kFALSE),RooFit.Save(kTRUE))

        canv = TCanvas("c1","c1",1500,1000)
        frame = massZ.frame()

        dataset_GenZ.plotOn(frame)

    
    model.plotOn(frame,RooFit.LineColor(kRed))
    model.plotOn(frame,RooFit.Components(singleCB),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlack))
    model.plotOn(frame,RooFit.Components(gaussShape1),RooFit.LineStyle(kDashed),RooFit.LineColor(kGreen))
    model.plotOn(frame,RooFit.Components(gaussShape2),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue))
    model.plotOn(frame,RooFit.Components(gaussShape3),RooFit.LineStyle(kDashed),RooFit.LineColor(kOrange))

    fit_results = {}
    fit_results_err = {}
    file_out = open("{}/GenZ_truelineshape{}.txt".format(outpath,args.channel), 'w')
    for par in parameters:
        fit_results[par] = fitResult.floatParsFinal().find(par).getVal()
        fit_results_err[par] = fitResult.floatParsFinal().find(par).getError()

        file_out.write("{}\t{:.5f}\n".format(par, fit_results[par]))
    
    file_out.close()

    frame.Draw()
    latex1 = TLatex()
    latex1.SetNDC()
    latex1.SetTextSize(0.04)
    latex1.SetTextFont(42)
    latex1.SetTextAlign(23)
    latex_i = 0.
    for par in parameters:
        latex1.DrawLatex(0.75, 0.9-latex_i, "{}: {:.3f} +/- {:.3f}".format(par, fit_results[par], fit_results_err[par]))
        latex_i = latex_i + 0.05

    # get the chi^2/dof
    npar = model.getParameters(dataset_GenZ).selectByAttrib("Constant",kFALSE).getSize()
    chi2ndf = frame.chiSquare(npar)
    #latex1.DrawLatex(0.25, 0.8, "chi^2/DOF = {:.3f}".format(chi2ndf))
    print("[[INFO]] chi^2: {}, DOF: {}".format(chi2ndf, npar))
    

    canv.SaveAs("{}/GenZ_lineshape{}.png".format(outpath,args.channel))
    canv.SaveAs("{}/GenZ_lineshape{}.pdf".format(outpath,args.channel))


def ZConstrainReFit(treeName, outpath, splitJobs):
    #cat = 'zero_to_one_jet'
    #cat = 'two_jet'
    cat = args.tree
    #analyzer_cfg = AC.Analyzer_Config(channel='inclusive', year=args.year, treename=cat)
    analyzer_cfg = AC.Analyzer_Config(channel='ggH', year=args.year, treename=cat)
    analyzer_cfg.Print_Config()

    analyzer_cfg.plot_output_path = "{}/ZRefit_plots".format(outpath)

    ntuples = LoadNtuples(analyzer_cfg)

    truelineshape_ele = Zrefit.getParameters("{}/GenZ_truelineshapeele.txt".format(args.truelineshape))
    #truelineshape_mu = Zrefit.getParameters("{}/GenZ_truelineshapemu.txt".format(args.truelineshape))
    truelineshape_mu = Zrefit.getParameters("{}/GenZ_truelineshapemu_anders.txt".format(args.truelineshape))

    #for sample in analyzer_cfg.samp_names:
    for sample in [args.sampleName]:
        print("[[INFO]] Start the Z contrain refit for sample: {}".format(sample))

        ntup = ntuples[sample] # just a short name
        nEvts = ntup.GetEntries()

        # Split jobs
        analyzer_cfg.out_dir = '{}/{}'.format(analyzer_cfg.sample_loc, sample)
        if splitJobs:
            nJobs = args.nJobs
            iJob = args.iJob
            analyzer_cfg.root_output_name = "{}_run2_refit_job{}.root".format(cat, iJob)

            if iJob < 1:
                print("[[ERROR]]: The job number can only start from 1")
                exit(0)

            splits = [nEvts//nJobs*(iJob-1), nEvts//nJobs*iJob]

            if iJob == nJobs:
                splits[1] = nEvts
            if iJob == 1:
                splits[0] = -1
        else:
            analyzer_cfg.root_output_name = "{}_run2_refit.root".format(cat)
            splits = [0, nEvts]
        
        # copy the output file name
        out_file = analyzer_cfg.Create_output_fiels()
        passedEvents = TTree(treeName,treeName)
        branch_names, branch = Zrefit.branch_setup(passedEvents)
        branch_names_new, branch_new = Zrefit.add_branch(passedEvents)
        if sample == 'data':
            branch_names = [b for b in branch_names if 'Gen' not in b]
            branch_names = [b for b in branch_names if 'n_iso_photons' not in b]

        #print("DEBUG: ", splits)
        print('\n\nOn sample: %s' %sample)
        print('total events: %d' %ntup.GetEntries())

        # loop through all the events
        for iEvt in range( nEvts ):
            ntup.GetEvent(iEvt)
            if ntup.event == args.ievent: 
                print(ntup.event, ntup.Z_lead_lepton_id)

            #if iEvt > 1: break
            #if iEvt > 5000: continue
            if args.single_check:
                if ntup.event != args.ievent: 
                    continue

            if (iEvt % 100000 == 1):
                print("looking at event %d" %iEvt)
            #if (iEvt == 5): break

            if (iEvt <= splits[0] or iEvt > splits[1]):
                continue

            
            for var in branch_names:
                branch[var][0] = getattr(ntup,var)

            l1 = TLorentzVector()
            l2 = TLorentzVector()
            l1_error = TLorentzVector()
            l2_error = TLorentzVector()
            gamma = TLorentzVector()
            gamma_fsr = TLorentzVector()
            Z = TLorentzVector()
            Z_fsr = TLorentzVector()
            Higgs = TLorentzVector()

            l1.SetPtEtaPhiM(ntup.Z_lead_lepton_pt, ntup.Z_lead_lepton_eta, ntup.Z_lead_lepton_phi, ntup.Z_lead_lepton_mass)
            l2.SetPtEtaPhiM(ntup.Z_sublead_lepton_pt, ntup.Z_sublead_lepton_eta, ntup.Z_sublead_lepton_phi, ntup.Z_sublead_lepton_mass)
            gamma.SetPtEtaPhiM(ntup.gamma_pt,ntup.gamma_eta,ntup.gamma_phi,0.0)
            gamma_fsr.SetPtEtaPhiM(ntup.gamma_fsr_pt, ntup.gamma_fsr_eta, ntup.gamma_fsr_phi, 0.0)

            if int(abs(ntup.Z_lead_lepton_id)) == 11:
                l1_f = (l1.E()+ntup.Z_lead_lepton_ptE_error)/l1.E()
                l2_f = (l2.E()+ntup.Z_sublead_lepton_ptE_error)/l2.E()
                l1_error.SetPxPyPzE(l1.Px()*l1_f,l1.Py()*l1_f,l1.Pz()*l1_f,l1.E()*l1_f)
                l2_error.SetPxPyPzE(l2.Px()*l2_f,l2.Py()*l2_f,l2.Pz()*l2_f,l2.E()*l2_f)

                l1_pterr = abs(l1.Pt()-l1_error.Pt())
                l2_pterr = abs(l2.Pt()-l2_error.Pt())
            elif int(abs(ntup.Z_lead_lepton_id)) == 13:
                l1_pterr = abs(ntup.Z_lead_lepton_ptE_error)
                l2_pterr = abs(ntup.Z_sublead_lepton_ptE_error)


            fsr_clean = True
            if ntup.gamma_fsr_pt < 0:
                fsr_clean = False
            if math.sqrt((gamma_fsr.Eta()-gamma.Eta())**2+(gamma_fsr.Phi()-gamma.Phi())**2) < 0.001 or math.sqrt((gamma_fsr.Eta()-l1.Eta())**2+(gamma_fsr.Phi()-l1.Phi())**2) < 0.001 or math.sqrt((gamma_fsr.Eta()-l2.Eta())**2+(gamma_fsr.Phi()-l2.Phi())**2) < 0.001:
                fsr_clean = False

            if fsr_clean:
                n_fsr = 1
                fsr_pterr = Zrefit.gamma_fsr_pterr(gamma_fsr)
            else:
                n_fsr = 0
                fsr_pterr = 0.0

            if int(abs(ntup.Z_lead_lepton_id)) == 11:
                l1_refit_pt,  l2_refit_pt, l1_refit_pt_err, l2_refit_pt_err = Zrefit.MakeZReFitModel(l1.Pt(), l1_pterr, l1.Theta(), l1.Phi(), l1.M(), l2.Pt(), l2_pterr, l2.Theta(), l2.Phi(), l2.M(), gamma_fsr.Pt(), fsr_pterr, gamma_fsr.Theta(), gamma_fsr.Phi(), n_fsr, truelineshape_ele)
            elif int(abs(ntup.Z_lead_lepton_id)) == 13:
                l1_refit_pt,  l2_refit_pt, l1_refit_pt_err, l2_refit_pt_err = Zrefit.MakeZReFitModel(l1.Pt(), l1_pterr, l1.Theta(), l1.Phi(), l1.M(), l2.Pt(), l2_pterr, l2.Theta(), l2.Phi(), l2.M(), gamma_fsr.Pt(), fsr_pterr, gamma_fsr.Theta(), gamma_fsr.Phi(), n_fsr, truelineshape_mu)

            l1_refit = TLorentzVector()
            l2_refit = TLorentzVector()
            Z_refit = TLorentzVector()
            Higgs_refit = TLorentzVector()

            l1_refit.SetPtEtaPhiM(l1_refit_pt, ntup.Z_lead_lepton_eta, ntup.Z_lead_lepton_phi, ntup.Z_lead_lepton_mass)
            l2_refit.SetPtEtaPhiM(l2_refit_pt, ntup.Z_sublead_lepton_eta, ntup.Z_sublead_lepton_phi, ntup.Z_sublead_lepton_mass)
            
            if fsr_clean:
                Z_refit = l1_refit + l2_refit + gamma_fsr
                Higgs_refit = Z_refit + gamma
                Z = l1 + l2
                Z_fsr = l1 + l2 + gamma_fsr
            else:
                Z_refit = l1_refit + l2_refit
                Higgs_refit = Z_refit + gamma
                Z = l1 + l2
                Z_fsr = l1 + l2

            Higgs = Z_fsr + gamma

            #print("[[INFO]]: ", l1.Pt(), l1_pterr, l1.Eta(), l1.Phi(), l1.M(), l2.Pt(), l2_pterr, l2.Eta(), l2.Phi(), l2.M())
            #print("[[INFO]]: ", ntup.Z_lead_lepton_ptE_error, ntup.Z_sublead_lepton_ptE_error)
            #print("[[INFO]]: id: {}, l1_pt: {}, l1_pt_err: {}, l2_pt: {}, l2_pt_err: {}, l1_refit_pt: {}, l2_refit_pt: {}, Z_mass1: {}, Z_mass2: {}, Z_mass_refit: {}".format(ntup.Z_lead_lepton_id, l1.Pt(), l1_pterr, l2.Pt(), l2_pterr, l1_refit_pt, l2_refit_pt, ntup.Z_mass, Z.M(), Z_refit.M()))
            #if n_fsr>0:
            #    print("[[INFO]] fsr event:", ntup.event)
            #    sys.exit(0)
            if args.single_check:
                print("[[INFO]] lepton1: l1_id: {}, l1_pt: {}, l1_eta: {}, l1_phi: {}, l1_pterr: {}, l1_pt_refit: {}, l1_eta_refit: {}, l1_phi_refit: {}".format(ntup.Z_lead_lepton_id, l1.Pt(), l1.Eta(), l1.Phi(), l1_pterr, l1_refit.Pt(), l1_refit.Eta(), l1_refit.Phi()))
                print("[[INFO]] lepton2: l2_id: {}, l2_pt: {}, l2_eta: {}, l2_phi: {}, l2_pterr: {}, l2_pt_refit: {}, l2_eta_refit: {}, l2_phi_refit: {}".format(ntup.Z_sublead_lepton_id, l2.Pt(), l2.Eta(), l2.Phi(), l2_pterr, l2_refit.Pt(), l2_refit.Eta(), l2_refit.Phi()))
                print("[[INFO]] photon: photon_pt: {}, photon_eta: {}, photon_phi: {}".format(gamma.Pt(), gamma.Eta(), gamma.Phi()))
                print("[[INFO]] photon_fsr: photon_fsr_pt: {}, photon_fsr_eta: {}, photon_fsr_phi: {}, photon_fsr_pterr: {}".format(gamma_fsr.Pt(), gamma_fsr.Eta(), gamma_fsr.Phi(), fsr_pterr))
                print("[[INFO]] Z: Z_mass: {}, Z_pt: {}, Z_eta: {}, Z_phi: {}".format(Z.M(), Z.Pt(), Z.Eta(), Z.Phi()))
                print("[[INFO]] Z_fsr: Z_fsr_mass: {}, Z_fsr_pt: {}, Z_fsr_eta: {}, Z_fsr_phi: {}".format(Z_fsr.M(), Z_fsr.Pt(), Z_fsr.Eta(), Z_fsr.Phi()))
                print("[[INFO]] Z_refit: Z_refit_mass: {}, Z_refit_pt: {}, Z_refit_eta: {}, Z_refit_phi: {}".format(Z_refit.M(), Z_refit.Pt(), Z_refit.Eta(), Z_refit.Phi()))
                print("[[INFO]] Z: H_mass: {}, H_pt: {}, H_eta: {}, H_phi: {}".format(Higgs.M(), Higgs.Pt(), Higgs.Eta(), Higgs.Phi()))
                print("[[INFO]] Z: H_refit_mass: {}, H_refit_pt: {}, H_refit_eta: {}, H_refit_phi: {}".format(Higgs_refit.M(), Higgs_refit.Pt(), Higgs_refit.Eta(), Higgs_refit.Phi()))
                sys.exit(0)

            branch_new["Z_lead_lepton_E"][0] = l1.E()
            branch_new["Z_sublead_lepton_E"][0] = l2.E()
            branch_new["Z_lead_lepton_pt_refit"][0] = l1_refit_pt
            branch_new["Z_sublead_lepton_pt_refit"][0] = l2_refit_pt
            branch_new["Z_pt_refit"][0] = Z_refit.Pt()
            branch_new["Z_eta_refit"][0] = Z_refit.Eta()
            branch_new["Z_phi_refit"][0] = Z_refit.Phi()
            branch_new["Z_mass_refit"][0] = Z_refit.M()
            branch_new["H_pt_refit"][0] = Higgs_refit.Pt()
            branch_new["H_eta_refit"][0] = Higgs_refit.Eta()
            branch_new["H_phi_refit"][0] = Higgs_refit.Phi()
            branch_new["H_mass_refit"][0] = Higgs_refit.M()


            passedEvents.Fill()
            
        # end loop
        
        out_file.Write()
        out_file.Close()


def ZConstrainReFit_plot(outpath):
    cat = args.tree
    analyzer_cfg = AC.Analyzer_Config(channel='inclusive', year=args.year, treename=cat)

    analyzer_cfg.Print_Config()
    analyzer_cfg.plot_output_path = "{}/ZRefit_plots".format(outpath)
    ntuples = LoadNtuples(analyzer_cfg)

    #for sample in analyzer_cfg.samp_names:
    for sample in analyzer_cfg.sig_names:
        print("[[INFO]] Start plot the Z contrain refit for sample: {}".format(sample))

        ntup = ntuples[sample] # just a short name

        var_names = ["H_mass", "H_mass_refit", "Z_mass", "Z_mass_refit"]
        hists = {}
        for var in var_names:
            if args.channel == 'ele':
                truelineshape = Zrefit.getParameters("{}/GenZ_truelineshapeele.txt".format(args.truelineshape))
                if "H" in var:
                    hists[var] = Zrefit.getVariableHists(ntup, var, 100, 110., 140.,"abs(Z_lead_lepton_id)==11")
                else:
                    hists[var] = Zrefit.getVariableHists(ntup, var, 100, 50., 130.,"abs(Z_lead_lepton_id)==11")
            elif args.channel == 'mu':
                truelineshape = Zrefit.getParameters("{}/GenZ_truelineshapemu.txt".format(args.truelineshape))
                if "H" in var:
                    hists[var] = Zrefit.getVariableHists(ntup, var, 100, 110., 140.,"abs(Z_lead_lepton_id)==13")
                else:
                    hists[var] = Zrefit.getVariableHists(ntup, var, 100, 50., 130.,"abs(Z_lead_lepton_id)==13")
            else:
                if "H" in var:
                    hists[var] = Zrefit.getVariableHists(ntup, var, 100, 110., 140.,"1,0")
                else:
                    hists[var] = Zrefit.getVariableHists(ntup, var, 100, 50., 130.,"1.0")


        Zrefit.makeComparisonPlot(hists["H_mass_refit"], hists["H_mass"], "H_mass", analyzer_cfg.plot_output_path, args.channel)
        Zrefit.makeComparisonPlot(hists["Z_mass_refit"], hists["Z_mass"], "Z_mass", analyzer_cfg.plot_output_path, args.channel)

        Zrefit.makeTrueLineshapeComparisonPlot(hists["Z_mass_refit"],truelineshape, "Z_mass_refit", analyzer_cfg.plot_output_path, args.channel)
        Zrefit.makeTrueLineshapeComparisonPlot(hists["Z_mass"],truelineshape, "Z_mass", analyzer_cfg.plot_output_path, args.channel)

def main():


    sw = ROOT.TStopwatch()
    sw.Start()

    print("[[INFO]] Start running Z refit")

    if args.trueLineshapeFit:
        print("[[INFO]] Start to get the true lineshape of Z by fitting Gen Z distribution")
        get_GenZ_lineshape(args.genFile, args.tree, args.output)

    if args.doReFit:
        print("[[INFO]] Start the Z contrain refit")
        ZConstrainReFit(args.tree, args.output, args.splitJobs)
        print("[[INFO]] Z contrain refit ended")

    if args.plotReFit:
        print("[[INFO]] Start plot the Z contrain refit")
        ZConstrainReFit_plot(args.output)
        print("[[INFO]] Z contrain refit plot ended")

    #log_file.close()
    sw.Stop()
    print('Real time: ' + str(round(sw.RealTime() / 60.0,2)) + ' minutes')
    print('CPU time:  ' + str(round(sw.CpuTime() / 60.0,2)) + ' minutes')

main()

