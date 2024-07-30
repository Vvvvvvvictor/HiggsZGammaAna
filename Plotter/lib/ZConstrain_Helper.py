from ROOT import *
from array import array
import math

'''
variables_name = [
    "H_pt",
    "H_eta",
    "H_phi",
    "H_mass",
    "Z_pt",
    "Z_eta",
    "Z_phi",
    "Z_mass",
    "Z_lead_lepton_pt",
    "Z_lead_lepton_eta",
    "Z_lead_lepton_phi",
    "Z_lead_lepton_mass",
    "Z_lead_lepton_charge",
    "Z_lead_lepton_id",
    "Z_lead_lepton_ptE_error",
    "Z_sublead_lepton_pt",
    "Z_sublead_lepton_eta",
    "Z_sublead_lepton_phi",
    "Z_sublead_lepton_mass",
    "Z_sublead_lepton_charge",
    "Z_sublead_lepton_id",
    "Z_sublead_lepton_ptE_error",
    "electron_1_energyErr",
    "electron_2_energyErr",
    "electron_3_energyErr",
    "electron_4_energyErr",
    "muon_1_ptErr",
    "muon_2_ptErr",
    "muon_3_ptErr",
    "muon_4_ptErr",
    "gamma_pt",
    "gamma_eta",
    "gamma_phi",
    "gamma_mass",
    "gamma_mvaID",
    "gamma_energyErr",
    "gamma_sieie",
    "gamma_hoe",
    "gamma_r9",
    "gamma_chiso",
    "gamma_alliso",
    "gamma_mvaID_WP80",
    "gamma_mvaID_WPL",
    "gamma_fsr_pt",
    "gamma_fsr_eta",
    "gamma_fsr_phi",
    "gamma_fsr_mass",
    "gamma_fsr_relIso03",
    "gamma_fsr_dROverEt2",
    "jet_1_pt",
    "jet_1_eta",
    "jet_1_phi",
    "jet_1_mass",
    "jet_1_btagDeepFlavB",
    "jet_2_pt",
    "jet_2_eta",
    "jet_2_phi",
    "jet_2_mass",
    "jet_2_btagDeepFlavB",
    "jet_3_pt",
    "jet_3_eta",
    "jet_3_phi",
    "jet_3_mass",
    "jet_3_btagDeepFlavB",
    "jet_4_pt",
    "jet_4_eta",
    "jet_4_phi",
    "jet_4_mass",
    "jet_4_btagDeepFlavB",
    "additional_lepton_1_pt",
    "additional_lepton_1_eta",
    "additional_lepton_1_phi",
    "additional_lepton_1_mass",
    "additional_lepton_1_charge",
    "additional_lepton_1_id",
    "additional_lepton_2_pt",
    "additional_lepton_2_eta",
    "additional_lepton_2_phi",
    "additional_lepton_2_mass",
    "additional_lepton_2_charge",
    "additional_lepton_2_id",
    "GenHzgHiggs_pt",
    "GenHzgHiggs_eta",
    "GenHzgHiggs_phi",
    "GenHzgHiggs_mass",
    "GenHzgLeadGenChild_pt",
    "GenHzgLeadGenChild_eta",
    "GenHzgLeadGenChild_phi",
    "GenHzgLeadGenChild_mass",
    "GenHzgLeadGenChild_pdgId",
    "GenHzgSubleadGenChild_pt",
    "GenHzgSubleadGenChild_eta",
    "GenHzgSubleadGenChild_phi",
    "GenHzgSubleadGenChild_mass",
    "GenHzgSubleadGenChild_pdgId",
    "GenHzgLeadGenChildChild1_pt",
    "GenHzgLeadGenChildChild1_eta",
    "GenHzgLeadGenChildChild1_phi",
    "GenHzgLeadGenChildChild1_mass",
    "GenHzgLeadGenChildChild1_pdgId",
    "GenHzgLeadGenChildChild1_pt",
    "GenHzgLeadGenChildChild2_eta",
    "GenHzgLeadGenChildChild2_phi",
    "GenHzgLeadGenChildChild2_mass",
    "GenHzgLeadGenChildChild2_pdgId",
    "n_jets",
    "n_leptons",
    "n_electrons",
    "n_muons",
    "n_iso_photons",
    "MET_pt",
    "MET_phi",
    "weight_central",
    "run",
    "event",
    "luminosityBlock",
    'HZ_relM',
    'H_relpt',
    'Z_relpt',
    'Z_lead_lepton_relpt',
    'Z_sublead_lepton_relpt',
    'gamma_relpt',
    'jet_1_relpt',
    'jet_2_relpt',
    'MET_relpt',
    'gamma_ptRelErr',
    'G_ECM',
    'Z_ECM',
    'Z_rapCM',
    'l_rapCM',
    'HZ_deltaRap',
    'l_cosProdAngle',
    'Z_cosProdAngle',
    'll_deltaR',
    'leadLG_deltaR',
    'ZG_deltaR',
    'subleadLG_deltaR',
    'H_deltaphi',
    'Z_deltaphi',
    'Z_lead_lepton_deltaphi',
    'Z_sublead_lepton_deltaphi',
    'additional_lepton_1_deltaphi',
    'additional_lepton_2_deltaphi',
    'MET_deltaphi',
    'weight',
    'mass_jj',
    'H_ptt',
    'H_al',
    'H_bt',
    'Z_cos_theta',
    'lep_cos_theta',
    'lep_phi',
    'l1g_deltaR',
    'l2g_deltaR',
    'delta_eta_jj',
    'delta_phi_jj',
    'delta_phi_zgjj',
    'photon_zeppenfeld',
    'H_zeppenfeld',
    'pt_balance',
    'pt_balance_0j',
    'pt_balance_1j',
    'is_center',
    'jet_1_deltaphi',
    'jet_2_deltaphi',
    'jet_3_deltaphi',
    'jet_4_deltaphi',
    'jet1G_deltaR',
    'jet2G_deltaR',
    'jet3G_deltaR',
    'jet4G_deltaR',
    'reweight',
    'BDT_score'
    ]
'''

def getParameters(file):
    params = {}
    
    f = open(file)
    for line in f.readlines():
        line = line.strip('\n')
        params[line.split()[0]] = float(line.split()[1])

    return params


def branch_setup(tree):
    global variables_branch, variables_name

    variables_branch = {}
    variables_name = [
    "H_pt",
    "H_eta",
    "H_phi",
    "H_mass",
    "Z_pt",
    "Z_eta",
    "Z_phi",
    "Z_mass",
    "Z_lead_lepton_pt",
    "Z_lead_lepton_eta",
    "Z_lead_lepton_phi",
    "Z_lead_lepton_mass",
    "Z_lead_lepton_charge",
    "Z_lead_lepton_id",
    "Z_lead_lepton_ptE_error",
    "Z_sublead_lepton_pt",
    "Z_sublead_lepton_eta",
    "Z_sublead_lepton_phi",
    "Z_sublead_lepton_mass",
    "Z_sublead_lepton_charge",
    "Z_sublead_lepton_id",
    "Z_sublead_lepton_ptE_error",
    "electron_1_energyErr",
    "electron_2_energyErr",
    "electron_3_energyErr",
    "electron_4_energyErr",
    "muon_1_ptErr",
    "muon_2_ptErr",
    "muon_3_ptErr",
    "muon_4_ptErr",
    "gamma_pt",
    "gamma_eta",
    "gamma_phi",
    "gamma_mass",
    "gamma_mvaID",
    "gamma_energyErr",
    "gamma_sieie",
    "gamma_hoe",
    "gamma_r9",
    "gamma_chiso",
    "gamma_alliso",
    "gamma_mvaID_WP80",
    "gamma_mvaID_WPL",
    "gamma_fsr_pt",
    "gamma_fsr_eta",
    "gamma_fsr_phi",
    "gamma_fsr_mass",
    "gamma_fsr_relIso03",
    "gamma_fsr_dROverEt2",
    "jet_1_pt",
    "jet_1_eta",
    "jet_1_phi",
    "jet_1_mass",
    "jet_1_btagDeepFlavB",
    "jet_2_pt",
    "jet_2_eta",
    "jet_2_phi",
    "jet_2_mass",
    "jet_2_btagDeepFlavB",
    "jet_3_pt",
    "jet_3_eta",
    "jet_3_phi",
    "jet_3_mass",
    "jet_3_btagDeepFlavB",
    "jet_4_pt",
    "jet_4_eta",
    "jet_4_phi",
    "jet_4_mass",
    "jet_4_btagDeepFlavB",
    "additional_lepton_1_pt",
    "additional_lepton_1_eta",
    "additional_lepton_1_phi",
    "additional_lepton_1_mass",
    "additional_lepton_1_charge",
    "additional_lepton_1_id",
    "additional_lepton_2_pt",
    "additional_lepton_2_eta",
    "additional_lepton_2_phi",
    "additional_lepton_2_mass",
    "additional_lepton_2_charge",
    "additional_lepton_2_id",
    "GenHzgHiggs_pt",
    "GenHzgHiggs_eta",
    "GenHzgHiggs_phi",
    "GenHzgHiggs_mass",
    "GenHzgLeadGenChild_pt",
    "GenHzgLeadGenChild_eta",
    "GenHzgLeadGenChild_phi",
    "GenHzgLeadGenChild_mass",
    "GenHzgLeadGenChild_pdgId",
    "GenHzgSubleadGenChild_pt",
    "GenHzgSubleadGenChild_eta",
    "GenHzgSubleadGenChild_phi",
    "GenHzgSubleadGenChild_mass",
    "GenHzgSubleadGenChild_pdgId",
    "GenHzgLeadGenChildChild1_pt",
    "GenHzgLeadGenChildChild1_eta",
    "GenHzgLeadGenChildChild1_phi",
    "GenHzgLeadGenChildChild1_mass",
    "GenHzgLeadGenChildChild1_pdgId",
    "GenHzgLeadGenChildChild1_pt",
    "GenHzgLeadGenChildChild2_eta",
    "GenHzgLeadGenChildChild2_phi",
    "GenHzgLeadGenChildChild2_mass",
    "GenHzgLeadGenChildChild2_pdgId",
    "n_jets",
    "n_leptons",
    "n_electrons",
    "n_muons",
    "n_iso_photons",
    "MET_pt",
    "MET_phi",
    "weight_central",
    "run",
    "event",
    "luminosityBlock"
    ]
    
    # setup the branchs
    for var in variables_name:
        variables_branch[var] = array('f',[0.])

        tree.Branch(var, variables_branch[var],"{}/F".format(var))

    return variables_name, variables_branch

def add_branch(tree):
    variables_branch_new = {}
    variables_name_new = ["Z_lead_lepton_E","Z_sublead_lepton_E","Z_lead_lepton_pt_refit","Z_sublead_lepton_pt_refit", "Z_pt_refit", "Z_eta_refit", "Z_phi_refit", "Z_mass_refit", "H_pt_refit", "H_eta_refit", "H_phi_refit", "H_mass_refit"]

    # setup the branchs
    for var in variables_name_new:
        variables_branch_new[var] = array('f',[0.])

        tree.Branch(var, variables_branch_new[var],"{}/F".format(var))

    return variables_name_new, variables_branch_new


def getEnergyResolutionEm(CorrectedEnergy, eta):
    if abs(eta) < 1.48:
        C = 0.35 / 100.
        S = 5.51 / 100.
        N = 98. / 1000.
    else:
        C = 0.
        S = 12.8 / 100.
        N = 440. / 1000.

    result = math.sqrt(C * C * CorrectedEnergy * CorrectedEnergy + S * S * CorrectedEnergy + N * N)

    return result

def gamma_fsr_pterr(gamma_fsr):
    perr = getEnergyResolutionEm(gamma_fsr.E(), gamma_fsr.Eta())

    pterr = perr*gamma_fsr.Pt()/gamma_fsr.P()

    return pterr


def MakeZReFitModel(PTRECO1_lep, PTErr1_lep, Theta1_lep, Phi1_lep, M1, PTRECO2_lep, PTErr2_lep, Theta2_lep, Phi2_lep, M2, PTRECO1_gamma, PTErr1_gamma, Theta1_gamma, Phi1_gamma, n_fsr, truelineshape_params):

    global pTRECO1_lep, pTRECO2_lep, pTMean1_lep, pTMean2_lep, pTSigma1_lep, pTSigma2_lep, theta1_lep, theta2_lep, phi1_lep, phi2_lep, m1, m2
    global pTRECO1_gamma, pTMean1_gamma, pTSigma1_gamma, theta1_gamma, phi1_gamma
    global gauss1_lep, gauss2_lep, gauss1_gamma, makeE_lep, E1_lep, E2_lep, makeE_gamma, E1_gamma, dotProduct_3d, p1v3D2, p1v3Dph1, p2v3Dph1, dotProduct_4d, p1D2, p1Dph1, p2Dph1
    global mZ, meanCB, sigmaCB, alphaCB, nCB, meanGauss1, sigmaGauss1, f1, meanGauss2, sigmaGauss2, f2, meanGauss3, sigmaGauss3, f3, singleCB, gaussShape1, CBplusGauss, gaussShape2, CBplusGaussplusGauss, gaussShape3, CBplusGaussplusGaussplusGauss, model, rastmp, pTs, r, covMatrix, finalPars, pT1_lep_refit, pT2_lep_refit, pTErr1_lep_refit, pTErr2_lep_refit


    #lep
    pTRECO1_lep = RooRealVar("pTRECO1_lep", "pTRECO1_lep", PTRECO1_lep, 5., 500.)
    pTRECO2_lep = RooRealVar("pTRECO2_lep", "pTRECO2_lep", PTRECO2_lep, 5., 500.)
    pTMean1_lep = RooRealVar("pTMean1_lep", "pTMean1_lep", PTRECO1_lep, max(5.0, PTRECO1_lep-2*PTErr1_lep), PTRECO1_lep+2*PTErr1_lep)
    pTMean2_lep = RooRealVar("pTMean2_lep", "pTMean2_lep", PTRECO2_lep, max(5.0, PTRECO2_lep-2*PTErr2_lep), PTRECO2_lep+2*PTErr2_lep)
    pTSigma1_lep = RooRealVar("pTSigma1_lep", "pTSigma1_lep", PTErr1_lep)
    pTSigma2_lep = RooRealVar("pTSigma2_lep", "pTSigma2_lep", PTErr2_lep)
    theta1_lep = RooRealVar("theta1_lep", "theta1_lep", Theta1_lep)
    theta2_lep = RooRealVar("theta2_lep", "theta2_lep", Theta2_lep)
    phi1_lep = RooRealVar("phi1_lep", "phi1_lep", Phi1_lep)
    phi2_lep = RooRealVar("phi2_lep", "phi2_lep", Phi2_lep)
    m1 = RooRealVar("m1", "m1", M1)
    m2 = RooRealVar("m2", "m2", M2)


    # gamma
    pTRECO1_gamma = RooRealVar("pTRECO1_gamma", "pTRECO1_gamma", PTRECO1_gamma, 5., 500.)
    pTMean1_gamma = RooRealVar("pTMean1_gamma", "pTMean1_gamma", PTRECO1_gamma, max(0.5, PTRECO1_gamma-2*PTErr1_gamma), PTRECO1_gamma+2*PTErr1_gamma)
    pTSigma1_gamma = RooRealVar("pTSigma1_gamma", "pTSigma1_gamma", PTErr1_gamma)
    theta1_gamma = RooRealVar("theta1_gamma", "theta1_gamma", Theta1_gamma)
    phi1_gamma = RooRealVar("phi1_gamma", "phi1_gamma", Phi1_gamma)
    
    #pTSigma1_lep.setConstant(kTRUE)
    #pTSigma2_lep.setConstant(kTRUE)
    #theta1_lep.setConstant(kTRUE)
    #theta2_lep.setConstant(kTRUE)
    #phi1_lep.setConstant(kTRUE)
    #phi2_lep.setConstant(kTRUE)
    #m1.setConstant(kTRUE)
    #m2.setConstant(kTRUE)
    #theta1_gamma.setConstant(kTRUE)
    #phi1_gamma.setConstant(kTRUE)
    #pTSigma1_gamma.setConstant(kTRUE)

    # gauss
    gauss1_lep = RooGaussian("gauss1_lep", "gauss1_lep", pTRECO1_lep, pTMean1_lep, pTSigma1_lep)
    gauss2_lep = RooGaussian("gauss2_lep", "gauss2_lep", pTRECO2_lep, pTMean2_lep, pTSigma2_lep)
    gauss1_gamma = RooGaussian("gauss1_gamma", "gauss1_gamma", pTRECO1_gamma, pTMean1_gamma, pTSigma1_gamma)

    makeE_lep = "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)"
    E1_lep = RooFormulaVar("E1_lep", makeE_lep, RooArgList(pTMean1_lep, theta1_lep, m1))
    E2_lep = RooFormulaVar("E2_lep", makeE_lep, RooArgList(pTMean2_lep, theta2_lep, m2))

    makeE_gamma = "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))"
    E1_gamma = RooFormulaVar("E1_gamma", makeE_gamma, RooArgList(pTMean1_gamma, theta1_gamma))

    # dotProduct 3d
    dotProduct_3d = "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3)))+(TMath::Cos(@4-@5)))"
    p1v3D2 = RooFormulaVar("p1v3D2", dotProduct_3d, RooArgList(pTMean1_lep, pTMean2_lep, theta1_lep, theta2_lep, phi1_lep, phi2_lep))
    p1v3Dph1 = RooFormulaVar("p1v3Dph1", dotProduct_3d, RooArgList(pTMean1_lep, pTMean1_gamma, theta1_lep, theta1_gamma, phi1_lep, phi1_gamma))
    p2v3Dph1 = RooFormulaVar("p2v3Dph1", dotProduct_3d, RooArgList(pTMean2_lep, pTMean1_gamma, theta2_lep, theta1_gamma, phi2_lep, phi1_gamma))

    dotProduct_4d = "@0*@1-@2"
    p1D2 = RooFormulaVar("p1D2", dotProduct_4d, RooArgList(E1_lep, E2_lep, p1v3D2))
    p1Dph1 = RooFormulaVar("p1Dph1", dotProduct_4d, RooArgList(E1_lep, E1_gamma, p1v3Dph1))
    p2Dph1 = RooFormulaVar("p2Dph1", dotProduct_4d, RooArgList(E2_lep, E1_gamma, p2v3Dph1))

    if n_fsr == 1:
        mZ = RooFormulaVar("mZ", "TMath::Sqrt(2*@0+2*@1+2*@2+@3*@3+@4*@4)", RooArgList(p1D2, p1Dph1, p2Dph1, m1, m2))
        #mZ = RooRealVar("mZ", "TMath::Sqrt(2*p1D2+2*p1Dph1+2*p2Dph1+m1*m1+m2*m2)", 91.0, 50., 130.)
    else:
        mZ = RooFormulaVar("mZ", "TMath::Sqrt(2*@0+@1*@1+@2*@2)", RooArgList(p1D2, m1, m2))
        #mZ = RooRealVar("mZ", "TMath::Sqrt(2*p1D2+m1*m1+m2*m2)", 91.0, 50., 130.)
    #mZ.setVal(91.0)

    # true shape
    meanCB = RooRealVar("meanCB","",truelineshape_params["meanCB"])
    sigmaCB = RooRealVar("sigmaCB","",truelineshape_params["sigmaCB"])
    alphaCB = RooRealVar("alphaCB","",truelineshape_params["alphaCB"])
    nCB = RooRealVar("nCB","",truelineshape_params["nCB"])
    meanGauss1 = RooRealVar("meanGauss1","",truelineshape_params["meanGauss1"])
    sigmaGauss1 = RooRealVar("sigmaGauss1","",truelineshape_params["sigmaGauss1"], 0.0001, 1000.)
    f1 = RooRealVar("f1","",truelineshape_params["f1"])
    meanGauss2 = RooRealVar("meanGauss2","",truelineshape_params["meanGauss2"])
    sigmaGauss2 = RooRealVar("sigmaGauss2","",truelineshape_params["sigmaGauss2"], 0.0001, 1000.)
    f2 = RooRealVar("f2","",truelineshape_params["f2"])
    meanGauss3 = RooRealVar("meanGauss3","",truelineshape_params["meanGauss3"])
    sigmaGauss3 = RooRealVar("sigmaGauss3","",truelineshape_params["sigmaGauss3"], 0.0001, 1000.)
    f3 = RooRealVar("f3","",truelineshape_params["f3"])

    singleCB = RooCBShape("singleCB", "", mZ, meanCB, sigmaCB, alphaCB, nCB)
    gaussShape1 = RooGaussian("gaussShape1", "", mZ, meanGauss1, sigmaGauss1)
    CBplusGauss = RooAddPdf("CBplusGauss", "", RooArgList(singleCB, gaussShape1), f1)
    gaussShape2 = RooGaussian("gaussShape2", "", mZ, meanGauss2, sigmaGauss2)
    CBplusGaussplusGauss = RooAddPdf("CBplusGaussplusGauss", "", RooArgList(CBplusGauss, gaussShape2), f2)
    gaussShape3 = RooGaussian("gaussShape3", "", mZ, meanGauss3, sigmaGauss3)
    CBplusGaussplusGaussplusGauss = RooAddPdf("CBplusGaussplusGaussplusGauss", "", RooArgList(CBplusGaussplusGauss, gaussShape3), f3)

    CBplusGauss.fixCoefNormalization(kTRUE)
    CBplusGaussplusGauss.fixCoefNormalization(kTRUE)
    CBplusGaussplusGaussplusGauss.fixCoefNormalization(kTRUE)


    meanCB.setConstant(kTRUE)
    sigmaCB.setConstant(kTRUE)
    alphaCB.setConstant(kTRUE)
    nCB.setConstant(kTRUE)
    meanGauss1.setConstant(kTRUE)
    sigmaGauss1.setConstant(kTRUE)
    f1.setConstant(kTRUE)
    meanGauss2.setConstant(kTRUE)
    sigmaGauss2.setConstant(kTRUE)
    f2.setConstant(kTRUE)
    meanGauss3.setConstant(kTRUE)
    sigmaGauss3.setConstant(kTRUE)
    f3.setConstant(kTRUE)

    model = RooProdPdf("model","",RooArgList(gauss1_lep, gauss2_lep, CBplusGaussplusGaussplusGauss))

    #print("DEBUG: MakeZReFitModel: model float parameters: ", E1_lep, E2_lep, p1v3D2, p1D2)
    #print("DEBUG: MakeZReFitModel: model float parameters: ", sigmaGauss1, pTRECO1_lep, pTMean1_lep, pTSigma1_lep, pTMean2_lep, pTSigma2_lep, mZ)

    # make fit
    rastmp = RooArgSet(pTRECO1_lep, pTRECO2_lep)

    pTs = RooDataSet("pTs","pTs", rastmp)
    pTs.add(rastmp)

    #r = model.fitTo(pTs, RooFit.Save(),RooFit.Constrain(mZ),RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))
    if n_fsr == 1:
        r = model.fitTo(pTs, RooFit.Save(),RooFit.Constrain(RooArgSet(pTMean1_lep, pTMean2_lep, pTMean1_gamma)),RooFit.PrintLevel(-1),RooFit.PrintEvalErrors(-1))
    else:
        r = model.fitTo(pTs, RooFit.Save(),RooFit.Constrain(RooArgSet(pTMean1_lep, pTMean2_lep)),RooFit.PrintLevel(-1),RooFit.PrintEvalErrors(-1))

    #print("DEBUG: MakeZReFitModel: model float parameters: ", sigmaGauss1, pTRECO1_lep, pTMean1_lep, pTSigma1_lep, pTMean2_lep, pTSigma2_lep, mZ)

    covMatrix = r.covarianceMatrix()
    finalPars = r.floatParsFinal()

    pT1_lep_refit = pTMean1_lep.getVal()
    pT2_lep_refit = pTMean2_lep.getVal()
    pTErr1_lep_refit = pTMean1_lep.getError()
    pTErr2_lep_refit = pTMean2_lep.getError()

    return pT1_lep_refit, pT2_lep_refit, pTErr1_lep_refit, pTErr2_lep_refit


def getVariableHists(Tree,varName,nBins,xMin,xMax,cut="1.0"):

    #Make a canvas to do the work on
    canvas = TCanvas(varName,varName,1000,800)

    #Extract the relevant variable from the trees
    #Tree.Draw("{0}>>hist{0}({1},{2},{3})".format(varName,nBins,xMin,xMax),"weight_central*({0})".format(cut))
    Tree.Draw("{0}>>hist{0}({1},{2},{3})".format(varName,nBins,xMin,xMax),"1.0*({0})".format(cut))
    #Tree.Draw("{0}>>hist{0}({1},{2},{3})".format(varName,nBins,xMin,xMax))
    Hist = gDirectory.Get("hist{0}".format(varName))

    canvas.Clear()

    return Hist

def makeLegend():
    #Make a legend in the corner of the canvas
    legend = TLegend(0.65,0.7,0.9,0.9)
    legend.SetFillStyle(1001)
    legend.SetBorderSize(1)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetShadowColor(0)
    legend.SetFillColor(kWhite)

    return legend

def makeComparisonPlot(sigHist1,bkgHist1,varName, outpath, channel):

    #We will be making changes to these histograms, but we don't want them done universally, so here we'll copy them
    (sigHist,bkgHist) = (sigHist1,bkgHist1)

    #Make a canvas to do the work on
    canvas = TCanvas(varName+"_comp",varName+"_comp",1000,800)
    
    #Mke a legend
    legend = makeLegend()

    #Normalise the distributions and set their appropriate colours
    sigHist.SetFillColor(0)
    sigHist.SetLineColor(kBlue)
    sigHist.SetLineWidth(3)
    sigHist.Draw("HIST")
    sigHist.GetXaxis().SetTitle(varName)
    legend.AddEntry(sigHist,"refit","l")

    bkgHist.SetFillColor(0)
    bkgHist.SetLineColor(kRed)
    bkgHist.SetLineWidth(3)
    bkgHist.Draw("sameHIST")
    legend.AddEntry(bkgHist,"w/o refit","l")

    #Do a quick trick to make sure we see the entire distribution
    if bkgHist.GetMaximum() > sigHist.GetMaximum() : sigHist.SetMaximum(1.2*bkgHist.GetMaximum())

    legend.Draw("same")

    canvas.SaveAs("{}/{}{}.png".format(outpath, varName, channel))
    canvas.SaveAs("{}/{}{}.pdf".format(outpath, varName, channel))


def makeTrueLineshapeComparisonPlot(hist,truelineshape_params,name,outpath, channel):

    mZ = RooRealVar("mZ","mZ", 50., 130.)
    histZ = RooDataHist( "histZ", "histZ", mZ, hist )


    # true shape
    meanCB = RooRealVar("meanCB","",truelineshape_params["meanCB"])
    sigmaCB = RooRealVar("sigmaCB","",truelineshape_params["sigmaCB"])
    alphaCB = RooRealVar("alphaCB","",truelineshape_params["alphaCB"])
    nCB = RooRealVar("nCB","",truelineshape_params["nCB"])
    meanGauss1 = RooRealVar("meanGauss1","",truelineshape_params["meanGauss1"])
    sigmaGauss1 = RooRealVar("sigmaGauss1","",truelineshape_params["sigmaGauss1"], 0.0001, 1000.)
    f1 = RooRealVar("f1","",truelineshape_params["f1"])
    meanGauss2 = RooRealVar("meanGauss2","",truelineshape_params["meanGauss2"])
    sigmaGauss2 = RooRealVar("sigmaGauss2","",truelineshape_params["sigmaGauss2"], 0.0001, 1000.)
    f2 = RooRealVar("f2","",truelineshape_params["f2"])
    meanGauss3 = RooRealVar("meanGauss3","",truelineshape_params["meanGauss3"])
    sigmaGauss3 = RooRealVar("sigmaGauss3","",truelineshape_params["sigmaGauss3"], 0.0001, 1000.)
    f3 = RooRealVar("f3","",truelineshape_params["f3"])

    singleCB = RooCBShape("singleCB", "", mZ, meanCB, sigmaCB, alphaCB, nCB)
    gaussShape1 = RooGaussian("gaussShape1", "", mZ, meanGauss1, sigmaGauss1)
    CBplusGauss = RooAddPdf("CBplusGauss", "", RooArgList(singleCB, gaussShape1), f1)
    gaussShape2 = RooGaussian("gaussShape2", "", mZ, meanGauss2, sigmaGauss2)
    CBplusGaussplusGauss = RooAddPdf("CBplusGaussplusGauss", "", RooArgList(CBplusGauss, gaussShape2), f2)
    gaussShape3 = RooGaussian("gaussShape3", "", mZ, meanGauss3, sigmaGauss3)
    CBplusGaussplusGaussplusGauss = RooAddPdf("CBplusGaussplusGaussplusGauss", "", RooArgList(CBplusGaussplusGauss, gaussShape3), f3)

    nbins = 50
    canv = TCanvas("c1","c1",1500,1000)
    mZ.setBins(nbins)
    frame = mZ.frame(nbins)
    histZ.plotOn(frame,RooFit.Binning(nbins))

    CBplusGaussplusGaussplusGauss.plotOn(frame, RooFit.LineColor(kRed))

    frame.Draw()

    canv.SaveAs("{}/TureLineCompare_{}_{}.png".format(outpath, name, channel))
    canv.SaveAs("{}/TureLineCompare_{}_{}.pdf".format(outpath, name, channel))
