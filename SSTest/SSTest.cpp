#include <fstream>
#include "RooExponential.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
//#include "boost/program_options.hpp"                                                                             
#include "boost/lexical_cast.hpp"
#include <RooGaussModel.h>
#include <RooTruthModel.h>
#include <RooDecay.h>
#include "RooAddModel.h"
#include <RooNumConvPdf.h>
#include <RooFFTConvPdf.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TFrame.h"
#include "RooFitResult.h"
// #include "Maketree.h"
// #include "ModGaus11.h"
using namespace RooFit;
using namespace std;

void SSTest(int cat = 0, int sig = 0, TString channel = "two_jet", TString bkg_fun = "bern2")
{
    double mgg_low=100, mgg_high=170, bin_size=280;
    if (channel.EqualTo("two_jet") & (cat == 0))
    {
        mgg_low = 100;
        mgg_high = 180;
        bin_size = 320;
    }
    // if (channel.EqualTo("two_jet") & (cat == 3))
    // {
    //     mgg_low = 100;
    //     mgg_high = 150;
    //     bin_size = 200;
    // }
    //background MC template
    TH1F* hbkg;
    TFile* fbkg = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
    if (fbkg->GetListOfKeys()->Contains((Form("bkg_%s_cat%d",channel.Data(), cat))))
        hbkg = (TH1F*)fbkg->Get(Form("bkg_%s_cat%d",channel.Data(), cat));
    else
        abort();
    // TFile* fbkg = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_template_v3_cut2/bkg/bkg_0sig_cat%d.root", cat));
    // hbkg = (TH1F*)fbkg->Get(Form("mass_cat%d", cat));
    double dataevents = hbkg->Integral();
    double mcsbevents = hbkg->Integral(0,122-mgg_low)+hbkg->Integral(bin_size-(mgg_high-128),bin_size);
	
    //signal MC
    TFile* fsig = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
    TH1F* hsig = (TH1F*)( fsig->Get(Form("sig_%s_cat%d", channel.Data(), cat)) );
    // TFile* fsig = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_template_v3_cut2/sig/FullSig_cat%d.root", cat));
    // TH1F* hsig = (TH1F*)fsig->Get(Form("mass_cat%d_FullSig", cat));
	double sigevents = hsig->Integral();

    //data full range
    TFile* ffr = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
    TH1F* hfr = (TH1F*)( ffr->Get(Form("data_full_%s_cat%d", channel.Data(), cat)) );
    // TFile* ffr = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_template_v3_cut2/data_sid/data_0sig_cat%d.root", cat));
    // TH1F* hfr = (TH1F*)ffr->Get(Form("mass_cat%d_data", cat));
	double frevents = hfr->Integral();

    //data side band
    TFile* fsb = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
    TH1F* hsb = (TH1F*)( fsb->Get(Form("data_%s_cat%d", channel.Data(), cat)) );
    // TFile* fsb = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_template_v3_cut2/data_sid/data_0sig_cat%d.root", cat));
    // TH1F* hsb = (TH1F*)fsb->Get(Form("mass_cat%d_data", cat));
	double sbevents = hsb->Integral();

    // hbkg->Scale(sbevents/mcsbevents);

    cout << "\n\tFinish preparing data!!!\n" << endl;

    //initializing
    RooRealVar mH("mH", "mH", 125.38, mgg_low, mgg_high);
    mH.setMin( mgg_low );
    mH.setMax( mgg_high );
    // mH.setMax( hbkg->GetBinCenter(hbkg->GetNbinsX()) + hbkg->GetBinWidth(hbkg->GetNbinsX()) );

    //background function fit
    RooDataHist* dbkg = new RooDataHist("bkg_mc","dataset with x", mH, hbkg);
    RooDataHist* dsb = new RooDataHist("data_sb","dataset with x", mH, hsb);
    RooDataHist* dfr = new RooDataHist("data_fr","dataset with x", mH, hfr);
    // cout<<"bkg hist integral "<<hbkg->Integral()<<" "<<dbkg->sumEntries()<<endl;

    RooRealVar nsig("nsig","nsig",sigevents,-100*sigevents,100*sigevents);
    //nsig.setConstant(kTRUE);
	RooRealVar nbkg("nbkg","nbkg",dataevents, 0.01*dataevents, 2*dataevents);

    RooRealVar mean(Form("mean_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("mean_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0);
    
    RooRealVar sigma_b1(Form("sigma_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.,1.,15.);
    RooRealVar sigma_b2(Form("sigma_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.,1.,15.);
    RooRealVar sigma_b3(Form("sigma_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.,1.,15.);
    RooRealVar sigma_b4(Form("sigma_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.,1.,15.);
    RooRealVar sigma_b5(Form("sigma_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.,1.,15.);
    
    RooRealVar step_b1(Form("step_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),110.,95,115.);
    RooRealVar step_b2(Form("step_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),110.,95.,115.);
    RooRealVar step_b3(Form("step_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),110.,95,115.);
    RooRealVar step_b4(Form("step_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),110.,95,115.);
    RooRealVar step_b5(Form("step_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("step_b5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),110.,95,115.);

    RooRealVar p0(Form("p0_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p0_env_pdf_ele_mu_cat%d_2020_13TeV",cat),15);
    RooRealVar b1p1(Form("b1p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b1p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b2p1(Form("b2p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b2p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b2p2(Form("b2p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b2p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1, -15.,15.);
    RooRealVar b3p1(Form("b3p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b3p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b3p2(Form("b3p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b3p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b3p3(Form("b3p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b3p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b4p1(Form("b4p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b4p2(Form("b4p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b4p3(Form("b4p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);//untag
    RooRealVar b4p4(Form("b4p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b4p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);//VBF&lepton
    RooRealVar b5p1(Form("b5p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);//untag
    RooRealVar b5p2(Form("b5p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b5p3(Form("b5p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b5p4(Form("b5p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);
    RooRealVar b5p5(Form("b5p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("b5p5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1,-15.,15.);

    RooRealVar sigma_pow1(Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,1.,8.);
    RooRealVar turnon_pow1(Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.,95,115);
    RooRealVar p1_pow1(Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.3453,-7,-2.);
    RooRealVar cp1_pow1(Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.7131e-05 ,0.0,1.);
    
    RooRealVar sigma_pow3(Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.5,1.,8.);
    RooRealVar turnon_pow3(Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107,95,115);
    RooRealVar p1_pow3(Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.3671,-10,-2.);
    RooRealVar cp1_pow3(Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.6543e-01,0.,1.);
    RooRealVar p3_pow3(Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.267,-8.,-2.);
    RooRealVar cp3_pow3(Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.1001e-06,0,1.);

    RooRealVar sigma_pow5(Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.0999,1.,8.);
    RooRealVar turnon_pow5(Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.21,95,115);
    RooRealVar p1_pow5(Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.8,-8,-2);
    RooRealVar cp1_pow5(Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.9425e-06,0.,1.);
    RooRealVar p3_pow5(Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.4921,-8.,-2.);
    RooRealVar cp3_pow5(Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.751e-02 ,0.,1.0);
    RooRealVar p5_pow5(Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.3125,-8.,-2.);
    RooRealVar cp5_pow5(Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.8914e-03,0.0,1.);

    RooRealVar sigma_lau1(Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.6496,1.,8.);
    RooRealVar turnon_lau1(Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.05,95,115);
    RooRealVar cl1_lau1(Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.e-12,0.0,1.);//untag
    RooRealVar cl2_lau1(Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.7 ,0.6,1.);

    RooRealVar sigma_lau2(Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.1145,1.,8.);
    RooRealVar turnon_lau2(Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.89,95,115);
    RooRealVar cl1_lau2(Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2209e-05,0.,1.);
    RooRealVar cl2_lau2(Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.2780e-01,0.5,1.0);
    RooRealVar cl3_lau2(Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2371e-07,0,1.);

    RooRealVar sigma_lau3(Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.53,1.,8.);
    RooRealVar turnon_lau3(Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),104.2,95,115);
    RooRealVar cl1_lau3(Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.3514e-10,0.,1.);
    RooRealVar cl2_lau3(Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5466e-08,0.,1.);
    RooRealVar cl3_lau3(Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.7939e-11,0,1.);
    RooRealVar cl4_lau3(Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.0697e-01,0.9,1.);

    RooRealVar sigma_lau4(Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.53,1.,8.);
    RooRealVar turnon_lau4(Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),102.2,95,115);
    RooRealVar cl1_lau4(Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.3514e-1,0.,1.);
    RooRealVar cl2_lau4(Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),9.5466e-1,0.,1.);
    RooRealVar cl3_lau4(Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),8.7939e-1,0,1.);
    RooRealVar cl4_lau4(Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.99,0.9,1.);
    RooRealVar cl5_lau4(Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.1,0,1.);

    RooRealVar sigma_lau5(Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.3253,1.,8.);
    RooRealVar turnon_lau5(Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),105.41,95,115);
    RooRealVar cl1_lau5(Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl1_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.8677e-05,0.,.1);
    RooRealVar cl2_lau5(Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl2_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.9106e-01,0.,.9);
    RooRealVar cl3_lau5(Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl3_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.7235e-09,0,.0001);
    RooRealVar cl4_lau5(Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl4_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),6.1470e-01 ,0.1,1.);
    RooRealVar cl5_lau5(Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl5_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.5457e-09,0,1e-4);
    RooRealVar cl6_lau5(Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cl6_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),2.7473e-02,0.0,0.5);

    RooRealVar sigma_exp1(Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.56,1.,7.);
    RooRealVar turnon_exp1(Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),106.97,95,115);
    RooRealVar p1_exp1(Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4e-02,-0.1,0.);
    RooRealVar cp1_exp1(Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.01,0,1.);

    RooRealVar sigma_exp3(Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.,1.,7.);
    RooRealVar turnon_exp3(Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.5,95.,115);
    RooRealVar p1_exp3(Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-3.4437e-02,-0.1,0.);
    RooRealVar cp1_exp3(Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),1.2335e-02,0.,1.);
    RooRealVar p3_exp3(Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-6.5070e-02,-0.5,0.);
    RooRealVar cp3_exp3(Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),5.9216e-01,0,1.);

    RooRealVar sigma_exp5(Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("sigma_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),4.12,1.,8);
    RooRealVar turnon_exp5(Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("turnon_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),107.21,95.,115);
    RooRealVar p1_exp5(Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-5.1816e-02,-0.1,0.);
    RooRealVar cp1_exp5(Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp1_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.5,0.0,1.);
    RooRealVar p3_exp5(Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-4.3062e-02,-0.5,0.);
    RooRealVar cp3_exp5(Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp3_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),0.05,0,1.);
    RooRealVar p5_exp5(Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("p5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),-9.2008e-04,-0.7,0.);
    RooRealVar cp5_exp5(Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("cp5_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),3.5533e-05,0,1.);
    
    RooGaussStepBernstein bern1(Form("bern1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_b1,step_b1, RooArgList(p0,b1p1));
    RooGaussStepBernstein bern2(Form("bern2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_b2,step_b2, RooArgList(p0,b2p1,b2p2));
    RooGaussStepBernstein bern3(Form("bern3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_b3,step_b3, RooArgList(p0,b3p1,b3p2,b3p3));
    RooGaussStepBernstein bern4(Form("bern4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_b4,step_b4, RooArgList(p0,b4p1,b4p2,b4p3,b4p4));
    RooGaussStepBernstein bern5(Form("bern5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("bern5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_b5,step_b5, RooArgList(p0,b5p1,b5p2,b5p3,b5p4,b5p5));
       
    
    RooGenericPdf step_pow1(Form("step_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*(@0)^(@2))", RooArgList(mH,turnon_pow1,p1_pow1,cp1_pow1));//step*(ax^b)
    RooGenericPdf step_pow3(Form("step_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4))", RooArgList(mH,turnon_pow3,p1_pow3,cp1_pow3,p3_pow3,cp3_pow3));//step*(ax^b+cx^d)
    RooGenericPdf step_pow5(Form("step_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(mH,turnon_pow5,p1_pow5,cp1_pow5,p3_pow5,cp3_pow5,p5_pow5,cp5_pow5));//step*(ax^b+cx^d+fx^g)
    RooGaussModel gau_pow1(Form("gau_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_pow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_pow1);
    RooGaussModel gau_pow3(Form("gau_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_pow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_pow3);
    RooGaussModel gau_pow5(Form("gau_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_pow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_pow5);
    RooFFTConvPdf gauxpow1(Form("gauxpow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxpow1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_pow1,gau_pow1);
    RooFFTConvPdf gauxpow3(Form("gauxpow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxpow3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_pow3,gau_pow3);
    RooFFTConvPdf gauxpow5(Form("gauxpow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxpow5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_pow5,gau_pow5);
    gauxpow1.setBufferFraction(0.25);
    gauxpow3.setBufferFraction(0.25);
    gauxpow5.setBufferFraction(0.25);
    
    
    RooGenericPdf step_lau1(Form("step_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(mH,turnon_lau1,cl1_lau1,cl2_lau1));//step*(ax^b)
    RooGenericPdf step_lau2(Form("step_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(mH,turnon_lau2,cl1_lau2,cl2_lau2,cl3_lau2));//step*(ax^b+cx^d+fx^g) 
    RooGenericPdf step_lau3(Form("step_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(mH,turnon_lau3,cl1_lau3,cl2_lau3,cl3_lau3,cl4_lau3));//step*(ax^b+cx^d)
    RooGenericPdf step_lau4(Form("step_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2))", RooArgList(mH,turnon_lau4,cl1_lau4,cl2_lau4,cl3_lau4,cl4_lau4,cl5_lau4));//step*(ax^b+cx^d)
    RooGenericPdf step_lau5(Form("step_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2)+@7*(@0)^(-7))", RooArgList(mH,turnon_lau5,cl1_lau5,cl2_lau5,cl3_lau5,cl4_lau5,cl5_lau5,cl6_lau5));//step*(ax^b+cx^d)

    RooGaussModel gau_lau1(Form("gau_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_lau1);
    RooGaussModel gau_lau2(Form("gau_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_lau2);
    RooGaussModel gau_lau3(Form("gau_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_lau3);
    RooGaussModel gau_lau4(Form("gau_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_lau4);
    RooGaussModel gau_lau5(Form("gau_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_lau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_lau5);

    RooFFTConvPdf gauxlau1(Form("gauxlau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_lau1,gau_lau1);
    RooFFTConvPdf gauxlau2(Form("gauxlau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau2_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_lau2,gau_lau2);
    RooFFTConvPdf gauxlau3(Form("gauxlau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_lau3,gau_lau3);
    RooFFTConvPdf gauxlau4(Form("gauxlau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau4_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_lau4,gau_lau4);
    RooFFTConvPdf gauxlau5(Form("gauxlau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxlau5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_lau5,gau_lau5);
    gauxlau1.setBufferFraction(0.25);
    gauxlau2.setBufferFraction(0.25);
    gauxlau3.setBufferFraction(0.25);
    gauxlau4.setBufferFraction(0.25);
    gauxlau5.setBufferFraction(0.25);
        
        
    RooGenericPdf step_exp1(Form("step_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(mH,turnon_exp1,p1_exp1,cp1_exp1));//step*(ax^b)
    RooGenericPdf step_exp3(Form("step_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(mH,turnon_exp3,p1_exp3,cp1_exp3,p3_exp3,cp3_exp3));//step*(ax^b+cx^d)
    RooGenericPdf step_exp5(Form("step_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), Form("step_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(mH,turnon_exp5,p1_exp5,cp1_exp5,p3_exp5,cp3_exp5,p5_exp5,cp5_exp5));//step*(ax^b+cx^d+fx^g)
    RooGaussModel gau_exp1(Form("gau_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_exp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_exp1);
    RooGaussModel gau_exp3(Form("gau_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_exp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_exp3);
    RooGaussModel gau_exp5(Form("gau_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gau_exp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,mean,sigma_exp5);
    RooFFTConvPdf gauxexp1(Form("gauxexp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxexp1_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_exp1,gau_exp1);
    RooFFTConvPdf gauxexp3(Form("gauxexp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxexp3_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_exp3,gau_exp3);
    RooFFTConvPdf gauxexp5(Form("gauxexp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),Form("gauxexp5_env_pdf_ele_mu_cat%d_2020_13TeV",cat),mH,step_exp5,gau_exp5);
    gauxexp1.setBufferFraction(0.25);
    gauxexp3.setBufferFraction(0.25);
    gauxexp5.setBufferFraction(0.25);

    cout << "\n\tFinish initializing fitting variables\n" << endl;

    RooAbsPdf* bkg_model; RooFitResult* bkg_model_fit;
    if(bkg_fun.Contains("gauxpow1", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxpow1.Clone();
    if(bkg_fun.Contains("gauxpow3", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxpow3.Clone();
    if(bkg_fun.Contains("gauxpow5", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxpow5.Clone();
    if(bkg_fun.Contains("gauxexp1", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxexp1.Clone();
    if(bkg_fun.Contains("gauxexp3", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxexp3.Clone();
    if(bkg_fun.Contains("gauxexp5", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxexp5.Clone();
    if(bkg_fun.Contains("gauxlau1", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxlau1.Clone();
    if(bkg_fun.Contains("gauxlau2", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxlau2.Clone();
    if(bkg_fun.Contains("gauxlau3", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxlau3.Clone();
    if(bkg_fun.Contains("gauxlau4", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxlau4.Clone();
    if(bkg_fun.Contains("gauxlau5", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) gauxlau5.Clone();
    if(bkg_fun.Contains("bern2", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) bern2.Clone();
    if(bkg_fun.Contains("bern3", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) bern3.Clone();
    if(bkg_fun.Contains("bern4", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) bern4.Clone();
    if(bkg_fun.Contains("bern5", TString::kIgnoreCase)) bkg_model = (RooAbsPdf*) bern5.Clone();
    // RooArgSet * params_bkg = bkg_model.getParameters((const RooArgSet*)(0));

    ofstream output(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/outputs/%s_%d_%dxsig.txt", channel.Data(), cat, sig), ofstream::app);

    TFile *f = new TFile(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/outputs/%s_%d_%dxsig.root", channel.Data(), cat, sig),"UPDATE"); 
    hbkg->Write(hbkg->GetName(), TObject::kOverwrite);
    hsb->Write(hsb->GetName(), TObject::kOverwrite);
    hfr->Write(hfr->GetName(), TObject::kOverwrite);
    hsig->Write(hsig->GetName(), TObject::kOverwrite);


    int flag; TString status;
    flag = 1; status = "Pass";
    double dmc, dss, ss, ss_mc, tot_err, ss_cor, delta, chi2=0., prob, nll;
    RooAbsCollection *floatPars; TIterator *iter;

    int bkg_npars, bkg_ndof;
    RooPlot *frame_bkg;

    mH.setRange("range_low",mgg_low,122);
    mH.setRange("signal",122,128);
    mH.setRange("range_high",128,mgg_high);

    //background function fit
    bkg_model_fit = bkg_model->fitTo(*dfr,Range("range"),SplitRange(true),Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE),EvalErrorWall(false)); 
    bkg_npars = bkg_model_fit->floatParsFinal().getSize();
    frame_bkg = mH.frame(Title(Form("Data side band with %s pdf", bkg_fun.Data())));
    bkg_ndof = bin_size-24-bkg_npars;
    dfr->plotOn(frame_bkg, Cut("mH>128 || mH<122"));
    bkg_model->plotOn(frame_bkg,NormRange("range_low,signal,range_high"));
    chi2 = frame_bkg->chiSquare(bkg_npars);
    bkg_model->SetName(bkg_fun);
    bkg_model->Write(bkg_model->GetName(), TObject::kOverwrite);
    nll = bkg_model_fit->minNll();
    prob = TMath::Prob(chi2*bkg_ndof, bkg_ndof);
    // if(prob<0.01) status = "Fail";
    output << "\t" << bkg_fun.Data() << "\tsb:\tnpars = " << bkg_npars << " \tchi^2 = " << chi2 << "\tprob = " << prob << "\tnll: " << nll << endl;
    bkg_model->paramOn(frame_bkg, Layout(0.34,0.96,0.89),Format("NEA",AutoPrecision(1)));
    frame_bkg->getAttText()->SetTextSize(0.03);
    frame_bkg->Draw();
    gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/data_sb_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

    cout << "\t=================================" << endl;
    cout << "\n\t Finish data side band fit\n" << endl;

    //background function fit
    bkg_model_fit = bkg_model->fitTo(*dbkg,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE),EvalErrorWall(false));
    bkg_npars = bkg_model_fit->floatParsFinal().getSize();
    bkg_ndof = bin_size-bkg_npars;
    frame_bkg = mH.frame(Title(Form("Background with %s pdf", bkg_fun.Data())));
    dbkg->plotOn(frame_bkg);
    bkg_model->plotOn(frame_bkg);
    bkg_model->paramOn(frame_bkg, Layout(0.34,0.96,0.89),Format("NEA",AutoPrecision(1)));
    frame_bkg->getAttText()->SetTextSize(0.03);
    bkg_model->SetName(bkg_fun);
    bkg_model->Write(bkg_model->GetName(), TObject::kOverwrite);
    nll = bkg_model_fit->minNll();
    output << "\t" << bkg_fun.Data() << "\tbkg:\tnpars = " << bkg_npars << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars)*bkg_ndof, bkg_ndof) << "\tnll: " << nll << endl;
    frame_bkg->Draw();
    // RooHist *hpull = frame_bkg->pullHist();
    // RooPlot *frame3 = mH.frame(Title("Pull Distribution"));
    // frame3->addPlotable(hpull, "P");
    // frame3->Draw();
    gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/mc_bkg_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

    // RooAbsCollection *m_bkgParameters = bkg_model->getParameters(RooArgSet())->selectByAttrib("Constant", false);
    // TIterator *it = m_bkgParameters->createIterator();
    // for (RooRealVar *p = (RooRealVar *)it->Next(); p != 0; p = (RooRealVar *)it->Next()) {
    //     p->setConstant(kTRUE);
    // }


    cout << "\t=================================\n";
    cout << "\n\tFinish background function fit\n" << endl;

    //signal function fit
    RooRealVar sigma("sigma","sigma",3.,0.1,5.); 
    RooRealVar MH("MH","MH",125., 124., 126.); 
    RooGaussian sig_gau("sig_gau","sig_gau",mH,MH,sigma);
    
    RooRealVar sigma_CB("sigma_CB","sigma_CB",0.5, 0.1, 4.); 
    RooRealVar alpha("alpha","alpha",0.5, 0., 1.0); 
    RooRealVar n_CB("n_CB","n_CB",12.,1.,50.); 
    RooRealVar fracG1("fracG1","fracG1",0.1,0.,1.0); 
    RooCBShape CBshape("CBShape", "CBShape", mH, MH, sigma_CB, alpha, n_CB);
    RooAddPdf* signal = new RooAddPdf("signal","signal",RooArgList(sig_gau, CBshape),fracG1);
    sigma.setConstant(false);MH.setConstant(false);sigma_CB.setConstant(false);alpha.setConstant(false);n_CB.setConstant(false);fracG1.setConstant(false);

    RooDataHist* dsig = new RooDataHist("data_bin","dataset with x", mH, hsig);
    RooFitResult *signal_fit;
    signal_fit = signal->fitTo(*dsig,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE),EvalErrorWall(false)); //FIXME kTRUE or kFALSE
    int sig_npars = signal_fit->floatParsFinal().getSize();
    int sig_ndof = bin_size-sig_npars;
    RooPlot *frame_sig = mH.frame(Title("Signal with distorted Gaussian pdf"));
    dsig->plotOn(frame_sig, DataError(RooAbsData::SumW2));
    signal->plotOn(frame_sig);
    signal->paramOn(frame_sig, Layout(0.4,0.96,0.89),Format("NEA",AutoPrecision(1)));
    signal->Write(signal->GetName(), TObject::kOverwrite);
    output << "\t" << bkg_fun.Data() << "\tsig:\tnpars = " << sig_npars << "\tchi^2 = " << frame_sig->chiSquare(sig_npars) << "\tprob = " << TMath::Prob(frame_sig->chiSquare(sig_npars)*sig_ndof, sig_ndof) << endl;
    frame_sig->Draw();
    gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/signal_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));


    sigma.setConstant(true);MH.setConstant(true);sigma_CB.setConstant(true);alpha.setConstant(true);n_CB.setConstant(true);fracG1.setConstant(true);

    cout << "\t=================================" << endl;
    cout << "\n\t Finish signal function fit\n" << endl;

    //inject signal
    TH1F* hdata = (TH1F*)hbkg->Clone();
    hdata->Add(hsig, sig);
    hdata->SetName("asimov data");
    hdata->Write(hdata->GetName(), TObject::kOverwrite);
    RooDataHist* ddata = new RooDataHist("data_bin","dataset with x", mH, hdata);

    // asimov dataset fit
    RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*signal, *bkg_model),RooArgList(nsig,nbkg));
    RooFitResult *model_fit;
    int data_npars, data_ndof, fit_status, tries;

    TCanvas *canv = new TCanvas();
    RooPlot *frame_data, *frame_data_trash;
    frame_data = mH.frame();
    frame_data_trash = mH.frame();
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
    pad1->SetBottomMargin(0.18);
    // pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.25);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    model_fit = model->fitTo(*ddata,Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE
    fit_status = model_fit->status();

    // fix the float parameters in model except nsig
    floatPars = model->getParameters(*ddata)->selectByAttrib("Constant", false);
    iter = floatPars->createIterator();
    for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
        if (var->GetName() != TString("nsig")) {
        var->setConstant(true);
        }
    }

    fit_status = -1;
    tries = 0;
    while((fit_status != 0) && (tries < 10)){
        // model_fit = model->chi2FitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),DataError(RooAbsData::SumW2)); //FIXME kTRUE or kFALSE
        model_fit = model->fitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE
        fit_status = model_fit->status();
        tries++;
    }

    ss_mc = nsig.getVal();
    dmc = nsig.getError();
    data_ndof = bin_size-data_npars;

    // if(ss < 0) frame_data_trash->SetMinimum(ss);
    ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
    RooHist *plotdata = (RooHist*)frame_data->getObject(frame_data->numItems()-1);
    ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
    model->plotOn(frame_data_trash, Name("fit"));
    chi2 = frame_data_trash->chiSquare(data_npars);
    prob = TMath::Prob(chi2*data_ndof, data_ndof);
    output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << ": " << chi2 << "\tprob = " << prob << "\tfitting status = " << fit_status << endl;

    // unfix the float parameters in model
    iter = floatPars->createIterator();
    for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
        var->setConstant(false);
    }

    // dmc = nsig.getError();
    // data_npars = model_fit->floatParsFinal().getSize();
    // data_ndof = bin_size-data_npars;
    // ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
    // RooHist *plotdata = (RooHist*)frame_data_trash->getObject(frame_data_trash->numItems()-1);
    // model->plotOn(frame_data_trash, Name("fitmcweight"));
    // chi2 = frame_data_trash->chiSquare(data_npars);
    // prob = TMath::Prob(chi2*data_ndof, data_ndof);
    // frame_data_trash->remove("fitmcweight");
    // output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << "\tprob = " << prob << "\tstatus = " << fit_status << endl;

    model_fit = model->fitTo(*ddata,Save(1),Minimizer("Minuit2","minimize"),AsymptoticError(kFALSE)); //FIXME kTRUE or kFALSE

    // fix the float parameters in model except nsig
    floatPars = model->getParameters(*ddata)->selectByAttrib("Constant", false);
    iter = floatPars->createIterator();
    for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
        if (var->GetName() != TString("nsig")) {
        var->setConstant(true);
        }
    }

    while((fit_status != 0) && (tries < 10)){
        model_fit = model->fitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kFALSE)); //FIXME kTRUE or kFALSE
        fit_status = model_fit->status();
        tries++;
    }

    ss = nsig.getVal();
    dss = nsig.getError();
    tot_err = sqrt(dss * dss + ss * ss);
    delta = fabs(ss) - 2 * dmc;
    ss_cor = (delta < 0) ? 0 : ((ss > 0) ? delta : -delta);
    if (delta > 0.2 * dss) status = "Fail";

    data_npars = model_fit->floatParsFinal().getSize();
    data_ndof = bin_size-data_npars;

    // if(ss < 0) frame_data->SetMinimum(ss);
    // ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
    model->plotOn(frame_data, Name("fit"));
    chi2 = frame_data->chiSquare(data_npars);
    prob = TMath::Prob(chi2*data_ndof, data_ndof);

    // RooHist *plotdata = (RooHist*)frame_data->getObject(frame_data->numItems()-1);
    // model->plotOn(frame_data, Name("signal"), Components(signal->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
    model->plotOn(frame_data, Name("background"), Components(bkg_model->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
    RooCurve* nomBkgCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);

    // if(prob<0.01) status="Fail";
    model->SetName(Form("%s_model", bkg_fun.Data()));
    frame_data->SetTitle(Form("Pesudo data with with x%d signal, prob: %.3f", sig, prob));
    frame_data->SetXTitle("");
    frame_data->SetLabelSize(0.042, "XY");
    frame_data->SetTitleSize(0.056, "Y");
    frame_data->SetTitleOffset(0.75, "Y");
    // model->Write(model->GetName(), TObject::kOverwrite);
    output << "\t" << bkg_fun.Data() << "\tdata(Psu):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << "\tprob = " << prob << "\tstatus = " << fit_status << endl;
    output << "\t" << bkg_fun.Data() << "\tSS:\tnsig = " << ss_mc << ":" << ss << "\tdmc = " << dmc << "\tss_cor = " << ss_cor << "\tdss = " << dss << "\ttot_err = " << tot_err << "\tstatus = " << status.Data() << "\n" << endl;
    // output << "\tnbkg = " << nbkg.getVal() << "\tnbkg_err = " << nbkg.getError() << "\n" << endl;

    TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->AddEntry(frame_data->findObject("data"),"MC","ep");
    leg->AddEntry(frame_data->findObject("fit"),"Bkg + Sig","l");
    // leg->AddEntry(frame_data->findObject("signal"),"Sig","l");
    leg->AddEntry(frame_data->findObject("background"),"Bkg","l");
    frame_data->Draw();
    leg->Draw("same");

    signal->plotOn(frame_data, RooFit::Name("signal"), RooFit::Normalization(ss,RooAbsReal::NumEvent), LineColor(kRed),LineWidth(4));
    // model->plotOn(frame_data, Name("signal"), Components(signal->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
    RooCurve* nomSigCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);

    pad2->cd();
    int npoints = plotdata->GetN();
    double xtmp,ytmp;//
    int point =0;
    TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
    for (int ipoint=0; ipoint<npoints; ++ipoint) {
        plotdata->GetPoint(ipoint, xtmp,ytmp);
        double bkgval = nomBkgCurve->interpolate(xtmp);
        // if ((xtmp > 122 ) && ( xtmp < 128) ) continue;
        double errhi = plotdata->GetErrorYhigh(ipoint);
        double errlow = plotdata->GetErrorYlow(ipoint);

        std::cout << "[INFO] Channel  " << channel.Data() << " setting point " << point <<" : xtmp "<< xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp-bkgval << std::endl;
        // if(fabs(ytmp)<1e-5) continue;
        hdatasub->SetPoint(point,xtmp,ytmp-bkgval);
        hdatasub->SetPointError(point,0.,0.,errlow,errhi );
        point++;
    }

    TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
    hdummy->SetStats(0);
    hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
    hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
    hdummy->GetYaxis()->SetTitle("data - bkg PDF");
    hdummy->GetYaxis()->SetTitleOffset(0.35);
    hdummy->GetYaxis()->SetTitleSize(0.12);
    hdummy->GetYaxis()->SetLabelSize(0.09);
    hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");//bing
    hdummy->GetXaxis()->SetTitleSize(0.12);
    hdummy->GetXaxis()->SetLabelSize(0.09);
    hdummy->Draw("HIST");
    hdummy->GetYaxis()->SetNdivisions(808);

    // TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
    // line3->SetLineColor(kRed);
    // //line3->SetLineStyle(kDashed);
    // line3->SetLineWidth(5.0);
    // line3->Draw();

    hdatasub->SetMarkerStyle(8);
    hdatasub->Draw("PESAME");
    nomSigCurve->Draw("L SAME");
    canv->SaveAs(Form("./test/pesudo_data_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

    // gPad->Print(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/test/pesudo_data_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

    // RooHist *hpull = frame_data->pullHist();
    // RooPlot *frame3 = mH.frame(Title("Pull Distribution"));
    // frame3->addPlotable(hpull, "P");
    // frame3->Draw();
    // gPad->Print(Form("test/%s_cat%d_%s_%dsig.pdf",channel.Data(),cat,bkg_fun.Data(),sig));

    cout << "\t=================================" << endl;
    cout << "\n\t Finish asimov data fit\n" << endl;
}