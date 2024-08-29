#include "TCanvas.h"

#include "RooPlot.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
// #include "../interface/RooPowerLaw.h"
// #include "../interface/RooPowerLawSum.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "RooRandom.h"
//bing
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooProdPdf.h"
#include "RooNumConvPdf.h"
#include "RooGaussModel.h"
#include "RooEffProd.h"
#include "RooFormulaVar.h"
#include "RooNLLVar.h"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

using namespace std;
using namespace RooFit;
using namespace boost;

RooAbsPdf* getsigPdfPdf(RooRealVar* obs_var){
    RooRealVar* sigma = new RooRealVar("sigma","sigma",2.,0.1,5.); 
    RooRealVar* MH = new RooRealVar("MH","MH",124.5, 123., 126.); 
    RooGaussian* sig_gau = new RooGaussian("sig_gau","sig_gau",*obs_var,*MH,*sigma);
    
    RooRealVar* sigma_CB = new RooRealVar("sigma_CB","sigma_CB",0.1, 0.1, 4.); 
    RooRealVar* alpha = new RooRealVar("alpha","alpha",0.1, 0., 1.0); 
    RooRealVar* n_CB = new RooRealVar("n_CB","n_CB",5.,1.,50.); 
    RooRealVar* fracG1 = new RooRealVar("fracG1","fracG1",0.1,0.,1.0); 
    RooCBShape* CBshape = new RooCBShape("CBShape", "CBShape", *obs_var, *MH, *sigma_CB, *alpha, *n_CB);
    RooAddPdf* sigPdf = new RooAddPdf("sigPdf","sigPdf",RooArgList(*sig_gau, *CBshape),*fracG1);
    return sigPdf;
}

RooAbsPdf* getBernsteinxZGMCShape(string prefix, int cat, int order, RooRealVar* obs_var){
  //bing add ZGMCShape
  // TFile *ZGMC_file = new TFile("./ZGCoreShape_01jet.root");
  // RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  // RooAbsPdf *ZGMCShape = w->pdf(Form("CoreShape_ZG_NAF_cat%d",cat));

  TFile *ZGMC_file = new TFile("./ZGCoreShape_01jet_histSmooth.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf(Form("CoreShape_ZG_NAF_cat%d",cat));

  // TFile *ZGMC_file = new TFile("./ZGCoreShape_fromMC_01jet_v4.root");
  // RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  // RooAbsPdf *ZGMCShape = w->pdf(Form("CoreShape_MC_cat%d",cat));

  // TFile *ZGMC_file = new TFile("./ZGCoreShape_fromFake_01jet.root");
  // RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  // RooAbsPdf *ZGMCShape = w->pdf(Form("CoreShape_MC_cat%d",cat));
  w->var("CMS_hzg_mass")->setRange(105, 170);
  ZGMC_file->Close();

  RooArgList *coeffList = new RooArgList();
  map<string, RooRealVar*> params;
  map<string, RooFormulaVar*> prods;
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),.1*(i+1)*pow(-1,i),-5.,5.);
    RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*prods[name]);
  }
  if (order==1) {
	  RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern1xZG",prefix.c_str()),Form("%s_bern1xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
    //return bern;
  } else if (order==2) {
	  RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern2xZG",prefix.c_str()),Form("%s_bern2xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
    //return bern;
  } else if (order==3) {
	  RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern3xZG",prefix.c_str()),Form("%s_bern3xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
    //return bern;
  } else if (order==4) {
	  RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern4xZG",prefix.c_str()),Form("%s_bern4xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
    //return bern;
  } else if (order==5) {
	  RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern5xZG",prefix.c_str()),Form("%s_bern5xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
    //return bern;
  } else if (order==6) {
    RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern6xZG",prefix.c_str()),Form("%s_bern6xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
    // return bern;
   } else if (order==7) {
  	RooBernsteinFast<7> *bern = new RooBernsteinFast<7>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern7xZG",prefix.c_str()),Form("%s_bern7xZG",prefix.c_str()), *ZGMCShape, *bern);
    return bernZGMC;
  	// return bern;
  // } else if (order==8) {
  // 	RooBernsteinFast<8> *bern = new RooBernsteinFast<8>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  //   RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern8xZG",prefix.c_str()),Form("%s_bern8xZG",prefix.c_str()), *ZGMCShape, *bern);
  //   return bernZGMC;
  // 	// return bern;
  // } else if (order==9) {
  // 	RooBernsteinFast<9> *bern = new RooBernsteinFast<9>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  //   RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern9xZG",prefix.c_str()),Form("%s_bern9xZG",prefix.c_str()), *ZGMCShape, *bern);
  //   return bernZGMC;
  // 	// return bern;
  // } else if (order==10) {
  // 	RooBernsteinFast<10> *bern = new RooBernsteinFast<10>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  //   RooEffProd *bernZGMC = new RooEffProd(Form("%s_bern10xZG",prefix.c_str()),Form("%s_bern10xZG",prefix.c_str()), *ZGMCShape, *bern);
  //   return bernZGMC;
  // 	// return bern;
  } else {
	  return NULL;
  }

  //return bern;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}

RooAbsPdf* getExponentialZGMCShape(string prefix, int order, RooRealVar* obs_var){
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("./ZGMCShape_fromMC_01J_v4.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");
  ZGMC_file->Close();

  if(order%2==0) return NULL;
  
  double par1_exp1, par1_exp3, par3_exp3, par1_exp5, par3_exp5, par5_exp5, par1_exp7, par3_exp7, par5_exp7, par7_exp7;
  double par1_hexp1, par1_hexp3, par3_hexp3, par1_hexp5, par3_hexp5, par5_hexp5, par1_hexp7, par3_hexp7, par5_hexp7, par7_hexp7;
  double par1_lexp1, par1_lexp3, par3_lexp3, par1_lexp5, par3_lexp5, par5_lexp5, par1_lexp7, par3_lexp7, par5_lexp7, par7_lexp7;
  double coeff1_exp1, coeff1_exp3, coeff3_exp3, coeff1_exp5, coeff3_exp5, coeff5_exp5, coeff1_exp7, coeff3_exp7, coeff5_exp7, coeff7_exp7;
  double coeff1_hexp1, coeff1_hexp3, coeff3_hexp3, coeff1_hexp5, coeff3_hexp5, coeff5_hexp5, coeff1_hexp7, coeff3_hexp7, coeff5_hexp7, coeff7_hexp7;
  double coeff1_lexp1, coeff1_lexp3, coeff3_lexp3, coeff1_lexp5, coeff3_lexp5, coeff5_lexp5, coeff1_lexp7, coeff3_lexp7, coeff5_lexp7, coeff7_lexp7;

  if(order==1){
  par1_exp1 = -0.1;          par1_lexp1 = -0.5;    par1_hexp1 = 0.;
  coeff1_exp1 = 0.9;          coeff1_lexp1 = -100.;    coeff1_hexp1 = 100.;
  }
  else if(order==3){
  par1_exp3 = -0.1;      par1_lexp3 = -0.5;    par1_hexp3 = 0.;
  coeff1_exp3 = 0.9;      coeff1_lexp3 = -100.;    coeff1_hexp3 = 100.;
  par3_exp3 = -0.001;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
  coeff3_exp3 = 0.4;      coeff3_lexp3 = -100.;    coeff3_hexp3 = 100.;
  } 
  else if(order==5)
  {
  par1_exp5 = -0.1;          par1_lexp5 = -0.5;    par1_hexp5 = 0.;
  coeff1_exp5 = 0.9;       coeff1_lexp5 = -100.;    coeff1_hexp5 = 100.;
  par3_exp5 = -0.001;          par3_lexp5 = -0.5;    par3_hexp5 = 0.;
  coeff3_exp5 = 0.4; coeff3_lexp5 = -100.;    coeff3_hexp5 = 100.;
  par5_exp5 = -0.005;      par5_lexp5 = -0.7;      par5_hexp5 = 0.;
  coeff5_exp5 = 0.4; coeff5_lexp5 = -100.;    coeff5_hexp5 = 100.;
  }
  else if(order==7)
  {
  par1_exp7 = -0.1;          par1_lexp7 = -0.5;    par1_hexp7 = 0.;
  coeff1_exp7 = 0.9;       coeff1_lexp7 = -100.;    coeff1_hexp7 = 100.;
  par3_exp7 = -0.001;          par3_lexp7 = -0.5;    par3_hexp7 = 0.;
  coeff3_exp7 = 0.4; coeff3_lexp7 = -100.;    coeff3_hexp7 = 100.;
  par5_exp7 = -0.005;      par5_lexp7 = -0.7;      par5_hexp7 = 0.;
  coeff5_exp7 = 0.4; coeff5_lexp7 = -100.;    coeff5_hexp7 = 100.;
  par7_exp7 = -0.005;      par7_lexp7 = -1.0;      par7_hexp7 = 0.;
  coeff7_exp7 = 0.4; coeff7_lexp7 = -100.;    coeff7_hexp7 = 100.;
  }
    if (order==1) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp1",prefix.c_str()),Form("%s_p1_exp1",prefix.c_str()),par1_exp1,par1_lexp1,par1_hexp1);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp1",prefix.c_str()),Form("%s_cp1_exp1",prefix.c_str()),coeff1_exp1,coeff1_lexp1,coeff1_hexp1);
      RooFormulaVar *exp = new RooFormulaVar(Form("%s_exp1",prefix.c_str()),Form("%s_exp1",prefix.c_str()), "@2*TMath::Exp(@0*@1)", RooArgList(*obs_var,*p1,*cp1));
      
      RooEffProd *expZGMC = new RooEffProd(Form("%s_exp1xZG",prefix.c_str()),Form("%s_exp1xZG",prefix.c_str()), *ZGMCShape, *exp);
      return expZGMC;
  } else if (order==3) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp3",prefix.c_str()),Form("%s_p1_exp3",prefix.c_str()),par1_exp3,par1_lexp3, par1_hexp3);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp3",prefix.c_str()),Form("%s_cp1_exp3",prefix.c_str()),coeff1_exp3,coeff1_lexp3,coeff1_hexp3);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp3",prefix.c_str()),Form("%s_p3_exp3",prefix.c_str()),par3_exp3,par3_lexp3, par3_hexp3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp3",prefix.c_str()),Form("%s_cp3_exp3",prefix.c_str()),coeff3_exp3,coeff3_lexp3,coeff3_hexp3);
      RooFormulaVar *exp = new RooFormulaVar(Form("%s_exp3",prefix.c_str()),Form("%s_exp3",prefix.c_str()), "@2*TMath::Exp(@0*@1)+@4*TMath::Exp(@0*@3)", RooArgList(*obs_var,*p1,*cp1,*p3,*cp3));

      RooEffProd *expZGMC = new RooEffProd(Form("%s_exp3xZG",prefix.c_str()),Form("%s_exp3xZG",prefix.c_str()), *ZGMCShape, *exp);
      return expZGMC;
  } else if (order==5) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp5",prefix.c_str()),Form("%s_p1_exp5",prefix.c_str()),par1_exp5,par1_lexp5, par1_hexp5);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp5",prefix.c_str()),Form("%s_cp1_exp5",prefix.c_str()),coeff1_exp5,coeff1_lexp5,coeff1_hexp5);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp5",prefix.c_str()),Form("%s_p3_exp5",prefix.c_str()),par3_exp5,par3_lexp5, par3_hexp5);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp5",prefix.c_str()),Form("%s_cp3_exp5",prefix.c_str()),coeff3_exp5,coeff3_lexp5,coeff3_hexp5);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_exp5",prefix.c_str()),Form("%s_p5_exp5",prefix.c_str()),par5_exp5,par5_lexp5, par5_hexp5);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_exp5",prefix.c_str()),Form("%s_cp5_exp5",prefix.c_str()),coeff5_exp5,coeff5_lexp5,coeff5_hexp5);
      RooFormulaVar *exp = new RooFormulaVar(Form("%s_exp5",prefix.c_str()),Form("%s_exp5",prefix.c_str()), "@2*TMath::Exp(@0*@1)+@4*TMath::Exp(@0*@3)+@6*TMath::Exp(@0*@5)", RooArgList(*obs_var,*p1,*cp1,*p3,*cp3,*p5,*cp5));
      
      RooEffProd *expZGMC = new RooEffProd(Form("%s_exp5xZG",prefix.c_str()),Form("%s_exp5xZG",prefix.c_str()), *ZGMCShape, *exp);
      return expZGMC;
  } else if (order==7) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp7",prefix.c_str()),Form("%s_p1_exp7",prefix.c_str()),par1_exp7,par1_lexp7, par1_hexp7);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp7",prefix.c_str()),Form("%s_cp1_exp7",prefix.c_str()),coeff1_exp7,coeff1_lexp7,coeff1_hexp7);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp7",prefix.c_str()),Form("%s_p3_exp7",prefix.c_str()),par3_exp7,par3_lexp7, par3_hexp7);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp7",prefix.c_str()),Form("%s_cp3_exp7",prefix.c_str()),coeff3_exp7,coeff3_lexp7,coeff3_hexp7);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_exp7",prefix.c_str()),Form("%s_p5_exp7",prefix.c_str()),par5_exp7,par5_lexp7, par5_hexp7);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_exp7",prefix.c_str()),Form("%s_cp5_exp7",prefix.c_str()),coeff5_exp7,coeff5_lexp7,coeff5_hexp7);
      RooRealVar *p7 = new RooRealVar(Form("%s_p7_exp7",prefix.c_str()),Form("%s_p7_exp7",prefix.c_str()),par7_exp7,par7_lexp7, par7_hexp7);
      RooRealVar *cp7 = new RooRealVar(Form("%s_cp7_exp7",prefix.c_str()),Form("%s_cp7_exp7",prefix.c_str()),coeff7_exp7,coeff7_lexp7,coeff7_hexp7);
      RooFormulaVar *exp = new RooFormulaVar(Form("%s_exp7",prefix.c_str()),Form("%s_exp7",prefix.c_str()), "@2*TMath::Exp(@0*@1)+@4*TMath::Exp(@0*@3)+@6*TMath::Exp(@0*@5)+@8*TMath::Exp(@0*@7)", RooArgList(*obs_var,*p1,*cp1,*p3,*cp3,*p5,*cp5,*p7,*cp7));
      
      RooEffProd *expZGMC = new RooEffProd(Form("%s_exp7xZG",prefix.c_str()),Form("%s_exp7xZG",prefix.c_str()), *ZGMCShape, *exp);
      return expZGMC;
  }
  else {
    return NULL;
  }
}

RooAbsPdf* getPowerLawZGMCShape(string prefix, int order, RooRealVar* obs_var){
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("./ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");
  ZGMC_file->Close();

  if(order%2==0) return NULL;
  double par1_pow1, par1_pow3, par3_pow3, par1_pow5, par3_pow5, par5_pow5, par1_pow7, par3_pow7, par5_pow7, par7_pow7;
  double par1_hpow1, par1_hpow3, par3_hpow3, par1_hpow5, par3_hpow5, par5_hpow5, par1_hpow7, par3_hpow7, par5_hpow7, par7_hpow7;
  double par1_lpow1, par1_lpow3, par3_lpow3, par1_lpow5, par3_lpow5, par5_lpow5, par1_lpow7, par3_lpow7, par5_lpow7, par7_lpow7;
  double coeff1_pow1, coeff1_pow3, coeff3_pow3, coeff1_pow5, coeff3_pow5, coeff5_pow5, coeff1_pow7, coeff3_pow7, coeff5_pow7, coeff7_pow7;
  double coeff1_hpow1, coeff1_hpow3, coeff3_hpow3, coeff1_hpow5, coeff3_hpow5, coeff5_hpow5, coeff1_hpow7, coeff3_hpow7, coeff5_hpow7, coeff7_hpow7;
  double coeff1_lpow1, coeff1_lpow3, coeff3_lpow3, coeff1_lpow5, coeff3_lpow5, coeff5_lpow5, coeff1_lpow7, coeff3_lpow7, coeff5_lpow7, coeff7_lpow7;

  if(order==1){
  par1_pow1 = -1.;          par1_lpow1 = -10.;    par1_hpow1 = 5.;
  coeff1_pow1 = 1.; coeff1_lpow1 = -1.;    coeff1_hpow1 = 1.;
  }
  else if(order==3){
  par1_pow3 = 2.7;      par1_lpow3 = -10.;    par1_hpow3 = 5.;
  coeff1_pow3 = 0.0001; coeff1_lpow3 = -1.;    coeff1_hpow3 = 1.;
  par3_pow3 = -3.;      par3_lpow3 = -10.;      par3_hpow3 = 5.;
  coeff3_pow3 = 0.99; coeff3_lpow3 = -1;    coeff3_hpow3 = 1.;
  } 
  else if(order==5)
  {
  par1_pow5 = -1.;          par1_lpow5 = -10.;    par1_hpow5 = 5.;
  coeff1_pow5 = 1.;       coeff1_lpow5 = -1.;    coeff1_hpow5 = 1.;
  par3_pow5 = -1.;          par3_lpow5 = -10.;    par3_hpow5 = 5.;
  coeff3_pow5 = 1.; coeff3_lpow5 = -1.;    coeff3_hpow5 = 1.;
  par5_pow5 = -1.;      par5_lpow5 = -10.;      par5_hpow5 = 5.;
  coeff5_pow5 = 1.; coeff5_lpow5 = -1.;    coeff5_hpow5 = 1.;
  }
  else if(order==7)
  {
  par1_pow7 = -1.;          par1_lpow7 = -10.;    par1_hpow7 = 5.;
  coeff1_pow7 = 1.;       coeff1_lpow7 = -10.;    coeff1_hpow7 = 10.;
  par3_pow7 = -1.;          par3_lpow7 = -10.;    par3_hpow7 = 5.;
  coeff3_pow7 = 1.; coeff3_lpow7 = -10.;    coeff3_hpow7 = 10.;
  par5_pow7 = -1.;      par5_lpow7 = -10.;      par5_hpow7 = 5.;
  coeff5_pow7 = 1.; coeff5_lpow7 = -10.;    coeff5_hpow7 = 10.;
  par7_pow7 = -1.;      par7_lpow7 = -10.;      par7_hpow7 = 5.;
  coeff7_pow7 = 1.; coeff7_lpow7 = -10.;    coeff7_hpow7 = 10.;
  }
  
    if (order==1) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow1",prefix.c_str()),Form("%s_p1_pow1",prefix.c_str()),par1_pow1,par1_lpow1,par1_hpow1);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow1",prefix.c_str()),Form("%s_cp1_pow1",prefix.c_str()),coeff1_pow1,coeff1_lpow1,coeff1_hpow1);
    	RooGenericPdf *pows = new RooGenericPdf(Form("%s_pow1",prefix.c_str()),Form("%s_pow1",prefix.c_str()), "@2*(@0)^(@1)", RooArgList(*obs_var,*p1,*cp1));

      RooEffProd *powZGMC = new RooEffProd(Form("%s_pow1xZG",prefix.c_str()),Form("%s_pow1xZG",prefix.c_str()), *ZGMCShape, *pows);
      return powZGMC;
  } else if (order==3) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow3",prefix.c_str()),Form("%s_p1_pow3",prefix.c_str()),par1_pow3,par1_lpow3,par1_hpow3);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow3",prefix.c_str()),Form("%s_cp1_pow3",prefix.c_str()),coeff1_pow3,coeff1_lpow3,coeff1_hpow3);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_pow3",prefix.c_str()),Form("%s_p3_pow3",prefix.c_str()),par3_pow3,par3_lpow3,par3_hpow3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_pow3",prefix.c_str()),Form("%s_cp3_pow3",prefix.c_str()),coeff3_pow3,coeff3_lpow3,coeff3_hpow3);
    	RooGenericPdf *pows = new RooGenericPdf(Form("%s_pow3",prefix.c_str()),Form("%s_pow3",prefix.c_str()), "@2*(@0)^(@1)+@4*(@0)^(@3)", RooArgList(*obs_var,*p1,*cp1,*p3,*cp3));

      RooEffProd *powZGMC = new RooEffProd(Form("%s_pow3xZG",prefix.c_str()),Form("%s_pow3xZG",prefix.c_str()), *ZGMCShape, *pows);
      return powZGMC;
  } else if (order==5) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow5",prefix.c_str()),Form("%s_p1_pow5",prefix.c_str()),par1_pow5,par1_lpow5,par1_hpow5);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow5",prefix.c_str()),Form("%s_cp1_pow5",prefix.c_str()),coeff1_pow5,coeff1_lpow5,coeff1_hpow5);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_pow5",prefix.c_str()),Form("%s_p3_pow5",prefix.c_str()),par3_pow5,par3_lpow5,par3_hpow5);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_pow5",prefix.c_str()),Form("%s_cp3_pow5",prefix.c_str()),coeff3_pow5,coeff3_lpow5,coeff3_hpow5);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_pow5",prefix.c_str()),Form("%s_p5_pow5",prefix.c_str()),par5_pow5,par5_lpow5,par5_hpow5);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_pow5",prefix.c_str()),Form("%s_cp5_pow5",prefix.c_str()),coeff5_pow5,coeff5_lpow5,coeff5_hpow5);
    	RooGenericPdf *pows = new RooGenericPdf(Form("%s_pow5",prefix.c_str()),Form("%s_pow5",prefix.c_str()), "@2*(@0)^(@1)+@4*(@0)^(@3)+@6*(@0)^(@5)", RooArgList(*obs_var,*p1,*cp1,*p3,*cp3,*p5,*cp5));
      
      RooEffProd *powZGMC = new RooEffProd(Form("%s_pow5xZG",prefix.c_str()),Form("%s_pow5xZG",prefix.c_str()), *ZGMCShape, *pows);
      return powZGMC;
  } else if (order==7) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow7",prefix.c_str()),Form("%s_p1_pow7",prefix.c_str()),par1_pow7,par1_lpow7,par1_hpow7);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow7",prefix.c_str()),Form("%s_cp1_pow7",prefix.c_str()),coeff1_pow7,coeff1_lpow7,coeff1_hpow7);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_pow7",prefix.c_str()),Form("%s_p3_pow7",prefix.c_str()),par3_pow7,par3_lpow7,par3_hpow7);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_pow7",prefix.c_str()),Form("%s_cp3_pow7",prefix.c_str()),coeff3_pow7,coeff3_lpow7,coeff3_hpow7);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_pow7",prefix.c_str()),Form("%s_p5_pow7",prefix.c_str()),par5_pow7,par5_lpow7,par5_hpow7);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_pow7",prefix.c_str()),Form("%s_cp5_pow7",prefix.c_str()),coeff5_pow7,coeff5_lpow7,coeff5_hpow7);
      RooRealVar *p7 = new RooRealVar(Form("%s_p7_pow7",prefix.c_str()),Form("%s_p7_pow7",prefix.c_str()),par7_pow7,par7_lpow7,par7_hpow7);
      RooRealVar *cp7 = new RooRealVar(Form("%s_cp7_pow7",prefix.c_str()),Form("%s_cp7_pow7",prefix.c_str()),coeff7_pow7,coeff7_lpow7,coeff7_hpow7);
    	RooGenericPdf *pows = new RooGenericPdf(Form("%s_pow7",prefix.c_str()),Form("%s_pow7",prefix.c_str()), "@2*(@0)^(@1)+@4*(@0)^(@3)+@6*(@0)^(@5)+@8*(@0)^(@7)", RooArgList(*obs_var,*p1,*cp1,*p3,*cp3,*p5,*cp5,*p7,*cp7));
      
      RooEffProd *powZGMC = new RooEffProd(Form("%s_pow7xZG",prefix.c_str()),Form("%s_pow7xZG",prefix.c_str()), *ZGMCShape, *pows);
      return powZGMC;
  } 
   else {
	  return NULL;
  }
}

RooAbsPdf* getLaurentZGMCShape(string prefix, int order, RooRealVar* obs_var){
 
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("./ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");
  ZGMC_file->Close();

  double coeff1_lau1, coeff2_lau1, coeff1_lau2, coeff2_lau2, coeff3_lau2, coeff1_lau3, coeff2_lau3, coeff3_lau3, coeff4_lau3, coeff1_lau4, coeff2_lau4, coeff3_lau4, coeff4_lau4, coeff5_lau4, coeff1_lau5, coeff2_lau5, coeff3_lau5, coeff4_lau5, coeff5_lau5, coeff6_lau5; 
  double coeff1_hlau1, coeff2_hlau1, coeff1_hlau2, coeff2_hlau2, coeff3_hlau2, coeff1_hlau3, coeff2_hlau3, coeff3_hlau3, coeff4_hlau3, coeff1_hlau4, coeff2_hlau4, coeff3_hlau4, coeff4_hlau4, coeff5_hlau4, coeff1_hlau5, coeff2_hlau5, coeff3_hlau5, coeff4_hlau5, coeff5_hlau5, coeff6_hlau5; 
  double coeff1_llau1, coeff2_llau1, coeff1_llau2, coeff2_llau2, coeff3_llau2, coeff1_llau3, coeff2_llau3, coeff3_llau3, coeff4_llau3, coeff1_llau4, coeff2_llau4, coeff3_llau4, coeff4_llau4, coeff5_llau4, coeff1_llau5, coeff2_llau5, coeff3_llau5, coeff4_llau5, coeff5_llau5, coeff6_llau5;  

  if(order==1){
  coeff1_lau1 = 0.1; coeff1_llau1 = -100.;    coeff1_hlau1 = 100.;
  coeff2_lau1 = 0.5; coeff2_llau1 = -100.;    coeff2_hlau1 = 100.;
  }
  else if(order==2){
  coeff1_lau2 = 0.1; coeff1_llau2 = -100.;    coeff1_hlau2 = 100.;
  coeff2_lau2 = 0.5; coeff2_llau2 = -100.;    coeff2_hlau2 = 100.;
  coeff3_lau2 = 0.1; coeff3_llau2 = -100.;    coeff3_hlau2 = 100.;
  } 
  else if(order==3)
  {
  coeff1_lau3 = 0.1; coeff1_llau3 = -100.;    coeff1_hlau3 = 100.;
  coeff2_lau3 = 0.5; coeff2_llau3 = -100.;    coeff2_hlau3 = 100.;
  coeff3_lau3 = 0.1; coeff3_llau3 = -100.;    coeff3_hlau3 = 100.;
  coeff4_lau3 = 0.5; coeff4_llau3 = -100.;    coeff4_hlau3 =100.;
  }
  else if(order==4)
  {
  coeff1_lau4 = 0.1; coeff1_llau4 = -100.;    coeff1_hlau4 = 100.;
  coeff2_lau4 = 0.5; coeff2_llau4 = -100.;    coeff2_hlau4 = 100.;
  coeff3_lau4 = 0.1; coeff3_llau4 = -100.;    coeff3_hlau4 = 100.;
  coeff4_lau4 = 0.5; coeff4_llau4 = -100.;    coeff4_hlau4 =100.;
  coeff5_lau4 = 0.1; coeff5_llau4 = -100.;    coeff5_hlau4 =100.;
  }
  else if(order==5)
  {
  coeff1_lau5 = 0.1; coeff1_llau5 = -100.;    coeff1_hlau5 = 100.;
  coeff2_lau5 = 0.5; coeff2_llau5 = -100.;    coeff2_hlau5 = 100.;
  coeff3_lau5 = 0.1; coeff3_llau5 = -100.;    coeff3_hlau5 = 100.;
  coeff4_lau5 = 0.5; coeff4_llau5 = -100.;    coeff4_hlau5 =100.;
  coeff5_lau5 = 0.1; coeff5_llau5 = -100.;    coeff5_hlau5 =100.;
  coeff6_lau5 = 0.5; coeff6_llau5 = -100.;    coeff6_hlau5 =100.;
  }

    if (order==1) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau1",prefix.c_str()),Form("%s_cp1_lau1",prefix.c_str()),coeff1_lau1,coeff1_llau1,coeff1_hlau1);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau1",prefix.c_str()),Form("%s_cp2_lau1",prefix.c_str()),coeff2_lau1,coeff2_llau1,coeff2_hlau1);
      RooGenericPdf *lau= new RooGenericPdf(Form("%s_lau1",prefix.c_str()),Form("%s_lau1",prefix.c_str()), "@1*(@0)^(-4)+@2*(@0)^(-5)", RooArgList(*obs_var,*cp1,*cp2));

      RooEffProd *lauZGMC = new RooEffProd(Form("%s_lau1xZG",prefix.c_str()),Form("%s_lau1xZG",prefix.c_str()), *ZGMCShape, *lau);
      return lauZGMC;
  } else if (order==2) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau2",prefix.c_str()),Form("%s_cp1_lau2",prefix.c_str()),coeff1_lau2,coeff1_llau2,coeff1_hlau2);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp1_lau2",prefix.c_str()),Form("%s_cp1_lau2",prefix.c_str()),coeff2_lau2,coeff2_llau2,coeff2_hlau2);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau2",prefix.c_str()),Form("%s_cp3_lau2",prefix.c_str()),coeff3_lau2,coeff3_llau2,coeff3_hlau2);
      RooGenericPdf *lau = new RooGenericPdf(Form("%s_lau2",prefix.c_str()),Form("%s_lau2",prefix.c_str()), "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)", RooArgList(*obs_var,*cp1,*cp2,*cp3));

      RooEffProd *lauZGMC = new RooEffProd(Form("%s_lau2xZG",prefix.c_str()),Form("%s_lau2xZG",prefix.c_str()), *ZGMCShape, *lau);
      return lauZGMC;
  } else if (order==3) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau3",prefix.c_str()),Form("%s_cp1_lau3",prefix.c_str()),coeff1_lau3,coeff1_llau3,coeff1_hlau3);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau3",prefix.c_str()),Form("%s_cp2_lau3",prefix.c_str()),coeff2_lau3,coeff2_llau3,coeff2_hlau3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau3",prefix.c_str()),Form("%s_cp3_lau3",prefix.c_str()),coeff3_lau3,coeff3_llau3,coeff3_hlau3);
      RooRealVar *cp4 = new RooRealVar(Form("%s_cp4_lau3",prefix.c_str()),Form("%s_cp4_lau3",prefix.c_str()),coeff4_lau3,coeff4_llau3,coeff4_hlau3);
      RooGenericPdf *lau = new RooGenericPdf(Form("%s_lau3",prefix.c_str()),Form("%s_lau3",prefix.c_str()), "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)+@4*(@0)^(-6)", RooArgList(*obs_var,*cp1,*cp2,*cp3,*cp4));

      RooEffProd *lauZGMC = new RooEffProd(Form("%s_lau3xZG",prefix.c_str()),Form("%s_lau3xZG",prefix.c_str()), *ZGMCShape, *lau);
      return lauZGMC;
  } else if (order==4) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau4",prefix.c_str()),Form("%s_cp1_lau4",prefix.c_str()),coeff1_lau4,coeff1_llau4,coeff1_hlau4);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau4",prefix.c_str()),Form("%s_cp2_lau4",prefix.c_str()),coeff2_lau4,coeff2_llau4,coeff2_hlau4);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau4",prefix.c_str()),Form("%s_cp3_lau4",prefix.c_str()),coeff3_lau4,coeff3_llau4,coeff3_hlau4);
      RooRealVar *cp4 = new RooRealVar(Form("%s_cp4_lau4",prefix.c_str()),Form("%s_cp4_lau4",prefix.c_str()),coeff4_lau4,coeff4_llau4,coeff4_hlau4);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_lau4",prefix.c_str()),Form("%s_cp5_lau4",prefix.c_str()),coeff5_lau4,coeff5_llau4,coeff5_hlau4);
      RooGenericPdf *lau = new RooGenericPdf(Form("%s_lau4",prefix.c_str()),Form("%s_lau4",prefix.c_str()), "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)+@4*(@0)^(-6)+@5*(@0)^(-2)", RooArgList(*obs_var,*cp1,*cp2,*cp3,*cp4,*cp5));

      RooEffProd *lauZGMC = new RooEffProd(Form("%s_lau4xZG",prefix.c_str()),Form("%s_lau4xZG",prefix.c_str()), *ZGMCShape, *lau);
      return lauZGMC;
  } else if (order==5) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau5",prefix.c_str()),Form("%s_cp1_lau5",prefix.c_str()),coeff1_lau5,coeff1_llau5,coeff1_hlau5);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau5",prefix.c_str()),Form("%s_cp2_lau5",prefix.c_str()),coeff2_lau5,coeff2_llau5,coeff2_hlau5);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau5",prefix.c_str()),Form("%s_cp3_lau5",prefix.c_str()),coeff3_lau5,coeff3_llau5,coeff3_hlau5);
      RooRealVar *cp4 = new RooRealVar(Form("%s_cp4_lau5",prefix.c_str()),Form("%s_cp4_lau5",prefix.c_str()),coeff4_lau5,coeff4_llau5,coeff4_hlau5);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_lau5",prefix.c_str()),Form("%s_cp5_lau5",prefix.c_str()),coeff5_lau5,coeff5_llau5,coeff5_hlau5);
      RooRealVar *cp6 = new RooRealVar(Form("%s_cp6_lau5",prefix.c_str()),Form("%s_cp6_lau5",prefix.c_str()),coeff6_lau5,coeff6_llau5,coeff6_hlau5);
      RooGenericPdf *lau = new RooGenericPdf(Form("%s_lau5",prefix.c_str()),Form("%s_lau5",prefix.c_str()), "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)+@4*(@0)^(-6)+@5*(@0)^(-2)+@6*(@0)^(-7)", RooArgList(*obs_var,*cp1,*cp2,*cp3,*cp4,*cp5,*cp6));

      RooEffProd *lauZGMC = new RooEffProd(Form("%s_lau5xZG",prefix.c_str()),Form("%s_lau5xZG",prefix.c_str()), *ZGMCShape, *lau);
      return lauZGMC;
  } 
   else {
	  return NULL;
  }
}

RooAbsPdf* getPdf(string type, int cat, int order, RooRealVar* obs_var, const char* ext=""){
  if (type=="Bernstein") return getBernsteinxZGMCShape(Form("bern%d",order),cat, order,obs_var);
  // else if (type=="Chebychev") return getChebychev(Form("cheb%d",order),order,*obs_var);
  else if (type=="Exponential") return getExponentialZGMCShape(Form("exp%d",order),order,obs_var);
  else if (type=="PowerLaw") return getPowerLawZGMCShape(Form("pow%d",order),order,obs_var);
  else if (type=="Laurent") return getLaurentZGMCShape(Form("lau%d",order),order,obs_var);
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void SSTest_core_function(int cat = 0, int sig = 0, TString channel = "zero_to_one_jet"){
  double mgg_low = 105, mgg_high = 170, bin_size = 4;

  //background MC template
  TH1F* hbkg;
  TFile* fbkg = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/NAF_SSTest_template.root"));
  if (fbkg->GetListOfKeys()->Contains((Form("bkg_%s_cat%d",channel.Data(), cat))))
      hbkg = (TH1F*)fbkg->Get(Form("bkg_%s_cat%d",channel.Data(), cat));
  else
      abort();
  // TFile* fbkg = TFile::Open(Form("./bkg_template_v3_cut2/bkg/bkg_0sig_cat%d.root", cat));
  // hbkg = (TH1F*)fbkg->Get(Form("mass_cat%d", cat));
  double dataevents = hbkg->Integral();
  double mcsbevents = hbkg->Integral(0,(122-mgg_low)*bin_size)+hbkg->Integral(bin_size*(mgg_high-128), bin_size*mgg_high);

  //sigPdf MC
  TFile* fsig = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  TH1F* hsig = (TH1F*)( fsig->Get(Form("sig_%s_cat%d", channel.Data(), cat)) );
  // TFile* fsig = TFile::Open(Form("./bkg_template_v3_cut2/sig/FullSig_cat%d.root", cat));
  // TH1F* hsig = (TH1F*)fsig->Get(Form("mass_cat%d_FullSig", cat));
  double sigevents = hsig->Integral();

  //data side band
  TFile* fsb = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  TH1F* hsb = (TH1F*)( fsb->Get(Form("data_%s_cat%d", channel.Data(), cat)) );
  // TFile* fsb = TFile::Open(Form("./bkg_template_v3_cut2/data_sid/data_0sig_cat%d.root", cat));
  // TH1F* hsb = (TH1F*)fsb->Get(Form("mass_cat%d_data", cat));
  double sbevents = hsb->Integral();

  hbkg->Scale(sbevents/mcsbevents);

  cout << "\n\tFinish preparing data!!!\n" << endl;

  vector<string> functionClasses;
  functionClasses.push_back("Bernstein");
  // functionClasses.push_back("Exponential");
  // functionClasses.push_back("PowerLaw");
  // functionClasses.push_back("Laurent");

  //initializing
  RooRealVar* CMS_hzg_mass = new RooRealVar("CMS_hzg_mass", "CMS_hzg_mass", 125.38, mgg_low, mgg_high);
  // CMS_hzg_mass->setMin( mgg_low. );
  // CMS_hzg_mass->setMax( hbkg->GetBinCenter(hbkg->GetNbinsX()) + hbkg->GetBinWidth(hbkg->GetNbinsX()) );

  //background function fit
  RooDataHist* dbkg = new RooDataHist("bkg","dataset with x", *CMS_hzg_mass, hbkg);
  RooDataHist* dsb = new RooDataHist("data_sb","dataset with x", *CMS_hzg_mass, hsb);
  // cout<<"bkg hist integral "<<hbkg->Integral()<<" "<<dbkg->sumEntries()<<endl;

  RooRealVar nsig("nsig","nsig",0,-50*sigevents,50*sigevents);
  RooRealVar nbkg("nbkg","nbkg",dataevents, 0.1*dataevents, 2*dataevents);

  for (vector<string>::iterator funcType=functionClasses.begin(); funcType!=functionClasses.end(); funcType++){
    int order = 1;
    TString bkg_fun;
    ofstream output(Form("./outputs/%s_%d_%dxsig.txt", channel.Data(), cat, sig), ofstream::app);
    // TFile *f = new TFile(Form("./outputs/%s_%d_%dxsig.root", channel.Data(), cat, sig),"UPDATE"); 
    while (order < 8){ //FIXME
      RooAbsPdf *bkgPdf = getPdf(*funcType,cat,order,CMS_hzg_mass,"");
      bkg_fun = Form("%s%d", funcType->c_str(), order);
      order++;
      int flag, fit_status, tries; TString status;
      double dmc, dss, ss, ss_mc, tot_err, ss_cor, delta, chi2, prob;
      RooFitResult *bkgPdf_fit, *bkg_model_fit;
      RooAbsPdf *bkg_model = bkgPdf;
      RooNLLVar nllt;
      RooAbsCollection *floatPars; TIterator *iter;
      if (bkgPdf){
        flag = 1; status = "Pass";
        // Data side band fitting
        // hbkg->Write(hbkg->GetName(), TObject::kOverwrite);
        // hsb->Write(hsb->GetName(), TObject::kOverwrite);
        // hsig->Write(hsig->GetName(), TObject::kOverwrite);

        int bkg_npars, bkg_ndof;
        double nll;
        RooPlot *frame_bkg;

        // cout << *bkgPdf << " " << *dsb << endl;

        CMS_hzg_mass->setRange("range_low",mgg_low,122);
        CMS_hzg_mass->setRange("signal",122,128);
        CMS_hzg_mass->setRange("range_high",128,mgg_high);
        CMS_hzg_mass->setRange("FULL",mgg_low,mgg_high); 

        // RooRealVar N("N", "Extended term", 0, 200000);
        // RooExtendPdf extmodel("extmodel", "Extended model", *bkgPdf, N, "FULL");

        //background function fit

        // bkg_model_fit = extmodel.fitTo(*dsb,Range("range"),SplitRange(true),Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kTRUE));

        bkg_model_fit = bkg_model->fitTo(*dsb,Range("range"),Save(1),Minimizer("Minuit2","minimize"),SumW2Error(kFALSE));

        fit_status = bkg_model_fit->status();
        bkg_npars = bkg_model_fit->floatParsFinal().getSize();
        frame_bkg = CMS_hzg_mass->frame(Title(Form("Data side band with %s pdf", bkg_fun.Data())));
        bkg_ndof = bin_size*(mgg_high-mgg_low-128+122)-bkg_npars;
        dsb->plotOn(frame_bkg, Cut("CMS_hzg_mass>128 | CMS_hzg_mass<122"));
        bkg_model->plotOn(frame_bkg,NormRange("range_low,range_high"));
        // bkg_model->plotOn(frame_bkg, NormRange("FULL"));
        chi2 = frame_bkg->chiSquare(bkg_npars);

        // extmodel.SetName(bkg_fun);
        // extmodel.Write(bkg_model->GetName(), TObject::kOverwrite);
        
        bkg_model->SetName(bkg_fun);
        bkg_model->Write(bkg_model->GetName(), TObject::kOverwrite);
        
        nll = bkg_model_fit->minNll();
        prob = TMath::Prob(chi2*bkg_ndof, bkg_ndof);
        // if(prob<0.05) status = "Fail";
        bkgPdf->Write(bkgPdf->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tsb:\tnpars = " << bkg_npars << " \tchi^2 = " << chi2 << "\tprob = " << prob << "\tfitting status = " << fit_status << endl;
        bkgPdf->paramOn(frame_bkg, RooFit::Layout(0.55,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        frame_bkg->getAttText()->SetTextSize(0.03);
        frame_bkg->Draw();
        gPad->Print(Form("./test/data_sb_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

        cout << "\t=================================" << endl;
        cout << "\n\t Finish data side band fit\n" << endl;

        // MC background fitting
        bkgPdf_fit = bkgPdf->fitTo(*dbkg,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));
        fit_status = bkgPdf_fit->status();

        bkg_npars = bkgPdf_fit->floatParsFinal().getSize();
        bkg_ndof = bin_size*(mgg_high-mgg_low)-bkg_npars;
        frame_bkg = CMS_hzg_mass->frame(Title(Form("Background with %s pdf", bkg_fun.Data())));
        
        TCanvas *canv; 
        TPad *pad1; TPad *pad2;
        canv = new TCanvas();
        pad1 = new TPad("pad1","pad1",0,0.25,1,1);
        pad2 = new TPad("pad2","pad2",0,0,1,0.35);
        pad1->SetBottomMargin(0.18);
        // pad2->SetTopMargin(0.00001);
        pad2->SetBottomMargin(0.25);
        pad1->Draw();
        pad2->Draw();
        pad1->cd();

        dbkg->plotOn(frame_bkg);
        RooHist *plotdata;
        plotdata = (RooHist*)frame_bkg->getObject(frame_bkg->numItems()-1);
        bkgPdf->plotOn(frame_bkg);
        RooCurve* nomBkgCurve;
        nomBkgCurve = (RooCurve*)frame_bkg->getObject(frame_bkg->numItems()-1);
        bkgPdf->paramOn(frame_bkg, RooFit::Layout(0.55,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        frame_bkg->getAttText()->SetTextSize(0.03);
        bkgPdf->SetName(Form("%s_model", bkg_fun.Data()));
        chi2 = frame_bkg->chiSquare(bkg_npars);
        prob = TMath::Prob(chi2*bkg_ndof, bkg_ndof);
        // if (prob < 0.05) status = "Fail";
        // bkgPdf->Write(bkgPdf->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tbkg:\tnpars = " << bkg_npars << " \tchi^2 = " << chi2 << "\tprob = " << prob << "\tfitting status = " << fit_status << endl;
        frame_bkg->Draw();

        pad2->cd();
        int npoints;
        npoints = plotdata->GetN();
        double xtmp,ytmp;//
        int point =0;
        TGraphAsymmErrors *hdatasub;
        hdatasub = new TGraphAsymmErrors(npoints);
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

        TH1 *hdummy;
        hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
        hdummy->SetStats(0);
        hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
        hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
        hdummy->GetYaxis()->SetTitle("data - fit PDF");
        hdummy->GetYaxis()->SetTitleOffset(0.35);
        hdummy->GetYaxis()->SetTitleSize(0.12);
        hdummy->GetYaxis()->SetLabelSize(0.09);
        hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");//bing
        hdummy->GetXaxis()->SetTitleSize(0.12);
        hdummy->GetXaxis()->SetLabelSize(0.09);
        hdummy->Draw("HIST");
        hdummy->GetYaxis()->SetNdivisions(808);

        TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
        line3->SetLineColor(kBlue);
        //line3->SetLineStyle(kDashed);
        line3->SetLineWidth(5.0);
        line3->Draw();
        hdatasub->SetMarkerStyle(8);
        hdatasub->Draw("PESAME");

        canv->SaveAs(Form("./test/mc_bkg_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));
        pad1->Close(); pad2->Close();

        cout << "\t=================================\n";
        cout << "\n\tFinish background function fit\n" << endl;

        RooAbsPdf* sigPdf = getsigPdfPdf(CMS_hzg_mass);
        RooFitResult *sigPdf_fit;
        RooDataHist* dsig = new RooDataHist("sig","dataset with x", *CMS_hzg_mass, hsig);
        sigPdf_fit = sigPdf->fitTo(*dsig,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE
        fit_status = sigPdf_fit->status();

        RooAbsCollection *m_sigParameters = sigPdf->getParameters(RooArgSet())->selectByAttrib("Constant", false);
        TIterator *sigit = m_sigParameters->createIterator();

        int sig_npars = sigPdf_fit->floatParsFinal().getSize();
        int sig_ndof = 2*bin_size*(mgg_high-mgg_low)-sig_npars;
        RooPlot *frame_sig = CMS_hzg_mass->frame(Title("sigPdf with distorted Gaussian pdf"));
        dsig->plotOn(frame_sig, DataError(RooAbsData::SumW2));
        sigPdf->plotOn(frame_sig);
        sigPdf->paramOn(frame_sig, RooFit::Layout(0.55,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        // sigPdf->Write(sigPdf->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tsig:\tnpars = " << sig_npars << "\tchi^2 = " << frame_sig->chiSquare(sig_npars) << "\tprob = " << TMath::Prob(frame_sig->chiSquare(sig_npars)*sig_ndof, sig_ndof) << "\tfitting status = " << fit_status << endl;
        frame_sig->Draw();
        gPad->Print(Form("./test/sigPdf_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));
        for (RooRealVar *p = (RooRealVar *)sigit->Next(); p != 0; p = (RooRealVar *)sigit->Next()) p->setConstant(kTRUE);

        cout << "\t=================================" << endl;
        cout << "\n\t Finish sigPdf function fit\n" << endl;

        //inject signal
        TH1F* hdata = (TH1F*)hbkg->Clone();
        hdata->Add(hsig, sig);
        hdata->SetName("asimov data");
        // hdata->Write(hdata->GetName(), TObject::kOverwrite);
        RooDataHist* ddata = new RooDataHist("data_bin","dataset with x", *CMS_hzg_mass, hdata);

        // pesudo dataset fit
        RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*sigPdf, *bkgPdf),RooArgList(nsig,nbkg));
        RooFitResult *model_fit;
        int data_npars, data_ndof;
        canv = new TCanvas();
        
        RooPlot *frame_data, *frame_data_trash;
        frame_data_trash = CMS_hzg_mass->frame();
        frame_data = CMS_hzg_mass->frame();
        pad1 = new TPad("pad1","pad1",0,0.25,1,1);
        pad2 = new TPad("pad2","pad2",0,0,1,0.35);
        pad1->SetBottomMargin(0.18);
        // pad2->SetTopMargin(0.00001);
        pad2->SetBottomMargin(0.25);
        pad1->Draw();
        pad2->Draw();
        pad1->cd();

        fit_status = -1;
        tries = 0;
        while((fit_status != 0) && (tries < 10)){
          model_fit = model->chi2FitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),DataError(RooAbsData::SumW2)); //FIXME kTRUE or kFALSE
          fit_status = model_fit->status();
          data_npars = model_fit->floatParsFinal().getSize();
          ss_mc = nsig.getVal();
          output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars << "\tss = " << ss_mc << "\tfitting status = " << fit_status << endl;
          tries++;
        }

        // fix the float parameters in model except nsig
        floatPars = model->getParameters(*ddata)->selectByAttrib("Constant", false);
        iter = floatPars->createIterator();
        for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
          if (var->GetName() != TString("nsig")) {
            var->setConstant(true);
          }
        }

        // RooChi2Var chi2Fit("chi2", "chi2", *model, *ddata);
        // RooMinuit minuit(chi2Fit);
        // minuit.migrad();
        // minuit.hesse();  // Calculate the hesse matrix
        // fit_status = minuit.save()->status();
        // data_npars = minuit.save()->floatParsFinal().getSize();

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
        data_ndof = bin_size*(mgg_high-mgg_low)-data_npars;

        // if(ss < 0) frame_data->SetMinimum(ss);
        ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
        plotdata = (RooHist*)frame_data->getObject(frame_data->numItems()-1);
        ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
        plotdata = (RooHist*)frame_data_trash->getObject(frame_data_trash->numItems()-1);
        model->plotOn(frame_data_trash, Name("fit"));
        chi2 = frame_data_trash->chiSquare(data_npars);
        prob = TMath::Prob(chi2*data_ndof, data_ndof);
        output << "\t" << bkg_fun.Data() << "\tdata(MC):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << ": " << chi2 << "\tprob = " << prob << "\tfitting status = " << fit_status << endl;
        
        // unfix the float parameters in model
        iter = floatPars->createIterator();
        for (RooRealVar* var = (RooRealVar*)iter->Next(); var != nullptr; var = (RooRealVar*)iter->Next()) {
          var->setConstant(false);
        }

        fit_status = -1;
        tries = 0;
        while((fit_status != 0) && (tries < 10)){
          model_fit = model->fitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kFALSE)); //FIXME kTRUE or kFALSE
          fit_status = model_fit->status();
          data_npars = model_fit->floatParsFinal().getSize();
          tries++;
        }

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
        data_ndof = bin_size*(mgg_high-mgg_low)-data_npars;
        tot_err = sqrt(dss*dss+ss*ss);
        delta = abs(ss)-2*dmc;
        if (delta<0) ss_cor = 0;
        else{
          if (ss>0) ss_cor = delta;
          else ss_cor = -1. * delta;
        }
        if (delta > 0.2 * dss) status = "Fail";

        // ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
        // RooHist *plotdata = (RooHist*)frame_data->getObject(frame_data->numItems()-1);
        model->plotOn(frame_data, Name("fit"));
        RooCurve* nomSumCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);
        chi2 = frame_data->chiSquare(data_npars);
        prob = TMath::Prob(chi2*data_ndof, data_ndof);
        model->plotOn(frame_data, Name("background"), Components(bkgPdf->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
        nomBkgCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);

        model->SetName(Form("%s_model", bkg_fun.Data()));
        frame_data->SetTitle(Form("Pesudo data with with x%d signal, prob: %.3f", sig, prob));
        frame_data->SetXTitle("");
        frame_data->SetLabelSize(0.042, "XY");
        frame_data->SetTitleSize(0.056, "Y");
        frame_data->SetTitleOffset(0.75, "Y");
        frame_data->Draw();
        
        sigPdf->plotOn(frame_data, RooFit::Name("signal"), RooFit::Normalization(ss,RooAbsReal::NumEvent), LineColor(kRed),LineWidth(4));
        // model->plotOn(frame_data, Name("signal"), Components(sigPdf->GetName()), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
        RooCurve* nomSigCurve = (RooCurve*)frame_data->getObject(frame_data->numItems()-1);

        // model->Write(model->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tdata(Psu):\tnpars = " << data_npars << "\tchi^2 = " << chi2 << "\tprob = " << prob << "\tfitting status = " << fit_status << endl;
        output << "\t" << bkg_fun.Data() << "\tSS:\tnsig = " << ss_mc << ":" << ss << "\tdmc = " << dmc << "\tss_cor = " << ss_cor << "\tdss = " << dss << "\ttot_err = " << tot_err << "\tstatus = " << status.Data() << "\n" << endl;
        // output << "\tnbkg = " << nbkg.getVal() << "\tnbkg_err = " << nbkg.getError() << "\n" << endl;
        TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
        leg->SetFillColor(0);
        leg->SetLineColor(0);
        leg->AddEntry(frame_data->findObject("data"),"MC","ep");
        leg->AddEntry(frame_data->findObject("fit"),"Bkg + Sig","l");
        // leg->AddEntry(frame_data->findObject("signal"),"Sig","l");
        leg->AddEntry(frame_data->findObject("background"),"Bkg","l");
        leg->Draw("same");

        pad2->cd();
        npoints = plotdata->GetN();
        point =0;
        hdatasub = new TGraphAsymmErrors(npoints);
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

        hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
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
        // line3->SetLineColor(kBlue);
        // //line3->SetLineStyle(kDashed);
        // line3->SetLineWidth(5.0);
        // line3->Draw();
        hdatasub->SetMarkerStyle(8);
        hdatasub->Draw("PESAME");
        // if (nomSumCurve) {
        //   TGraph *scaledSigCurve = new TGraph(nomSumCurve->GetN());
        //   for (int i = 0; i < nomSumCurve->GetN(); ++i) {
        //     double x, y, xb, yb;
        //     nomSumCurve->GetPoint(i, x, y);
        //     nomBkgCurve->GetPoint(i, xb, yb); 
        //     cout << x << " " << y << " " << xb << " " << yb << endl;
        //     scaledSigCurve->SetPoint(i, x, y - yb);
        //   }
        //   scaledSigCurve->SetLineColor(kRed); // Set the color to differentiate it
        //   scaledSigCurve->SetLineWidth(3);
        //   scaledSigCurve->SetLineStyle(ELineStyle::kDashed); // Keep the dashed line style
        //   scaledSigCurve->Draw("L SAME");
        // }
        nomSigCurve->Draw("L SAME");
        canv->SaveAs(Form("./test/pesudo_data_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));
        pad1->Close(); pad2->Close();
        // gPad->Print(Form("./test/pesudo_data_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

        // RooHist *hpull = frame_data->pullHist();
        // RooPlot *frame3 = CMS_hzg_mass->frame(Title("Pull Distribution"));
        // frame3->addPlotable(hpull, "P");
        // frame3->Draw();
        // gPad->Print(Form("test/%s_cat%d_%s_%dsig.pdf",channel.Data(),cat,bkg_fun.Data(),sig));

        cout << "\t=================================" << endl;
        cout << "\n\t Finish asimov data fit\n" << endl;
      }
    }
  }
}