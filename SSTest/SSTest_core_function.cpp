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
    RooRealVar* sigma = new RooRealVar("sigma","sigma",2.,0.,5.); 
    RooRealVar* MH = new RooRealVar("MH","MH",125., 124., 126.); 
    RooGaussian* sig_gau = new RooGaussian("sig_gau","sig_gau",*obs_var,*MH,*sigma);
    
    RooRealVar* sigma_CB = new RooRealVar("sigma_CB","sigma_CB",1.0, 0.01, 1.); 
    RooRealVar* alpha = new RooRealVar("alpha","alpha",0.5, 0., 1.0); 
    RooRealVar* n_CB = new RooRealVar("n_CB","n_CB",7.,1.,15.); 
    RooRealVar* fracG1 = new RooRealVar("fracG1","fracG1",0.3,0.,0.5); 
    RooCBShape* CBshape = new RooCBShape("CBShape", "CBShape", *obs_var, *MH, *sigma_CB, *alpha, *n_CB);
    RooAddPdf* sigPdf = new RooAddPdf("sigPdf","sigPdf",RooArgList(*sig_gau, *CBshape),*fracG1);
    return sigPdf;
}

RooAbsPdf* getBernsteinxZGMCShape(string prefix, int order, RooRealVar* obs_var){
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("./ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");
  ZGMC_file->Close();

  RooArgList *coeffList = new RooArgList();
  map<string, RooRealVar*> params;
  map<string, RooFormulaVar*> prods;
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),-5.,5.);
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
    //return bern;
  //  } else if (order==7) {
  //	RooBernsteinFast<7> *bern = new RooBernsteinFast<7>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  // 	return bern;
  } else {
	  return NULL;
  }

  //return bern;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}

RooAbsPdf* getExponentialZGMCShape(string prefix, int order, RooRealVar* obs_var){
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("./ZGMCShape.root");
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
  par1_exp1 = -0.1;          par1_lexp1 = -0.2;    par1_hexp1 = 0.;
  coeff1_exp1 = 0.9;          coeff1_lexp1 = -100.;    coeff1_hexp1 = 100.;
  }
  else if(order==3){
  par1_exp3 = -0.1;      par1_lexp3 = -0.2;    par1_hexp3 = 0.;
  coeff1_exp3 = 0.9;      coeff1_lexp3 = -100.;    coeff1_hexp3 = 100.;
  par3_exp3 = -0.001;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
  coeff3_exp3 = 0.4;      coeff3_lexp3 = -100.;    coeff3_hexp3 = 100.;
  } 
  else if(order==5)
  {
  par1_exp5 = -0.1;          par1_lexp5 = -0.2;    par1_hexp5 = 0.;
  coeff1_exp5 = 0.9;       coeff1_lexp5 = -100.;    coeff1_hexp5 = 100.;
  par3_exp5 = -0.001;          par3_lexp5 = -0.5;    par3_hexp5 = 0.;
  coeff3_exp5 = 0.4; coeff3_lexp5 = -100.;    coeff3_hexp5 = 100.;
  par5_exp5 = -0.005;      par5_lexp5 = -0.7;      par5_hexp5 = 0.;
  coeff5_exp5 = 0.4; coeff5_lexp5 = -100.;    coeff5_hexp5 = 100.;
  }
  else if(order==7)
  {
  par1_exp7 = -0.1;          par1_lexp7 = -0.2;    par1_hexp7 = 0.;
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
  coeff1_pow1 = 1.; coeff1_lpow1 = -10.;    coeff1_hpow1 = 10.;
  }
  else if(order==3){
  par1_pow3 = 2.7;      par1_lpow3 = 2.;    par1_hpow3 = 3.;
  coeff1_pow3 = 0.0001; coeff1_lpow3 = 0.;    coeff1_hpow3 = 0.001;
  par3_pow3 = -3.;      par3_lpow3 = -4;      par3_hpow3 = -2.;
  coeff3_pow3 = 0.99; coeff3_lpow3 = 0.9;    coeff3_hpow3 = 1.;
  } 
  else if(order==5)
  {
  par1_pow5 = -1.;          par1_lpow5 = -10.;    par1_hpow5 = 5.;
  coeff1_pow5 = 1.;       coeff1_lpow5 = -10.;    coeff1_hpow5 = 10.;
  par3_pow5 = -1.;          par3_lpow5 = -10.;    par3_hpow5 = 5.;
  coeff3_pow5 = 1.; coeff3_lpow5 = -10.;    coeff3_hpow5 = 10.;
  par5_pow5 = -1.;      par5_lpow5 = -10.;      par5_hpow5 = 5.;
  coeff5_pow5 = 1.; coeff5_lpow5 = -10.;    coeff5_hpow5 = 10.;
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

  double coeff1_lau1, coeff2_lau1, coeff1_lau2, coeff2_lau2, coeff3_lau2, coeff1_lau3, coeff2_lau3, coeff3_lau3, coeff4_lau3, coeff1_lau4, coeff2_lau4, coeff3_lau4, coeff4_lau4, coeff5_lau4; 
  double coeff1_hlau1, coeff2_hlau1, coeff1_hlau2, coeff2_hlau2, coeff3_hlau2, coeff1_hlau3, coeff2_hlau3, coeff3_hlau3, coeff4_hlau3, coeff1_hlau4, coeff2_hlau4, coeff3_hlau4, coeff4_hlau4, coeff5_hlau4; 
  double coeff1_llau1, coeff2_llau1, coeff1_llau2, coeff2_llau2, coeff3_llau2, coeff1_llau3, coeff2_llau3, coeff3_llau3, coeff4_llau3, coeff1_llau4, coeff2_llau4, coeff3_llau4, coeff4_llau4, coeff5_llau4; 

  if(order==1){
  coeff1_lau1 = 0.1; coeff1_llau1 = -100.;    coeff1_hlau1 = 100.;
  coeff2_lau1 = 0.5; coeff2_llau1 = -100.;    coeff2_hlau1 = 100.;
  }
  else if(order==2){
  coeff1_lau2 = 0.01; coeff1_llau2 = -100.;    coeff1_hlau2 = 100.;
  coeff2_lau2 = 0.5; coeff2_llau2 = -100.;    coeff2_hlau2 = 100.;
  coeff3_lau2 = 0.01; coeff3_llau2 = -100.;    coeff3_hlau2 = 100.;
  } 
  else if(order==3)
  {
  coeff1_lau3 = 0.1; coeff1_llau3 = -100.;    coeff1_hlau3 = 100.;
  coeff2_lau3 = 0.5; coeff2_llau3 = -100.;    coeff2_hlau3 = 100.;
  coeff3_lau3 = 0.01; coeff3_llau3 = -100.;    coeff3_hlau3 = 100.;
  coeff4_lau3 = 0.5; coeff4_llau3 = -100.;    coeff4_hlau3 =100.;
  }
  else if(order==4)
  {
  coeff1_lau4 = 0.1; coeff1_llau4 = -100.;    coeff1_hlau4 = 100.;
  coeff2_lau4 = 0.5; coeff2_llau4 = -100.;    coeff2_hlau4 = 100.;
  coeff3_lau4 = 0.01; coeff3_llau4 = -100.;    coeff3_hlau4 = 100.;
  coeff4_lau4 = 0.5; coeff4_llau4 = -100.;    coeff4_hlau4 =100.;
  coeff5_lau4 = 0.5; coeff5_llau4 = -100.;    coeff5_hlau4 =100.;
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
      RooGenericPdf *lau = new RooGenericPdf(Form("%s_lau4",prefix.c_str()),Form("%s_lau4",prefix.c_str()), "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)+@4*(@0)^(-6)+@5*(@0)^(-7)", RooArgList(*obs_var,*cp1,*cp2,*cp3,*cp4,*cp5));

      RooEffProd *lauZGMC = new RooEffProd(Form("%s_lau3xZG",prefix.c_str()),Form("%s_lau3xZG",prefix.c_str()), *ZGMCShape, *lau);
      return lauZGMC;
  } 
   else {
	  return NULL;
  }
}

RooAbsPdf* getPdf(string type, int order, RooRealVar* obs_var, const char* ext=""){
  if (type=="Bernstein") return getBernsteinxZGMCShape(Form("bern%d",order),order,obs_var);
  // else if (type=="Chebychev") return getChebychev(Form("cheb%d",order),order,*obs_var);
  else if (type=="Exponential") return getExponentialZGMCShape(Form("exp%d",order),order,obs_var);
  else if (type=="PowerLaw") return getPowerLawZGMCShape(Form("pow%d",order),order,obs_var);
  else if (type=="Laurent") return getLaurentZGMCShape(Form("lau%d",order),order,obs_var);
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void SSTest_core_function(int cat = 1, int sig = 0, TString channel = "untagged"){
  //background MC template
  TH1F* hbkg;
  TFile* fbkg = TFile::Open(Form("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/bkg_sig_template.root"));
  if (fbkg->GetListOfKeys()->Contains((Form("bkg_%s_cat%d",channel.Data(), cat))))
      hbkg = (TH1F*)fbkg->Get(Form("bkg_%s_cat%d",channel.Data(), cat));
  else
      abort();
  // TFile* fbkg = TFile::Open(Form("./bkg_template_v3_cut2/bkg/bkg_0sig_cat%d.root", cat));
  // hbkg = (TH1F*)fbkg->Get(Form("mass_cat%d", cat));
  double dataevents = hbkg->Integral();
  double mcsbevents = hbkg->Integral(0,22)+hbkg->Integral(28,80);

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
  functionClasses.push_back("Exponential");
  functionClasses.push_back("Laurent");
  functionClasses.push_back("PowerLaw");

  //initializing
  RooRealVar* CMS_hzg_mass = new RooRealVar("CMS_hzg_mass", "CMS_hzg_mass", 125.38, 105, 170);
  CMS_hzg_mass->setMin( 105. );
  CMS_hzg_mass->setMax( hbkg->GetBinCenter(hbkg->GetNbinsX()) + hbkg->GetBinWidth(hbkg->GetNbinsX()) );

  //background function fit
  RooDataHist* dbkg = new RooDataHist("data_bin","dataset with x", *CMS_hzg_mass, hbkg);
  RooDataHist* dsb = new RooDataHist("data_bin","dataset with x", *CMS_hzg_mass, hsb);
  // cout<<"bkg hist integral "<<hbkg->Integral()<<" "<<dbkg->sumEntries()<<endl;

  RooRealVar nsig("nsig","nsig",0,-100*sigevents,100*sigevents);
  RooRealVar nbkg("nbkg","nbkg",dataevents, 0.01*dataevents, 2*dataevents);

  for (vector<string>::iterator funcType=functionClasses.begin(); funcType!=functionClasses.end(); funcType++){
    int order = 1;
    TString bkg_fun;
    ofstream output(Form("./outputs/%s_%d_%dxsig.txt", channel.Data(), cat, sig), ofstream::app);
    // TFile *f = new TFile(Form("./outputs/%s_%d_%dxsig.root", channel.Data(), cat, sig),"UPDATE"); 
    while (order < 6){ //FIXME
      RooAbsPdf *bkgPdf = getPdf(*funcType,order,CMS_hzg_mass,"");
      bkg_fun = Form("%s%d", funcType->c_str(), order);
      order++;
      RooFitResult* bkgPdf_fit;
      if (bkgPdf){
        // Data side band fitting
        // hbkg->Write(hbkg->GetName(), TObject::kOverwrite);
        // hsb->Write(hsb->GetName(), TObject::kOverwrite);
        // hsig->Write(hsig->GetName(), TObject::kOverwrite);

        int bkg_npars, bkg_ndof;
        RooPlot *frame_bkg;

        CMS_hzg_mass->setRange("R1",105,122);
        CMS_hzg_mass->setRange("R3",122,128);
        CMS_hzg_mass->setRange("R2",128,170);

        //background function fit
        bkgPdf_fit = bkgPdf->fitTo(*dsb,Range("R1,R2"),RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); 
        bkg_npars = bkgPdf_fit->floatParsFinal().getSize();
        bkg_ndof = 74-bkg_npars;
        frame_bkg = CMS_hzg_mass->frame(Title(Form("Data side band with %s pdf", bkg_fun.Data())));
        dsb->plotOn(frame_bkg);
        bkgPdf->plotOn(frame_bkg,Range("R1"));
        bkgPdf->plotOn(frame_bkg,Range("R2"));
        bkgPdf->SetName(bkg_fun);
        cout << "start" << endl;
        // bkgPdf->Write(bkgPdf->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tsb:\tnpars = " << bkg_npars << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars)*bkg_ndof, bkg_ndof) << endl;
        bkgPdf->plotOn(frame_bkg,Range("R3"));
        frame_bkg->Draw();
        gPad->Print(Form("./test/data_sb_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

        cout << "\t=================================" << endl;
        cout << "\n\t Finish data side band fit\n" << endl;

        // MC background fitting
        bkgPdf_fit = bkgPdf->fitTo(*dbkg,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));

        RooAbsCollection *m_bkgParameters = bkgPdf->getParameters(RooArgSet())->selectByAttrib("Constant", false);
        TIterator *bkgit = m_bkgParameters->createIterator();
        for (RooRealVar *p = (RooRealVar *)bkgit->Next(); p != 0; p = (RooRealVar *)bkgit->Next()) p->setConstant(kTRUE);

        bkg_npars = bkgPdf_fit->floatParsFinal().getSize();
        bkg_ndof = 80-bkg_npars;
        frame_bkg = CMS_hzg_mass->frame(Title(Form("Background with %s pdf", bkg_fun.Data())));
        dbkg->plotOn(frame_bkg);
        bkgPdf->plotOn(frame_bkg);
        bkgPdf->SetName(Form("%s_model", bkg_fun.Data()));
        // bkgPdf->Write(bkgPdf->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tbkg:\tnpars = " << bkg_npars << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars)*bkg_ndof, bkg_ndof) << endl;
        frame_bkg->Draw();
        // RooHist *hpull = frame_bkg->pullHist();
        // RooPlot *frame3 = CMS_hzg_mass->frame(Title("Pull Distribution"));
        // frame3->addPlotable(hpull, "P");
        // frame3->Draw();
        gPad->Print(Form("./test/mc_bkg_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

        cout << "\t=================================\n";
        cout << "\n\tFinish background function fit\n" << endl;

        RooAbsPdf* sigPdf = getsigPdfPdf(CMS_hzg_mass);
        RooFitResult *sigPdf_fit;
        RooDataHist* dsig = new RooDataHist("data_bin","dataset with x", *CMS_hzg_mass, hsig);
        sigPdf_fit = sigPdf->fitTo(*dsig,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE

        RooAbsCollection *m_sigParameters = sigPdf->getParameters(RooArgSet())->selectByAttrib("Constant", false);
        TIterator *sigit = m_sigParameters->createIterator();
        for (RooRealVar *p = (RooRealVar *)sigit->Next(); p != 0; p = (RooRealVar *)sigit->Next()) p->setConstant(kTRUE);

        int sig_npars = sigPdf_fit->floatParsFinal().getSize();
        int sig_ndof = 80-sig_npars;
        RooPlot *frame_sig = CMS_hzg_mass->frame(Title("sigPdf with distorted Gaussian pdf"));
        dsig->plotOn(frame_sig, DataError(RooAbsData::SumW2));
        sigPdf->plotOn(frame_sig);
        // sigPdf->Write(sigPdf->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tsig:\tnpars = " << sig_npars << "\tchi^2 = " << frame_sig->chiSquare(sig_npars) << "\tprob = " << TMath::Prob(frame_sig->chiSquare(sig_npars)*sig_ndof, sig_ndof) << endl;
        frame_sig->Draw();
        gPad->Print(Form("./test/sigPdf_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

        cout << "\t=================================" << endl;
        cout << "\n\t Finish sigPdf function fit\n" << endl;

        //inject signal
        TH1F* hdata = (TH1F*)hbkg->Clone();
        hdata->Add(hsig, sig);
        hdata->SetName("asimov data");
        // hdata->Write(hdata->GetName(), TObject::kOverwrite);
        RooDataHist* ddata = new RooDataHist("data_bin","dataset with x", *CMS_hzg_mass, hdata);

        // asimov dataset fit
        RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*sigPdf, *bkgPdf),RooArgList(nsig,nbkg));
        RooFitResult *model_fit;
        model_fit = model->fitTo(*ddata,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE)); //FIXME kTRUE or kFALSE
        int data_npars = model_fit->floatParsFinal().getSize();
        int data_ndof = 80-data_npars;
        RooPlot *frame_data = CMS_hzg_mass->frame(Title(Form("Asimove data with with x%d signal", sig)));
        ddata->plotOn(frame_data, DataError(RooAbsData::SumW2));
        model->plotOn(frame_data);
        model->SetName(Form("%s_model", bkg_fun.Data()));
        // model->Write(model->GetName(), TObject::kOverwrite);
        output << "\t" << bkg_fun.Data() << "\tdata:\tnpars = " << data_npars << "\tchi^2 = " << frame_data->chiSquare(data_npars) << "\tprob = " << TMath::Prob(frame_data->chiSquare(data_npars)*data_ndof, data_ndof) << endl << "\t" << bkg_fun.Data() << "\tSS:\tnsig = " << nsig.getVal() << "\tdmc = " << nsig.getError() << "\tnbkg = " << nbkg.getVal() << "\tnbkg_err = " << nbkg.getError() << "\n" << endl;
        frame_data->Draw();
        gPad->Print(Form("./test/asimov_data_shape_%s_cat%d_%s.pdf",channel.Data(),cat,bkg_fun.Data()));

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