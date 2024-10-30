#include "TCanvas.h"

#include "RooPlot.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "../interface/RooPowerLaw.h"
#include "../interface/RooPowerLawSum.h"
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

#include "../interface/PdfModelBuilder.h"

#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

using namespace std;
using namespace RooFit;
using namespace boost;

PdfModelBuilder::PdfModelBuilder():
  obs_var_set(false),
  signal_modifier_set(false),
  signal_set(false),
  bkgHasFit(false),
  sbHasFit(false),
  keysPdfAttributesSet(false),
  verbosity(0)
{
  
  recognisedPdfTypes.push_back("Bernstein");
  recognisedPdfTypes.push_back("Exponential");
  recognisedPdfTypes.push_back("PowerLaw");
  recognisedPdfTypes.push_back("Laurent");
  recognisedPdfTypes.push_back("KeysPdf");
  recognisedPdfTypes.push_back("File");
  //bing
  recognisedPdfTypes.push_back("BernsteinStepxGau");
  recognisedPdfTypes.push_back("PowerLawStepxGau");
  recognisedPdfTypes.push_back("LaurentStepxGau");
  recognisedPdfTypes.push_back("ExponentialStepxGau");
  
  recognisedPdfTypes.push_back("BernsteinxZGMCShape");
  recognisedPdfTypes.push_back("ExponentialZGMCShape");
  recognisedPdfTypes.push_back("PowerLawZGMCShape");
  recognisedPdfTypes.push_back("LaurentZGMCShape");


  wsCache = new RooWorkspace("PdfModelBuilderCache");

};

PdfModelBuilder::~PdfModelBuilder(){};

void PdfModelBuilder::setObsVar(RooRealVar *var){
  obs_var=var;
  obs_var_set=true;
}

void PdfModelBuilder::setSignalModifier(RooRealVar *var){
  signalModifier=var;
  signal_modifier_set=true;
}

void PdfModelBuilder::setSignalModifierVal(float val){
  signalModifier->setVal(val);
}

void PdfModelBuilder::setSignalModifierConstant(bool val){
  signalModifier->setConstant(val);
}

RooAbsPdf* PdfModelBuilder::getChebychev(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    //RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    //prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*params[name]);
  }
  //RooChebychev *cheb = new RooChebychev(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}

RooAbsPdf* PdfModelBuilder::getBernstein(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
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
  //bing
  RooRealVar *fgaus = new RooRealVar("fgaus", "gaus fraction",0.5,0.,1.) ;
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),10,0.0,10.0) ;
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma",prefix.c_str()),Form("%s_sigma",prefix.c_str()),8,-10.,20.) ;
  RooGaussian *gaus = new RooGaussian(Form("%s_gaus",prefix.c_str()),Form("%s_gaus",prefix.c_str()),*obs_var,*mean,*sigma) ;

  RooRealVar *step_value = new RooRealVar(Form("%s_step",prefix.c_str()), Form("%s_step",prefix.c_str()), 105., 100., 110.);
    
  //RooBernstein *bern = new RooBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  if (order==1) {
	 RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	//bing
    //RooAddPdf *bern_gaus = new RooAddPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(),RooArgList(*bern,*gaus),*fgaus) ;

    //RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(), *obs_var, *bern, *gaus);
    //RooGenericPdf *step_func = new RooGenericPdf("step_func", "step_func", "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));

    RooGenericPdf *step_func = new RooGenericPdf(Form("%s_stepfunc",prefix.c_str()), Form("%s_stepfunc",prefix.c_str()), "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));
    obs_var->setRange(-200, 500);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *step_func);
    obs_var->setRange(95, 180);
    //return bern;
    return bern_gaus;
  } else if (order==2) {
	RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	//bing
    //RooAddPdf *bern_gaus = new RooAddPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(),RooArgList(*bern,*gaus),*fgaus) ;
    //RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(), *obs_var, *bern, *gaus);
    RooGenericPdf *step_func = new RooGenericPdf(Form("%s_stepfunc",prefix.c_str()), Form("%s_stepfunc",prefix.c_str()), "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));
    obs_var->setRange(-200, 500);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *step_func);
    obs_var->setRange(95, 180);
    //return bern;
    return bern_gaus;
  } else if (order==3) {
	RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	//bing
    //RooAddPdf *bern_gaus = new RooAddPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(),RooArgList(*bern,*gaus),*fgaus) ;
    //RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(), *obs_var, *bern, *gaus);
     RooGenericPdf *step_func = new RooGenericPdf(Form("%s_stepfunc",prefix.c_str()), Form("%s_stepfunc",prefix.c_str()), "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));
    obs_var->setRange(-200, 500);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *step_func);
    obs_var->setRange(95, 180);
    //return bern;
    return bern_gaus;
  } else if (order==4) {
	RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	//bing
    //RooAddPdf *bern_gaus = new RooAddPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(),RooArgList(*bern,*gaus),*fgaus) ;
    //RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(), *obs_var, *bern, *gaus);
     RooGenericPdf *step_func = new RooGenericPdf(Form("%s_stepfunc",prefix.c_str()), Form("%s_stepfunc",prefix.c_str()), "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));
    obs_var->setRange(-200, 500);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *step_func);
    obs_var->setRange(95, 180);
    //return bern;
    return bern_gaus;
  } else if (order==5) {
	RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	//bing
    //RooAddPdf *bern_gaus = new RooAddPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(),RooArgList(*bern,*gaus),*fgaus) ;
    //RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(), *obs_var, *bern, *gaus);
     RooGenericPdf *step_func = new RooGenericPdf(Form("%s_stepfunc",prefix.c_str()), Form("%s_stepfunc",prefix.c_str()), "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));
    obs_var->setRange(-200, 500);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *step_func);
    obs_var->setRange(95, 180);
    //return bern;
    return bern_gaus;
  } else if (order==6) {
	RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	//bing
    //RooAddPdf *bern_gaus = new RooAddPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(),RooArgList(*bern,*gaus),*fgaus) ;
    //RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(prefix.append("gaus").c_str(),prefix.append("gaus").c_str(), *obs_var, *bern, *gaus);
     RooGenericPdf *step_func = new RooGenericPdf(Form("%s_stepfunc",prefix.c_str()), Form("%s_stepfunc",prefix.c_str()), "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *bern));
    obs_var->setRange(-200, 500);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *step_func);
    obs_var->setRange(95, 180);
    //return bern;
    return bern_gaus;
//  } else if (order==7) {
//	RooBernsteinFast<7> *bern = new RooBernsteinFast<7>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
 // 	return bern;
  } else {
	return NULL;
  }

  //return bern;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}

RooAbsPdf* PdfModelBuilder::getBernsteinStepxGau(string prefix, int order){
  
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),5,2.,6);
  //RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),1.5,0.1,3);
  //RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),1.5,0.1,3);//6GeV
  //RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),2,0.1,3);//9GeV
  //RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),1.5,0.1,3);//5GeV
  //RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),10,1,20);//1GeV
  //RooRealVar *sigma = new RooRealVar(Form("%s_sigma_b%d",prefix.c_str(),order),Form("%s_sigma_b%d",prefix.c_str(),order),1.5,0.1,3);//16GeV
  //RooRealVar *step = new RooRealVar(Form("step_b%d",order),Form("step_b%d",order),115,110,125);//1GeV
  RooRealVar *step = new RooRealVar(Form("%s_step_b%d",prefix.c_str(),order),Form("%s_step_b%d",prefix.c_str(),order),105,100,110);

  RooArgList *coeffList = new RooArgList();
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  
  for (int i=0; i<order+1; i++){
    string name = Form("%s_b%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.03,-14.,15.);
    if (i==0){
      param->setVal(10.);
      param->setConstant(kTRUE);
    }
    //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.03,-14.,15.);
    //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.3,-14.,15.);//6GeV
    //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),6,-14.,15.);//5GeV
    params.insert(pair<string,RooRealVar*>(name,param));
    coeffList->add(*params[name]);
  }
  //RooBernstein *bern = new RooBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  if (order==1) {
  RooGaussStepBernstein *bern =new RooGaussStepBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*mean,*sigma,*step, *coeffList);
    	return bern;
  } else if (order==2) {
    RooGaussStepBernstein *bern =new RooGaussStepBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*mean,*sigma,*step, *coeffList);
  	return bern;
  } else if (order==3) {
    RooGaussStepBernstein *bern =new RooGaussStepBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*mean,*sigma,*step, *coeffList);
  	return bern;
  } else if (order==4) {
    RooGaussStepBernstein *bern =new RooGaussStepBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*mean,*sigma,*step, *coeffList);
  	return bern;
  } else if (order==5) {
    RooGaussStepBernstein *bern =new RooGaussStepBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*mean,*sigma,*step, *coeffList);
  	return bern;
  } 
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getPowerLawStepxGau(string prefix, int order, int cat){
  if(order%2==0) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  double sigma_pow,sigma_lpow,sigma_hpow;
  // double turnon_pow = 107.;
  double turnon_pow,turnon_lpow,turnon_hpow;
  double par1_pow1, par1_pow3, par3_pow3, par1_pow5, par3_pow5, par5_pow5;
  double par1_hpow1, par1_hpow3, par3_hpow3, par1_hpow5, par3_hpow5, par5_hpow5;
  double par1_lpow1, par1_lpow3, par3_lpow3, par1_lpow5, par3_lpow5, par5_lpow5;
  double coeff1_pow1, coeff1_pow3, coeff3_pow3, coeff1_pow5, coeff3_pow5, coeff5_pow5;
  double coeff1_hpow1, coeff1_hpow3, coeff3_hpow3, coeff1_hpow5, coeff3_hpow5, coeff5_hpow5;
  double coeff1_lpow1, coeff1_lpow3, coeff3_lpow3, coeff1_lpow5, coeff3_lpow5, coeff5_lpow5;

  if(cat==1)
  {
    if(order==1){
    turnon_pow = 108.32;      turnon_lpow = 106.;   turnon_hpow = 115.;
    sigma_pow = 6.7;          sigma_lpow = 2.;      sigma_hpow = 7.;
    par1_pow1 = -6.;          par1_lpow1 = -15.;    par1_hpow1 = -5.;
    coeff1_pow1 = 3.6576e-04; coeff1_lpow1 = 0.;    coeff1_hpow1 = 0.5;
    }
    else if(order==3){
    turnon_pow = 108.2;       turnon_lpow = 106.;   turnon_hpow = 115.;
    sigma_pow = 6.6;          sigma_lpow = 2.;      sigma_hpow = 8.;
    par1_pow3 = -6.0519;      par1_lpow3 = -10.;    par1_hpow3 = -2.;
    coeff1_pow3 = 4.3299e-05; coeff1_lpow3 = 0.;    coeff1_hpow3 = 0.5;
    par3_pow3 = -4.5113;      par3_lpow3 = -8;      par3_hpow3 = 2.;
    coeff3_pow3 = 4.0951e-13; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    turnon_pow = 108;         turnon_lpow = 106.;   turnon_hpow = 115.;
    sigma_pow = 6.;           sigma_lpow = 2.;      sigma_hpow = 8.;
    par1_pow5 = -9.;          par1_lpow5 = -15.;    par1_hpow5 = -5.;
    coeff1_pow5 = 0.04;       coeff1_lpow5 = 0.;    coeff1_hpow5 = 1.;
    par3_pow5 = -7.;          par3_lpow5 = -10.;    par3_hpow5 = -2.;
    coeff3_pow5 = 4.0951e-13; coeff3_lpow5 = 0.;    coeff3_hpow5 = 1.;
    par5_pow5 = -4.5113;      par5_lpow5 = -8;      par5_hpow5 = 2.;
    coeff5_pow5 = 4.0951e-13; coeff5_lpow5 = 0.;    coeff5_hpow5 = 1.;
    }
  }
  else if(cat==2)
  {
    if(order==1){
    sigma_pow = 4;          sigma_lpow = 3.;      sigma_hpow = 15.;
    turnon_pow = 108.;        turnon_lpow = 106.;   turnon_hpow = 110.;
    par1_pow1 = -5.6;         par1_lpow1 = -15.;    par1_hpow1 = 5.;
    coeff1_pow1 = 0.002;      coeff1_lpow1 = 0.;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 4.;           sigma_lpow = 3.;      sigma_hpow = 15.;
    turnon_pow = 108.5;       turnon_lpow = 105.;   turnon_hpow = 110.;
    par1_pow3 = -7;           par1_lpow3 = -10.;    par1_hpow3 = 5.;
    coeff1_pow3 = 7.7276e-01; coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -5.9269;      par3_lpow3 = -8;      par3_hpow3 = 5.;
    coeff3_pow3 = 6.3834e-05; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 4.;           sigma_lpow = 3.;      sigma_hpow = 15.;
    turnon_pow = 108.64;      turnon_lpow = 105.;   turnon_hpow = 110.;
    par1_pow5 = -10.4;        par1_lpow5 = -15.;    par1_hpow5 = 5.;
    coeff1_pow5 = 7.4192e-01; coeff1_lpow5 = 0.;    coeff1_hpow5 = 1.;
    par3_pow5 = -7.0937;      par3_lpow5 = -10.;    par3_hpow5 = 2.;
    coeff3_pow5 = 9.0955e-01; coeff3_lpow5 = 0.;    coeff3_hpow5 = 1.;
    par5_pow5 = -6.2482;      par5_lpow5 = -10;     par5_hpow5 = 1.;
    coeff5_pow5 =2.1590e-03;  coeff5_lpow5 = 0.001; coeff5_hpow5 = 1.;
    }
  }
  else if(cat==3)
  {
    if(order==1){
    sigma_pow = 3.5;          sigma_lpow = 3.;      sigma_hpow = 8.;
    turnon_pow = 106.;        turnon_lpow = 105.;   turnon_hpow = 108.;
    par1_pow1 = -6.5;         par1_lpow1 = -15.;    par1_hpow1 = -5.;
    coeff1_pow1 = 0.02;       coeff1_lpow1 = 0.;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 3.5;          sigma_lpow = 3.;      sigma_hpow = 5.;
    turnon_pow = 107.57;      turnon_lpow = 106.;   turnon_hpow = 108.;
    par1_pow3 = -8.8928;      par1_lpow3 = -10.;    par1_hpow3 = -5.;
    coeff1_pow3 = 0.89;       coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -6.5;         par3_lpow3 = -10;     par3_hpow3 = -2.;
    coeff3_pow3 = 0.001;      coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 3.8618;       sigma_lpow = 2.;      sigma_hpow = 5.;
    turnon_pow = 107.34;      turnon_lpow = 105.;   turnon_hpow = 108.;
    par1_pow5 = -7.0724;      par1_lpow5 = -10.;    par1_hpow5 = -5.;
    coeff1_pow5 = 9.1185e-01; par1_lpow5 = 0.;      par1_hpow5 = 1.;
    par3_pow5 = -7.4414;      par3_lpow5 = -10.;    par3_hpow5 = -2.;
    coeff3_pow5 = 5.7468e-06; par3_lpow5 = 0.;      par3_hpow5 = 1.;
    par5_pow5 = -8.4821;      par5_lpow5 = -10;     par5_hpow5 = -1.;
    coeff5_pow5 = 2.8846e-03; par5_lpow5 = 0.;      par5_hpow5 = 1.;
    }
  }
  else if(cat==4)
  {
    if(order==1){
    sigma_pow = 3.3;          sigma_lpow = 2.5;      sigma_hpow = 8.;
    turnon_pow = 106.9;       turnon_lpow = 105.;   turnon_hpow = 108.;
    par1_pow1 = -6.3453;      par1_lpow1 = -10.;    par1_hpow1 = -5.;
    coeff1_pow1 = 3.3018e-02; coeff1_lpow1 = 0.;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 3.5;          sigma_lpow = 3.;      sigma_hpow = 8.;
    turnon_pow = 107;         turnon_lpow = 105.;   turnon_hpow = 108.;
    par1_pow3 = -6.3671;      par1_lpow3 = -10.;    par1_hpow3 = -5.;
    coeff1_pow3 = 9.6543e-01; coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -4.267;       par3_lpow3 = -8;      par3_hpow3 = -2.;
    coeff3_pow3 = 1.1001e-06; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 3.5546;       sigma_lpow = 2.;      sigma_hpow = 5.;
    turnon_pow = 106.91;      turnon_lpow = 106.;   turnon_hpow = 109.;
    par1_pow5 = -8.8;         par1_lpow5 = -11.;    par1_hpow5 = -5.;
    coeff1_pow5 = 1.e-05;     par1_lpow5 = 0.;      par1_hpow5 = 1.;
    par3_pow5 = -6.5;         par3_lpow5 = -8.;     par3_hpow5 = -3.;
    coeff3_pow5 = 1.9518e-02; par3_lpow5 = 0.;      par3_hpow5 = 1.;
    par5_pow5 = -6.3;         par5_lpow5 = -8;      par5_hpow5 = -3.;
    coeff5_pow5 = 2.0522e-02; par5_lpow5 = 0.;      par5_hpow5 = 1.;
    }
  }
  else if(cat==501)
  {
    if(order==1){
    sigma_pow = 6.;           sigma_lpow = 0.1;     sigma_hpow = 10.;
    turnon_pow = 108.;        turnon_lpow = 106.;   turnon_hpow = 109.;
    par1_pow1 = -5.6;         par1_lpow1 = -10.;    par1_hpow1 = -5.;
    coeff1_pow1 = 0.0002;     coeff1_lpow1 = 0.;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 8.5;          sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 108;         turnon_lpow = 100.;   turnon_hpow = 109.;
    par1_pow3 = -6.4708;      par1_lpow3 = -10.;    par1_hpow3 = -2.;
    coeff1_pow3 = 8.4306e-01; coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -6.1246;      par3_lpow3 = -8;      par3_hpow3 = -2.;
    coeff3_pow3 = 1.4141e-06; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 6.;           sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 108;         turnon_lpow = 100.;   turnon_hpow = 109.;
    par1_pow5 = -9;           par1_lpow5 = -15.;    par1_hpow5 = -5.;
    coeff1_pow5 = 0.04;       par1_lpow5 = 0.;      par1_hpow5 = 1.;
    par3_pow5 = -7.;          par3_lpow5 = -10.;    par3_hpow5 = -2.;
    coeff3_pow5 = 0.1;        par3_lpow5 = 0.;      par3_hpow5 = 1.;
    par5_pow5 = -5.6;         par5_lpow5 = -10;     par5_hpow5 = -1.;
    coeff5_pow5 = 0.001;      par5_lpow5 = 0.001;   par5_hpow5 = 1.;
    }
  }
  else if(cat==502)
  {
    if(order==1){
    sigma_pow = 6.;           sigma_lpow = 0.1;     sigma_hpow = 10.;
    turnon_pow = 104.6;        turnon_lpow = 104.;   turnon_hpow = 105.;
    par1_pow1 = -5.6;         par1_lpow1 = -8.;    par1_hpow1 = -5.;
    coeff1_pow1 = 0.0006;     coeff1_lpow1 = 0.;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 5.;          sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 104.6;      turnon_lpow = 100.;   turnon_hpow = 110.;
    par1_pow3 = -6.1;         par1_lpow3 = -10.;   par1_hpow3 = -5.;
    coeff1_pow3 = 0.72;       coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -5.1;      par3_lpow3 = -6;      par3_hpow3 = -1.;
    coeff3_pow3 = 0.0003; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 6.;           sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 104.6;         turnon_lpow = 100.;   turnon_hpow = 110.;
    par1_pow5 = -9;           par1_lpow5 = -15.;    par1_hpow5 = -5.;
    coeff1_pow5 = 0.04;       par1_lpow5 = 0.;      par1_hpow5 = 1.;
    par3_pow5 = -7.;          par3_lpow5 = -10.;    par3_hpow5 = -2.;
    coeff3_pow5 = 0.1;        par3_lpow5 = 0.;      par3_hpow5 = 1.;
    par5_pow5 = -5.6;         par5_lpow5 = -10;     par5_hpow5 = -1.;
    coeff5_pow5 = 0.002;      par5_lpow5 = 0.001;   par5_hpow5 = 1.;
    }
  }
  else if(cat==503)
  {
    if(order==1){
    sigma_pow = 4.;           sigma_lpow = 0.1;     sigma_hpow = 5.;
    turnon_pow = 104.6;        turnon_lpow = 104.;   turnon_hpow = 105.;
    par1_pow1 = -5.2;         par1_lpow1 = -8.;    par1_hpow1 = -1.;
    coeff1_pow1 = 0.001;     coeff1_lpow1 = 0.000001;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 5.;          sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 104.6;      turnon_lpow = 100.;   turnon_hpow = 110.;
    par1_pow3 = -5;         par1_lpow3 = -8.;   par1_hpow3 = -3.;
    coeff1_pow3 = 0.3;       coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -4;      par3_lpow3 = -5;      par3_hpow3 = -1.;
    coeff3_pow3 = 0.01; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 4.8;           sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 104.6;         turnon_lpow = 100.;   turnon_hpow = 110.;
    par1_pow5 = -6;           par1_lpow5 = -12.;    par1_hpow5 = -5.;
    coeff1_pow5 = 0.07;       par1_lpow5 = 0.;      par1_hpow5 = 1.;
    par3_pow5 = -5.;          par3_lpow5 = -10.;    par3_hpow5 = -2.;
    coeff3_pow5 = 0.03;        par3_lpow5 = 0.;      par3_hpow5 = 1.;
    par5_pow5 = -4;         par5_lpow5 = -10;     par5_hpow5 = -1.;
    coeff5_pow5 = 0.0001;      par5_lpow5 = 0.00001;   par5_hpow5 = 1.;
    }
  }
  else if(cat==6789)
  {
    if(order==1){
    sigma_pow = 4.95;           sigma_lpow = 1.;     sigma_hpow = 8.;
    turnon_pow = 104.6;        turnon_lpow = 104.;   turnon_hpow = 105.;
    par1_pow1 = -5.2;         par1_lpow1 = -8.;    par1_hpow1 = -1.;
    coeff1_pow1 = 0.001;     coeff1_lpow1 = 0.000001;    coeff1_hpow1 = 1.;
    }
    else if(order==3){
    sigma_pow = 5.;          sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 104.6;      turnon_lpow = 100.;   turnon_hpow = 110.;
    par1_pow3 = -5;         par1_lpow3 = -8.;   par1_hpow3 = -3.;
    coeff1_pow3 = 0.3;       coeff1_lpow3 = 0.;    coeff1_hpow3 = 1.;
    par3_pow3 = -4;      par3_lpow3 = -5;      par3_hpow3 = -1.;
    coeff3_pow3 = 0.01; coeff3_lpow3 = 0.;    coeff3_hpow3 = 1.;
    } 
    else if(order==5)
    {
    sigma_pow = 4.8;           sigma_lpow = 3.;      sigma_hpow = 10.;
    turnon_pow = 104.6;         turnon_lpow = 100.;   turnon_hpow = 110.;
    par1_pow5 = -6;           par1_lpow5 = -12.;    par1_hpow5 = -5.;
    coeff1_pow5 = 0.07;       par1_lpow5 = 0.;      par1_hpow5 = 1.;
    par3_pow5 = -5.;          par3_lpow5 = -10.;    par3_hpow5 = -2.;
    coeff3_pow5 = 0.03;        par3_lpow5 = 0.;      par3_hpow5 = 1.;
    par5_pow5 = -4;         par5_lpow5 = -10;     par5_hpow5 = -1.;
    coeff5_pow5 = 0.0001;      par5_lpow5 = 0.00001;   par5_hpow5 = 1.;
    }
  }
  
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),order),Form("%s_sigma_p%d",prefix.c_str(),order),sigma_pow,sigma_lpow,sigma_hpow);
  RooRealVar *turnon = new RooRealVar(Form("%s_turnon_p%d",prefix.c_str(),order),Form("%s_turnon_p%d",prefix.c_str(),order),turnon_pow,turnon_lpow,turnon_hpow);
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  
    if (order==1) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow1",prefix.c_str()),Form("%s_p1_pow1",prefix.c_str()),par1_pow1,par1_lpow1,par1_hpow1);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow1",prefix.c_str()),Form("%s_cp1_pow1",prefix.c_str()),coeff1_pow1,coeff1_lpow1,coeff1_hpow1);
    	RooGenericPdf *step = new RooGenericPdf(Form("%s_step_pow1",prefix.c_str()),Form("%s_step_pow1",prefix.c_str()), "1e-20+(@0 > @1)*(@3*(@0)^(@2))", RooArgList(*obs_var,*turnon,*p1,*cp1));
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_pow1",prefix.c_str()),Form("%s_gau_pow1",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxpow = new RooFFTConvPdf(Form("%s_gauxpow1",prefix.c_str()),Form("%s_gauxpow1",prefix.c_str()),*obs_var,*step,*gau);
      return gauxpow;
  } else if (order==3) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow3",prefix.c_str()),Form("%s_p1_pow3",prefix.c_str()),par1_pow3,par1_lpow3,par1_hpow3);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow3",prefix.c_str()),Form("%s_cp1_pow3",prefix.c_str()),coeff1_pow3,coeff1_lpow3,coeff1_hpow3);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_pow3",prefix.c_str()),Form("%s_p3_pow3",prefix.c_str()),par3_pow3,par3_lpow3,par3_hpow3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_pow3",prefix.c_str()),Form("%s_cp3_pow3",prefix.c_str()),coeff3_pow3,coeff3_lpow3,coeff3_hpow3);
    	RooGenericPdf *step = new RooGenericPdf(Form("%s_step_pow3",prefix.c_str()),Form("%s_step_pow3",prefix.c_str()), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4))", RooArgList(*obs_var,*turnon,*p1,*cp1,*p3,*cp3));
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_pow3",prefix.c_str()),Form("%s_gau_pow3",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxpow = new RooFFTConvPdf(Form("%s_gauxpow3",prefix.c_str()),Form("%s_gauxpow3",prefix.c_str()),*obs_var,*step,*gau);
      return gauxpow;
  } else if (order==5) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_pow5",prefix.c_str()),Form("%s_p1_pow5",prefix.c_str()),par1_pow5,par1_lpow5,par1_hpow5);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_pow5",prefix.c_str()),Form("%s_cp1_pow5",prefix.c_str()),coeff1_pow5,coeff1_lpow5,coeff1_hpow5);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_pow5",prefix.c_str()),Form("%s_p3_pow5",prefix.c_str()),par3_pow5,par3_lpow5,par3_hpow5);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_pow5",prefix.c_str()),Form("%s_cp3_pow5",prefix.c_str()),coeff3_pow5,coeff3_lpow5,coeff3_hpow5);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_pow5",prefix.c_str()),Form("%s_p5_pow5",prefix.c_str()),par5_pow5,par5_lpow5,par5_hpow5);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_pow5",prefix.c_str()),Form("%s_cp5_pow5",prefix.c_str()),coeff5_pow5,coeff5_lpow5,coeff5_hpow5);
    	RooGenericPdf *step = new RooGenericPdf(Form("%s_step_pow5",prefix.c_str()),Form("%s_step_pow5",prefix.c_str()), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(*obs_var,*turnon,*p1,*cp1,*p3,*cp3,*p5,*cp5));
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_pow5",prefix.c_str()),Form("%s_gau_pow5",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxpow = new RooFFTConvPdf(Form("%s_gauxpow5",prefix.c_str()),Form("%s_gauxpow5",prefix.c_str()),*obs_var,*step,*gau);
      return gauxpow;
  } 
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getPowerLawGeneric(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    string formula="";
    RooArgList *dependents = new RooArgList();
    dependents->add(*obs_var);
    // first do recursive fraction
    if (order>1) {
      formula += "(1.-";
      for (int i=1; i<=nfracs; i++){
        if (i<nfracs) formula += Form("@%d-",i);
        else formula += Form("@%d)*",i);
        string name =  Form("%s_f%d",prefix.c_str(),i);
        params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.1,0.,1.)));
        dependents->add(*params[name]);
      }
    }
    for (int i=1; i<=npows; i++){
      string pname =  Form("%s_p%d",prefix.c_str(),i);
      string fname =  Form("%s_f%d",prefix.c_str(),i-1);
      params.insert(pair<string,RooRealVar*>(pname, new RooRealVar(pname.c_str(),pname.c_str(),TMath::Max(-10.,-2.*(i+1)),-10.,0.)));
      if (i==1) {
        formula += Form("TMath::Power(@0,@%d)",nfracs+i);
        dependents->add(*params[pname]);
      }
      else {
        formula += Form(" + @%d*TMath::Power(@0,@%d)",i-1,nfracs+i);
        dependents->add(*params[pname]);
      }
    }
    cout << "FORM -- " << formula << endl;
    dependents->Print("v");
    RooGenericPdf *pow = new RooGenericPdf(prefix.c_str(),prefix.c_str(),formula.c_str(),*dependents);
    pow->Print("v");
    return pow;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));

  }
}

RooAbsPdf* PdfModelBuilder::getPowerLaw(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-2.;
    double low=-10.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *pow = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return pow;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));

}

RooAbsPdf* PdfModelBuilder::getExponential(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-1.;
    double low=-2.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *exp = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return exp;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(exp->GetName(),exp));

}

RooAbsPdf* PdfModelBuilder::getExponentialStepxGau(string prefix, int order, int cat){
  if(order%2==0) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  double sigma_exp,sigma_lexp,sigma_hexp;
  double turnon_exp,turnon_lexp,turnon_hexp;
  double par1_exp1, par1_exp3, par3_exp3, par1_exp5, par3_exp5, par5_exp5;
  double par1_hexp1, par1_hexp3, par3_hexp3, par1_hexp5, par3_hexp5, par5_hexp5;
  double par1_lexp1, par1_lexp3, par3_lexp3, par1_lexp5, par3_lexp5, par5_lexp5;
  double coeff1_exp1, coeff1_exp3, coeff3_exp3, coeff1_exp5, coeff3_exp5, coeff5_exp5;
  double coeff1_hexp1, coeff1_hexp3, coeff3_hexp3, coeff1_hexp5, coeff3_hexp5, coeff5_hexp5;
  double coeff1_lexp1, coeff1_lexp3, coeff3_lexp3, coeff1_lexp5, coeff3_lexp5, coeff5_lexp5;
  if(cat==1)
  {
    if(order==1){
    turnon_exp = 108;      turnon_lexp = 106.;   turnon_hexp = 115.;
    sigma_exp = 5;          sigma_lexp = 3.;      sigma_hexp = 10.;
    par1_exp1 = -0.03;          par1_lexp1 = -0.7;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.5; coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    turnon_exp = 108;       turnon_lexp = 107.;   turnon_hexp = 115.;
    sigma_exp = 5;          sigma_lexp = 3.;      sigma_hexp = 10.;
    par1_exp3 = -0.04;      par1_lexp3 = -0.2;    par1_hexp3 = 0.;
    coeff1_exp3 = 0.3; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -0.03;      par3_lexp3 = -0.2;      par3_hexp3 = 0.;
    coeff3_exp3 = 0.1; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    turnon_exp = 108;         turnon_lexp = 105.;   turnon_hexp = 115.;
    sigma_exp = 5.5;           sigma_lexp = 3.;      sigma_hexp = 8.;
    par1_exp5 = -0.07;          par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.5;       coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;          par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;      par5_hexp5 = 0.;
    coeff5_exp5 = 0.002; coeff5_lexp5 = 0.;    coeff5_hexp5 = 1.;
    }
  }
  else if(cat==2)
  {
    if(order==1){
    sigma_exp = 6;          sigma_lexp = 2.5;      sigma_hexp = 8.;
    turnon_exp = 108.;        turnon_lexp = 105.;   turnon_hexp = 110.;
    par1_exp1 = -0.03;         par1_lexp1 = -0.5;    par1_hexp1 = 0.5;
    coeff1_exp1 = 0.5;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 6.;           sigma_lexp = 2;      sigma_hexp = 10.;
    turnon_exp = 108.5;       turnon_lexp = 105.;   turnon_hexp = 110.;
    par1_exp3 = -7.0678e-02;  par1_lexp3 = -0.5;    par1_hexp3 = 0.5;
    coeff1_exp3 = 8.0177e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.1603e-02;      par3_lexp3 = -0.5;      par3_hexp3 = 0.5;
    coeff3_exp3 = 2.6703e-02; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 6.;           sigma_lexp = 2;      sigma_hexp = 10.;
    turnon_exp = 108.64;      turnon_lexp = 105.;   turnon_hexp = 110.;
    par1_exp5 = -0.07;        par1_lexp5 = -0.5;    par1_hexp5 = 0.5;
    coeff1_exp5 = 0.05; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;      par3_lexp5 = -0.5;    par3_hexp5 = 0.5;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;     par5_hexp5 = 0.5;
    coeff5_exp5 =0.002;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  else if(cat==3)
  {
    if(order==1){
    sigma_exp = 3.85;          sigma_lexp = 2.;      sigma_hexp = 5.;
    turnon_exp = 100.;        turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp1 = -0.05;         par1_lexp1 = -0.5;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.05;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 3.85;           sigma_lexp = 2.;      sigma_hexp = 5.;
    turnon_exp = 100;       turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp3 = -5.4503e-02;  par1_lexp3 = -0.5;    par1_hexp3 = 0.;
    coeff1_exp3 = 9.6120e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.8149e-02;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
    coeff3_exp3 = 1.8388e-07; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 3.85;           sigma_lexp = 2.;      sigma_hexp = 5.;
    turnon_exp = 100;      turnon_lexp = 100.;   turnon_hexp = 109.;
    par1_exp5 = -0.175;        par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.8; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.05;      par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.2; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.00001;      par5_lexp5 = -0.5;     par5_hexp5 = 0.;
    coeff5_exp5 =0.000001;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  else if(cat==4)
  {
    if(order==1){
    sigma_exp = 5.;          sigma_lexp = 1.;      sigma_hexp = 8.;
    turnon_exp = 105.;        turnon_lexp = 100.;   turnon_hexp = 115.;
    par1_exp1 = -0.3;         par1_lexp1 = -1.5;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.5;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 5.;           sigma_lexp = 1.;      sigma_hexp = 10.;
    turnon_exp = 105.5;       turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp3 = -0.5;  par1_lexp3 = -1.5;    par1_hexp3 = 0.;
    coeff1_exp3 = 8.0177e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.1603e-02;      par3_lexp3 = -1.5;      par3_hexp3 = 0.;
    coeff3_exp3 = 2.6703e-02; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 105.64;      turnon_lexp = 100.;   turnon_hexp = 109.;
    par1_exp5 = -0.07;        par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.05; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;      par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;     par5_hexp5 = 0.;
    coeff5_exp5 =0.002;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  else if(cat==501)
  {
    if(order==1){
    sigma_exp = 5;          sigma_lexp = 3.;      sigma_hexp = 8.;
    turnon_exp = 108.;        turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp1 = -0.03;         par1_lexp1 = -0.5;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.5;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.5;       turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp3 = -7.0678e-02;  par1_lexp3 = -0.5;    par1_hexp3 = 0.;
    coeff1_exp3 = 8.0177e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.1603e-02;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
    coeff3_exp3 = 2.6703e-02; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.64;      turnon_lexp = 105.;   turnon_hexp = 109.;
    par1_exp5 = -0.07;        par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.05; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;      par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;     par5_hexp5 = 0.;
    coeff5_exp5 =0.002;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  else if(cat==502)
  {
    if(order==1){
    sigma_exp = 5;          sigma_lexp = 3.;      sigma_hexp = 8.;
    turnon_exp = 108.;        turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp1 = -0.03;         par1_lexp1 = -0.5;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.5;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.5;       turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp3 = -7.0678e-02;  par1_lexp3 = -0.5;    par1_hexp3 = 0.;
    coeff1_exp3 = 8.0177e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.1603e-02;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
    coeff3_exp3 = 2.6703e-02; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.64;      turnon_lexp = 105.;   turnon_hexp = 109.;
    par1_exp5 = -0.07;        par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.05; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;      par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;     par5_hexp5 = 0.;
    coeff5_exp5 =0.002;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  else if(cat==503)
  {
    if(order==1){
    sigma_exp = 5;          sigma_lexp = 3.;      sigma_hexp = 8.;
    turnon_exp = 108.;        turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp1 = -0.03;         par1_lexp1 = -0.5;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.5;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.5;       turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp3 = -7.0678e-02;  par1_lexp3 = -0.5;    par1_hexp3 = 0.;
    coeff1_exp3 = 8.0177e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.1603e-02;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
    coeff3_exp3 = 2.6703e-02; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.64;      turnon_lexp = 105.;   turnon_hexp = 109.;
    par1_exp5 = -0.07;        par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.05; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;      par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;     par5_hexp5 = 0.;
    coeff5_exp5 =0.002;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  else if(cat==6789)
  {
    if(order==1){
    sigma_exp = 5;          sigma_lexp = 3.;      sigma_hexp = 8.;
    turnon_exp = 108.;        turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp1 = -0.03;         par1_lexp1 = -0.5;    par1_hexp1 = 0.;
    coeff1_exp1 = 0.5;      coeff1_lexp1 = 0.;    coeff1_hexp1 = 1.;
    }
    else if(order==3){
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.5;       turnon_lexp = 100.;   turnon_hexp = 110.;
    par1_exp3 = -7.0678e-02;  par1_lexp3 = -0.5;    par1_hexp3 = 0.;
    coeff1_exp3 = 8.0177e-01; coeff1_lexp3 = 0.;    coeff1_hexp3 = 1.;
    par3_exp3 = -4.1603e-02;      par3_lexp3 = -0.5;      par3_hexp3 = 0.;
    coeff3_exp3 = 2.6703e-02; coeff3_lexp3 = 0.;    coeff3_hexp3 = 1.;
    } 
    else if(order==5)
    {
    sigma_exp = 5.;           sigma_lexp = 3.;      sigma_hexp = 10.;
    turnon_exp = 108.64;      turnon_lexp = 105.;   turnon_hexp = 109.;
    par1_exp5 = -0.07;        par1_lexp5 = -0.5;    par1_hexp5 = 0.;
    coeff1_exp5 = 0.05; coeff1_lexp5 = 0.;    coeff1_hexp5 = 1.;
    par3_exp5 = -0.04;      par3_lexp5 = -0.5;    par3_hexp5 = 0.;
    coeff3_exp5 = 0.0002; coeff3_lexp5 = 0.;    coeff3_hexp5 = 1.;
    par5_exp5 = -0.02;      par5_lexp5 = -0.5;     par5_hexp5 = 0.;
    coeff5_exp5 =0.002;  coeff5_lexp5 = 0.; coeff5_hexp5 = 1.;
    }
  }
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),order),Form("%s_sigma_p%d",prefix.c_str(),order),sigma_exp,sigma_lexp,sigma_hexp);
  RooRealVar *turnon = new RooRealVar(Form("%s_turnon_p%d",prefix.c_str(),order),Form("%s_turnon_p%d",prefix.c_str(),order),turnon_exp,turnon_lexp,turnon_hexp);
  
    if (order==1) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp1",prefix.c_str()),Form("%s_p1_exp1",prefix.c_str()),par1_exp1,par1_lexp1,par1_hexp1);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp1",prefix.c_str()),Form("%s_cp1_exp1",prefix.c_str()),coeff1_exp1,coeff1_lexp1,coeff1_hexp1);
    	RooGenericPdf *step = new RooGenericPdf(Form("%s_step_exp1",prefix.c_str()),Form("%s_step_exp1",prefix.c_str()), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(*obs_var,*turnon,*p1,*cp1));
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp1",prefix.c_str()),Form("%s_gau_exp1",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s_gauxexp1",prefix.c_str()),Form("%s_gauxexp1",prefix.c_str()),*obs_var,*step,*gau);
      return gauxexp;
  } else if (order==3) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp3",prefix.c_str()),Form("%s_p1_exp3",prefix.c_str()),par1_exp3,par1_lexp3, par1_hexp3);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp3",prefix.c_str()),Form("%s_cp1_exp3",prefix.c_str()),coeff1_exp3,coeff1_lexp3,coeff1_hexp3);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp3",prefix.c_str()),Form("%s_p3_exp3",prefix.c_str()),par3_exp3,par3_lexp3, par3_hexp3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp3",prefix.c_str()),Form("%s_cp3_exp3",prefix.c_str()),coeff3_exp3,coeff3_lexp3,coeff3_hexp3);
    	RooGenericPdf *step = new RooGenericPdf(Form("%s_step_exp3",prefix.c_str()),Form("%s_step_exp3",prefix.c_str()), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(*obs_var,*turnon,*p1,*cp1,*p3,*cp3));
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp3",prefix.c_str()),Form("%s_gau_exp3",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s_gauxexp3",prefix.c_str()),Form("%s_gauxexp3",prefix.c_str()),*obs_var,*step,*gau);
      return gauxexp;
  } else if (order==5) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp5",prefix.c_str()),Form("%s_p1_exp5",prefix.c_str()),par1_exp5,par1_lexp5, par1_hexp5);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp5",prefix.c_str()),Form("%s_cp1_exp5",prefix.c_str()),coeff1_exp5,coeff1_lexp5,coeff1_hexp5);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp5",prefix.c_str()),Form("%s_p3_exp5",prefix.c_str()),par3_exp5,par3_lexp5, par3_hexp5);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp5",prefix.c_str()),Form("%s_cp3_exp5",prefix.c_str()),coeff3_exp5,coeff3_lexp5,coeff3_hexp5);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_exp5",prefix.c_str()),Form("%s_p5_exp5",prefix.c_str()),par5_exp5,par5_lexp5, par5_hexp5);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_exp5",prefix.c_str()),Form("%s_cp5_exp5",prefix.c_str()),coeff5_exp5,coeff5_lexp5,coeff5_hexp5);
    	RooGenericPdf *step = new RooGenericPdf(Form("%s_step_exp5",prefix.c_str()),Form("%s_step_exp5",prefix.c_str()), "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(*obs_var,*turnon,*p1,*cp1,*p3,*cp3,*p5,*cp5));
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp5",prefix.c_str()),Form("%s_gau_exp5",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s_gauxexp5",prefix.c_str()),Form("%s_gauxexp5",prefix.c_str()),*obs_var,*step,*gau);
      return gauxexp;
  } 
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getPowerLawSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *pows = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.,1.)));
      //params[name]->removeRange();
      fracs->add(*params[name]);
    }
    for (int i=1; i<=npows; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-9.,-1.*(i+1)),-9.,1.)));
      //params[name]->removeRange();
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooPower(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      pows->add(*utilities[ename]);
    }
    //cout << "RooArgLists..." << endl;
    //fracs->Print("v");
    //pows->Print("v");
    //cout << "Function..." << endl;
    RooAbsPdf *pow = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*fracs,true); 
    //pow->Print("v");
    return pow;

    //bing
    RooRealVar *fgaus = new RooRealVar("fgaus", "gaus fraction",0.5,0.,1.) ;
    RooRealVar *mean = new RooRealVar("mean","mean",0,-10.0,10.0) ;
    RooRealVar *sigma = new RooRealVar("sigma","sigma",1,0.,10.) ;
    RooGaussian *gaus = new RooGaussian("gaus","gaus",*obs_var,*mean,*sigma) ;
    //RooAddPdf *pow_gaus = new RooAddPdf("pow_gaus","pow_gaus",RooArgList(*pow,*gaus),*fgaus) ;

    RooRealVar *step_value = new RooRealVar("step_value", "step value", 115., 110., 130.);
    RooGenericPdf *step_func = new RooGenericPdf("step_func", "step_func", "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *pow));
    obs_var->setRange(-400.0,500.0);
    RooFFTConvPdf *pow_gaus = new RooFFTConvPdf("pow_gaus", "pow_gaus", *obs_var, *gaus, *step_func);
    obs_var->setRange(110.0,180.0);
    //return pow_gaus;

    //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));
  }
}

RooAbsPdf* PdfModelBuilder::getLaurentSeries(string prefix, int order){
 
  int nlower=int(ceil(order/2.));
  int nhigher=order-nlower;
  // first do 0th order
  RooArgList *pows = new RooArgList();
  RooArgList *plist = new RooArgList();
  string pname =  Form("%s_pow0",prefix.c_str());
  utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.))));
  pows->add(*utilities[pname]);

  // even terms
  for (int i=1; i<=nlower; i++){
    string name = Form("%s_l%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powl%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.-i))));
    pows->add(*utilities[pname]);
  }
  // odd terms
  for (int i=1; i<=nhigher; i++){
    string name = Form("%s_h%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powh%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.+i))));
    pows->add(*utilities[pname]);
  }
  RooAddPdf *pdf = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*plist,true);
  return pdf;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));

  //bing
  RooRealVar *mean1 = new RooRealVar("mean1","mean1",0.0) ;
    RooRealVar *sigma1 = new RooRealVar("sigma1","sigma1",5,-10.,20.) ;
    RooGaussian *gaus1 = new RooGaussian("gaus1","gaus1",*obs_var,*mean1,*sigma1) ;
    RooRealVar *step_value1 = new RooRealVar("step_value1", "step_value1",115.,110.,130.) ;
    RooGenericPdf *step_func1 = new RooGenericPdf("step_func1","step_func1","(1e-20+( @0 > @1)) * @2",RooArgSet(*obs_var,*step_value1,*pdf));
    obs_var->setRange(-400.0,500.0);
    RooFFTConvPdf *pdf_gaus = new RooFFTConvPdf("pdf_gaus","pdf_gaus", *obs_var, *gaus1, *step_func1);
    obs_var->setRange(110.0,180.0);
  //return pdf_gaus;
}

RooAbsPdf* PdfModelBuilder::getLaurentStepxGau(string prefix, int order, int cat){
  //if(order>2) return NULL;
  
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  double sigma_lau,sigma_llau,sigma_hlau;
  double turnon_lau,turnon_llau,turnon_hlau;
  double coeff1_lau1, coeff2_lau1, coeff1_lau2, coeff2_lau2, coeff3_lau2, coeff1_lau3, coeff2_lau3, coeff3_lau3, coeff4_lau3, coeff1_lau4, coeff2_lau4, coeff3_lau4, coeff4_lau4, coeff5_lau4; 
  double coeff1_hlau1, coeff2_hlau1, coeff1_hlau2, coeff2_hlau2, coeff3_hlau2, coeff1_hlau3, coeff2_hlau3, coeff3_hlau3, coeff4_hlau3, coeff1_hlau4, coeff2_hlau4, coeff3_hlau4, coeff4_hlau4, coeff5_hlau4; 
  double coeff1_llau1, coeff2_llau1, coeff1_llau2, coeff2_llau2, coeff3_llau2, coeff1_llau3, coeff2_llau3, coeff3_llau3, coeff4_llau3, coeff1_llau4, coeff2_llau4, coeff3_llau4, coeff4_llau4, coeff5_llau4; 

 if(cat==1)
  {
    if(order==1){
    turnon_lau = 108;      turnon_llau = 104.;   turnon_hlau = 115.;
    sigma_lau = 4.;          sigma_llau = 2.;      sigma_hlau = 8.;
    coeff1_lau1 = 5.8468e-08; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.5;
    coeff2_lau1 = 0.25; coeff2_llau1 = 0.;    coeff2_hlau1 = 1.;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 106.;   turnon_hlau = 115.;
    sigma_lau = 4.;          sigma_llau = 2.;      sigma_hlau = 8.;
    coeff1_lau2 = 1.2503e-06; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.5;
    coeff2_lau2 = 9.7549e-01; coeff2_llau2 = 0.1;    coeff2_hlau2 = 1.;
    coeff3_lau2 = 0.; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 106;   turnon_hlau = 115.;
    sigma_lau = 4.;           sigma_llau = 2.;      sigma_hlau = 8.;
    coeff1_lau3 = 5.e-09; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.001;
    coeff2_lau3 = 4.e-08; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.1;
    coeff3_lau3 = 7e-11; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.0001;
    coeff4_lau3 = 0.5; coeff4_llau3 = 0.;    coeff4_hlau3 = 1.;
    }
    else if(order==4)
    {
    turnon_lau = 106.5;         turnon_llau = 104.;   turnon_hlau = 115.;
    sigma_lau = 4.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau4 = 1.8872e-15; coeff1_llau4 = 0.;    coeff1_hlau4 = 0.01;
    coeff2_lau4 = 3.2066e-13; coeff2_llau4 = 0.;    coeff2_hlau4 = 0.1;
    coeff3_lau4 = 9.5240e-17; coeff3_llau4 = 0.;    coeff3_hlau4 = 0.0001;
    coeff4_lau4 = 9.9935e-01; coeff4_llau4 = 0.1;    coeff4_hlau4 = 1.;
    coeff5_lau4 = 9.9935e-01; coeff5_llau4 = 0.9;    coeff5_hlau4 = 1.;
    }
  }
  else if(cat==2)
  {
    if(order==1){
    turnon_lau = 108;      turnon_llau = 106.;   turnon_hlau = 109.;
    sigma_lau = 4.;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau1 = 5.8468e-08; coeff1_llau1 = 0.;    coeff1_hlau1 = 1.;
    coeff2_lau1 = 0.25; coeff2_llau1 = 0.;    coeff2_hlau1 = 1.;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 106.;   turnon_hlau = 109.;
    sigma_lau = 4.;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau2 = 1.2503e-06; coeff1_llau2 = 0.;    coeff1_hlau2 = 1.;
    coeff2_lau2 = 9.7549e-01; coeff2_llau2 = 0.1;    coeff2_hlau2 = 1.;
    coeff3_lau2 = 0.; coeff3_llau2 = 0.;    coeff3_hlau2 = 1.;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 106;   turnon_hlau = 109.;
    sigma_lau = 4.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau3 = 5.e-09; coeff1_llau3 = 0.;    coeff1_hlau3 = 1.;
    coeff2_lau3 = 4.e-08; coeff2_llau3 = 0.;    coeff2_hlau3 = 1.;
    coeff3_lau3 = 7e-11; coeff3_llau3 = 0.;    coeff3_hlau3 = 1.;
    coeff4_lau3 = 0.5; coeff4_llau3 = 0.;    coeff4_hlau3 = 1.;
    }
    else if(order==4)
    {
    turnon_lau = 106.5;         turnon_llau = 104.;   turnon_hlau = 110.;
    sigma_lau = 4.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau4 = 1.8872e-15; coeff1_llau4 = 0.;    coeff1_hlau4 = 1.;
    coeff2_lau4 = 3.2066e-13; coeff2_llau4 = 0.;    coeff2_hlau4 = 1.;
    coeff3_lau4 = 9.5240e-17; coeff3_llau4 = 0.;    coeff3_hlau4 = 1.;
    coeff4_lau4 = 9.9935e-01; coeff4_llau4 = 0.1;    coeff4_hlau4 = 1.;
    coeff5_lau4 = 9.9935e-01; coeff5_llau4 = 0.9;    coeff5_hlau4 = 1.;
    }
  }
  else if(cat==3)
  {
    if(order==1){
    turnon_lau = 108.32;      turnon_llau = 104.;   turnon_hlau = 120.;
    sigma_lau = 6.7;          sigma_llau = 1.;      sigma_hlau = 8.;
    coeff1_lau1 = 3.6576e-04; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.05;
    coeff2_lau1 = 3.6576e-04; coeff2_llau1 = 0.;    coeff2_hlau1 = 0.5;
    }
    else if(order==2){
    turnon_lau = 106.2;       turnon_llau = 104.;   turnon_hlau = 120.;
    sigma_lau = 6.6;          sigma_llau = 1.;      sigma_hlau = 8.;
    coeff1_lau2 = 3.6576e-04; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.005;
    coeff2_lau2 = 3.6576e-04; coeff2_llau2 = 0.;    coeff2_hlau2 = 0.01;
    coeff3_lau2 = 3.6576e-04; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 106.5;         turnon_llau = 104.;   turnon_hlau = 120.;
    sigma_lau = 3.;           sigma_llau = 1.;      sigma_hlau = 8.;
    coeff1_lau3 = 1.8872e-15; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.01;
    coeff2_lau3 = 3.2066e-13; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.1;
    coeff3_lau3 = 9.5240e-17; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.0001;
    coeff4_lau3 = 9.9935e-01; coeff4_llau3 = 0.9;    coeff4_hlau3 = 1.;
    }
    else if(order==4)
    {
    turnon_lau = 106.5;         turnon_llau = 104.;   turnon_hlau = 120.;
    sigma_lau = 3.;           sigma_llau = 1.;      sigma_hlau = 8.;
    coeff1_lau4 = 1.8872e-15; coeff1_llau4 = 0.;    coeff1_hlau4 = 0.01;
    coeff2_lau4 = 3.2066e-13; coeff2_llau4 = 0.;    coeff2_hlau4 = 0.1;
    coeff3_lau4 = 9.5240e-17; coeff3_llau4 = 0.;    coeff3_hlau4 = 0.0001;
    coeff4_lau4 = 9.9935e-01; coeff4_llau4 = 0.1;    coeff4_hlau4 = 1.;
    coeff5_lau4 = 9.9935e-01; coeff5_llau4 = 0.9;    coeff5_hlau4 = 1.;
    }
  }
  else if(cat==4)
  {
    if(order==1){
    turnon_lau = 108.32;      turnon_llau = 100.;   turnon_hlau = 120.;
    sigma_lau = 6.7;          sigma_llau = 5.;      sigma_hlau = 8.;
    coeff1_lau1 = 3.6576e-04; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.5;
    coeff2_lau1 = 3.6576e-04; coeff2_llau1 = 0.;    coeff2_hlau1 = 0.5;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.6;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau2 = 3.6576e-04; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.5;
    coeff2_lau2 = 3.6576e-04; coeff2_llau2 = 0.;    coeff2_hlau2 = 0.5;
    coeff3_lau2 = 3.6576e-04; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau3 = 3.6576e-04; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.5;
    coeff2_lau3 = 3.6576e-04; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.5;
    coeff3_lau3 = 3.6576e-04; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.5;
    coeff4_lau3 = 3.6576e-04; coeff4_llau3 = 0.;    coeff4_hlau3 = 0.5;
    }
    else if(order==4)
    {
    turnon_lau = 106.5;         turnon_llau = 104.;   turnon_hlau = 120.;
    sigma_lau = 3.;           sigma_llau = 1.;      sigma_hlau = 8.;
    coeff1_lau4 = 1.8872e-15; coeff1_llau4 = 0.;    coeff1_hlau4 = 0.01;
    coeff2_lau4 = 3.2066e-13; coeff2_llau4 = 0.;    coeff2_hlau4 = 0.1;
    coeff3_lau4 = 9.5240e-17; coeff3_llau4 = 0.;    coeff3_hlau4 = 0.0001;
    coeff4_lau4 = 9.9935e-01; coeff4_llau4 = 0.1;    coeff4_hlau4 = 1.;
    coeff5_lau4 = 9.9935e-01; coeff5_llau4 = 0.9;    coeff5_hlau4 = 1.;
    }
  }
  else if(cat==501)
  {
    if(order==1){
    turnon_lau = 108.32;      turnon_llau = 100.;   turnon_hlau = 120.;
    sigma_lau = 6.7;          sigma_llau = 5.;      sigma_hlau = 8.;
    coeff1_lau1 = 3.6576e-04; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.5;
    coeff2_lau1 = 3.6576e-04; coeff2_llau1 = 0.;    coeff2_hlau1 = 0.5;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.6;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau2 = 3.6576e-04; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.5;
    coeff2_lau2 = 3.6576e-04; coeff2_llau2 = 0.;    coeff2_hlau2 = 0.5;
    coeff3_lau2 = 3.6576e-04; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau3 = 3.6576e-04; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.5;
    coeff2_lau3 = 3.6576e-04; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.5;
    coeff3_lau3 = 3.6576e-04; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.5;
    coeff4_lau3 = 3.6576e-04; coeff4_llau3 = 0.;    coeff4_hlau3 = 0.5;
    }
  }
  else if(cat==502)
  {
    if(order==1){
    turnon_lau = 108.32;      turnon_llau = 100.;   turnon_hlau = 120.;
    sigma_lau = 6.7;          sigma_llau = 5.;      sigma_hlau = 8.;
    coeff1_lau1 = 3.6576e-04; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.5;
    coeff2_lau1 = 3.6576e-04; coeff2_llau1 = 0.;    coeff2_hlau1 = 0.5;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.6;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau2 = 3.6576e-04; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.5;
    coeff2_lau2 = 3.6576e-04; coeff2_llau2 = 0.;    coeff2_hlau2 = 0.5;
    coeff3_lau2 = 3.6576e-04; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau3 = 3.6576e-04; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.5;
    coeff2_lau3 = 3.6576e-04; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.5;
    coeff3_lau3 = 3.6576e-04; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.5;
    coeff4_lau3 = 3.6576e-04; coeff4_llau3 = 0.;    coeff4_hlau3 = 0.5;
    }
  }
  else if(cat==503)
  {
    if(order==1){
    turnon_lau = 108.32;      turnon_llau = 100.;   turnon_hlau = 120.;
    sigma_lau = 6.7;          sigma_llau = 5.;      sigma_hlau = 8.;
    coeff1_lau1 = 3.6576e-04; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.5;
    coeff2_lau1 = 3.6576e-04; coeff2_llau1 = 0.;    coeff2_hlau1 = 0.5;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.6;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau2 = 3.6576e-04; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.5;
    coeff2_lau2 = 3.6576e-04; coeff2_llau2 = 0.;    coeff2_hlau2 = 0.5;
    coeff3_lau2 = 3.6576e-04; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau3 = 3.6576e-04; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.5;
    coeff2_lau3 = 3.6576e-04; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.5;
    coeff3_lau3 = 3.6576e-04; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.5;
    coeff4_lau3 = 3.6576e-04; coeff4_llau3 = 0.;    coeff4_hlau3 = 0.5;
    }
  }
  else if(cat==6789)
  {
    if(order==1){
    turnon_lau = 108.32;      turnon_llau = 100.;   turnon_hlau = 120.;
    sigma_lau = 6.7;          sigma_llau = 5.;      sigma_hlau = 8.;
    coeff1_lau1 = 3.6576e-04; coeff1_llau1 = 0.;    coeff1_hlau1 = 0.5;
    coeff2_lau1 = 3.6576e-04; coeff2_llau1 = 0.;    coeff2_hlau1 = 0.5;
    }
    else if(order==2){
    turnon_lau = 108.2;       turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.6;          sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau2 = 3.6576e-04; coeff1_llau2 = 0.;    coeff1_hlau2 = 0.5;
    coeff2_lau2 = 3.6576e-04; coeff2_llau2 = 0.;    coeff2_hlau2 = 0.5;
    coeff3_lau2 = 3.6576e-04; coeff3_llau2 = 0.;    coeff3_hlau2 = 0.5;
    } 
    else if(order==3)
    {
    turnon_lau = 108;         turnon_llau = 107.;   turnon_hlau = 120.;
    sigma_lau = 6.;           sigma_llau = 3.;      sigma_hlau = 8.;
    coeff1_lau3 = 3.6576e-04; coeff1_llau3 = 0.;    coeff1_hlau3 = 0.5;
    coeff2_lau3 = 3.6576e-04; coeff2_llau3 = 0.;    coeff2_hlau3 = 0.5;
    coeff3_lau3 = 3.6576e-04; coeff3_llau3 = 0.;    coeff3_hlau3 = 0.5;
    coeff4_lau3 = 3.6576e-04; coeff4_llau3 = 0.;    coeff4_hlau3 = 0.5;
    }
  }
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),order),Form("%s_sigma_p%d",prefix.c_str(),order),sigma_lau,sigma_llau,sigma_hlau);
  RooRealVar *turnon = new RooRealVar(Form("%s_turnon_p%d",prefix.c_str(),order),Form("%s_turnon_p%d",prefix.c_str(),order),turnon_lau,turnon_llau,turnon_hlau);
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  
    if (order==1) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau1",prefix.c_str()),Form("%s_cp1_lau1",prefix.c_str()),coeff1_lau1,coeff1_llau1,coeff1_hlau1);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau1",prefix.c_str()),Form("%s_cp2_lau1",prefix.c_str()),coeff2_lau1,coeff2_llau1,coeff2_hlau1);
      RooGenericPdf *step= new RooGenericPdf(Form("%s_step_lau1",prefix.c_str()),Form("%s_step_lau1",prefix.c_str()), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5))", RooArgList(*obs_var,*turnon,*cp1,*cp2));//step*(ax^b)
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau1",prefix.c_str()),Form("%s_gau_lau1",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s_gauxlau1",prefix.c_str()),Form("%s_gauxlau1",prefix.c_str()),*obs_var,*step,*gau);
      return gauxlau;
  } else if (order==2) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau2",prefix.c_str()),Form("%s_cp1_lau2",prefix.c_str()),coeff1_lau2,coeff1_llau2,coeff1_hlau2);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp1_lau2",prefix.c_str()),Form("%s_cp1_lau2",prefix.c_str()),coeff2_lau2,coeff2_llau2,coeff2_hlau2);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau2",prefix.c_str()),Form("%s_cp3_lau2",prefix.c_str()),coeff3_lau2,coeff3_llau2,coeff3_hlau2);
      RooGenericPdf *step = new RooGenericPdf(Form("%s_step_lau2",prefix.c_str()), Form("%s_step_lau2",prefix.c_str()), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3))", RooArgList(*obs_var,*turnon,*cp1,*cp2,*cp3));//step*(ax^b+cx^d+fx^g) 
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau2",prefix.c_str()),Form("%s_gau_lau2",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s_gauxlau2",prefix.c_str()),Form("%s_gauxlau2",prefix.c_str()),*obs_var,*step,*gau);
      return gauxlau;
  } else if (order==3) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau3",prefix.c_str()),Form("%s_cp1_lau3",prefix.c_str()),coeff1_lau3,coeff1_llau3,coeff1_hlau3);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau3",prefix.c_str()),Form("%s_cp2_lau3",prefix.c_str()),coeff2_lau3,coeff2_llau3,coeff2_hlau3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau3",prefix.c_str()),Form("%s_cp3_lau3",prefix.c_str()),coeff3_lau3,coeff3_llau3,coeff3_hlau3);
      RooRealVar *cp4 = new RooRealVar(Form("%s_cp4_lau3",prefix.c_str()),Form("%s_cp4_lau3",prefix.c_str()),coeff4_lau3,coeff4_llau3,coeff4_hlau3);
      RooGenericPdf *step = new RooGenericPdf(Form("%s_step_lau3",prefix.c_str()),Form("%s_step_lau3",prefix.c_str()), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))", RooArgList(*obs_var,*turnon,*cp1,*cp2,*cp3,*cp4));//step*(ax^b+cx^d)
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau3",prefix.c_str()),Form("%s_gau_lau3",prefix.c_str()),*obs_var,*mean,*sigma);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s_gauxlau3",prefix.c_str()),Form("%s_gauxlau3",prefix.c_str()),*obs_var,*step,*gau);
      return gauxlau;
  } else if (order==4) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau4",prefix.c_str()),Form("%s_cp1_lau4",prefix.c_str()),coeff1_lau4,coeff1_llau4,coeff1_hlau4);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau4",prefix.c_str()),Form("%s_cp2_lau4",prefix.c_str()),coeff2_lau4,coeff2_llau4,coeff2_hlau4);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau4",prefix.c_str()),Form("%s_cp3_lau4",prefix.c_str()),coeff3_lau4,coeff3_llau4,coeff3_hlau4);
      RooRealVar *cp4 = new RooRealVar(Form("%s_cp4_lau4",prefix.c_str()),Form("%s_cp4_lau4",prefix.c_str()),coeff4_lau4,coeff4_llau4,coeff4_hlau4);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_lau4",prefix.c_str()),Form("%s_cp5_lau4",prefix.c_str()),coeff5_lau4,coeff5_llau4,coeff5_hlau4);
      RooGenericPdf *step = new RooGenericPdf(Form("%s_step_lau4",prefix.c_str()),Form("%s_step_lau4",prefix.c_str()), "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-7))", RooArgList(*obs_var,*turnon,*cp1,*cp2,*cp3,*cp4,*cp5));//step*(ax^b+cx^d)
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau4",prefix.c_str()),Form("%s_gau_lau4",prefix.c_str()),*obs_var,*mean,*sigma);

      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s_gauxlau4",prefix.c_str()),Form("%s_gauxlau4",prefix.c_str()),*obs_var,*step,*gau);
      return gauxlau;
  }
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getBernsteinxZGMCShape(string prefix, int order){
  
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/ZGMCShape/ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");

  RooArgList *coeffList = new RooArgList();
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),-5.,5.);
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1)*pow(-1,i),-5.,5.);//bing
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

RooAbsPdf* PdfModelBuilder::getExponentialZGMCShape(string prefix, int order){
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/ZGMCShape/ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");

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
  par3_exp3 = -0.001;      par3_lexp3 = -0.01;      par3_hexp3 = 0.;
  coeff3_exp3 = 0.4;      coeff3_lexp3 = -100.;    coeff3_hexp3 = 100.;
  } 
  else if(order==5)
  {
  par1_exp5 = -0.1;          par1_lexp5 = -0.2;    par1_hexp5 = 0.;
  coeff1_exp5 = 0.9;       coeff1_lexp5 = -100.;    coeff1_hexp5 = 100.;
  par3_exp5 = -0.001;          par3_lexp5 = -0.5;    par3_hexp5 = 0.;
  coeff3_exp5 = 0.4; coeff3_lexp5 = -100.;    coeff3_hexp5 = 100.;
  par5_exp5 = -0.005;      par5_lexp5 = -0.5;      par5_hexp5 = 0.;
  coeff5_exp5 = 0.4; coeff5_lexp5 = -100.;    coeff5_hexp5 = 100.;
  }
  else if(order==7)
  {
  par1_exp7 = -0.1;          par1_lexp7 = -0.2;    par1_hexp7 = 0.;
  coeff1_exp7 = 0.9;       coeff1_lexp7 = -100.;    coeff1_hexp7 = 100.;
  par3_exp7 = -0.001;          par3_lexp7 = -0.5;    par3_hexp7 = 0.;
  coeff3_exp7 = 0.4; coeff3_lexp7 = -100.;    coeff3_hexp7 = 100.;
  par5_exp7 = -0.005;      par5_lexp7 = -0.5;      par5_hexp7 = 0.;
  coeff5_exp7 = 0.4; coeff5_lexp7 = -100.;    coeff5_hexp7 = 100.;
  par7_exp7 = -0.005;      par7_lexp7 = -0.5;      par7_hexp7 = 0.;
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


RooAbsPdf* PdfModelBuilder::getPowerLawZGMCShape(string prefix, int order){
  
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/ZGMCShape/ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");

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
  par3_pow3 = -3.;      par3_lpow3 = -4;      par3_hpow3 = 0.;
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

RooAbsPdf* PdfModelBuilder::getLaurentZGMCShape(string prefix, int order){
 
  //bing add ZGMCShape
  TFile *ZGMC_file = new TFile("/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/ZGMCShape/ZGMCShape.root");
  RooWorkspace *w = (RooWorkspace *)ZGMC_file->Get("w");
  RooAbsPdf *ZGMCShape = w->pdf("ZGMCShape");

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


RooAbsPdf* PdfModelBuilder::getKeysPdf(string prefix){
  if (!keysPdfAttributesSet){
    cerr << "ERROR -- keysPdf attributes not set" << endl;
    exit(1);
  }
  return new RooKeysPdf(prefix.c_str(),prefix.c_str(),*obs_var,*keysPdfData,RooKeysPdf::MirrorBoth,keysPdfRho);
}

RooAbsPdf* PdfModelBuilder::getPdfFromFile(string &prefix){
  vector<string> details;
  split(details,prefix,boost::is_any_of(","));

  string fname = details[2];
  string wsname = details[1];
  string pdfname = details[0];

  TFile *tempFile = TFile::Open(fname.c_str());
  if (!tempFile){
    cerr << "PdfModelBuilder::getPdfFromFile -- file not found " << fname << endl;
    assert(0);
  }
  RooWorkspace *tempWS = (RooWorkspace*)tempFile->Get(wsname.c_str());
  if (!tempWS){
    cerr << "PdfModelBuilder::getPdfFromFile -- workspace not found " << wsname << endl;
    assert(0);
  }
  RooAbsPdf *tempPdf = (RooAbsPdf*)tempWS->pdf(pdfname.c_str());
  if (!tempPdf){
    cerr << "PdfModelBuilder::getPdfFromFile -- pdf not found " << pdfname << endl;
    assert(0);
  }
  prefix = pdfname;
  RooAbsPdf *pdf = (RooAbsPdf*)tempPdf->Clone(prefix.c_str());
  tempFile->Close();
  delete tempFile;
  return pdf;
}

RooAbsPdf* PdfModelBuilder::getExponentialSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addExponential -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int nexps=order-nfracs;
    assert(nfracs==nexps-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *exps = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.0001,0.9999)));
      fracs->add(*params[name]);
    }
    for (int i=1; i<=nexps; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-1.,-0.04*(i+1)),-1.,0.)));
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooExponential(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      exps->add(*utilities[ename]);
    }
    //fracs->Print("v");
    //exps->Print("v");
    RooAbsPdf *exp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*exps,*fracs,true);

    //bing
    RooRealVar *mean = new RooRealVar("mean","mean",0.0) ;
    RooRealVar *sigma = new RooRealVar("sigma","sigma",5,-10.,20.) ;
    RooGaussian *gaus = new RooGaussian("gaus","gaus",*obs_var,*mean,*sigma) ;
    RooRealVar *step_value = new RooRealVar("step_value", "step value",115.,100.,130.) ;
    // RooGenericPdf *step_func = new RooGenericPdf("step_func","step_func","(abs(step_value-obs_var)/(step_value-obs_var)+1.0)/2.0",RooArgSet(*obs_var,*step_value));
    RooGenericPdf *step_func = new RooGenericPdf("step_func","step_func","(1e-20+( @0 > @1)) * @2",RooArgSet(*obs_var,*step_value,*exp));
    //RooAddPdf *exp_gaus = new RooAddPdf("exp_gaus","exp_gaus",RooArgList(*exp,*gaus),*fgaus) ;
    obs_var->setRange(-200.0,200.0);
    //cout << "-----------[[test1]]---------------" << endl;
    //obs_var->Print("v");
    //cout << "-----------[[test3]]---------------" << endl;
    RooFFTConvPdf *exp_gaus = new RooFFTConvPdf("exp_gaus","exp_gaus", *obs_var, *gaus, *step_func);
    obs_var->setRange(110.0,180.0);
    //cout << "-----------[[test2]]---------------" << endl;
    //obs_var->Print("v");
    //exp->Print("v");
    cout << "--------------------------" << endl;
    return exp;
    //bing
    //return exp_gaus;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(exp->GetName(),exp));

  }
}


void PdfModelBuilder::addBkgPdf(string type, int nParams, string name, bool cache){
 
  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf=0;// avoid uninitialised variable error in cmssw  

  if (type=="Bernstein") pdf = getBernstein(name,nParams);
  if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
  if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
  if (type=="Laurent") pdf = getLaurentSeries(name,nParams);
  if (type=="KeysPdf") pdf = getKeysPdf(name);
  if (type=="File") pdf = getPdfFromFile(name);

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }

}

void PdfModelBuilder::setKeysPdfAttributes(RooDataSet *data, double rho){
  keysPdfData = data;
  keysPdfRho = rho;
  keysPdfAttributesSet=true;
}

void PdfModelBuilder::setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm){
  sigPdf=pdf;
  sigNorm=norm;
  signal_set=true;
}

void PdfModelBuilder::setSignalPdfFromMC(RooDataSet *data){
  
  RooDataHist *sigMCBinned = new RooDataHist(Form("roohist_%s",data->GetName()),Form("roohist_%s",data->GetName()),RooArgSet(*obs_var),*data);
  sigPdf = new RooHistPdf(Form("pdf_%s",data->GetName()),Form("pdf_%s",data->GetName()),RooArgSet(*obs_var),*sigMCBinned);
  sigNorm = new RooConstVar(Form("sig_events_%s",data->GetName()),Form("sig_events_%s",data->GetName()),data->sumEntries());
  signal_set=true;
}

void PdfModelBuilder::makeSBPdfs(bool cache){
  
  if (!signal_set){
    cerr << "ERROR - no signal model set!" << endl;
    exit(1);
  }
  if (!signal_modifier_set){
    cerr << "ERROR - no signal modifier set!" << endl;
    exit(1);
  }
 
  if (sigNorm) {
    sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*signalModifier,*sigNorm));
  }
  else {
    sigYield = signalModifier;
  }
  bkgYield = new RooRealVar("bkg_yield","bkg_yield",1000.,0.,1.e6);

  for (map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++){
    RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.c_str()),Form("sb_%s",bkg->first.c_str()),RooArgList(*(bkg->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
    if (cache) {
      wsCache->import(*sbMod,RecycleConflictNodes());
      RooAbsPdf *cachePdf = (RooAbsPdf*)wsCache->pdf(sbMod->GetName());
      signalModifier = (RooRealVar*)wsCache->var(signalModifier->GetName());
      sbPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
    }
    else {
      sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
    }
  }
}

map<string,RooAbsPdf*> PdfModelBuilder::getBkgPdfs(){
  return bkgPdfs;
}

map<string,RooAbsPdf*> PdfModelBuilder::getSBPdfs(){
  return sbPdfs;
}

RooAbsPdf* PdfModelBuilder::getSigPdf(){
  return sigPdf;
}

void PdfModelBuilder::plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly,string specificPdfName){
  
  TCanvas *canv = new TCanvas();
  bool specPdf=false;
  if (specificPdfName!="") specPdf=true;

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (specPdf && it->first!=specificPdfName && specificPdfName!="NONE") continue;
    RooPlot *plot = obs_var->frame();
    data->plotOn(plot,Binning(binning));
    if (specificPdfName!="NONE") {
	 it->second->plotOn(plot);
	 it->second->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
    }	
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",name.c_str(),it->first.c_str()));
    canv->Print(Form("%s_%s.png",name.c_str(),it->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::fitToData(RooAbsData *data, bool bkgOnly, bool cache, bool print){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
    if (print){
      cout << "Fit Res Before: " << endl;
      fit->floatParsInit().Print("v");
      cout << "Fit Res After: " << endl;
      fit->floatParsFinal().Print("v");
    }
    if (cache) {
      RooArgSet *fitargs = (RooArgSet*)it->second->getParameters(*obs_var);
      // remove the signal strength since this will be set AFTER fitting the background 
      fitargs->remove(*signalModifier); 
      wsCache->defineSet(Form("%s_params",it->first.c_str()),*fitargs);
      wsCache->defineSet(Form("%s_observs",it->first.c_str()),*obs_var);
      wsCache->saveSnapshot(it->first.c_str(),*fitargs,true);
      if (print) {
        cout << "Cached values: " << endl;
        fitargs->Print("v");
      }
    }
  }
  if (bkgOnly) bkgHasFit=true;
  else sbHasFit=true;
}

void PdfModelBuilder::setSeed(int seed){
  RooRandom::randomGenerator()->SetSeed(seed);
}

RooDataSet* PdfModelBuilder::makeHybridDataset(vector<float> switchOverMasses, vector<RooDataSet*> dataForHybrid){
  
  assert(switchOverMasses.size()==dataForHybrid.size()-1);

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  obs_var->Print("v");
  assert(cut_strings.size()==dataForHybrid.size());
  
	RooDataSet *data=0;// avoid uninitialised variable error in cmssw  
  for (unsigned int i=0; i<dataForHybrid.size(); i++){
    RooDataSet *cutData = (RooDataSet*)dataForHybrid[i]->reduce(Name("hybridToy"),Title("hybridToy"),CutRange(cut_strings[i].c_str()));
    //RooDataSet *cutData = new RooDataSet("hybridToy","hybridToy",RooArgSet(*obs_var),Import(*dataForHybrid[i]),CutRange(cut_strings[i].c_str()));
    if (i==0) data=cutData;
    else data->append(*cutData);
  }
  return data;
}

void PdfModelBuilder::throwHybridToy(string postfix, int nEvents, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly, bool binned, bool poisson, bool cache){
  
  assert(switchOverMasses.size()==functions.size()-1);
  toyHybridData.clear();

  // have to throw unbinned for the hybrid
  throwToy(postfix,nEvents,bkgOnly,false,poisson,cache);

  vector<RooDataSet*> dataForHybrid;
  string hybridName = "hybrid";
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    hybridName += "_"+*func;
    for (map<string,RooDataSet*>::iterator it=toyDataSet.begin(); it!=toyDataSet.end(); it++){
      if (it->first.find(*func)!=string::npos){
        dataForHybrid.push_back(it->second);
      }
    }
  }
  if (dataForHybrid.size()!=functions.size()){
    cerr << "One of the requested hybrid functions has not been found" << endl;
    exit(1);
  }

  RooDataSet *hybridData = makeHybridDataset(switchOverMasses,dataForHybrid);
  toyHybridData.clear();
  if (binned) {
    RooDataHist *hybridDataHist = hybridData->binnedClone();
    hybridDataHist->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridDataHist->GetName(),hybridDataHist));
  }
  else {
    hybridData->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridData->GetName(),hybridData));
  }
}

void PdfModelBuilder::throwToy(string postfix, int nEvents, bool bkgOnly, bool binned, bool poisson, bool cache){

  toyData.clear();
  toyDataSet.clear();
  toyDataHist.clear();
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
    if (!bkgHasFit) cerr << "WARNING -- bkg has not been fit to data. Are you sure this is wise?" << endl; 
  }
  else {
    pdfSet = sbPdfs;
    if (!sbHasFit) cerr << "WARNING -- sb has not been fit to data. Are you sure this is wise?" << endl;
  }
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (cache) {
      wsCache->loadSnapshot(it->first.c_str());
      cout << "Loaded snapshot, params at.." << endl;
      it->second->getVariables()->Print("v");
    }
    RooAbsData *toy;
    if (binned){
      RooDataHist *toyHist;
      if (poisson) toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataHist.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));
      toy=toyHist;
    }
    else {
      RooDataSet *toySet;
      if (poisson) toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataSet.insert(pair<string,RooDataSet*>(toySet->GetName(),toySet));
      toy=toySet;
    }
    toyData.insert(pair<string,RooAbsData*>(toy->GetName(),toy));
  }
  
}

map<string,RooAbsData*> PdfModelBuilder::getToyData(){
  return toyData;
}

map<string,RooAbsData*> PdfModelBuilder::getHybridToyData(){
  return toyHybridData;
}

void PdfModelBuilder::plotHybridToy(string prefix, int binning, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly){

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  
  int tempColors[4] = {kBlue,kRed,kGreen+2,kMagenta};

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  RooPlot *plot = obs_var->frame();
  TCanvas *canv = new TCanvas();
  int i=0;
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
      // check if in list of hybrid functions
      if (pdfIt->first.find(*func)!=string::npos) {
        for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
          //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
          if (toyIt->first.find(pdfIt->first)!=string::npos){
            RooAbsData *data = toyIt->second->reduce(CutRange(cut_strings[i].c_str()));
            data->plotOn(plot,Binning(binning),MarkerColor(tempColors[i]),LineColor(tempColors[i]),CutRange(cut_strings[i].c_str()));
            pdfIt->second->plotOn(plot,LineColor(tempColors[i]),Range(cut_strings[i].c_str()));
            i++;
          }
        }
      }
    }
  }
  for (map<string,RooAbsData*>::iterator hybrid=toyHybridData.begin(); hybrid!=toyHybridData.end(); hybrid++){
    hybrid->second->plotOn(plot,Binning(binning),MarkerSize(0.8),MarkerStyle(kFullSquare));
    plot->SetMinimum(0.0001);
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",prefix.c_str(),hybrid->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::plotToysWithPdfs(string prefix, int binning, bool bkgOnly){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  TCanvas *canv = new TCanvas();
  for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
    for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
      //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
      if (toyIt->first.find(pdfIt->first)!=string::npos){
        RooPlot *plot = obs_var->frame();
        toyIt->second->plotOn(plot,Binning(binning));
        pdfIt->second->plotOn(plot,LineColor(kRed));
        pdfIt->second->paramOn(plot,LineColor(kRed),RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        plot->Draw();
        canv->Print(Form("%s_%s.pdf",prefix.c_str(),pdfIt->first.c_str()));
        canv->Print(Form("%s_%s.png",prefix.c_str(),pdfIt->first.c_str()));
      }
    }
  }
  delete canv;

}

void PdfModelBuilder::saveWorkspace(string filename){
  
  TFile *outFile = new TFile(filename.c_str(),"RECREATE");
  outFile->cd();
  wsCache->Write();
  outFile->Close();
  delete outFile;
}

void PdfModelBuilder::saveWorkspace(TFile *file){
  file->cd();
  wsCache->Write();
}
