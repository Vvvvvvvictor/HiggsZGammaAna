#include "Background/interface/PdfModelBuilder.h"
#include "Background/interface/WSTFileWrapper.h"
#include "RooChi2Var.h"
#include "RooRealIntegral.h"
#include "RooAbsPdf.h"

using namespace std;
using namespace RooFit;


void test(){
    TFile file = TFile("/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Background/HZGamma_BkgModel_UL_xingchen_chi2/fit_results_run2_cat3/CMS-HGG_mva_13TeV_multipdf_cat3.root");
    RooWorkspace *inWS = (RooWorkspace*)file.Get("multipdf");

    RooRealVar *mzg = (RooRealVar*)inWS->var("CMS_hzg_mass");//FIXED
    RooDataHist *hist_bkg = (RooDataHist*)inWS->data("roohist_data_mass_cat0");//FIXED
    RooAbsPdf *pdf_bkg = (RooAbsPdf*)inWS->pdf("CMS_hgg_untag_cat3_bkgshape");//FIXED

    pdf_bkg->Print();
    RooArgSet *params = pdf_bkg->getParameters(*mzg);
    RooFIter iter = params->fwdIterator();
    RooAbsArg* param;
    while ((param = iter.next())) {
        RooRealVar* realParam = dynamic_cast<RooRealVar*>(param);
        if (realParam) {
            // Modify parameter values as needed
            
            if (realParam->GetName() == TString("env_pdf_untag_cat3_bern4_b1")) {
                realParam->setVal(0.30002); // Modify mean parameter
            }
            if (realParam->GetName() == TString("env_pdf_untag_cat3_bern4_b2")) {
                realParam->setVal(1.1396); // Modify mean parameter
            }
            if (realParam->GetName() == TString("env_pdf_untag_cat3_bern4_b3")) {
                realParam->setVal(0.19042); // Modify mean parameter
            }
            if (realParam->GetName() == TString("env_pdf_untag_cat3_bern4_b4")) {
                realParam->setVal(0.42523); // Modify mean parameter
            }
            if (realParam->GetName() == TString("env_pdf_untag_cat3_bern4_sigma_b4")) {
                realParam->setVal(5.0937); // Modify mean parameter
            }
            if (realParam->GetName() == TString("step_b4")) {
                realParam->setVal(103.57); // Modify mean parameter
            }
            cout<<realParam->GetName()<<":"<<realParam->getVal()<<endl;
            
        }
    }

    RooChi2Var *stat = new RooChi2Var("a", "a", *pdf_bkg, *hist_bkg, RooFit::DataError(RooAbsData::Poisson));
    cout<<stat->getVal()<<":"<<params<<endl;
    //file.multipdf.Print();
}