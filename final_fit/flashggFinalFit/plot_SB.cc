#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooRealVar.h>
#include <RooBinning.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooSimultaneous.h>

void plot_SB() {
    // TFile *f = new TFile("higgsCombine.cat1_bestfit.MultiDimFit.mH125.root");
    // RooWorkspace *w = dynamic_cast<RooWorkspace*>(f->Get("w"));
    TFile *f = new TFile("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/Combine_results_cf/CMS-HGG_mva_13TeV_multipdf_cat0.root");
    RooWorkspace *w = dynamic_cast<RooWorkspace*>(f->Get("multipdf"));
    if (!w) {
        delete f;
        return;
    }

    int n_bins = 65;
    RooBinning binning1(17, 105, 122);
    RooBinning binning2(42, 128, 170);

    TCanvas *can = new TCanvas();
    RooRealVar *MH = w->var("CMS_hzg_mass");
    RooPlot *plot = MH->frame();
    // RooAbsData *data_obs = w->data("data_obs");
    // data_obs->plotOn(plot, RooFit::Binning(binning1));
    // data_obs->plotOn(plot, RooFit::Binning(binning2));

    // Load the S+B model
    // RooSimultaneous *sb_model = dynamic_cast<RooSimultaneous*>(w->pdf("model_s"));
    //RooSimultaneous *b_model = dynamic_cast<RooSimultaneous*>(w->pdf("model_b"));
    // if (!sb_model) {
    //     delete f;
    //     delete can;
    //     delete plot;
    //     return;
    // }
    RooAbsPdf *pdf_cat1 = w->pdf("CMS_hgg_cat0_13TeV_bkgshape");
    // RooAbsPdf *pdf_cat1 = sb_model->getPdf("cat1");
    //RooAbsPdf *bkgpdf_cat1 = b_model->getPdf("cat1");
    if (!pdf_cat1) {
        delete f;
        delete can;
        delete plot;
        return;
    }
    //RooAbsReal *integral = pdf_cat1->createIntegral(*MH, RooFit::NormSet(*MH));
    //std::cout << integral->getVal() << std::endl;
    // w->var("r")->setVal(1);
    Double_t N_SB = pdf_cat1->expectedEvents(*MH);
    //Double_t N_B = bkgpdf_cat1->expectedEvents(*MH);
    std::cout << "N_SB = " << N_SB << std::endl;
    //std::cout << "N_B = " << N_B << std::endl;
    //std::cout << "N_S = " << N_SB - N_B << std::endl;

    // Prefit
    pdf_cat1->plotOn(plot, RooFit::Components(*pdf_cat1), RooFit::LineColor(2), RooFit::Name("prefit"));

    // Postfit
    w->loadSnapshot("multipdf");
    pdf_cat1->plotOn(plot, RooFit::Components(*pdf_cat1), RooFit::LineColor(4), RooFit::Name("postfit"));
    // double r_bestfit = w->var("r")->getVal();
    //w->var("r")->setVal(0);
    Double_t N_B = pdf_cat1->expectedEvents(*MH);
    std::cout << "N_B = " << N_B << std::endl;
    std::cout << "N_S = " << N_SB - N_B << std::endl;

    // plot->Draw();

    // TLegend *leg = new TLegend(0.55, 0.6, 0.85, 0.85);
    // leg->AddEntry(plot->findObject("prefit"), "Prefit S+B model (r=1.00)", "L");
    // leg->AddEntry(plot->findObject("postfit"), Form("Postfit S+B model (r=%.2f)", r_bestfit), "L");
    // leg->Draw("Same");

    // can->Update();
    // can->SaveAs("cat1_sb_model.png");

    //delete f;
    //delete can;
    //delete plot;
    //delete leg;
}
