#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <algorithm>
#include <THStack.h>
#include <TStyle.h>
#include <TColor.h>
#include "/afs/cern.ch/user/j/jiehan/private/boot.h"

using namespace std;

void draw_2D(TString mode = "test", TString dataset = "bkgmc")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    double X_axis[3] = {100, 0, 1};
    double Y_axis[3] = {100, 0, 1};
    // double X_axis[3] = {120, 0.000000000001, 0.03};

    TString y_name = "bdt_score_t";
    Float_t Y;
    TString x_name = "bdt_score_VBF";
    Float_t X;
    TString file_dir = "/eos/home-j/jiehan/root/outputs/two_jet";
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml";

    int isbarrel = 1;

    cout << "\n\n\tFinish basic setting" << endl;

    // TCanvas *can[3];
    setStyle();
    gStyle->SetTitleXSize(0.07);
    gStyle->SetTitleYSize(0.06);
    TCanvas *can = new TCanvas("can", "can", 800, 600);
    can->SetFillColorAlpha(0, 0);
    can->SetFrameFillColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetOptStat(0);
    can->SetLeftMargin(0.15);
    can->SetTopMargin(0.05);
    can->SetRightMargin(0.15);
    can->SetBottomMargin(0.17);
    // can->SetLogy(1);
    can->SetLogz(1);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Double_t w, pt, eta, bit, match, isb, isp, sieie, chiso, subl, mva;
    Long64_t n;
    double max = 0.0;

    TH2D *h2;

    // h2 = new TH2D("h2", Form("%s; %s:%s; %s", mode.Data(), mode.Data(), x_name.Data(), y_name.Data()), X_axis[0], X_axis[1], X_axis[2], Y_axis[0], Y_axis[1], Y_axis[2]);
    h2 = new TH2D("h2", Form("%s; VBF BDT score; Dijet BDT score", mode.Data()), X_axis[0], X_axis[1], X_axis[2], Y_axis[0], Y_axis[1], Y_axis[2]);
    // h2->SetLineStyle(1);
    // h2->SetLineWidth(3);
    cout << "\n\tLoading file" << Form("%s/%s.root", file_dir.Data(), dataset.Data()) << "!......" << endl;

    TFile *f = new TFile(Form("%s/%s.root", file_dir.Data(), dataset.Data()));

    TTree *t = (TTree *)f->Get(mode.Data());
    Long64_t nentries = t->GetEntries();
    t->SetBranchAddress(x_name, &X);
    t->SetBranchAddress(y_name, &Y);
    // t->SetBranchAddress("gamma_pt", &pt);
    // t->SetBranchAddress("gamma_eta", &eta);
    // t->SetBranchAddress("gamma_photon_match", &match);
    // t->SetBranchAddress("photon_selection", &bit);
    // t->SetBranchAddress("photon_is_barrel", &isb);
    // t->SetBranchAddress("gamma_genPartFlav", &isp);
    // t->SetBranchAddress("gamma_chiso", &chiso);
    t->SetBranchAddress("weight", &w);
    // t->SetBranchAddress("gamma_mvaID", &mva);
    // t->SetBranchAddress("n_jets", &n);
    // t->SetBranchAddress("Z_sublead_lepton_pt", &subl);

    cout << "\n\tStart drawing " << x_name.Data() << "\n\tentries are : " << nentries << " ! oh my god ! too much !"
            << "\n\tloading...loading...zzzzzz\n"
            << endl;

    double sum = 0.0;
    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);
        // if ((pt<20) | (pt>30)) continue;
        // if (isp!=0) continue;
        // if (chiso*pt>2) continue;
        // if (abs(eta)<1.5) if ((mva<0.42)|(mva>1)) continue;
        // if (abs(eta)>1.5) if ((mva<-0.14)|(mva>1)) continue;
        // if (abs(eta)>1.5) continue;
        // if (abs(eta)<1.5) if (chiso*pt<1.141) continue;
        // if (abs(eta)>1.5) if (chiso*pt<1.051) continue;
        // if (subl<15) continue;
        // if (n!=0) continue;
        // if ((match==1) & ((int(bit/(1<<5))%2==1) | (int(bit/(1<<5))%2==1)) & (isb==isbarrel))
        // if (isb==isbarrel)
        {
            h2->Fill(X, Y, w);
            sum+=w;
            // h2->Fill(X, Y);
            // h2->Fill(X, Y*pt);
            // h2->Fill(X*pt, Y);
            // h2->Fill(X*pt, Y/pt);
        }
    }

    // f = new TFile(Form("%s/overlap/zg.root", str_dir.Data()));

    // t = (TTree *)f->Get(mode.Data());
    // nentries = t->GetEntries();
    // t->SetBranchAddress(x_name, &X);
    // t->SetBranchAddress(y_name, &Y);
    // // t->SetBranchAddress("gamma_pt", &pt);
    // // t->SetBranchAddress("gamma_eta", &eta);
    // // t->SetBranchAddress("gamma_photon_match", &match);
    // // t->SetBranchAddress("photon_selection", &bit);
    // // t->SetBranchAddress("photon_is_barrel", &isb);
    // // t->SetBranchAddress("gamma_genPartFlav", &isp);
    // // t->SetBranchAddress("gamma_chiso", &chiso);
    // t->SetBranchAddress("weight", &w);
    // // t->SetBranchAddress("gamma_mvaID", &mva);
    // // t->SetBranchAddress("n_iso_photons", &n);
    // // t->SetBranchAddress("Z_sublead_lepton_pt", &subl);

    // cout << "\n\tStart drawing " << x_name.Data() << "\n\tentries are : " << nentries << " ! oh my god ! too much !"
    //         << "\n\tloading...loading...zzzzzz\n"
    //         << endl;

    // // double sum = 0.0;
    // for (Long64_t i = 0; i < nentries; i++)
    // {
    //     t->GetEntry(i);
    //     // if ((pt<30) | (pt>30)) continue;
    //     // if (chiso*pt>2) continue;
    //     // if (abs(eta)<1.5) if ((mva<-0.97)|(mva>-0.4)) continue;
    //     // if (abs(eta)>1.5) if ((mva<-0.97)|(mva>-0.59)) continue;
    //     // if (subl<15) continue;
    //     // if (isp!=1) continue;
    //     // if ((match==1) & ((int(bit/(1<<5))%2==1) | (int(bit/(1<<5))%2==1)))
    //     // if (isb==isbarrel)
    //     {
    //         h2->Fill(X, Y, w);
    //         // h2->Fill(X, Y);
    //         // h2->Fill(X, Y*pt);
    //         // h2->Fill(X, Y/pt);
    //         // h2->Fill(X*pt, Y/pt);
    //     }
    // }

    cout << "\tFinished loading File: " << sum << " events in plot." << endl;
    // h2->SetAxisRange(X_axis[1], X_axis[2], "X");
    // h2->SetAxisRange(Y_axis[1], Y_axis[2], "Y");
    h2->Draw("COLZ");

    double vb = 0.41;

    TLine *l0 = new TLine(vb,Y_axis[1],vb,Y_axis[2]);
    // TLine *l0 = new TLine(-0.97,5,X_axis[2],5);
    l0->SetLineColor(kBlue);
    l0->SetLineWidth(3);
    l0->Draw();

    // TLine *l10 = new TLine(X_axis[1],0.02,vb,0.02);
    // // TLine *l1 = new TLine(-0.97,5,X_axis[2],5);
    // l10->SetLineColor(kBlue);
    // l10->SetLineWidth(3);
    // l10->Draw();

    TLine *l1 = new TLine(X_axis[1],0.27,vb,0.27);
    // TLine *l1 = new TLine(-0.97,5,X_axis[2],5);
    l1->SetLineColor(kBlue);
    l1->SetLineWidth(3);
    l1->Draw();

    TLine *l2 = new TLine(X_axis[1],0.52,vb,0.52);
    // TLine *l2 = new TLine(X_axis[1],1.5,-0.97,1.5);
    l2->SetLineColor(kBlue);
    l2->SetLineWidth(3);
    l2->Draw();

    TLine *l3 = new TLine(X_axis[1],0.64,vb,0.64);
    // TLine *l3 = new TLine(-0.97, 0, -0.97, 5);
    l3->SetLineColor(kBlue);
    l3->SetLineWidth(3);
    l3->Draw();

    TLine *l4 = new TLine(vb,0.39,X_axis[2],0.39);
    // TLine *l4 = new TLine(-0.97,5,X_axis[2],5);
    l4->SetLineColor(kBlue);
    l4->SetLineWidth(3);
    l4->Draw();

    TLine *l5 = new TLine(vb,0.62,X_axis[2],0.62);
    // TLine *l5 = new TLine(X_axis[1],1.5,-0.97,1.5);
    l5->SetLineColor(kBlue);
    l5->SetLineWidth(3);
    l5->Draw();

    TLine *l6 = new TLine(vb,0.76,X_axis[2],0.76);
    // TLine *l6 = new TLine(-0.97, 0, -0.97, 5);
    l6->SetLineColor(kBlue);
    l6->SetLineWidth(3);
    l6->Draw();

    // TLatex *text = new TLatex();
    // text->SetNDC();
    // text->SetTextColor(kRed);
    // text->SetTextSize(0.08);
    // text->SetTextAlign(22);
    // text->DrawLatex(0.5, 0.4, "Reserved");
    // text->DrawLatex(0.3, 0.47, "A");
    // text->DrawLatex(0.6, 0.47, "B");
    // text->DrawLatex(0.3, 0.85, "C");
    // text->DrawLatex(0.6, 0.85, "D");

    can->Update();

    // can->SetGrid(0, 2);
    can->SaveAs(Form("%s/figures/2d_%s_%s_%s.pdf", str_dir.Data(), x_name.Data(), y_name.Data(), dataset.Data()));

    // for (int typei = 0; typei < ltypes; typei++)
    //     h2[typei]->Reset("ICES");

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/figures/2d_%s_%s_%s.pdf", str_dir.Data(), x_name.Data(), y_name.Data(), dataset.Data())) << endl;
}
