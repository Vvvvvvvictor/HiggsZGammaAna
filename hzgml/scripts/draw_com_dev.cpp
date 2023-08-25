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

void draw_com_dev(TString mode = "inclusive")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    // TString s_LG[] = {"SM ZG", "DYJets", "TT", "VBS ZG"};//, "Diboson"};
    // TString l_LG[] = {"data"};
    double Y_axis[3] = {0, 0.02, 100};
    int stypes = 4, ltypes = 1;

    TString x_name = "n_jets";
    Long64_t data;
    double x_axis[3] = {6, 0, 6};
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    TString stype[] = {"zg", "dy0", "tt", "llajj"};//, {"ww", "wz", "zz", "0"}};
    TString ltype[] = {"data"};

    cout << "\n\n\tFinish basic setting" << endl;

    setStyle();
    gStyle->SetTitleXSize(0.1);
    gStyle->SetTitleYSize(0.1);
    TCanvas *can = new TCanvas("can", "can", 800, 300);
    can->SetFillColorAlpha(0, 0);
    can->SetFrameFillColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetOptStat(0);
    can->SetLeftMargin(0.1);
    can->SetTopMargin(0.05);
    can->SetRightMargin(0.05);
    can->SetBottomMargin(0.22);
    can->SetLogy(0);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Double_t w, cut, n;
    double max = 0.0, num[ltypes + stypes];

    // THStack *hs = new THStack("hs", "test stacked histograms");
    TH1D *h1 = new TH1D("h1", Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);
    h1->SetLineColor(1);
    h1->SetFillColor(1);
    // h1[typei]->SetFillStyle(3004);
    h1->SetLineWidth(2);

    for (int typei = 0; typei < stypes; typei++)
    {
        cout << "\n\tLoading file" << Form("%s/overlap/%s.root", str_dir.Data(), stype[typei].Data()) << "!......" << endl;

        TFile *f = new TFile(Form("%s/overlap/%s.root", str_dir.Data(), stype[typei].Data()));

        TTree *t = (TTree *)f->Get(mode.Data());
        Long64_t nentries = t->GetEntries();
        t->SetBranchAddress(x_name.Data(), &data);
        t->SetBranchAddress("weight_central", &w);
        if (typei == 1)
            t->SetBranchAddress("n_iso_photons", &n);
        t->SetBranchAddress("H_mass", &cut);

        cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
            << "\n\tloading...loading...zzzzzz\n"
            << endl;

        num[typei] = 0.0;
        for (Long64_t i = 0; i < nentries; i++)
        {
            t->GetEntry(i);
            if (typei == 1 && n != 0) continue;
            if (cut < 120 || cut > 130)
            {
                h1->Fill(data, w);
                num[typei] += w;
            }
        }
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei] << " items are in it!!";
    }

    TH1D *h2 = new TH1D("h2", Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);

    for (int typei = 0; typei < ltypes; typei++)
    {
        cout << "\n\tLoading file" << Form("%s/overlap/%s.root", str_dir.Data(), ltype[typei].Data()) << "!......" << endl;

        TFile *f = new TFile(Form("%s/overlap/%s.root", str_dir.Data(), ltype[typei].Data()));

        TTree *t = (TTree *)f->Get(mode.Data());
        Long64_t nentries = t->GetEntries();
        t->SetBranchAddress(x_name.Data(), &data);
        t->SetBranchAddress("weight_central", &w);
        t->SetBranchAddress("H_mass", &cut);

        cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                << "\n\tloading...loading...zzzzzz\n"
                << endl;

        num[typei + stypes] = 0.0;
        for (Long64_t i = 0; i < nentries; i++)
        {
            t->GetEntry(i);
            if (cut < 120 || cut > 130)
            {
                h2->Fill(data, w);
                num[stypes+typei] += w;
            }
        }
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei + stypes] << " items are in it!!";
    }

    h1->Divide(h2);
    h1->Draw("E1");
    h1->GetXaxis()->SetTitle(Form("with DY(madgraph): %s", x_name.Data()));
    h1->GetYaxis()->SetTitle("MC/data");
    h1->GetYaxis()->SetTitleOffset(0.5);
    can->Update();

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry("h1", "MC/data", "lep");
    legend->SetBorderSize(0);
    // legend->SetFillColor(0);
    legend->SetFillColorAlpha(0, 0);
    legend->Draw();
    cout << "\n\tFinish Create TLegend" << endl;

    // can->SetGrid(0, 1);
    can->SaveAs(Form("%s/overlap/figures/com/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data()));

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/com/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data())) << endl;
}
