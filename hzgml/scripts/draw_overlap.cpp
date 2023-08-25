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

void draw_overlap(TString mode = "inclusive", int id = 11)
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString LG[]={"Run3(2022E)", "Run2(2017, Scaled)"};
    TString F[]={"datae", "data"}, ID;
    double Y_axis[3] = {0, 0.02, 100};
    double pos[2] = {0.6, 1.00};
    double ratio = 1.3; // the ratio of top limit to the highest bin
    int ltypes = 2;

    //TODO: check things below!
    Double_t data;
    if (id==11) ID = "ele";
    if (id==13) ID = "mu";
    TString x_name = "H_phi";
    TString x_title = "#phi_{ll}";
    double X_axis[3] = {32, -3.2, 3.2};
    TString y_title = Form("Events/%.1f" , (X_axis[2]-X_axis[1])/X_axis[0]);
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    
    cout << "\n\n\tFinish basic setting" << endl;

    setStyle();
    gStyle->SetTitleXSize(0.07);
    gStyle->SetTitleYSize(0.06);
    TCanvas *can = new TCanvas("can", "can", 600, 600);
    can->SetFillColorAlpha(0, 0);
    can->SetFrameFillColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetOptStat(0);
    can->SetLeftMargin(0.05);
    can->SetTopMargin(0);
    can->SetRightMargin(0.05);
    can->SetBottomMargin(0);
    can->SetLogy(0);
    TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.30);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Double_t w, n, pt, wp, hmass, zmass;
    Long64_t flav;
    Double_t pass1, pass2;
    double max = 0.0, num[ltypes], sf;
    int sum_num[ltypes];

    TH1D *h2[ltypes];
    can->cd(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    pad1->SetFillColorAlpha(0, 0);
    pad1->SetFillStyle(0);

    for (int typei = 0; typei < ltypes; typei++)
    {
        h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; %s: %s; events", mode.Data(), ID.Data(), x_title.Data()), X_axis[0], X_axis[1], X_axis[2]);
        h2[typei]->SetLineColor(1);
        if (typei==1) {
            h2[typei]->SetLineColor(TColor::GetColorDark(typei + 3));
            h2[typei]->SetFillColorAlpha(TColor::GetColorDark(typei + 3), 0.5);
            h2[typei]->SetFillStyle(1001);
        }
        h2[typei]->SetLineStyle(1);
        h2[typei]->SetLineWidth(3);
        cout << "\n\tLoading file" << Form("%s/overlap/%s.root", str_dir.Data(), F[typei].Data()) << "!......" << endl;

        TFile *f = new TFile(Form("%s/overlap/%s.root", str_dir.Data(), F[typei].Data()));

        TTree *t = (TTree *)f->Get(mode.Data());
        Long64_t nentries = t->GetEntries();
        t->SetBranchAddress(x_name, &data);
        t->SetBranchAddress("weight_central", &w);
        t->SetBranchAddress("gamma_mvaID_WP80", &wp);
        t->SetBranchAddress("Z_lead_lepton_id", &flav);
        t->SetBranchAddress("H_mass", &hmass);
        t->SetBranchAddress("Z_mass", &zmass);
        // t->SetBranchAddress("n_iso_photons", &n);
        // t->SetBranchAddress("gamma_genPartFlav", &pass1);
        // t->SetBranchAddress("Z_sublead_lepton_pt", &pass2);
        //  t->SetBranchAddress("gamma_pt", &pt);

        cout << "\n\tStart drawing " << x_name.Data() << "\n\tentries are : " << nentries << " ! oh my god ! too much !"
             << "\n\tloading...loading...zzzzzz\n"
             << endl;

        num[typei] = 0.0;
        sum_num[typei] = 0;
        double sum = 0.0;
        for (Long64_t i = 0; i < nentries; i++)
        {
            t->GetEntry(i);
            sum+=w;
            // if (pass2 < 15) continue;
            // if (n >= 1 & typei == 1)
            //     {
            //         h2[typei]->Fill(data, w);
            //         num[typei] += w;
            //         sum_num[typei]++;
            //     }
            // if (n >= 1 & typei == 0)
            // {
            //     h2[typei]->Fill(data, w);
            //     num[typei] += w;
            //     sum_num[typei]++;
            // }
            if (flav!=id) continue;
            if (wp!=1) continue;
            if ((zmass<80) | (zmass>100)) continue;
            if ((hmass>120) & (hmass<130)) continue;
            if (x_name.Contains("H_mass", TString::kIgnoreCase)) data = hmass;
            if (x_name.Contains("Z_mass", TString::kIgnoreCase)) data = zmass;
            // if (typei==1) w = w;
            h2[typei]->Fill(data, w);
            num[typei] += w;
            sum_num[typei]++;
        }
        // h2[typei]->Scale(1. / h2[typei]->Integral(), "width");
        if (typei==1) {
            sf = num[0]/num[typei];
            h2[typei]->Scale(sf);
            num[typei] = num[0];
        }
        if (max < h2[typei]->GetBinContent(h2[typei]->GetMaximumBin()))
            max = h2[typei]->GetBinContent(h2[typei]->GetMaximumBin());
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << sum_num[typei] << "/" << nentries  << " items are in it!!";
        cout << "\n\t" << num[typei] << "/" << sum  << " items are in it!!";
    }
    for (int typei = 0; typei < ltypes; typei++)
    {
        if (typei==0){
            h2[typei]->Draw("E1, same");
            continue;
        }
        h2[typei]->Draw("HIST, same");
    }
    h2[0]->SetMaximum(max * ratio);
    // h2[0]->GetYaxis()->SetLabelFont(22);
    // h2[0]->GetYaxis()->SetLabelSize(20);
    // h2[0]->GetYaxis()->SetTitleFont(22);
    // h2[0]->GetYaxis()->SetLabelOffset(0.1);
    // h2[0]->GetYaxis()->SetTitleSize(24);
    // h2[0]->GetXaxis()->SetLabelFont(22);
    // h2[0]->GetXaxis()->SetLabelSize(20);
    h2[0]->GetXaxis()->SetTitleSize(0);
    h2[0]->GetYaxis()->SetTitle(y_title);
    // h2[0]->GetYaxis()->SetTitleOffset(1.25);
    // h2[0]->GetXaxis()->SetTitleOffset(1);

    pad1->SetLeftMargin(0.15);
    pad1->SetTopMargin(0.05);
    // pad1->SetRightMargin(0.15);
    pad1->SetBottomMargin(0.02);

    // for (int typei = 0; typei < stypes; typei++)
    // {
    //     h1[typei]->Reset("ICES");
    // }
    can->Update();

    TLegend *legend = new TLegend(pos[0], 0.825, pos[1], 0.95);
    legend->SetTextFont(22);
    for (int typei = 0; typei < ltypes; typei++)
        legend->AddEntry(Form("h2_%d", typei), Form("%s(%.0f)", LG[typei].Data(), num[typei]), "l");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0);
    legend->Draw();
    can->Update();
    cout << "\n\tFinish Create TLegend" << endl;

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(22);
    text->DrawLatex(0.65, 0.775, Form("L=5.67fb^{-1}, SF: %.2f%%", sf*100));

    // TLatex *text = new TLatex();
    // text->SetNDC();
    // text->SetTextSize(0.05);
    // // if (pt_bin[1]!=5000)
    // //     text->DrawLatex(0.175, 0.875, Form("%.0fGeV<P_{T}^{#gamma}<%.0fGeV", pt_bin[0], pt_bin[1]));
    // // else
    // text->DrawLatex(0.4, 0.9, Form("L=", pt_bin[0]));

    // TLine *l1 = new TLine(20, 0, 20, 2350);
    // l1->SetLineColor(kRed);
    // l1->SetLineWidth(3);
    // l1->Draw();

    can->Update();

    can->cd(0);
    pad2->cd();
    pad2->SetGrid(0, 3);
    pad2->SetLeftMargin(0.15);
    pad2->SetBottomMargin(0.35);
    pad2->SetFillColorAlpha(0, 0);
    pad2->SetFillStyle(0);
    setStyle();
    gStyle->SetTitleXSize(0.28);
    gStyle->SetTitleYSize(0.32);
    TH1D *div_h, *numerator_h;
    div_h = (TH1D*) h2[0]->Clone();
    numerator_h = (TH1D*) h2[1]->Clone();
    div_h->SetLineColor(1);
    div_h->SetLineWidth(2);
    div_h->Divide(numerator_h);

    double min;
    int max_bin = div_h->GetMaximumBin();
    max = div_h->GetBinContent(max_bin);
    double max_err = div_h->GetBinError(max_bin);
    min = div_h->GetMinimum(0);
    int min_bin = 0;
    div_h->GetBinWithContent(min, min_bin, 0);
    double min_err = div_h->GetBinError(min_bin);
    div_h->SetMaximum(max + max_err);
    div_h->SetMinimum(min - min_err*2);

    // div_h->SetMaximum(3);
    // div_h->SetMinimum(0);

    div_h->Draw("E1");
    // div_h->GetXaxis()->SetLabelFont(22);
    div_h->GetXaxis()->SetLabelSize(0.125);
    // div_h->GetXaxis()->SetLabelOffset(0.02);
    // div_h->GetXaxis()->SetTitleFont(22);
    div_h->GetXaxis()->SetTitleSize(0.175);
    // div_h->GetXaxis()->SetTitleOffset(0.5);
    // div_h->GetYaxis()->SetLabelFont(22);
    div_h->GetYaxis()->SetLabelSize(0.100);
    // div_h->GetYaxis()->SetLabelOffset(0.02);
    // div_h->GetYaxis()->SetTitleFont(22);
    div_h->GetYaxis()->SetTitleSize(0.150);
    div_h->GetYaxis()->SetTitle("#frac{2022E}{2017}");
    div_h->GetYaxis()->SetTitleOffset(0.4); 
    // div_h->GetXaxis()->SetTickLength(0.1);
    div_h->GetXaxis()->SetTitle(x_title);
    // div_h->GetYaxis()->SetNdivisions(205);

    pad2->SetLeftMargin(0.15);
    pad2->SetTopMargin(0.00);
    // pad2->SetRightMargin(0.15);
    pad2->SetBottomMargin(0.40);

    can->Update();
    cout << "\tFinish drawing subpad!!!" << endl;

    // can->SetGrid(0, 2);
    can->SaveAs(Form("%s/overlap/figures/over/%s_%s_%s.png", str_dir.Data(), x_name.Data(), mode.Data(), ID.Data()));

    // for (int typei = 0; typei < ltypes; typei++)
    //     h2[typei]->Reset("ICES");

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/over/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data())) << endl;
}
