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

void draw_bdtDis(TString mode = "zero_jet", int region = 0, int button = 0)
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString LG[][6] = {{"Fake photon", "SM ZG", "EWK ZG", "TTGJets", "TGJets", "Data"}, {"DY+fake", "SM ZG", "EWK ZG", "TT", "Diboson", "Data"}};
    double Y_axis[3] = {1, 0.2, 100};

    TString x_name = "bdt_score_t";
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    TString stype[][6][4] = {{{"data_fake", "mc_true", "0"}, {"mc_med", "0"}, {"mc_med", "0"}, {"mc_med", "0"}, {"mc_med", "0"}}, {{"DYJetsToLL", "0"}, {"ZGToLLG", "0"}, {"LLAJJ", "0"}, {"TT", "0"}, {"WW", "WZ", "ZZ", "0"}}};//, {"mc_med", "0"}};//, "sig"};
    TString ltype[][4] = {{"data", "0"}};//{"sig", "0"}, {"ggH", "0"}, {"VBF", "0"}, {"WminusH", "WplusH", "ZH", "0"}, {"ttH", "0"}};

    int bins = 10;
    int stypes = 2, ltypes = 1;
    if (button)
        stypes = 5;
    else
        stypes = 5;
    // string tag = "zero_jet";
    // if (mode.Data() == tag)
    //     ltypes = 4;
    cout << "\n\tThere are " << ltypes << " kinds of Higgs generation mode!!!" << endl;
    Double_t num[stypes + ltypes];
    // double res, reserr;
    // ifstream file(Form("%s/outputs/bin_binaries_%s.txt", str_dir.Data(), mode.Data()));
    // cout << "\n\tThe binaries of bin are: ";
    // file >> bins >> setw(1);
    // cout << bins << ", ";
    // Float_t edges[bins + 1];
    // for (int n = 0; n < bins + 1; n++)
    // {
    //     file >> edges[n] >> setw(1);
    //     cout << edges[n] << " ";
    // }
    // file >> res >> setw(1) >> reserr;

    cout << "\n\n\tFinish basic setting" << endl;

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
    can->SetRightMargin(0.05);
    can->SetBottomMargin(0.13);
    can->SetLogy(1);
    TPad *pad1 = new TPad("pad1", "", 0, 0.23, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Float_t data;
    Double_t w, hmass, wp;
    Long64_t tag;
    double max=0.0, min=0.0;

    THStack *hs = new THStack("hs", "test stacked histograms");
    // TH1D *mc_h = new TH1D("mc", Form("%s; %s; Events", mode.Data(), x_name.Data()), bins, edges);
    // TH1D *data_h =  new TH1D("data", Form("%s; %s; Events", mode.Data(), x_name.Data()), bins, edges);
    TH1D *mc_h = new TH1D("mc", Form("%s; %s; Events", mode.Data(), x_name.Data()), bins, 0, 1);
    TH1D *data_h =  new TH1D("data", Form("%s; %s; Events", mode.Data(), x_name.Data()), bins, 0, 1);
    TH1D *h1[stypes];

    can->cd(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    double sum=0.;
    for (int typei = 0; typei < stypes; typei++)
    {
        // h1[typei] = new TH1D(Form("h1_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), bins, edges);
        h1[typei] = new TH1D(Form("h1_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), bins, 0, 1);
        h1[typei]->SetLineColor(TColor::GetColorDark(typei + 2));
        h1[typei]->SetFillColor(TColor::GetColorDark(typei + 2));
        h1[typei]->SetLineWidth(1);

        int len = 0;
        num[typei] = 0.0;
        while (stype[button][typei][len] != "0")
        {
        cout << "\n\tLoading file" << Form("%s/outputs/%s/%s.root", str_dir.Data(), mode.Data(), stype[button][typei][len].Data()) << "!......" << endl;
        TFile *f = new TFile(Form("%s/outputs/%s/%s.root", str_dir.Data(), mode.Data(), stype[button][typei][len].Data()));

        TTree *t = (TTree *)f->Get("test");
        Long64_t nentries = t->GetEntries();
        t->SetBranchAddress(x_name.Data(), &data);
        t->SetBranchAddress("weight", &w);
        t->SetBranchAddress("H_mass", &hmass);
        if (!button)
            {if (typei) t->SetBranchAddress("tagger", &tag);}

        cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
             << "\n\tloading...loading...zzzzzz\n"
             << endl;

        for (Long64_t i = 0; i < nentries; i++)
        {
            t->GetEntry(i);
            if (button)
                {if (wp==0) continue;}
            else
                {if (typei>0) if (tag!=typei) continue;}

            if (region)
            {
                if (hmass > 120 && hmass < 130)
                {
                    h1[typei]->Fill(data, w);
                    mc_h->Fill(data, w);
                    num[typei] += w;
                }
            }
            else
            {
                if (hmass < 120 || hmass > 130)
                {
                    h1[typei]->Fill(data, w);
                    mc_h->Fill(data, w);
                    num[typei] += w;
                    sum += w;
                }
            }
        }
        f->Close();
        len++;
        }
        hs->Add(h1[typei]);
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei] << " items are in it!!";
    }
    max = mc_h->GetBinContent(mc_h->GetMaximumBin());
    min = mc_h->GetBinContent(mc_h->GetMinimumBin());
    hs->Draw("HIST");
    // hs->GetXaxis()->SetTitle(Form("%s:%s", mode.Data(), x_name.Data()));
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.05);
    pad1->SetLogy(1);
    hs->GetYaxis()->SetLabelFont(43);
    hs->GetYaxis()->SetLabelSize(20);
    hs->GetYaxis()->SetTitleFont(43);
    hs->GetYaxis()->SetTitleSize(24);
    hs->GetXaxis()->SetLabelFont(43);
    hs->GetXaxis()->SetLabelSize(20);
    hs->GetXaxis()->SetTitleFont(43);
    hs->GetXaxis()->SetTitleSize(24);
    hs->GetYaxis()->SetTitle(Form("Events/%.2f", 1./bins));
    hs->GetYaxis()->SetTitleOffset(1.2);
    if (Y_axis[0] > 0)
    {
        hs->SetMinimum(TMath::Power(10, TMath::Log10(min)-0.5));
        hs->SetMaximum(TMath::Power(10, TMath::Log10(max)+0.2));
    }

    TH1D *h2[ltypes];

    for (int typei = 0; typei < ltypes; typei++)
    {
        // h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), bins, edges);
        h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), bins, 0, 1);
        h2[typei]->SetLineColor(typei + 1);
        h2[typei]->SetLineStyle(1);
        h2[typei]->SetLineWidth(2);
        int len = 0;
        num[typei + stypes] = 0.0;
        while (ltype[typei][len] != "0")
        {
            cout << "\n\tLoading file" << Form("%s/outputs/%s/%s.root", str_dir.Data(), mode.Data(), ltype[typei][len].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/outputs/%s/%s.root", str_dir.Data(), mode.Data(), ltype[typei][len].Data()));

            TTree *t = (TTree *)f->Get("test");
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name.Data(), &data);
            t->SetBranchAddress("weight", &w);
            t->SetBranchAddress("H_mass", &hmass);

            cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                 << "\n\tloading...loading...zzzzzz\n"
                 << endl;

            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                if (region)
                {
                    if (hmass > 120 && hmass < 130)
                    {
                        h2[typei]->Fill(data, w);
                        data_h->Fill(data, w);
                        num[typei + stypes] += w;
                    }
                }
                else
                {
                    if (hmass < 120 || hmass > 130)
                    {
                        h2[typei]->Fill(data, w);
                        data_h->Fill(data, w);
                        num[typei + stypes] += w;
                    }
                }
            }
            f->Close();
            len++;
        }
        h2[typei]->Draw("E1, same");
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei + stypes] << " items are in it!!";
    }

    double factor = num[stypes] / sum;
    for (int typei=0; typei<stypes; typei++){
        h1[typei]->Scale(factor);
        num[typei] = num[typei]*factor;
    }
    if (max == mc_h->GetBinContent(mc_h->GetMaximumBin()))
        max = max * factor;
    mc_h->Scale(factor);
    // hs->SetMaximum(max*1.25);

    data_h->Divide(mc_h);
    data_h->Sumw2();
    double kai2 = 0.;
    int good_bin = 0;
    for (int bini=0; bini<bins; bini++){
        double value = data_h->GetBinContent(bini);
        double error = data_h->GetBinError(bini);
        if ((value>0) & (error>0)){
            kai2 += (value-1)/error*(value-1)/error;
            good_bin += 1;
        }
    }

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextColor(1);
    text->SetTextSize(0.04);
    text->SetTextAlign(22);
    text->DrawLatex(0.3, 0.9, Form("Scale factor: %.4f", factor));
    text->DrawLatex(0.5, 0.9, Form("#chi^{2}/n: %.4f", kai2/good_bin));
    // text->DrawLatex(0.425, 0.85, Form("Significance: %.4f #pm %.4f", res, reserr));

    TLegend *legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    for (int typei = 0; typei < stypes; typei++)
        legend->AddEntry(Form("h1_%d", typei), Form("%s(%.2f)", LG[button][typei].Data(), num[typei]), "f");
    for (int typei = 0; typei < ltypes; typei++)
        legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", LG[button][typei + stypes].Data(), num[typei + stypes]), "lep");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0);
    legend->Draw();
    cout << "\n\tFinish Create TLegend" << endl;

    can->cd(0);
    pad2->cd();
    pad2->SetGrid(0, 3);
    pad2->SetLeftMargin(0.12);
    pad2->SetBottomMargin(0.34);
    setStyle();
    gStyle->SetTitleXSize(0.28);
    gStyle->SetTitleYSize(0.32);
    data_h->SetLineColor(1);
    data_h->SetLineWidth(1);
    data_h->Draw("E1");
    data_h->GetXaxis()->SetLabelFont(43);
    data_h->GetXaxis()->SetLabelSize(20);
    data_h->GetXaxis()->SetLabelOffset(0.02);
    data_h->GetXaxis()->SetTitleFont(43);
    data_h->GetXaxis()->SetTitleSize(30);
    data_h->GetXaxis()->SetTitleOffset(3);
    data_h->GetYaxis()->SetLabelFont(43);
    data_h->GetYaxis()->SetLabelSize(20);
    data_h->GetYaxis()->SetLabelOffset(0.02);
    data_h->GetYaxis()->SetTitleFont(43);
    data_h->GetYaxis()->SetTitleSize(24);
    data_h->GetYaxis()->SetTitle("Data/MC");
    data_h->GetYaxis()->SetTitleOffset(1.2);
    data_h->GetXaxis()->SetTickLength(0.1);
    data_h->GetXaxis()->SetTitle(Form("%s: %s", mode.Data(), x_name.Data()));
    can->Update();
    cout << "\tFinish drawing subpad!!!" << endl;

    can->SaveAs(Form("%s/plots/bdt_dis/%s_%d.pdf", str_dir.Data(), mode.Data(), region));

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/plots/bdt_dis/%s_%d.pdf", str_dir.Data(), mode.Data(), region)) << endl;
}
