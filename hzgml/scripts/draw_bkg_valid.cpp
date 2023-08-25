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

void draw_bkg_valid(TString mode = "zero_jet", int button = 1)
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    int sig_ratio = 1000;
    TString LG[][7] = {{"Fake photon", "SM ZG", "EWK ZG", "TTGJets", "TGJets", "Data", Form("Signal(x%d)", sig_ratio)}, {"DY+fake", "SM ZG", "EWK ZG", "TT", "Diboson", "Data", Form("Signal(x%d)", sig_ratio)}};
    double Y_axis[3] = {0, 0.2, 100};

    TString x_name = "Z_cos_theta";
    int X_axis[]={40, -1, 1};
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    TString stype[][7][4] = {{{"data_fake", "mc_true", "0"}, {"mc_med", "0"}, {"mc_med", "0"}, {"mc_med", "0"}, {"mc_med", "0"}}, {{"DYJetsToLL", "0"}, {"ZGToLLG", "0"}, {"LLAJJ", "0"}, {"TT", "0"}, {"WW", "WZ", "ZZ", "0"}}};//, {"mc_med", "0"}};//, "sig"};
    TString ltype[][7] = {{"data", "0"}, {"ggH", "VBF", "WminusH", "WplusH", "ttH", "0"}};//, {"ggH", "0"}, {"VBF", "0"}, {"WminusH", "WplusH", "ZH", "0"}, {"ttH", "0"}};

    int bins;
    int stypes = 2, ltypes = 2;
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
    // ifstream file(Form("%s/significances/bin_binaries_%s.txt", str_dir.Data(), mode.Data()));
    // cout << "\n\tThe binaries of bin are: ";
    // file >> bins >> setw(1);
    // cout << bins << ", ";
    // Float_t edges[bins + 1];
    // Float_t Significance[bins];
    // for (int n = 0; n < bins + 1; n++)
    // {
    //     file >> edges[n] >> setw(1);
    //     cout << edges[n] << " ";
    // }
    // for (int n = 0; n < bins; n++)
    // {
    //     file >> Significance[n] >> setw(1);
    //     cout << Significance[n] << " ";
    // }
    // file >> res >> setw(1) >> reserr;

    cout << "\n\n\tFinish basic setting" << endl;

    Double_t w, hmass, wp, data, bdt;
    Long64_t tag;
    double max=0.0, min=0.0;

    // if (mode.Contains("VH_ttH", TString::kIgnoreCase)) X_axis[0]=16;

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
    // can->SetLogy(1);
    TPad *pad1 = new TPad("pad1", "", 0, 0.23, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    THStack *hs = new THStack("hs", "test stacked histograms");
    // TH1D *mc_h = new TH1D("mc", Form("%s; %s; Events", mode.Data(), x_name.Data()), bins, edges);
    // TH1D *data_h =  new TH1D("data", Form("%s; %s; Events", mode.Data(), x_name.Data()), bins, edges);
    TH1D *mc_h = new TH1D("mc", Form("%s; %s; Events", mode.Data(), x_name.Data()), X_axis[0], X_axis[1], X_axis[2]);
    TH1D *data_h =  new TH1D("data", Form("%s; %s; Events", mode.Data(), x_name.Data()), X_axis[0], X_axis[1], X_axis[2]);
    TH1D *h1[stypes];

    can->cd(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    double sum = 0.;

    for (int typei = 0; typei < stypes; typei++)
    {
        // h1[typei] = new TH1D(Form("h1_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), bins, edges);
        h1[typei] = new TH1D(Form("h1_%d", typei), Form("%s; bdt_score_t; num/1GeV", mode.Data()), X_axis[0], X_axis[1], X_axis[2]);
        h1[typei]->SetLineColor(TColor::GetColorDark(typei + 2));
        h1[typei]->SetFillColor(TColor::GetColorDark(typei + 2));
        h1[typei]->SetLineWidth(1);

        int len = 0;
        num[typei] = 0.0;
        while (stype[button][typei][len] != "0")
        {
        cout << "\n\tLoading file" << Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), stype[button][typei][len].Data()) << "!......" << endl;
        TFile *f = new TFile(Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), stype[button][typei][len].Data()));

        TTree *t = (TTree *)f->Get(mode);
        Long64_t nentries = t->GetEntries();
        t->SetBranchAddress(x_name.Data(), &data);
        t->SetBranchAddress("weight", &w);
        t->SetBranchAddress("H_mass", &hmass);
        t->SetBranchAddress("gamma_mvaID_WP80", &bdt);
        if (!button) if (typei) t->SetBranchAddress("tagger", &tag);

        cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
             << "\n\tloading...loading...zzzzzz\n"
             << endl;

        for (Long64_t i = 0; i < nentries; i++)
        {
            t->GetEntry(i);
            // if (bdt<edges[n] | bdt>edges[n+1]) continue;
            // cout<<i<<endl;
            if (!bdt) continue;
            if (!button) if (typei>0) if (tag!=typei) continue;
            if (!x_name.Contains("H_mass", TString::kIgnoreCase)) if (hmass>120 & hmass<130) continue;
            if (x_name.Contains("H_mass", TString::kIgnoreCase)) data = hmass;
            h1[typei]->Fill(data, w);
            mc_h->Fill(data, w);
            num[typei] += w;
            if (hmass<120 | hmass>130) sum += w;
        }
        // f->Close();
        len++;
        }
        hs->Add(h1[typei]);
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei] << " items are in it!!";
    }
    max = mc_h->GetBinContent(mc_h->GetMaximumBin());
    min = mc_h->GetBinContent(mc_h->GetMinimumBin());

    TH1D *h2[ltypes];

    for (int typei = 0; typei < ltypes; typei++)
    {
        // h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), bins, edges);
        h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; bdt_score_t; num", mode.Data()), X_axis[0], X_axis[1], X_axis[2]);
        h2[typei]->SetLineColor(typei + 1);
        if (typei){
            h2[typei]->SetLineColor(kAzure+1);
            h2[typei]->SetFillColor(kAzure+1);
            h2[typei]->SetFillStyle(3325);
        }
        h2[typei]->SetLineStyle(1);
        h2[typei]->SetLineWidth(2);
        int len = 0;
        num[typei + stypes] = 0.0;
        while (ltype[typei][len] != "0")
        {
            cout << "\n\tLoading file" << Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), ltype[typei][len].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), ltype[typei][len].Data()));

            TTree *t = (TTree *)f->Get(mode);
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name.Data(), &data);
            t->SetBranchAddress("weight", &w);
            t->SetBranchAddress("H_mass", &hmass);
            t->SetBranchAddress("gamma_mvaID_WP80", &bdt);

            cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                 << "\n\tloading...loading...zzzzzz\n"
                 << endl;

            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                // if (bdt<edges[n] | bdt>edges[n+1]) continue;
                if (!bdt) continue;
                if (x_name.Contains("H_mass", TString::kIgnoreCase)) data = hmass;
                if ((typei==0) & (hmass < 120 | hmass > 130))
                {
                    h2[typei]->Fill(data, w);
                    data_h->Fill(data, w);
                    num[typei + stypes] += w;
                }
                if (typei)
                {
                    w = w*sig_ratio;
                    h2[typei]->Fill(data, w);
                    num[typei + stypes] += w;
                }
            }
            // f->Close();
            len++;
        }
        if (h2[typei]->GetBinContent(h2[typei]->GetMaximumBin()) > max)
            max = h2[typei]->GetBinContent(h2[typei]->GetMaximumBin());
        // hs->SetMaximum(TMath::Power(10, TMath::Log10(max)+0.1));
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei + stypes] << " items are in it!!";
    }

    // double factor = num[stypes] / sum;
    // for (int typei=0; typei<stypes; typei++){
    //     h1[typei]->Scale(factor);
    //     num[typei] = num[typei]*factor;
    // }
    // if (max == mc_h->GetBinContent(mc_h->GetMaximumBin()))
    //     max = max * factor;
    // mc_h->Scale(factor);
    hs->SetMaximum(max*1.4);
    hs->Draw("HIST");
    // hs->GetXaxis()->SetTitle(Form("%s:%s", mode.Data(), x_name.Data()));
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.05);
    // pad1->SetLogy(1);
    hs->GetYaxis()->SetLabelFont(43);
    hs->GetYaxis()->SetLabelSize(20);
    hs->GetYaxis()->SetTitleFont(43);
    hs->GetYaxis()->SetTitleSize(24);
    hs->GetXaxis()->SetLabelFont(43);
    hs->GetXaxis()->SetLabelSize(20);
    hs->GetXaxis()->SetTitleFont(43);
    hs->GetXaxis()->SetTitleSize(24);
    hs->GetYaxis()->SetTitle(Form("Events/%dGeV", (X_axis[2]-X_axis[1])/X_axis[0]));
    hs->GetYaxis()->SetTitleOffset(1.2);
    if (Y_axis[0] > 0)
    {
        hs->SetMinimum(TMath::Power(10, TMath::Log10(min)-0.5));
        hs->SetMaximum(TMath::Power(10, TMath::Log10(max)+0.2));
    }
    h2[0]->Draw("E1, same");
    for (int typei=1; typei<ltypes; typei++){
        h2[typei]->Draw("HIST, same");
    }
    data_h->Divide(mc_h);
    data_h->Sumw2();
    double kai2 = 0.;
    int good_bin = 0;
    for (int bini=0; bini<X_axis[0]; bini++){
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
    // text->DrawLatex(0.3, 0.9, "Binary: (0.00, 1.00)");
    // text->DrawLatex(0.5, 0.9, Form("Scale factor: %.4f", factor));
    text->DrawLatex(0.4, 0.9, Form("#chi^{2}/n: %.4f", kai2/good_bin));

    TLegend *legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    for (int typei = 0; typei < stypes; typei++)
        legend->AddEntry(Form("h1_%d", typei), Form("%s(%.2f)", LG[button][typei].Data(), num[typei]), "f");
    for (int typei = 0; typei < ltypes; typei++){
        if (typei==0) legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", LG[button][typei + stypes].Data(), num[typei + stypes]), "lep");
        else legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", LG[button][typei + stypes].Data(), num[typei + stypes]), "f");
    }
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
    data_h->GetXaxis()->SetTitle(Form("%s: M_{ll#gamma}[GeV]", mode.Data()));
    can->Update();
    cout << "\tFinish drawing subpad!!!" << endl;

    can->SaveAs(Form("%s/plots/variables/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data()));

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/plots/variables/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data())) << endl;
}
