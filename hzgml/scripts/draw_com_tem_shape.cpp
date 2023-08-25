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

void draw_com_tem_shape(TString mode = "inclusive", int lep_id = 13, int isbarrel = 1)
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString s_LG[] = {"Fake tem.", "Non-prompt CR", "True tem."};//"SM ZG", "VBS ZG", "TTG", "TG"};
    TString l_LG[] = {"Ref. tem."};
    double Y_axis[3] = {0, 0.02, 100};
    double pos[2] = {0.55, 0.95};
    double ratio = 1.3; // the ratio of top limit to the highest bin
    int stypes = 3, ltypes = 1;

    int ptbins_temp1[] = {15,17,20,25,30,35,40,50,5000};
    int ptbins_temp0[] = {15,17,20,25,30,40,5000};

    TString lepton, region;
    int *ptbins, num_ptbin;

    if (lep_id == 11) lepton = "ele";
    if (lep_id == 13) lepton = "mu";

    if (isbarrel == 1){
        region = "Barrel";
        ptbins = ptbins_temp1;
        num_ptbin = 8;
    }
    if (isbarrel == 0){
        region = "Endcap";
        ptbins = ptbins_temp0;
        num_ptbin = 6;
    }

    TString x_name = "l1g_deltaR";
    TString unit = "";
    Double_t data;
    double x_axis[3] = {30, 0, 4.8};
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    TString ftem[][6] = {{"data", "ZGToLLG", "TTGJets", "TGJets", "ZG2JToG2L2J", "0"}, {"data", "ZGToLLG", "TTGJets", "TGJets", "ZG2JToG2L2J", "0"}, {"ZGToLLG", "TTGJets", "TGJets", "ZG2JToG2L2J", "0"}};//{"ZGToLLG", "ZG2JToG2L2J", "TTGJets", "TGJets", "0"}};//, {"ZG2JToG2L2J", "0"}, {"TTGJets", "0"}, {"TGJets", "0"}, {"dy", "0"}, {"tt", "0"}, {"ZG2JToG2L2J", "0"}, {"ww", "wz", "zz", "0"}};
    TString fref[][6] = {{"data", "ZGToLLG", "TTGJets", "TGJets", "ZG2JToG2L2J", "0"}};

    cout << "\n\n\tFinish basic setting" << endl;

    for (int bini=0; bini<num_ptbin; bini++){

    // cout << "\n\tTarget: " << Form("%s/overlap/figures/com/check_fake_template_data_%s_%s_%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data(), region.Data(), lepton.Data()) << endl;
    cout << "\n\tTarget: " << Form("%s/overlap/figures/com/check_fake_template_data_%s_%s_%s_%s_pt%d_%d.pdf", str_dir.Data(), x_name.Data(), mode.Data(), region.Data(), lepton.Data(), ptbins[bini], ptbins[bini+1]) << endl;

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
    TPad *pad1 = new TPad("pad1", "", 0, 0.2, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.28);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Double_t w, fw, tw, iso, sieie, cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9;
    Long64_t n, id;
    double max = 0.0, min = 0.0, num[ltypes + stypes];

    TH1D *h1[stypes];
    // TH1D *ftem_h = new TH1D("ftem", Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);
    TH1D *fref_h =  new TH1D("fref", Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);
    // ftem_h->SetLineColor(TColor::GetColorDark(2));
    // ftem_h->SetFillColorAlpha(2, 0.5);
    // ftem_h->SetFillStyle(1001);;
    // ftem_h->SetLineWidth(3);
    fref_h->SetLineColor(TColor::GetColorDark(2));
    fref_h->SetFillColorAlpha(2, 0.5);
    fref_h->SetFillStyle(1001);;
    fref_h->SetLineWidth(3);

    can->cd(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    for (int typei = 0; typei < stypes; typei++)
    {
        h1[typei] = new TH1D(Form("h1_%d", typei), Form("%s; %s; Events/%2f", mode.Data(), x_name.Data(), (x_axis[2]-x_axis[1])/x_axis[0]), x_axis[0], x_axis[1], x_axis[2]);
        h1[typei]->SetLineColor(TColor::GetColorDark(typei + 3));
        // h1[typei]->SetFillColor(typei + 2);
        // h1[typei]->SetFillStyle(3004);
        h1[typei]->SetLineWidth(3);

        int len = 0;
        num[typei] = 0.0;
        while(ftem[typei][len] != "0"){
            cout << "\n\tLoading file" << Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), ftem[typei][len].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), ftem[typei][len].Data()));

            TTree *t = (TTree *)f->Get(mode.Data());
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name.Data(), &data);
            t->SetBranchAddress("weight", &w);
            t->SetBranchAddress("H_mass", &cut1);
            t->SetBranchAddress("n_jets", &n);
            t->SetBranchAddress("Z_lead_lepton_id", &id);
            t->SetBranchAddress("gamma_pt", &cut3);
            t->SetBranchAddress("gamma_chiso", &iso);
            // t->SetBranchAddress("gamma_sieie", &sieie);
            if ((len>0 & typei<2) | (typei==2))
            t->SetBranchAddress("gamma_genPartFlav", &cut4);
            t->SetBranchAddress("gamma_eta", &cut5);
            t->SetBranchAddress("gamma_mvaID", &cut6);

            cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                << "\n\tloading...loading...zzzzzz\n"
                << endl;

            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                // if (n!=0) continue;
                // if (typei == 1 && n != 0) continue;
                // if (cut5/cut4 > 0.48) continue;
                // if (cut2<15) continue;
                // if (typei==1) w = w * tw;
                // if (typei==0) w = w * fw;
                if ((cut3<ptbins[bini]) | (cut3>ptbins[bini+1])) continue;
                if ((len>0 & typei<2) | (typei==2)) if (cut4!=1) continue;
                if (id!=lep_id) continue;
                // if (typei==0) if (iso*cut3>3) continue;
                // if (typei==1) if (iso*cut3>5) continue;
                // if (cut6<-0.97) if (iso*cut3>1.5) continue;
                if (isbarrel==1){
                    if (abs(cut5)>1.5) continue;
                    if (typei==0) {
                        if (cut6<-0.97 | cut6>-0.4) continue;
                        if (iso*cut3>1.141) continue;
                        // if (iso*cut3>1.694) continue;
                        }

                    if (typei==1) {
                        if (cut6<-0.4 | cut6>0.42) continue;
                        if (iso*cut3>1.141) continue;
                        }
                    if (iso*cut3>1.141) continue;
                    if (typei==2) if (cut6<0.42) continue;
                    // if (sieie>0.01015) continue;
                }
                if (isbarrel==0){
                    if ((abs(cut5)>2.6) | (abs(cut5)<1.5)) continue;
                    if (typei==0) {
                        if (cut6<-0.97 | cut6>-0.59) continue;
                        if (iso*cut3>1.051) continue;
                        // if (iso*cut3>2.089) continue;
                        }
                    if (typei==1) {
                        if (cut6<-0.59 | cut6>0.14) continue;
                        if (iso*cut3>1.051) continue;
                        }
                    if (iso*cut3>1.051) continue;
                    if (typei==2) if (cut6<0.14) continue;
                    // if (sieie>0.0272) continue;
                }
                if ((len>0) & (typei<2)) w = -w;
                // if ((abs(cut5)<1.5 & cut6>-0.4) | (abs(cut5)>1.5 & cut6>-0.58)) continue;
                if ( (cut1 < 120 | cut1 > 130))//&& cut3 > 0.22 && cut4 < 10 && cut5 > 0.7)
                {
                    h1[typei]->Fill(data, w);
                    num[typei] += w;
                    // ftem_h->Fill(data, w);
                }
            }
            len++ ;
        }
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei] << " items are in it!!";
    }
    // ftem_h->Scale(1. / ftem_h->Integral(), "width");
    // max = ftem_h->GetBinContent(ftem_h->GetMaximumBin());
    // ftem_h->Draw("HIST");
    pad1->SetLeftMargin(0.15);
    pad1->SetTopMargin(0.05);

    // for (int typei = 0; typei < stypes; typei++)
    // {
    //     h1[typei]->Reset("ICES");
    // }
    can->Update();

    // TH1F *h2[ltypes];

    for (int typei = 0; typei < ltypes; typei++)
    {
        // h2[typei] = new TH1F(Form("h2_%d", typei), Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);
        // h2[typei]->SetLineColor(TColor::GetColorDark(typei+1));
        // h2[typei]->SetLineStyle(1);
        // h2[typei]->SetLineWidth(1);
        
        int len = 0;
        num[typei + stypes] = 0.0;
        while (fref[typei][len] != "0")
        {
            cout << "\n\tLoading file" << Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), fref[typei][len].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), fref[typei][len].Data()));

            TTree *t = (TTree *)f->Get(mode.Data());
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name.Data(), &data);
            t->SetBranchAddress("weight", &w);
            t->SetBranchAddress("H_mass", &cut1);
            t->SetBranchAddress("n_jets", &n);
            t->SetBranchAddress("Z_lead_lepton_id", &id);
            t->SetBranchAddress("gamma_pt", &cut3);
            t->SetBranchAddress("gamma_chiso", &iso);
            if (len) t->SetBranchAddress("gamma_genPartFlav", &cut4);
            t->SetBranchAddress("gamma_mvaID", &cut6);
            t->SetBranchAddress("gamma_eta", &cut5);

            cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                 << "\n\tloading...loading...zzzzzz\n"
                 << endl;

            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                // if (n!=1) continue;
                // if (cut5/cut4 > 0.48) continue;
                // if (cut4<20) continue;
                // if (cut2<15) continue;
                // if (cut3<0.4) continue;
                // if (cut6<-0.8) continue;
                if ((cut3<ptbins[bini]) | (cut3>ptbins[bini+1])) continue;
                if (len) if (cut4!=1) continue;
                if (id!=lep_id) continue;
                if (isbarrel==1){
                    if (abs(cut5)>1.5) continue;
                    if (cut6<0.42) continue;
                    // if (iso*cut3<1.141) continue;
                }
                if (isbarrel==0){
                    if (abs(cut5)<1.5) continue;
                    if (cut6<0.14) continue;
                    // if (iso*cut3<1.051) continue;
                }
                if (len) w = -w;
                // if ((abs(cut5)<1.5 & cut6<0.42) | (abs(cut5)>1.5 & cut6<0.14)) continue;
                if ( (cut1 < 120 | cut1 > 130) )//&& cut3 > 0.22 && cut4 < 10 && cut5 > 0.7)
                {
                    // h2[typei]->Fill(cut1);
                    num[stypes+typei] += w;
                    fref_h->Fill(data, w);
                }
            }
            cout << "\n\tNumbers in HIST: " << Form("%.3f", num[stypes+typei]) << endl;
            len++;
        }
        fref_h->Scale(1. / fref_h->Integral(), "width");
        if (max < fref_h->GetBinContent(fref_h->GetMaximumBin()))
            max = fref_h->GetBinContent(fref_h->GetMaximumBin());
        fref_h->Draw("HIST, same");
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei + stypes] << " items are in it!!";
    }
    for (int typei = 0; typei < stypes; typei++)
    {
        h1[typei]->Scale(1. / h1[typei]->Integral(), "width");
        if (max < h1[typei]->GetBinContent(h1[typei]->GetMaximumBin()))
            max = h1[typei]->GetBinContent(h1[typei]->GetMaximumBin());
        h1[typei]->Draw("HIST, same");
    }
    fref_h->GetYaxis()->SetLabelFont(43);
    fref_h->GetYaxis()->SetLabelSize(20);
    fref_h->GetYaxis()->SetTitleFont(43);
    fref_h->GetYaxis()->SetTitleSize(24);
    fref_h->GetXaxis()->SetLabelFont(43);
    fref_h->GetXaxis()->SetLabelSize(20);
    fref_h->GetXaxis()->SetTitleFont(43);
    fref_h->GetXaxis()->SetTitleSize(24);
    // fref_h->GetXaxis()->SetTitle(Form("%s:%s", mode.Data(), x_name.Data()));
    fref_h->GetXaxis()->SetTitle("");
    fref_h->GetYaxis()->SetTitle(Form("Rel. value/%.2g%s" , (x_axis[2]-x_axis[1])/x_axis[0], unit.Data()));
    fref_h->GetYaxis()->SetTitleOffset(1.6);
    fref_h->SetMaximum(max * ratio);
    can->Update();
    cout << "\n\tFinish drawing lines" << endl;

    TLegend *legend = new TLegend(pos[0], 0.775, pos[1], 0.95);
    for (int typei = 0; typei < stypes; typei++)
        legend->AddEntry(Form("h1_%d", typei), Form("%s(%.2f)", s_LG[typei].Data(), num[typei]), "l");
    // for (int typei = 0; typei < ltypes; typei++)
    //     legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", l_LG[typei].Data(), num[typei + stypes]), "l");
    legend->AddEntry("fref", Form("%s(%.2f)", l_LG[0].Data(), num[stypes]), "l");
    // legend->AddEntry("fref", Form("%s(%.2f)", s_LG[0].Data(), num[0]), "l");
    legend->SetBorderSize(0);
    // legend->SetFillColor(0);
    legend->SetFillColorAlpha(0, 0);
    legend->Draw();

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextColor(1);
    text->SetTextSize(0.04);
    text->SetTextAlign(22);
    if (bini==num_ptbin-1) text->DrawLatex(0.35, 0.875, Form("%s, %s, P_{T}^{#gamma} > %d", region.Data(), lepton.Data(), ptbins[bini])); // TODO: check it before run
    // else if (i==0) text->DrawLatex(0.5, 0.5, Form("%s < %1.1f", y_lg.Data(), max)); // TODO: check it before run
    else // TODO: check it before run
        text->DrawLatex(0.35, 0.875, Form("%s, %s, %d < P_{T}^{#gamma} < %d", region.Data(), lepton.Data(), ptbins[bini], ptbins[bini+1]));
    // text->DrawLatex(0.35, 0.875, Form("%s, %s", region.Data(), lepton.Data()));

    can->Update();
    cout << "\n\tFinish Create TLegend" << endl;

    can->cd(0);
    pad2->cd();
    pad2->SetGrid(0, 3);
    pad2->SetLeftMargin(0.15);
    pad2->SetBottomMargin(0.35);
    setStyle();
    gStyle->SetTitleXSize(0.28);
    gStyle->SetTitleYSize(0.32);
    TH1D *div_h;
    div_h = (TH1D*) h1[1]->Clone();
    div_h->SetLineColor(1);
    div_h->SetLineWidth(1);
    div_h->Divide(h1[0]);
    int max_bin = div_h->GetMaximumBin();
    max = div_h->GetBinContent(max_bin);
    double max_err = div_h->GetBinError(max_bin);
    min = div_h->GetMinimum(0);
    int min_bin = 0;
    div_h->GetBinWithContent(min, min_bin, 0);
    double min_err = div_h->GetBinError(min_bin);
    div_h->SetMaximum(max + max_err*1.5);
    div_h->SetMinimum(min - min_err*1.5);
    div_h->SetMaximum(3);
    div_h->SetMinimum(0);
    div_h->Draw("E1");
    div_h->GetXaxis()->SetLabelFont(43);
    div_h->GetXaxis()->SetLabelSize(20);
    div_h->GetXaxis()->SetLabelOffset(0.02);
    div_h->GetXaxis()->SetTitleFont(43);
    div_h->GetXaxis()->SetTitleSize(30);
    div_h->GetXaxis()->SetTitleOffset(3);
    div_h->GetYaxis()->SetLabelFont(43);
    div_h->GetYaxis()->SetLabelSize(20);
    div_h->GetYaxis()->SetLabelOffset(0.02);
    div_h->GetYaxis()->SetTitleFont(43);
    div_h->GetYaxis()->SetTitleSize(24);
    div_h->GetYaxis()->SetTitle("Non./Fake");
    div_h->GetYaxis()->SetTitleOffset(1.6); 
    div_h->GetXaxis()->SetTickLength(0.1);
    div_h->GetXaxis()->SetTitle(Form("%s%s", x_name.Data(), unit.Data())); //TODO: check the unit of variables
    can->Update();
    cout << "\tFinish drawing subpad!!!" << endl;

    // can->SetGrid(2, 1);
    can->SaveAs(Form("%s/overlap/figures/com/check_fake_template_data_%s_%s_%s_%s_pt%d_%d.pdf", str_dir.Data(), x_name.Data(), mode.Data(), region.Data(), lepton.Data(), ptbins[bini], ptbins[bini+1]));
    // can->SaveAs(Form("%s/overlap/figures/com/check_fake_template_data_%s_%s_%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data(), region.Data(), lepton.Data()));

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/com/check_fake_template_data_%s_%s_%s_%s_pt%d_pt%d.pdf", str_dir.Data(), x_name.Data(), mode.Data(), region.Data(), lepton.Data(), ptbins[bini], ptbins[bini+1])) << endl;

    fref_h->Reset("ICES");
    for (int typei = 0; typei < stypes; typei++)
    {
        h1[typei]->Reset("ICES");
    }
    // ftem_h->Reset("ICES");
    div_h->Reset("ICES");
    can->Clear();
    }
}
