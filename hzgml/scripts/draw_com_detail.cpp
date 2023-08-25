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

void draw_com_detail(TString mode = "inclusive")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString s_LG[] = {"Fake photons1", "Fake photons2", "Fake photons3", "Fake photons4", "SM ZG", "VBS ZG", "TTG", "TG"};//, "Diboson"};
    TString l_LG[] = {"data"};
    double Y_axis[3] = {0, 0.02, 100};
    double pos[2] = {0.75, 0.95};
    double ratio = 1.1; // the ratio of top limit to the highest bin
    int stypes = 8, ltypes = 1;

    TString x_name = "Z_lead_lepton_eta";
    TString X_name = "Inclusive: #eta_{lead l}";
    Double_t data;
    double x_axis[3] = {50, -2.5, 2.5};
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    TString stype[][4] = {{"data_fake", "mc_true", "0"}, {"data_fake", "mc_true", "0"}, {"data_fake", "mc_true", "0"}, {"data_fake", "mc_true", "0"}, {"zg", "0"}, {"llajj", "0"}, {"ttg", "0"}, {"tg", "0"}};//, {"dy", "0"}, {"tt", "0"}, {"llajj", "0"}, {"ww", "wz", "zz", "0"}};
    TString ltype[][4] = {{"data", "0"}};

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
    can->SetLeftMargin(0.05);
    can->SetTopMargin(0);
    can->SetRightMargin(0.05);
    can->SetBottomMargin(0);
    can->SetLogy(0);
    TPad *pad1 = new TPad("pad1", "", 0, 0.19, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Double_t w, gam_w, ele_w, mu_w, fw, cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9, n;
    Long64_t ismu, nj;
    double max = 0.0, num[ltypes + stypes];

    THStack *hs = new THStack("hs", "test stacked histograms");
    TH1D *h1[stypes];
    TH1D *mc_h = new TH1D("mc", Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);
    TH1D *data_h =  new TH1D("data", Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);

    can->cd(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    for (int typei = 0; typei < stypes; typei++)
    {
        h1[typei] = new TH1D(Form("h1_%d", typei), Form("%s; %s; Events/%2f", mode.Data(), x_name.Data(), (x_axis[2]-x_axis[1])/x_axis[0]), x_axis[0], x_axis[1], x_axis[2]);
        h1[typei]->SetLineColor(typei + 2);
        h1[typei]->SetFillColor(typei + 2);
        // h1[typei]->SetFillStyle(3004);
        h1[typei]->SetLineWidth(1);

        int len = 0;
        num[typei] = 0.0;
        while(stype[typei][len] != "0"){
            cout << "\n\tLoading file" << Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), stype[typei][len].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), stype[typei][len].Data()));

            TTree *t = (TTree *)f->Get(mode.Data());
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name.Data(), &data);
            t->SetBranchAddress("weight_central", &w);
            // t->SetBranchAddress("weight_photon_id_sf_Photon_central", &gam_w);
            // t->SetBranchAddress("weight_electron_id_sf_SelectedElectron_central", &ele_w);
            // t->SetBranchAddress("weight_muon_id_sf_SelectedMuon_central", &mu_w);
            if (typei < 4)
                t->SetBranchAddress("fake_weight", &fw);
            t->SetBranchAddress("H_mass", &cut1);
            // t->SetBranchAddress("H_pt", &cut2);
            // t->SetBranchAddress("gamma_eta", &cut3);
            t->SetBranchAddress("Z_sublead_lepton_pt", &cut4);
            t->SetBranchAddress("gamma_pt", &cut5);
            t->SetBranchAddress("photon_is_barrel", &cut6);
            t->SetBranchAddress("gamma_photon_match", &cut7);
            t->SetBranchAddress("photon_selection", &cut8);
            t->SetBranchAddress("Z_sublead_lepton_id", &ismu);
            // t->SetBranchAddress("n_jets", &nj);
            if (typei > 3)
                t->SetBranchAddress("photon_isprompt", &cut9);

            cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                << "\n\tloading...loading...zzzzzz\n"
                << endl;

            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                // if (typei == 1 && n != 0) continue;
                // if (cut5/cut4 > 0.48) continue;
                // if (cut2>10) continue;
                if (cut4<15) continue;
                if (cut5<20) continue;
                // if (abs(cut3)>1.4442 && abs(cut3)<1.566) continue;
                // if (abs(ismu)!=11) continue;
                // if (nj>0) continue;
                if ((typei < 4) && (int(cut8/(1<<(typei+1)))%2 != 1)) continue;
                // if (typei < 4) continue;
                if ((typei > 3) && ((cut7 != 1) || (int(cut8/(1<<0))%2 != 1) || (cut9 != 1))) continue;
                if ( (cut1 < 120 | cut1 > 130) && (typei > 3))//&& cut3 > 0.22 && cut4 < 10 && cut5 > 0.7)
                {
                    h1[typei]->Fill(data, w);//*ele_w*mu_w);//*gam_w);
                    num[typei] += w;//*ele_w*mu_w;//*gam_w;
                    mc_h->Fill(data, w);//*ele_w*mu_w);//*gam_w);
                }
                if ( (cut1 < 120 | cut1 > 130) && (typei < 4))
                {
                    h1[typei]->Fill(data, w*fw);
                    num[typei] += w*fw;
                    mc_h->Fill(data, w*fw);
                }
            }
            len++ ;
        }
        hs->Add(h1[typei]);
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei] << " items are in it!!";
    }
    max = mc_h->GetBinContent(mc_h->GetMaximumBin());
    hs->Draw("HIST");
    pad1->SetLeftMargin(0.15);
    pad1->SetTopMargin(0.05);
    hs->GetYaxis()->SetLabelFont(43);
    hs->GetYaxis()->SetLabelSize(20);
    hs->GetYaxis()->SetTitleFont(43);
    hs->GetYaxis()->SetTitleSize(28);
    hs->GetXaxis()->SetLabelFont(43);
    hs->GetXaxis()->SetLabelSize(20);
    hs->GetXaxis()->SetTitleFont(43);
    hs->GetXaxis()->SetTitleSize(28);
    hs->GetXaxis()->SetTitle(" ");
    hs->GetYaxis()->SetTitle(Form("Events/%.2g" , (x_axis[2]-x_axis[1])/x_axis[0]));
    hs->GetYaxis()->SetTitleOffset(1.5);

    // for (int typei = 0; typei < stypes; typei++)
    // {
    //     h1[typei]->Reset("ICES");
    // }
    can->Update();

    TH1D *h2[ltypes];

    for (int typei = 0; typei < ltypes; typei++)
    {
        h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; %s; num", mode.Data(), x_name.Data()), x_axis[0], x_axis[1], x_axis[2]);
        h2[typei]->SetLineColor(TColor::GetColorDark(typei+1));
        h2[typei]->SetLineStyle(1);
        h2[typei]->SetLineWidth(1);
        
        int len = 0;
        while (ltype[typei][len] != "0")
        {
            cout << "\n\tLoading file" << Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), ltype[typei][len].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/skimmed_ntuples/%s/2017.root", str_dir.Data(), ltype[typei][len].Data()));

            TTree *t = (TTree *)f->Get(mode.Data());
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name.Data(), &data);
            t->SetBranchAddress("weight_central", &w);
            //t->SetBranchAddress("weight_photon_id_sf_Photon_central", gam_w);
            //t->SetBranchAddress("weight_electron_id_sf_SelectedElectron_central", ele_w);
            //t->SetBranchAddress("weight_muon_id_sf_SelectedMuon_central", mu_w);
            t->SetBranchAddress("H_mass", &cut1);
            // t->SetBranchAddress("H_pt", &cut2);
            // t->SetBranchAddress("gamma_eta", &cut3);
            t->SetBranchAddress("Z_sublead_lepton_pt", &cut4);
            t->SetBranchAddress("gamma_pt", &cut5);
            t->SetBranchAddress("photon_is_barrel", &cut6);
            t->SetBranchAddress("gamma_photon_match", &cut7);
            t->SetBranchAddress("photon_selection", &cut8);
            t->SetBranchAddress("Z_sublead_lepton_id", &ismu);
            // t->SetBranchAddress("n_jets", &nj);

            cout << "\n\tStart drawing\n\tentries are : " << nentries << " ! oh my god ! too much !"
                 << "\n\tloading...loading...zzzzzz\n"
                 << endl;

            num[typei + stypes] = 0.0;
            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                // if (cut5/cut4 > 0.48) continue;
                // if (cut2>10) continue;
                if (cut4<15) continue;
                if (cut5<20) continue;
                // if (abs(cut3)>1.4442 && abs(cut3)<1.566) continue;
                // if (abs(ismu)!=11) continue;
                // if (nj>0) continue;
                if ((cut7 != 1) || (int(cut8/(1<<0))%2 != 1)) continue;
                if ( (cut1 < 120 | cut1 > 130) )//&& cut3 > 0.22 && cut4 < 10 && cut5 > 0.7)
                {
                    h2[typei]->Fill(data);
                    num[stypes+typei] += 1;
                    data_h->Fill(data);
                }
            }
            len++;
        }
        if (max < data_h->GetBinContent(data_h->GetMaximumBin()))
            max = data_h->GetBinContent(data_h->GetMaximumBin());
        // h2->GetXaxis()->SetTitle(" ");
        h2[typei]->Draw("E1, same");
        cout << "\tFinished loading File" << endl;
        cout << "\n\t" << num[typei + stypes] << " items are in it!!";
    }

    hs->SetMaximum(max * ratio);
    can->Update();
    TLegend *legend = new TLegend(pos[0], 0.7, pos[1], 0.95);
    for (int typei = 0; typei < stypes; typei++)
        legend->AddEntry(Form("h1_%d", typei), Form("%s(%.2f)", s_LG[typei].Data(), num[typei]), "f");
    for (int typei = 0; typei < ltypes; typei++)
        legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", l_LG[typei].Data(), num[typei + stypes]), "lep");
    legend->SetBorderSize(0);
    // legend->SetFillColor(0);
    legend->SetFillColorAlpha(0, 0);
    legend->Draw();
    can->Update();
    cout << "\n\tFinish Create TLegend" << endl;

    can->cd(0);
    pad2->cd();
    pad2->SetGrid(0, 3);
    pad2->SetLeftMargin(0.15);
    pad2->SetBottomMargin(0.4);
    setStyle();
    gStyle->SetTitleXSize(0.28);
    gStyle->SetTitleYSize(0.32);
    data_h->SetLineColor(1);
    data_h->SetLineWidth(1);
    data_h->Divide(mc_h);
    data_h->Draw("E1");
    data_h->GetXaxis()->SetLabelFont(43);
    data_h->GetXaxis()->SetLabelSize(20);
    data_h->GetXaxis()->SetLabelOffset(0.02);
    data_h->GetXaxis()->SetTitleFont(43);
    data_h->GetXaxis()->SetTitleSize(28);
    data_h->GetXaxis()->SetTitleOffset(3.5);
    data_h->GetYaxis()->SetLabelFont(43);
    data_h->GetYaxis()->SetLabelSize(20);
    data_h->GetYaxis()->SetLabelOffset(0.02);
    data_h->GetYaxis()->SetTitleFont(43);
    data_h->GetYaxis()->SetTitleSize(28);
    data_h->GetYaxis()->SetTitle("data/MC");
    data_h->GetYaxis()->SetTitleOffset(1.5); 
    data_h->GetYaxis()->SetRangeUser(data_h->GetMinimum(0)*0.5,data_h->GetBinContent(data_h->GetMaximumBin())*1.2);
    data_h->GetXaxis()->SetTickLength(0.1);
    data_h->GetXaxis()->SetTitle(X_name);

    // TH1D h_mc_stat = mc_h->Clone(mc_h->GetName()+'_mcstat');
    // h_mc_stat->Divide(mc_h);
    // h_mc_stat->SetFillColor(kBlue);
    // h_mc_stat->SetFillStyle(3244);
    // h_mc_stat->Draw();

    can->Update();
    cout << "\tFinish drawing subpad!!!" << endl;

    // can->SetGrid(2, 1);
    can->SaveAs(Form("%s/overlap/figures/com/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data()));

    cout << Form("\tPath is %s", str_dir.Data()) << endl;
    cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/com/%s_%s.pdf", str_dir.Data(), x_name.Data(), mode.Data())) << endl;
}
