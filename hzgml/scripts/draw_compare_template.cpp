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

void draw_compare_template(TString mode = "inclusive")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString LG[]{"Data template", "Fake template", "True template(zg)", "Fake(dy,medium,non-prompt)", "Fake(dy,chiso SB)"};
    TString F[]{"data", "data", "zg", "dy", "dy"};
    double Y_axis[3] = {0, 0.02, 100};

    double X_axis[3] = {60, 0.00515, 0.02015};
    int isBarrel = 1;
    int pcut[] = {15, 17, 20, 25, 30, 5000};
    int n_cut = 4;
    
    // double X_axis[3] = {40, 0.0172, 0.0572};
    // int isBarrel = 0;
    // int pcut[] = {15, 20, 25, 5000};
    // int n_cut = 2;

    int ltypes = 3; //TODO: check the number of type
    int sum_num[ltypes];

    Double_t data;
    TString x_name = "photon_sieie";
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    // TString ltype[] = {"SM ZG", "DY+fake"};

    cout << "\n\tThere are " << ltypes << " kinds of Higgs to ZGamma background!!!" << endl;
    Double_t num[ltypes];

    cout << "\n\n\tFinish basic setting" << endl;

    Double_t w, n;
    Double_t pass1, pass2, isprompt, isbarrel, chiso, pt, err;
    double max;

    TH1D *h2[ltypes];

    for (int n_bin = n_cut; n_bin >= 0 ; n_bin--)
    {
        max = 0.0;
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
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);
        can->SetLogy(0);
        cout << "\n\tFinish Set Style and Canvas" << endl;

        for (int typei = 0; typei < ltypes; typei++)
        {
            h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; %s:%s; events", mode.Data(), mode.Data(), x_name.Data()), X_axis[0], X_axis[1], X_axis[2]);
            h2[typei]->SetLineColor(TColor::GetColorDark(typei + 3));
            h2[typei]->SetLineStyle(1);
            h2[typei]->SetLineWidth(2);
            cout << "\n\tLoading file" << Form("%s/overlap/%s.root", str_dir.Data(), F[typei].Data()) << "!......" << endl;

            TFile *f = new TFile(Form("%s/overlap/%s.root", str_dir.Data(), F[typei].Data()));

            TTree *t = (TTree *)f->Get(mode.Data());
            Long64_t nentries = t->GetEntries();
            t->SetBranchAddress(x_name, &data);
            if (typei>=2){
                t->SetBranchAddress("weight_central", &w);
                t->SetBranchAddress("photon_genPartFlav", &isprompt);
            }
            // t->SetBranchAddress("n_iso_photons", &n);
            t->SetBranchAddress("photon_chiso", &chiso);
            t->SetBranchAddress("photon_pt", &pt);
            t->SetBranchAddress("photon_is_barrel", &isbarrel);
            t->SetBranchAddress("gamma_photon_match", &pass1);
            t->SetBranchAddress("photon_selection", &pass2);
            t->SetBranchAddress("gamma_energyErr", &err);

            cout << "\n\tStart drawing " << x_name.Data() << "\n\tentries are : " << nentries << " ! oh my god ! too much !"
                << "\n\tloading...loading...zzzzzz\n"
                << endl;

            num[typei] = 0.0;
            sum_num[typei] = 0;
            double sum = 0.0;
            for (Long64_t i = 0; i < nentries; i++)
            {
                t->GetEntry(i);
                if (typei >= 2)
                    sum += w;
                else sum += 1;
                if (pass1!=1) continue;
                if (isbarrel != isBarrel) continue;
                if ((pt < pcut[n_bin]) | (pt > pcut[n_bin+1])) continue;
                if (err/pt>0.48) continue;
                if ((typei==4) & (int(pass2/(1<<6))%2==1) & (chiso*pt>5) & (chiso*pt<10))
                    {
                        h2[typei]->Fill(data, w);
                        num[typei] += w;
                        sum_num[typei]++;
                    }
                if ((typei==3) & (isprompt==0) & (int(pass2/(1<<5))%2==1))
                    {
                        h2[typei]->Fill(data, w);
                        num[typei] += w;
                        sum_num[typei]++;
                    }
                if ((typei==2) & (isprompt==1) & (int(pass2/(1<<5))%2==1))
                    {
                        h2[typei]->Fill(data, w);
                        num[typei] += w;
                        sum_num[typei]++;
                    }
                if ((typei==1) & (int(pass2/(1<<6))%2==1) & (chiso*pt>5) & (chiso*pt<10))
                    {
                        h2[typei]->Fill(data);
                        num[typei] += 1;
                        sum_num[typei]++;
                    }
                if ((typei==0) & (int(pass2/(1<<5))%2==1))
                    {
                        h2[typei]->Fill(data);
                        num[typei] += 1;
                        sum_num[typei]++;
                    }
            }
            h2[typei]->Scale(1. / h2[typei]->Integral(), "width");
            if (max < h2[typei]->GetBinContent(h2[typei]->GetMaximumBin()))
                max = h2[typei]->GetBinContent(h2[typei]->GetMaximumBin());
            cout << "\tFinished loading File" << endl;
            // cout << "\n\t" << sum_num[typei] << "/" << nentries  << " items are in it!!";
            cout << "\n\t" << num[typei] << "/" << sum  << " items are in it!!";
        }
        h2[0]->SetMaximum(max * 1.15);
        h2[0]->GetYaxis()->SetTitle(Form("Relative Value/%1.1e" , (X_axis[2]-X_axis[1])/X_axis[0]));
        h2[0]->GetYaxis()->SetTitleOffset(1.25);
        for (int typei = 0; typei < ltypes; typei++)
        {
            h2[typei]->Draw("HIST, same");
        }
        can->Update();

        TLegend *legend = new TLegend(0.6, 0.8, 0.95, 0.95);
        for (int typei = 0; typei < ltypes; typei++)
            legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", LG[typei].Data(), num[typei]), "l");
        legend->SetBorderSize(0);
        legend->SetFillColorAlpha(0, 0);
        legend->Draw();

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextColor(1);
        text->SetTextSize(0.03);
        text->SetTextAlign(22);
        text->DrawLatex(0.8, 0.75, Form("%dGeV < P_{T}^{#gamma} < %dGeV", pcut[n_bin], pcut[n_bin+1]));
        text->DrawLatex(0.8, 0.7, Form("chiso SB: 5GeV < chiso > 10GeV"));

        can->Update();
        cout << "\n\tFinish Create TLegend" << endl;

        // can->SetGrid(0, 2);
        can->SaveAs(Form("%s/overlap/figures/template/%s_%d_%d_to_%d.pdf", str_dir.Data(), mode.Data(), isBarrel, pcut[n_bin], pcut[n_bin+1]));

        for (int typei = 0; typei < ltypes; typei++)
            h2[typei]->Reset("ICES");
        can->Clear();

        cout << Form("\tPath is %s", str_dir.Data()) << endl;
        cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/template/%s_%d_%d_to_%d.pdf", str_dir.Data(), mode.Data(), isBarrel, pcut[n_bin], pcut[n_bin+1])) << endl;
    }
}
