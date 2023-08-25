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

void draw_variable(TString mode = "inclusive")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString LG[]={"Barrel", "Early Endcaps"};
    TString F[]= {"data", "zg", "0"};// "zg", "llajj", "ttg", "tg", "0"};
    double Y_axis[3] = {0, 0.02, 100};
    int ltypes = 2, sum_num[ltypes][5];
    Double_t num[ltypes];

    // int isbarrel = 1;
    // TString region = "Barrel";
    // int ptbins[] = {15,17,20,25,30,35,40,50,60,5000};
    double X_axis[3] = {100, 0, 10};
    // int num_ptbin = 9;
    
    // int isbarrel = 0;
    // TString region = "Endcap";
    // int ptbins[] = {15,20,25,35,50,5000};
    // double X_axis[3] = {50, 0, 2};
    // int num_ptbin = 5;

    Double_t data;
    TString x_name = "gamma_chiso";
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    // TString ltype[] = {"Barrel", "Early Endcaps"};

    // for (int bini=0; bini<num_ptbin; bini++){
        cout << "\n\tThere are " << ltypes << " kinds of shape to plot!!!" << endl;

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
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);
        can->SetLogy(1);
        cout << "\n\tFinish Set Style and Canvas" << endl;

        Double_t w, n, minus;
        Double_t sel, chiso, pt, eta, mva;
        double max = 0.0;

        TH1D *h2[ltypes];

        for (int typei = 0; typei < ltypes; typei++)
        {
            h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; %s:%s; rel. value", mode.Data(), mode.Data(), x_name.Data()), X_axis[0], X_axis[1], X_axis[2]);
            h2[typei]->SetLineColor(TColor::GetColorDark(typei + 3));
            h2[typei]->SetLineStyle(1);
            h2[typei]->SetLineWidth(3);

            num[typei] = 0.0;
            double sum = 0.0;

            int filei = 0;
            while (F[filei] != "0"){
                sum_num[typei][filei] = 0;
                cout << "\n\tLoading file" << Form("%s/overlap/%s.root", str_dir.Data(), F[filei].Data()) << "!......" << endl;
                TFile *f = new TFile(Form("%s/overlap/%s.root", str_dir.Data(), F[filei].Data()));

                TTree *t = (TTree *)f->Get(mode.Data());
                Long64_t nentries = t->GetEntries();
                t->SetBranchAddress(x_name, &data);
                t->SetBranchAddress("weight", &w);
                // t->SetBranchAddress("n_iso_photons", &n);
                // t->SetBranchAddress("photon_selection", &sel);
                t->SetBranchAddress("gamma_mvaID", &mva);
                t->SetBranchAddress("gamma_eta", &eta);
                t->SetBranchAddress("gamma_pt", &pt);
                // t->SetBranchAddress("photon_is_barrel", &isb);
                // t->SetBranchAddress("gamma_energyErr", &eerr);
                if (filei)
                    minus = -1.0;
                else
                    minus = 1.0;

                cout << "\n\tStart drawing " << x_name.Data() << "\n\tentries are : " << nentries << " ! oh my god ! too much !"
                    << "\n\tloading...loading...zzzzzz\n"
                    << endl;

                for (Long64_t i = 0; i < nentries; i++)
                {
                    t->GetEntry(i);
                    sum+=w;
                    // if ((pt<ptbins[bini]) | (pt>ptbins[bini+1]) | (isb!=isbarrel) | (match!=1)) continue;
                    // if (eerr/pt>0.48) continue;
                    if ((typei==0) & (abs(eta)<1.5)) continue;
                    if ((typei==1) & (abs(eta)>1.5)) continue;
                    if (mva>-0.97 | mva<-1) continue;
                    w = minus*w;
                    h2[typei]->Fill(data*pt, w);
                    num[typei] += w;
                    sum_num[typei][filei]++;
                }
                cout << "\tFinished loading File" << endl;
                cout << "\n\t" << sum_num[typei][filei] << "/" << nentries  << " items are in it!!";
                cout << "\n\t" << num[typei] << "/" << sum  << " items are in it!!";

                filei++;
            }
            h2[typei]->Scale(1. / h2[typei]->Integral(), "width");
            if (max < h2[typei]->GetBinContent(h2[typei]->GetMaximumBin()))
                max = h2[typei]->GetBinContent(h2[typei]->GetMaximumBin());
        }
        h2[0]->SetMaximum(max * 1.1);
        h2[0]->GetYaxis()->SetTitle(Form("Rel. Value/%.2g" , (X_axis[2]-X_axis[1])/X_axis[0]));
        h2[0]->GetXaxis()->SetTitle("ISO_{ch}^{#gamma}");
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
        can->Update();
        cout << "\n\tFinish Create TLegend" << endl;

        // can->SetGrid(0, 2);
        can->SaveAs(Form("%s/overlap/figures/over/chiso.pdf", str_dir.Data()));

        for (int typei = 0; typei < ltypes; typei++)
                h2[typei]->Reset("ICES");
        can->Clear();

        cout << Form("\tPath is %s", str_dir.Data()) << endl;
        cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/over/energyErr_pt.pdf", str_dir.Data())) << endl;
    // }
}
