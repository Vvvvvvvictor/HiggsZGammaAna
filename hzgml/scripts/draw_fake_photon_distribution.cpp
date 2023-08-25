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

void draw_fake_photon_distribution(TString mode = "inclusive")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    // TString LG[]{"Reference(data-true, pass ID)", "Fake(data-true, inversed ID)", "Fake(dy, pass ID)", "Fake(dy, inversed ID)"};
    // TString F[][6]= {{"data", "zg", "llajj", "ttg", "tg", "0"}, {"data", "zg", "llajj", "ttg", "tg", "0"}, {"dy", "0"}, {"dy", "0"}};
    TString LG[]{"True(dy, >0.4)", "Fake(dy, >0.4)", "Fake(dy, <-0.5)", "Fake(dy, -0.5~0.4)"};
    TString F[][6]= {{"dy", "0"}, {"dy", "0"}, {"dy", "0"}, {"dy", "0"}};
    double Y_axis[3] = {0, 0.02, 100};
    int ltypes = 4, sum_num[ltypes][5];
    Double_t num[ltypes];

    int lep_id = 11;
    TString lepton = "ele";

    // int lep_id = 13;
    // TString lepton = "mu";

    // int isbarrel = 1;
    // TString region = "Barrel";
    // int ptbins[] = {15,17,20,25,30,35,40,50,5000};
    // int num_ptbin = 8;
    
    int isbarrel = 0;
    TString region = "Endcap";
    int ptbins[] = {15,20,25,35,40,5000};
    int num_ptbin = 5;

    Double_t data;
    TString x_name = "Z_eta";
    double X_axis[3] = {35, -7, 7};
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    // TString ltype[] = {"chiso medium", "chiso side band"};

    for (int bini=0; bini<num_ptbin; bini++){
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
        can->SetRightMargin(0.1);
        can->SetBottomMargin(0.15);
        can->SetLogy(0);
        cout << "\n\tFinish Set Style and Canvas" << endl;

        Double_t w, minus;
        Long64_t n, islep;
        Double_t sel, chiso, pt, eta, isp, subl, mva;
        double max = 0.0;

        TH1D *h2[ltypes];

        for (int typei = 0; typei < ltypes; typei++)
        {
            // if (typei<2) continue;
            h2[typei] = new TH1D(Form("h2_%d", typei), Form("%s; %s:%s; rel. value", mode.Data(), mode.Data(), x_name.Data()), X_axis[0], X_axis[1], X_axis[2]);
            h2[typei]->SetLineColor(TColor::GetColorDark(typei + 3));
            h2[typei]->SetLineStyle(1);
            h2[typei]->SetLineWidth(3);

            num[typei] = 0.0;
            double sum = 0.0;

            int filei = 0;
            while (F[typei][filei] != "0"){
                sum_num[typei][filei] = 0;
                cout << "\n\tLoading file" << Form("%s/overlap/%s.root", str_dir.Data(), F[typei][filei].Data()) << "!......" << endl;
                TFile *f = new TFile(Form("%s/overlap/%s.root", str_dir.Data(), F[typei][filei].Data()));

                TTree *t = (TTree *)f->Get(mode.Data());
                Long64_t nentries = t->GetEntries();
                t->SetBranchAddress(x_name, &data);
                t->SetBranchAddress("weight", &w);
                // t->SetBranchAddress("n_jets", &n);
                t->SetBranchAddress("gamma_chiso", &chiso);
                t->SetBranchAddress("gamma_pt", &pt);
                t->SetBranchAddress("gamma_eta", &eta);
                t->SetBranchAddress("Z_sublead_lepton_pt", &subl);
                t->SetBranchAddress("gamma_mvaID", &mva);
                // if ((typei<2 & filei>0) | (typei>1))
                t->SetBranchAddress("gamma_genPartFlav", &isp);
                t->SetBranchAddress("Z_lead_lepton_id", &islep);
                // if ((typei<2) & (filei>0))
                //     minus = -1.0;
                // else
                minus = 1.0;

                cout << "\n\tStart drawing " << x_name.Data() << "\n\tentries are : " << nentries << " ! oh my god ! too much !"
                    << "\n\tloading...loading...zzzzzz\n"
                    << endl;

                for (Long64_t i = 0; i < nentries; i++)
                {
                    t->GetEntry(i);
                    sum+=w;
                    if (isbarrel==1){
                        if ((eta<-1.4442) | (eta>1.4442)) continue;
                    }
                    if (isbarrel==0){
                        if ((abs(eta)>2.5) | (abs(eta)<1.566)) continue;
                    }
                    if ((pt<ptbins[bini]) | (pt>ptbins[bini+1])) continue;
                    if (chiso*pt>5) continue;
                    // if (n!=1) continue;
                    if (subl<15) continue;
                    // if (data>0.012) continue;
                    if (abs(islep)!=lep_id) continue;
                    if ((typei==0) & (isp!=1)) continue;
                    if ((typei>0) & (isp!=0)) continue;
                    if ((typei<2) & ((isbarrel==1 & mva<0.42) | (isbarrel==0 & mva<0.14))) continue;
                    if ((typei==2) & (mva>-0.5)) continue;  
                    if ((typei==3) & (mva<-0.5 | mva>0.4)) continue;
                    w = minus*w;
                    // cout << i << endl;
                    h2[typei]->Fill(data, w);
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
        h2[0]->SetMaximum(max * 1.2);
        h2[0]->GetYaxis()->SetTitle(Form("Events/%.2g" , (X_axis[2]-X_axis[1])/X_axis[0]));
        h2[0]->GetYaxis()->SetTitleOffset(1.25);
        for (int typei = 0; typei < ltypes; typei++)
        {
            h2[typei]->Draw("HIST, same");
        }
        can->Update();

        TLegend *legend = new TLegend(0.55, 0.8, 0.90, 0.95);
        for (int typei = 0; typei < ltypes; typei++)
            legend->AddEntry(Form("h2_%d", typei), Form("%s(%.2f)", LG[typei].Data(), num[typei]), "l");
        legend->SetBorderSize(0);
        legend->SetFillColorAlpha(0, 0);
        legend->Draw();
        can->Update();
        cout << "\n\tFinish Create TLegend" << endl;

        // can->SetGrid(0, 2);
        can->SaveAs(Form("%s/overlap/figures/over/fake_%s_%s_%s_pt%d_pt%d.pdf", str_dir.Data(), x_name.Data(), region.Data(), lepton.Data(), ptbins[bini], ptbins[bini+1]));

        for (int typei = 0; typei < ltypes; typei++)
                h2[typei]->Reset("ICES");
        can->Clear();

        cout << Form("\tPath is %s", str_dir.Data()) << endl;
        cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/overlap/figures/over/fake_%s_%s_pt%d_pt%d.pdf", str_dir.Data(), x_name.Data(), region.Data(), ptbins[bini], ptbins[bini+1])) << endl;
    }
}
