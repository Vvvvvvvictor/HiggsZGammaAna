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

void draw_binaff(TString mode = "inclusive")
{

    cout << "\tNOW!!! We will draw the most beautiful PICTURE in the world" << endl;
    cout << "\tIt worth coutless $$$$" << endl;
    cout << "\tLET'S GO!!!!!!" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;

    TString x_name = "gamma_mvaID";
    TString y_name = "H_mass"; //TODO: check it before run
    TString y_lg = "M_{ll#gamma}"; // TODO: check it before run
    double width = 10, zero = 100, min, max; // check it before run
    TString str_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml";
    TString bkg[] = {"data"};
    // TString sig[] = {"ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"};
    TString sig[] = {"zg"};

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
    can->SetBottomMargin(0.15);
    can->SetLogy(0);
    cout << "\n\tFinish Set Style and Canvas" << endl;

    Double_t sigw, bkgw, sigx, sigy, bkgx, bkgy, top;

    cout << "\n\tLoading " << Form("%s file from", mode.Data()) << Form("%s/overlap", str_dir.Data()) << "!......" << endl;

    TChain *sigt = new TChain(Form("%s", mode.Data()));
    TChain *bkgt = new TChain(Form("%s", mode.Data()));

    for (int i = 0; i < 1; i++)
        bkgt->Add(Form("%s/overlap/%s.root", str_dir.Data(), bkg[i].Data()));
    for (int i = 0; i < 1; i++)
        sigt->Add(Form("%s/overlap/%s.root", str_dir.Data(), sig[i].Data()));

    Long64_t signentries = sigt->GetEntries();
    sigt->SetBranchAddress(x_name.Data(), &sigx);
    sigt->SetBranchAddress("weight", &sigw);
    sigt->SetBranchAddress(y_name.Data(), &sigy);

    cout << "\n\tStart drawing\n\t Signal entries are : " << signentries << " ! oh my god ! too much !"
         << "\n\tloading...loading...zzzzzz\n"
         << endl;

    TH1D *sigh = new TH1D("Prompt photon", Form("%s; #gamma mvaID; Ralative value/0.025", mode.Data()), 80, -1, 1);
    sigh->SetLineColor(kRed);
    sigh->SetLineWidth(3);

    Long64_t bkgnentries = bkgt->GetEntries();
    bkgt->SetBranchAddress(x_name.Data(), &bkgx);
    bkgt->SetBranchAddress("weight", &bkgw);
    bkgt->SetBranchAddress(y_name.Data(), &bkgy);

    cout << "\n\tStart drawing\n\t Background entries are : " << bkgnentries << " ! oh my god ! too much !"
         << "\n\tloading...loading...zzzzzz\n"
         << endl;

    TH1D *bkgh = new TH1D("Non-prompt photon", Form("%s; #gamma mvaID; Ralative value/0.025", mode.Data()), 80, -1, 1);
    bkgh->SetLineColor(kBlack);
    bkgh->SetLineWidth(3);

    double sigsum, bkgsum, sigrat, bkgrat;
    for (int i = 0; i < 8; i++)
    {
        sigsum = 0; bkgsum = 0; top=0.;
        min = i * width + zero;
        max = min + width;
        // if (i == 7) max = 10000; //TODO: check it before run
        // else if (i == 0) min = -10000; //TODO: check it before run
        cout << "\n\tThis picture is from " << min << " to " << max << endl;
        cout << "\n\tNow is processing signal part!";

        for (Long64_t i = 0; i < signentries; i++)
        {
            sigt->GetEntry(i);
            sigsum += sigw;
            if (sigy > min && sigy < max)
                sigh->Fill(sigx, sigw);
        }
        cout << "\tFinished loading signal data" << endl;

        cout << "\n\tNow is processing background part!";

        for (Long64_t i = 0; i < bkgnentries; i++)
        {
            bkgt->GetEntry(i);
            bkgsum += bkgw;
            if (bkgy > min && bkgy < max)
                bkgh->Fill(bkgx, bkgw);
        }
        cout << "\n\tFinished loading signal data" << endl;

        bkgh->Add(sigh, -1);
        sigrat = sigh->Integral(54,78) / sigh->Integral();
        sigh->Scale(1./sigh->Integral(), "width");
        top = sigh->GetBinContent(sigh->GetMaximumBin());
        bkgrat = bkgh->Integral(54,78) / bkgh->Integral();
        bkgh->Scale(1./bkgh->Integral(), "width");
        if (top < bkgh->GetBinContent(bkgh->GetMaximumBin()))
            sigh->SetMaximum(bkgh->GetBinContent(bkgh->GetMaximumBin())*1.05);
        sigh->Draw("HIST");
        bkgh->Draw("HIST, same");
        bkgh->GetYaxis()->SetTitleOffset(1.2);

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextColor(1);
        text->SetTextSize(0.08);
        text->SetTextAlign(22);
        // if (i==7) text->DrawLatex(0.5, 0.5, Form("%s > %.0f", y_lg.Data(), min)); // TODO: check it before run
        // else if (i==0) text->DrawLatex(0.5, 0.5, Form("%s < %1.1f", y_lg.Data(), max)); // TODO: check it before run
        // else // TODO: check it before run
        text->DrawLatex(0.5, 0.5, Form("%.0f < %s < %.0f", min, y_lg.Data(), max));// TODO: check it before run

        TLegend *legend = new TLegend(0.6, 0.75, 0.9, 0.95);
        legend->AddEntry("Prompt photon", Form("Prompt photon(%2.1f%%)", 100*sigrat), "l");
        legend->AddEntry("Non-prompt photon", Form("Non-prompt photon(%2.2f%%)", 100*bkgrat), "l");
        legend->SetBorderSize(0);
        legend->SetFillColorAlpha(0, 0);
        legend->SetTextSize(0.03);
        legend->Draw();
        cout << "\n\tFinish Create TLegend" << endl;
        can->SaveAs(Form("%s/plots/bin_cut/%s_%s_%d.pdf", str_dir.Data(), mode.Data(), y_name.Data(), i));
        can->Clear();

        cout << Form("\tPath is %s", str_dir.Data()) << endl;
        cout << Form("\tWe get the PICTURE:%s!!!!!!!", Form("%s/plots/bin_cut/%s_%s_%d.pdf", str_dir.Data(), x_name.Data(), y_name.Data(), i)) << endl;

        sigh -> Reset("ICES");
        bkgh -> Reset("ICES");
    }
}
