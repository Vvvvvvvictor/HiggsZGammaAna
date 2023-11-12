//#include "setTDRStyle.C"
#include "TMath.h"
//#include "TROOT.h"
#include <TH1D.h>
#include <TH1D.h>
#include <TH1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeafF.h>
#include <TChain.h>
#include <TFile.h>
#include "TSystem.h"
#include <TChain.h>
#include "TSystem.h"
#include <TString.h>
#include <iostream>
#include <vector>
#include <TPostScript.h>
#include <iostream>
#include <iomanip>  //for precision

//==============
const int debug=1;
//const TString TreeName = "ZmmgTree";
const TString TreeName = "passedEvents";//bing
//const TString PrintInfor="CMS Preliminary   #sqrt{s} = 13 TeV   L = 2.7 fb^{-1}";
//const TString PrintInfor="CMS Preliminary                      #sqrt{s} = 13 TeV";
const TString PrintInfor1="#bf{CMS} #it{} #it{Preliminary}";
//const TString PrintInfor2="12.9 fb^{-1} (13TeV)";
//const TString PrintInfor2="35.9 fb^{-1} (13TeV)"; //"2017 data (13TeV)"; //"35.9 fb^{-1} (13TeV)";
//const TString PrintInfor2="41.5 fb^{-1} (13TeV)";
const TString PrintInfor2="138 fb^{-1} (13TeV)";
//const double MCXSweight =  6025.2*2136.361*1.0/28747969;
// const double MCXSweight =  6025.2*2488.*1.0/28747969;
//const double MCXSweight =  6025.2*2534.565*1.0/28747969;

const TString PrintEtaEB="|#eta|<1.4442";
const TString PrintEtaEE="|#eta|>1.566";
///////////////////////////////////////

void DrawMyPlots(string Object, string Selections,  string XTitle, string YUnit, string PlotName, int BinTotal, double BinXLow, double BinXHig, double mySys, int LegendLR, int IfLogY, int IfEB=-1){


  TCanvas *c1 = new TCanvas("reconstruction1","reconstruction1");
  c1->SetFillColor(0);

  TLegend *legend;
  if(LegendLR==0) legend = new TLegend(0.69,0.63,0.89,0.9);
  // else legend = new TLegend(0.25,0.63,0.45,0.9);
  else legend = new TLegend(0.2,0.63,0.4,0.9);
   legend->SetFillColor(0);

  //====add root file
  //TChain *Data_Tree=new TChain(TreeName);
  //TChain *MC_Tree=new TChain(TreeName);

  TChain *Data_Tree=new TChain("passedEvents");
  TChain *MC_Tree=new TChain("passedEvents");

  //====================
  //Data_Tree->Add("./data2016_all_all.root");
  //MC_Tree->Add("./DY_all_all.root");

  Data_Tree->Add("/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/run2/ALP_data.root");//bing
  MC_Tree->Add("/publicfs/cms/user/wangzebing/ALP/Analysis_out/UL/run2/ALP_DYJetsToLL.root");//bing

  //Data_Tree->Add("/publicfs/cms/user/wangzebing/ALP/Analysis_out/18/massInde/ALP_data.root");//bing
  //MC_Tree->Add("/publicfs/cms/user/wangzebing/ALP/Analysis_out/18/massInde/ALP_DYJetsToLL.root");//bing

  //=========entries================
  int entries_Data = Data_Tree->GetEntries();
  if(debug==1) cout <<"JTao: nEntries_Data = "<<entries_Data<<endl;

  int entries_MC = MC_Tree->GetEntries();
  if(debug==1) cout <<"JTao: nEntries_MC = "<<entries_MC<<endl;

  c1->cd();

  //=================
  char *myLimits= new char[100];
  sprintf(myLimits,"(%d,%f,%f)",BinTotal,BinXLow,BinXHig);
  TString Taolimits(myLimits);

  //====data=======
  TString variable_Data = Object + ">>Histo_Data_temp" + Taolimits;
  Data_Tree->Draw(variable_Data, Selections.c_str());
  TH1D *h_data = (TH1D*)gDirectory->Get("Histo_Data_temp");
  h_data->SetTitle("");
  c1->Clear();


  double Ntot_Data=h_data->Integral();
  if( debug==1 ) cout<<"JTao: N_Data= "<<Ntot_Data<<endl;

  Double_t scale_Data = 1.0/Ntot_Data;
  h_data->Sumw2();
  //h_data->Scale(scale_Data);

  TString variable_MC = Object + ">>Histo_MC_temp" + Taolimits;
  MC_Tree->Draw(variable_MC, Selections.c_str());
  TH1D *h_MC = (TH1D*)gDirectory->Get("Histo_MC_temp");
  c1->Clear();
  double Ntot_MC=h_MC->Integral();
  if( debug==1 ) cout<<"JTao: N_MC= "<<Ntot_MC<<endl;
  Double_t scale_MC = Ntot_Data*1.0/Ntot_MC;
  h_MC->Sumw2();
  h_MC->Scale(scale_MC);
  //h_MC->Scale(MCXSweight);
  double Ntot_MCXSW=h_MC->Integral();
  if( debug==1 ) cout<<"JTao: Weighted N_MC= "<<Ntot_MCXSW<<endl;

  //======

  TString ObjectPlusSys = Form("%f*", 1.0+mySys) + Object;
  //TString ObjectPlusSys = "pho_MVAOutputCorr_Down1Sigma";

  TString variable_MCP = ObjectPlusSys + ">>Histo_MCP_temp" + Taolimits;

  MC_Tree->Draw(variable_MCP, Selections.c_str());
  TH1D *h_MCP = (TH1D*)gDirectory->Get("Histo_MCP_temp");
  c1->Clear();
  h_MCP->Scale(scale_MC);

  TString ObjectMinusSys = Form("%f*", 1.0-mySys) + Object;
  //TString ObjectMinusSys = "pho_MVAOutputCorr_Up1Sigma";
  TString variable_MCM = ObjectMinusSys + ">>Histo_MCM_temp" + Taolimits;
  MC_Tree->Draw(variable_MCM, Selections.c_str());
  TH1D *h_MCM = (TH1D*)gDirectory->Get("Histo_MCM_temp");
  c1->Clear();
  h_MCM->Scale(scale_MC);

  //======
  //if( debug==1 ) cout<<"JTao: Weighted N_MC= "<<Ntot_MCXSW<<endl;
  double YMean[BinTotal], YMeanNorm[BinTotal], XMean[BinTotal];
  double X_ErrH[BinTotal], X_ErrL[BinTotal];
  double MCSys_ErrH[BinTotal], MCSys_ErrL[BinTotal];
  double MCSys_RelErrH[BinTotal], MCSys_RelErrL[BinTotal];
  TH1D *h_MCUP = new TH1D("h_MCUP","",BinTotal,BinXLow,BinXHig);
  TH1D *h_MCDN = new TH1D("h_MCDN","",BinTotal,BinXLow,BinXHig);
  TH1D *h_MCUPNorm = new TH1D("h_MCUPNorm","",BinTotal,BinXLow,BinXHig);
  TH1D *h_MCDNNorm = new TH1D("h_MCDNNorm","",BinTotal,BinXLow,BinXHig);
  double WidthBin=(BinXHig-BinXLow)/BinTotal;
  double Chi2=0.;
  for(int ibin=0; ibin<BinTotal; ibin++){
    double Nd = h_data->GetBinContent(ibin+1);
    double Nm = h_MC->GetBinContent(ibin+1);
    double NmErr = h_MC->GetBinError(ibin+1);

    Chi2 += fabs(NmErr)>1e-9?(Nm-Nd)*(Nm-Nd)*1.0/(NmErr*NmErr):0.0;

    //===
    YMean[ibin] = Nm;  YMeanNorm[ibin] = 1.0; XMean[ibin] = BinXLow+ (ibin+0.5)*WidthBin;
    X_ErrH[ibin] = 0.5*WidthBin;   X_ErrL[ibin] = 0.5*WidthBin;
    double dNmp = h_MCP->GetBinContent(ibin+1) - Nm;
    double dNmm = h_MCM->GetBinContent(ibin+1) - Nm;
    double dplus= dNmp>dNmm?dNmp:dNmm; if(dplus<0) dplus *= -1;  //dplus = 0.;
    double dminus=dNmp<dNmm?dNmp:dNmm; if(dminus>0) dminus *= -1; dminus *= -1; //dminus = 0;  dminus *= -1;
    MCSys_ErrH[ibin] = dplus; MCSys_ErrL[ibin] = dminus;
    h_MCUP->SetBinContent(ibin+1, Nm+dplus);
    h_MCDN->SetBinContent(ibin+1, Nm-dminus);
    MCSys_RelErrH[ibin] = Nm>0?dplus/Nm:0.;
    MCSys_RelErrL[ibin] = Nm>0?dminus/Nm:0.;
    h_MCUPNorm->SetBinContent(ibin+1, 1+MCSys_RelErrH[ibin]);
    h_MCDNNorm->SetBinContent(ibin+1, 1-MCSys_RelErrL[ibin]);
  }
  TGraphAsymmErrors * gr_MCSys = new TGraphAsymmErrors(BinTotal, XMean, YMean, X_ErrL, X_ErrH, MCSys_ErrL, MCSys_ErrH);
  TGraphAsymmErrors * gr_MCSysNorm = new TGraphAsymmErrors(BinTotal, XMean, YMeanNorm, X_ErrL, X_ErrH, MCSys_RelErrL, MCSys_RelErrH);
  gr_MCSys->SetFillColor(kRed-10); //2
  gr_MCSys->SetFillStyle(3001); //3001
  gr_MCSysNorm->SetFillColor(kRed-10); //2
  gr_MCSysNorm->SetFillStyle(3001);

  h_MCUP->SetLineColor(2);
  h_MCUP->SetLineStyle(1);
  h_MCUP->SetLineWidth(2);
  h_MCDN->SetLineColor(2);
  h_MCDN->SetLineStyle(1);
  h_MCDN->SetLineWidth(2);
  h_MCUPNorm->SetLineColor(2);
  h_MCUPNorm->SetLineStyle(1);
  h_MCUPNorm->SetLineWidth(2);
  h_MCDNNorm->SetLineColor(2);
  h_MCDNNorm->SetLineStyle(1);
  h_MCDNNorm->SetLineWidth(2);
  cout<<"JTao: chi2 = "<<Chi2<<endl;

  double maxY=max(h_data->GetMaximum(),h_MC->GetMaximum());
  maxY = max(maxY, h_MCUP->GetMaximum());
  maxY = max(maxY, h_MCDN->GetMaximum());
  double minY=min(h_data->GetMinimum(),h_MC->GetMinimum());
  h_data->GetYaxis()->SetRangeUser(0.95*minY, 1.05*maxY);
  minY = min(minY, h_MCUP->GetMinimum());
  minY = min(minY, h_MCDN->GetMinimum());
  //h_data->SetMaximum(1.2*maxY);
  if(IfLogY==1) h_data->GetYaxis()->SetRangeUser(1.0, 1.5*maxY);

  h_data->SetLineColor(1);
  h_data->SetFillColor(0);
  h_data->SetLineStyle(1);
  h_data->SetLineWidth(2);
  //h_data->GetXaxis()->SetTitle(XTitle.c_str());
  //double WidthBin=(BinXHig-BinXLow)/BinTotal;
  //TString TitleY( Form("A.U. / %.2g GeV",WidthBin) );
  //TString TitleY( Form("No. of Entries in data / %.2g GeV",WidthBin) );
  //TString TitleY = "A.U";
  //string PreTitleY( Form("No. of Entries / %.2g ",WidthBin) );
  string PreTitleY( Form("Events / %.2g ",WidthBin) );
  string TitleY =  PreTitleY + YUnit;
  h_data->GetYaxis()->SetTitle(TitleY.c_str());

  //h_data->SetTitleSize(0.05,"X");
  h_data->SetTitleSize(0.05,"X");
  h_data->SetTitleSize(0.05,"Y");
  //h_data->SetTitleOffset(1.3, "Y");
  h_data->SetTitleOffset(1.1, "Y");

  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerSize(1.0);
  //h_data->SetMarkerSize(1.2);
  h_data->SetMarkerStyle(20);

  h_MC->SetFillColor(8);
  h_MC->SetMarkerStyle(0);
  h_MC->SetLineColor(8);
  h_MC->SetLineStyle(1);
  h_MC->SetLineWidth(2);

  //legend->AddEntry(h_data,"data","pe");
  legend->AddEntry(h_data,"Data","pe");
  //legend->AddEntry(h_MC,"MC","f");
  //legend->AddEntry(h_MC,"Zll#gamma","f");
  //legend->AddEntry(h_MC,"Z#rightarrow#mu#mu#gamma","f");
  legend->AddEntry(h_MC,"DYJets","f");//bing
  legend->AddEntry(gr_MCSys,"MC syst.","f");

  //prepare 2 pads
  const Int_t nx=1;
  const Int_t ny=2;
  const Double_t x1[2] = {0.0,0.0};
  const Double_t x2[2] = {1.0,1.0};
  //const Double_t y1[] = {1.0,0.3};
  //const Double_t y2[] = {0.3,0.00};
  const Double_t y1[2] = {0.3,0.0};
  const Double_t y2[2] = {1.0,0.3};
  Double_t psize[2];
  TPad *pad;
  const char *myname = "c";
  char *name2 = new char [strlen(myname)+6];
  Int_t n = 0;
  for (int iy=0;iy<ny;iy++) {
    for (int ix=0;ix<nx;ix++) {
      n++;
      sprintf(name2,"%s_%d",myname,n);
      if(ix==0){
        gStyle->SetPadLeftMargin(.166);
      }else{
        gStyle->SetPadLeftMargin(.002);
        gStyle->SetPadTopMargin(.002);
      }

      if(iy==0){//upper
        gStyle->SetPadTopMargin(0.05*(1./0.7)); // 0.05
        gStyle->SetPadBottomMargin(.02);
      }
      if(iy==(ny-1)){//lower pad
        gStyle->SetPadTopMargin(.05);
        //gStyle->SetPadBottomMargin(.13*(1./0.3));
        gStyle->SetPadBottomMargin(.40);


      }
      pad = new TPad(name2,name2,x1[ix],y1[iy],x2[ix],y2[iy]);
      pad->SetNumber(n);
      pad->Draw();
      psize[iy]=y1[iy]-y2[iy];
      //if(iy>0 )pad->SetGrid(kTRUE);
    }// end of loop over x
  }// end of loop over y
  delete [] name2;

  //===Drawing====
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);

  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  h_data->GetXaxis()->SetLabelColor(0);
  h_data->SetNdivisions(510 ,"X");

  c1->cd(1);
  gPad->SetLogy(IfLogY);
  //=========
  h_data->Draw("PE1");
//  h_data->Draw("PEX0");
  h_MC->Draw("hist,same");
  gr_MCSys->Draw("same2");
  h_MCUP->Draw("same");
  h_MCDN->Draw("same");
  legend->Draw("same");
  h_data->Draw("samePE1");
  h_data->Draw("Axissame");

  /*
  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.05);
  a.DrawLatex(0.2,0.94, PrintInfor);
  */
  TLatex * tex = new TLatex(0.16,0.94, PrintInfor1);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.70,0.94, PrintInfor2);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  //TLatex * texEta = new TLatex(0.72,0.50, PrintEtaEB);
  //if(IfEB==0) texEta = new TLatex(0.72,0.50, PrintEtaEE);

  /*
  TLatex *texEta = new TLatex(0.22,0.5, PrintEtaEB);
  if(IfEB==0) texEta = new TLatex(0.22,0.5, PrintEtaEE);
  if(LegendLR==0) {
   texEta = new TLatex(0.72,0.5, PrintEtaEB);
   if(IfEB==0) texEta = new TLatex(0.72,0.5, PrintEtaEE);
  }


  texEta->SetNDC();
  texEta->SetTextFont(42);
  texEta->SetTextSize(0.045);
  texEta->SetLineWidth(2);
  if(IfEB>-1) texEta->Draw();
  */

  ///====
  TLine *Line1 = new TLine(h_data->GetBinLowEdge(1),1,h_data->GetBinLowEdge(h_data->GetNbinsX())+ h_data->GetBinWidth(h_data->GetNbinsX()),1);
  Line1->SetLineColor(1);
  Line1->SetLineWidth(2);
  Line1->SetLineStyle(4);

  TH1D *histoRatio = new TH1D(*h_data);
  histoRatio->Divide(h_data, h_MC, 1., 1.);
  histoRatio->SetLineColor(1);
  histoRatio->SetLineStyle(1);
  histoRatio->SetLineWidth(2);
  histoRatio->SetMarkerColor(1);
  histoRatio->SetMarkerStyle(20);

  c1->cd(2);
  gPad->SetLogy(0);
  histoRatio->SetTitleOffset(1,"X");
  histoRatio->SetTitleSize(0.12,"X");
  histoRatio->SetLabelSize(0.1,"X");
  histoRatio->GetXaxis()->SetTitle(XTitle.c_str());
  //histoRatio->GetYaxis()->SetTitle("data/MC");
  histoRatio->GetYaxis()->SetTitle("Data / MC");
  //histoRatio->SetTitleOffset(0.5,"Y");
  histoRatio->SetTitleOffset(0.4,"Y");
  //histoRatio->SetTitleSize(0.12,"Y");
  histoRatio->SetTitleSize(0.14,"Y");
  histoRatio->SetLabelSize(0.1,"Y");
  histoRatio->SetLabelColor(1,"X");
  histoRatio->GetYaxis()->CenterTitle();
  //histoRatio->SetNdivisions(505 ,"Y");
  histoRatio->SetNdivisions(510 ,"X");

  gPad->SetTickx(1);
  gPad->SetTicky(1);

  histoRatio->GetXaxis()->SetTickLength(0.08);
  histoRatio->GetYaxis()->SetTickLength(0.06);
  histoRatio->GetYaxis()->SetNdivisions(503);

  //histoRatio->SetMinimum(0.); //0.);
  //histoRatio->SetMaximum(3.); //3.0);
  //if(IfEB==0) histoRatio->SetMaximum(5.);
//  histoRatio->GetYaxis()->SetRangeUser(0., 3.);
  histoRatio->GetYaxis()->SetRangeUser(0.2, 3.);
  //if(IfEB==0) histoRatio->GetYaxis()->SetRangeUser(0., 2.);
  histoRatio->Draw();
  gr_MCSysNorm->Draw("same2");
  h_MCUPNorm->Draw("same");
  h_MCDNNorm->Draw("same");
  histoRatio->Draw("same");
  Line1->Draw("same");

 //===================================
  //string nameplots="DataMCsys_"+PlotName+"_75.png";
  string nameplots="./plots_sys/DataMCsys_"+PlotName+".png";
  //string nameplots="DataMCsys_"+PlotName+"_range.pdf";
  c1->Print(nameplots.c_str());

  string nameplotspdf="./plots_sys/DataMCsys_"+PlotName+".pdf";
  c1->Print(nameplotspdf.c_str());
  //c1->Clear();
  //legend->Clear();

}



void DrawSys(){

  //gROOT->ProcessLine(".x hggPaperStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);

  //DrawMyPlots(string var, string selections,  string xtitle, string plotname, int nbin, double xlow, double xhig, double mySys, int legendLR, int IfLogY)
  cout<<endl;
  cout<<"===========================Preselected========================"<<endl;
  /*
  string PrePreselections="event_PassHLT_DiMu>0 && ( (fabs(pho_sceta) < 1.5 && pho_full5x5_r9Corr > 0.85 && pho_full5x5_r9Corr > 0.5) || (fabs(pho_sceta) < 1.5  && pho_full5x5_r9Corr < 0.85 && pho_full5x5_r9Corr > 0.5 && pho_full5x5_sieieCorr < 0.015 && pho_pfPhoIso03Corr_effA < 4. && pho_TrackIsoMuCorr < 6.) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9Corr > 0.9 && pho_full5x5_r9Corr > 0.85) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9Corr < 0.9 && pho_full5x5_r9Corr > 0.8 && pho_full5x5_sieieCorr < 0.035 &&  pho_pfPhoIso03Corr_effA < 4. && pho_TrackIsoMuCorr < 6.) ) && pho_hoe < 0.08 && pho_MVAOutputCorr>-0.9 && pho_PassEleVeto>0";

  string PrePreselections25="event_PassHLT_DiMu>0 && ( (fabs(pho_sceta) < 1.5 && pho_full5x5_r9Corr > 0.85 && pho_full5x5_r9Corr > 0.5) || (fabs(pho_sceta) < 1.5  && pho_full5x5_r9Corr < 0.85 && pho_full5x5_r9Corr > 0.5 && pho_full5x5_sieieCorr < 0.015 && pho_pfPhoIso03Corr_effA < 4. && pho_TrackIsoMuCorr < 6.) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9Corr > 0.9 && pho_full5x5_r9Corr > 0.85) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9Corr < 0.9 && pho_full5x5_r9Corr > 0.8 && pho_full5x5_sieieCorr < 0.035 &&  pho_pfPhoIso03Corr_effA < 4. && pho_TrackIsoMuCorr < 6.) ) && pho_hoe < 0.08 && pho_MVAOutputCorr>-0.9 && pho_PassEleVeto>0 && pho_pt>25.";


//  string PrePreselections="event_PassHLT_DiMu>0 && ( (fabs(pho_sceta) < 1.5 && pho_full5x5_r9 > 0.85 && pho_full5x5_r9 > 0.5) || (fabs(pho_sceta) < 1.5  && pho_full5x5_r9 < 0.85 && pho_full5x5_r9 > 0.5 && pho_full5x5_sieie < 0.015 && pho_pfPhoIso03_effA < 4. && pho_TrackIsoMuCorr < 6.) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9 > 0.9 && pho_full5x5_r9 > 0.85) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9 < 0.9 && pho_full5x5_r9 > 0.8 && pho_full5x5_sieie < 0.035 && pho_pfPhoIso03_effA < 4. && pho_TrackIsoMuCorr < 6.) ) && pho_hoe < 0.08 && pho_MVAOutput>-0.9 && pho_PassEleVeto>0 ";

//  string PrePreselections25="event_PassHLT_DiMu>0 && ( (fabs(pho_sceta) < 1.5 && pho_full5x5_r9 > 0.85 && pho_full5x5_r9 > 0.5) || (fabs(pho_sceta) < 1.5  && pho_full5x5_r9 < 0.85 && pho_full5x5_r9 > 0.5 && pho_full5x5_sieie < 0.015 && pho_pfPhoIso03_effA < 4. && pho_TrackIsoMuCorr < 6.) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9 > 0.9 && pho_full5x5_r9 > 0.85) || (fabs(pho_sceta) > 1.5  && pho_full5x5_r9 < 0.9 && pho_full5x5_r9 > 0.8 && pho_full5x5_sieie < 0.035 && pho_pfPhoIso03_effA < 4. && pho_TrackIsoMuCorr < 6.) ) && pho_hoe < 0.08 && pho_MVAOutput>-0.9 && pho_PassEleVeto>0 && pho_pt>25.";
//
  string MorePreselections = " && ( pho_full5x5_r9Corr > 0.8 || pho_egChargedHadronIso < 20. || pho_egChargedHadronIso/pho_pt<0.3 ) ";
//  string MorePreselections = " && ( pho_full5x5_r9 > 0.8 || pho_egChargedHadronIso < 20. || pho_egChargedHadronIso/pho_pt<0.3 ) ";
string Preselections= PrePreselections + MorePreselections;
string Preselections25= PrePreselections25 + MorePreselections;
//
  cout<<"===========================Preselected: EB========================"<<endl;
  string  EBSelection = " fabs(pho_sceta) < 1.5 ";
  string FinalEBSelections = "PUWeight*XSWeight*MuIDWeight*(" + EBSelection + " && " + Preselections + ")";
  //string FinalEBSelections = "RhoWeight*PUWeight*XSWeight*MuIDWeight*(" + EBSelection + " && " + Preselections + ")";
  string FinalEBSelections25 = "PUWeight*XSWeight*MuIDWeight*(" + EBSelection + " && " + Preselections25 + ")";
  //string FinalEBSelections25 = "RhoWeight*PUWeight*XSWeight*MuIDWeight*(" + EBSelection + " && " + Preselections25 + ")";

  DrawMyPlots("pho_sigEOverERegSmear",FinalEBSelections,"Photon #sigma_{E}/E", "", "pho_sigEOverERegSmear_PreEB", 40, 0, 0.08, 0.05, 0, 0, 1);
  DrawMyPlots("pho_sigEOverERegSmear",FinalEBSelections25,"Photon #sigma_{E}/E", "", "pho_sigEOverERegSmear_PreEB_PT25", 40, 0, 0.08, 0.05, 0, 0, 1);

  DrawMyPlots("pho_sigEOverEReg",FinalEBSelections,"Photon #sigma_{E}/E", "", "pho_sigEOverEReg_PreEB", 40, 0, 0.08, 0.05, 0, 0, 1);
  DrawMyPlots("pho_sigEOverEReg",FinalEBSelections25,"Photon #sigma_{E}/E", "", "pho_sigEOverEReg_PreEB_PT25", 40, 0, 0.08, 0.05, 0, 0, 1);

  //=== EE ===
  cout<<endl;
  cout<<"===========================Preselected: EE========================"<<endl;
  string  EESelection = " fabs(pho_sceta) > 1.5 ";
  string FinalEESelections = "PUWeight*XSWeight*MuIDWeight*(" + EESelection + " && " + Preselections + ")";
  //string FinalEESelections = "RhoWeight*PUWeight*XSWeight*MuIDWeight*(" + EESelection + " && " + Preselections + ")";
  string FinalEESelections25 = "PUWeight*XSWeight*MuIDWeight*(" + EESelection + " && " + Preselections25 + ")";
  //string FinalEESelections25 = "RhoWeight*PUWeight*XSWeight*MuIDWeight*(" + EESelection + " && " + Preselections25 + ")";

  DrawMyPlots("pho_sigEOverERegSmear",FinalEESelections,"Photon #sigma_{E}/E", "", "pho_sigEOverERegSmear_PreEE", 40, 0., 0.08, 0.05, 1, 0, 0);
  DrawMyPlots("pho_sigEOverERegSmear",FinalEESelections25,"Photon #sigma_{E}/E", "", "pho_sigEOverERegSmear_PreEE_PT25", 40, 0., 0.08, 0.05, 1, 0, 0);

  DrawMyPlots("pho_sigEOverEReg",FinalEESelections,"Photon #sigma_{E}/E", "", "pho_sigEOverEReg_PreEE", 40, 0., 0.08, 0.05, 1, 0, 0);
  DrawMyPlots("pho_sigEOverEReg",FinalEESelections25, "Photon #sigma_{E}/E", "", "pho_sigEOverEReg_PreEE_PT25", 40, 0., 0.08, 0.05, 1, 0, 0);
  */

  //string selection = "factor*(passChaHadIso && passNeuHadIso && passHOverE && passdR_gl && (!passH_m) && passBDT)";
  //DrawMyPlots("Z_m",selection, "Z_m", "", "Z_m", 100, 50., 130., 0.05, 1, 0, 0);
  //DrawMyPlots("ALP_m",selection, "ALP_m", "", "ALP_m_2017", 50, 0., 45., 0.1, 1, 0, 0);

  //DrawMyPlots("pho1Pt",selection, "pho1Pt", "", "pho1Pt_2017", 50, 5., 50., 0.1, 0, 0, 0);
  //DrawMyPlots("pho1eta",selection, "pho1eta", "", "pho1eta_2017", 50, -3., 3., 0.1, 1, 0, 0);
  //DrawMyPlots("pho1phi",selection, "pho1phi", "", "pho1phi_2017", 50, -4., 4., 0.1, 1, 0, 0);
  //DrawMyPlots("pho1R9",selection, "pho1R9", "", "pho1R9_2017", 50, 0., 1., 0.1, 1, 0, 0);
  //DrawMyPlots("pho1IetaIeta55",selection, "pho1IetaIeta55", "", "pho1IetaIeta55_2017", 50, 0., 0.07, 0.1, 0, 0, 0);

  //DrawMyPlots("pho2Pt",selection, "pho2Pt", "", "pho2Pt_2017", 50, 5., 50., 0.1, 0, 0, 0);
  //DrawMyPlots("pho2eta",selection, "pho2eta", "", "pho2eta_2017", 50, -3., 3., 0.1, 1, 0, 0);
  //DrawMyPlots("pho2phi",selection, "pho2phi", "", "pho2phi_2017", 50, -4., 4., 0.1, 1, 0, 0);
  //DrawMyPlots("pho2R9",selection, "pho2R9", "", "pho2R9_2017", 50, 0., 1., 0.1, 1, 0, 0);
  //DrawMyPlots("pho2IetaIeta55",selection, "pho2IetaIeta55", "", "pho2IetaIeta55_2017", 50, 0., 0.07, 0.1, 0, 0, 0);

  //DrawMyPlots("ALP_calculatedPhotonIso",selection, "ALP_calculatedPhotonIso", "", "ALP_calculatedPhotonIso_2017", 50, 0., 5., 0.1, 0, 1, 0);
  //DrawMyPlots("var_dR_g1Z",selection, "var_dR_g1Z", "", "var_dR_g1Z_2017", 50, 0., 8., 0.1, 0, 0, 0);
  //DrawMyPlots("var_Pta",selection, "var_Pta", "", "var_Pta_2017", 50, 0., 60., 0.1, 1, 0, 0);
  //DrawMyPlots("var_MhMZ",selection, "var_MhMZ", "", "var_MhMZ_2017", 50, 120., 600., 0.1, 0, 0, 0);
  //DrawMyPlots("H_pt",selection, "H_pt", "", "H_pt_2017", 50, 0., 200., 0.1, 0, 0, 0);

  string selection = "factor*pho1SFs*pho2SFs*(passChaHadIso && (!passNeuHadIso) && passHOverE && passdR_gl && H_m>95 && H_m<180)";
  //DrawMyPlots("Val_BDT",selection, "Val_BDT", "", "Val_BDT_2017", 50, 0., 1.0, 0.05, 0, 0, 0);
  //DrawMyPlots("pho2R9",selection, "pho2R9", "", "pho2R9_2018UL", 50, 0., 1.1, 0.01, 1, 0, 0);
  //DrawMyPlots("pho1R9",selection, "pho1R9", "", "pho1R9_UL", 50, 0., 1.1, 0.005, 1, 0, 0);
  //DrawMyPlots("pho1IetaIeta55",selection, "pho1IetaIeta55", "", "pho1IetaIeta55_UL", 50, 0., 0.07, 0.01, 0, 0, 0);
  //DrawMyPlots("pho1PIso_noCorr",selection, "pho1PIso_noCorr", "", "pho1PIso_noCorr_UL", 50, 0., 40.0, 0.01, 0, 1, 0);

  //DrawMyPlots("pho2R9",selection, "pho2R9", "", "pho2R9_UL", 50, 0., 1.1, 0.005, 1, 0, 0);
  //DrawMyPlots("pho2IetaIeta55",selection, "pho2IetaIeta55", "", "pho2IetaIeta55_UL", 50, 0., 0.07, 0.02, 0, 0, 0);
  //DrawMyPlots("pho2PIso_noCorr",selection, "pho2PIso_noCorr", "", "pho2Iso_noCorr_UL", 50, 0., 40.0, 0.01, 0, 1, 0);

  DrawMyPlots("ALP_calculatedPhotonIso",selection, "ALP_calculatedPhotonIso", "", "ALP_calculatedPhotonIso_UL", 50, 0., 5, 0.05, 0, 1, 0);

  // #############2018
  //DrawMyPlots("pho1Pt",selection, "pho1Pt", "", "pho1Pt_2018", 50, 5., 50., 0.1, 0, 0, 0);
  //DrawMyPlots("pho1eta",selection, "pho1eta", "", "pho1eta_2018", 50, -3., 3., 0.1, 1, 0, 0);
  //DrawMyPlots("pho1phi",selection, "pho1phi", "", "pho1phi_2018", 50, -4., 4., 0.1, 1, 0, 0);
  //DrawMyPlots("pho1R9",selection, "pho1R9", "", "pho1R9_2018", 50, 0., 1., 0.1, 1, 0, 0);
  //DrawMyPlots("pho1IetaIeta55",selection, "pho1IetaIeta55", "", "pho1IetaIeta55_2018", 50, 0., 0.07, 0.1, 0, 0, 0);

  //DrawMyPlots("pho2Pt",selection, "pho2Pt", "", "pho2Pt_2018", 50, 5., 50., 0.1, 0, 0, 0);
  //DrawMyPlots("pho2eta",selection, "pho2eta", "", "pho2eta_2018", 50, -3., 3., 0.1, 1, 0, 0);
  //DrawMyPlots("pho2phi",selection, "pho2phi", "", "pho2phi_2018", 50, -4., 4., 0.1, 1, 0, 0);
  //DrawMyPlots("pho2R9",selection, "pho2R9", "", "pho2R9_2018", 50, 0., 1., 0.1, 1, 0, 0);
  //DrawMyPlots("pho2IetaIeta55",selection, "pho2IetaIeta55", "", "pho2IetaIeta55_2018", 50, 0., 0.07, 0.1, 0, 0, 0);

  //DrawMyPlots("ALP_calculatedPhotonIso",selection, "ALP_calculatedPhotonIso", "", "ALP_calculatedPhotonIso_2018", 50, 0., 5., 0.1, 0, 1, 0);
  //DrawMyPlots("var_dR_g1Z",selection, "var_dR_g1Z", "", "var_dR_g1Z_2018", 50, 0., 8., 0.1, 0, 0, 0);
  //DrawMyPlots("var_Pta",selection, "var_Pta", "", "var_Pta_2018", 50, 0., 60., 0.1, 1, 0, 0);
  //DrawMyPlots("var_MhMZ",selection, "var_MhMZ", "", "var_MhMZ_2018", 50, 120., 600., 0.1, 0, 0, 0);
  //DrawMyPlots("H_pt",selection, "H_pt", "", "H_pt_2018", 50, 0., 200., 0.1, 0, 0, 0);
}
