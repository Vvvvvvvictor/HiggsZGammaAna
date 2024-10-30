#include <iostream>
#include "TDirectory.h"
#include "RooAddPdf.h"
using namespace std;

//int main(){
void RenameWorkspace_DataTagger(){
  TString fn("CMS-HGG_sigfit.root");
  //TString dsetName("ggF_X250_WWgg_qqlnugg_13TeV_HHWWggTag_0"); 
  TString dsetName("ggF_250_13TeV_HHWWggTag_0_125");
  //TString pdfName("hggpdfsmrel_13TeV_ggF_HHWWggTag_0");
  TString newdsetName("sig_ggF_mass_m125_HHWWggTag_0");
  TString wsName("wsig_13TeV");
  TFile *_file0 = TFile::Open(fn,"UPDATE");

  //cout << "_file0 = " << _file0 << endl;
  //TDirectory *direc = _file0->GetDirectory("HHWWggCandidateDumper");
  //RooWorkspace *ws = (RooWorkspace*)direc->Get("cms_HHWWgg_13TeV");
  RooWorkspace *ws = (RooWorkspace*)_file0->Get(wsName);
  //RooWorkspace *ws_clone = ws->Clone();
  //std::cout << "ws = " << ws << std::endl; 
  //direc->cd();
  //TList *l = direc->GetListofKeys();
  //cout << "l = " << l << endl;
  //RooWorkspace *ws = cms_HHWWgg_13TeV;
  //RooAbsData *d = ws->data("ggF_X250_WWgg_qqlnugg_13TeV_SL");
  RooAbsData *d = ws->data(dsetName);
  //RooAbsPdf *p =ws->pdf(pdfName);
  d->SetName(newdsetName);
  //ws->Write();
  //p->SetName("sigpdfrelHHWWggTag_0_allProcs");
  //ws->writeToFile("UpdatedwsName.root");
  
  // Create new file 
  TFile *newFile = TFile::Open("Updated.root","UPDATE");
  
  ws->Write();
}

