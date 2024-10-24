#include <iostream>
#include "TDirectory.h"
#include <string>
#include <cstdlib>
#include <sstream>
#include <cstring> 
using namespace std;

//int main(){
void RenameWorkspace_SignalTagger(){
int mass[19]={250,260,270,280,300,320,350,400,500,550,600,650,700,800,850,900,1000,1250};
//750
//
  for (int i : mass)
{
  std::string NAME="Siganl_X";
  std::string NAMEend="_.root";
  ostringstream oss;
  oss << NAME << i<<NAMEend;
  std::cout<< oss.str() << std::endl;
  ostringstream ossdset;
  std::string dsetNAME="ggF_X";
  std::string dsetNAMEend="_WWgg_qqlnugg_13TeV_SL";
   ossdset<< dsetNAME << i <<dsetNAMEend;
std::cout<< ossdset.str() << std::endl;
  TString fn(oss.str());
  std::cout<< fn << std::endl;
  //TString dsetName("ggF_X250_WWgg_qqlnugg_13TeV_HHWWggTag_0"); 
  TString dsetName(ossdset.str());
 std::cout<< dsetName << std::endl;
 TString newdsetName("ggF_125_13TeV_HHWWggTag_0");
  TString wsName("cms_HHWWgg_13TeV");
  TString newwsName("cms_hgg_13TeV");
TString direcName("HHWWggCandidateDumper");
  TString newdirecName("tagsDumper");
  TFile *_file0 = TFile::Open(fn,"UPDATE");
  cout << "_file0 = " << _file0 << endl;
  //TDirectory *direc = _file0->GetDirectory("HHWWggCandidateDumper");
  TDirectory *direc = _file0->GetDirectory(direcName);
  //RooWorkspace *ws = (RooWorkspace*)direc->Get("cms_HHWWgg_13TeV");
  RooWorkspace *ws = (RooWorkspace*)direc->Get(wsName);
  //RooWorkspace *ws_clone = ws->Clone();
  //std::cout << "ws = " << ws << std::endl; 
  //direc->cd();
  //TList *l = direc->GetListofKeys();
  //cout << "l = " << l << endl;
  //RooWorkspace *ws = cms_HHWWgg_13TeV;
  //RooAbsData *d = ws->data("ggF_X250_WWgg_qqlnugg_13TeV_SL");
  RooAbsData *d = ws->data(dsetName);
  d->SetName(newdsetName);
  direc->cd();
  //ws->Write();
  //ws->writeToFile("UpdatedwsName.root");
  ws->SetName(newwsName);
  
  // Create new file
  ostringstream ossout;
std::string outhead="output_ggF_X"; 
std::string outend = "_WWgg_qqlnugg.root";
  ossout<< outhead << i << outend;
  TString Outfn(ossout.str()) ;
  TFile *newFile = TFile::Open(Outfn,"UPDATE");
  newFile->mkdir(newdirecName);
  TDirectory *newDirec = newFile->GetDirectory(newdirecName);
  newDirec->cd();
  
  ws->Write();
}
}


