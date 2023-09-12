#include "TStyle.h"
#include "TPad.h"
#include "iostream"

using namespace std;

void setStyle(){

        TStyle *myStyle = new TStyle("myStyle","my Style");
        myStyle->SetPalette(0);
        //myStyle->SetPadTickX(1);
        //myStyle->SetPadTickY(1);
        //myStyle->SetFillColor(0);
        myStyle->SetOptFit(0);             //1= display fitted result, 0 = no
        myStyle->SetOptTitle(0);             //1= display fitted result, 0 = no
        myStyle->SetOptStat(0);            //1= display the statistics, 0 = no
        myStyle->SetTitleFont(22);         // Font for X, Y, and Top title
        myStyle->SetNdivisions(508,"X");
        //myStyle->SetNdivisions(202,"X");
        myStyle->SetNdivisions(508,"Y");
        myStyle->SetNdivisions(508,"Z");

        myStyle->SetPadLeftMargin(0.2);         // pad left margin  for writing Y title
        myStyle->SetPadBottomMargin(0.15);       // pad bottom margin for writing X title
        myStyle->SetPadRightMargin(0.04);       // 0.04 pad right margin  for writing Y title
        myStyle->SetPadTopMargin(0.12);         // 0.04 pad top margin for writing X title

        myStyle->SetMarkerSize(0.8);
        myStyle->SetStatFont(22);          // Font for statistics
        myStyle->SetLabelFont(22,"X");     // Font for X label
        myStyle->SetLabelFont(22,"Y");     // Font for Y label
        myStyle->SetLabelFont(22,"Z");
        myStyle->SetTitleFont(22,"X");     // Font for X title 
        myStyle->SetTitleFont(22,"Y");     // Font for X title 
        myStyle->SetTitleFont(22,"Z");
        //myStyle->SetLabelSize(0.08,"X");   // Size for X label
        //myStyle->SetLabelSize(0.08,"Y");   // Size for Y label
        myStyle->SetLabelSize(0.05,"X");   // Size for X label
        //myStyle->SetLabelSize(0.00,"X");   // Size for X label
        myStyle->SetLabelSize(0.05,"Y");   // Size for Y label
        myStyle->SetLabelSize(0.05,"Z");
        myStyle->SetLabelOffset(0.015,"X");   // Size for X label
        myStyle->SetLabelOffset(0.015,"Y");   // Size for Y label
        myStyle->SetLabelOffset(0.015,"Z");
        myStyle->SetTitleXOffset(1.);     // X title offset
        myStyle->SetTitleYOffset(1.);     // Y title offset
        myStyle->SetTitleXSize(0.08);      // X title size
        myStyle->SetTitleYSize(0.08);      // Y title size
        myStyle->SetHistLineWidth(2);      // Histogram line width
        myStyle->SetStatX(0.95);           // Centroid X position of statistics box
        myStyle->SetStatY(0.95);           // Centroid Y position of statistics box
        myStyle->SetTitleX(0.2);          // Centroid X position of title box
        myStyle->SetTitleY(0.2);          // Centroid Y position of title box
            gROOT->SetStyle("myStyle");
}

void setPad()
{
        //gPad->SetGrid();
        gPad->SetBorderMode(0);
        gPad->SetBorderSize(0);
        gPad->SetFrameFillColor(0);
        gPad->SetFrameBorderMode(0);
        gPad->SetFillColor(0);
        gPad->SetLeftMargin(0.12);         // pad left margin  for writing Y title
        gPad->SetBottomMargin(0.18);       // pad bottom margin for writing X title
        gPad->SetRightMargin(0.13);       // pad right margin  for writing Y title
        gPad->SetTopMargin(0.06);         // pad top margin for writing X title
        //gPad->SetLogy(1);  
}

//#endif