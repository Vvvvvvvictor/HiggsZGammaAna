import numpy as np
import math
from array import array
import ROOT
from ROOT import RooFit

import numpy as np
from scipy import stats
from scipy.optimize import brenth

def getPdfYield(File_name, WS_name, Obj_name, mass, bkg=0):
    f = ROOT.TFile(File_name)
    w = f.Get(WS_name)
    pdf = w.pdf(Obj_name)
    # w.argSet("CMS_hzg_mass").writeToFile("bkg.txt")
    # n = pdf.expectedEvents( w.argSet(mass) ) * w.var(Obj_name+"_norm").getVal()
    w.var(mass).setVal(125)
    w.var("CMS_hzg_mass").setVal(125)
    # w.var(mass).setRange(110, 140)
    # w.allVars().writeToFile("bkg.txt")
    n_tot = pdf.getNormIntegral( w.argSet(mass) ).getVal()
    # print(n_tot)
    w.var(mass).setRange(122, 128)
    n_cen = pdf.getNormIntegral( w.argSet(mass) ).getVal()
    # print(n_cen)
    if bkg:
        n = n_cen / n_tot * w.var(Obj_name+"_norm").getVal()
    else:
        n = n_cen
    # print(n)

    mass_ = w.var("CMS_hzg_mass")
    mass_.setRange(105.0, 170.0)
    mass_.setBins(65)
    canv = ROOT.TCanvas("c1","c1",1250,1000)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    frame = mass_.frame()
    pdf.plotOn(frame)
    frame.Draw()
    canv.SaveAs("bkg.pdf")
    # NBg = pdf.expectedEvents( w.obj('multipdf:roohist_data_mass_cat0').get() )
    return n

if __name__ == "__main__":
    for i in range(4):
        nbkg = getPdfYield("Combine_results_cf/CMS-HGG_mva_13TeV_multipdf_cat{}.root".format(i), "multipdf", "CMS_hgg_cat0_13TeV_bkgshape", "CMS_hzg_mass", 1)
    #     # nsig = getPdfYield("Combine_results_cf/CMS-HGG_sigfit_data_ggh_cat{}.root".format(i), "wsig_13TeV", "hggpdfsmrel_13TeV_ggh_cat0", "MH", 0)
    #     # nbkg = getPdfYield("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_10_2_13/src/flashggFinalFit/ZGMCShape/ZGMCShape.root", "w", "ZGMCShape")
        print "bkg:", nbkg
    #     print "sig:", nsig, "; bkg:", nbkg 
    #     print i, ":", nsig/(nbkg+nsig)**0.5 
    #     # print(nsig)
    # nbkg = getPdfYield("Combine_results_tf/CMS-HGG_mva_13TeV_multipdf_cat0.root", "multipdf", "CMS_hgg_cat0_13TeV_bkgshape", "CMS_hzg_mass", 1)