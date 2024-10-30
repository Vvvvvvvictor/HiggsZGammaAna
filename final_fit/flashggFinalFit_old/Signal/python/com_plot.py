import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, kBlack, TLatex, gPad
import CMS_lumi
import tdrstyle

import argparse
import math

import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-SM", "--SM_Radion", action="store_true",
                    default=False, help="Display SM limits", required=False)
parser.add_argument("-pb", "--pb", action="store_true", default=False,
                    help="Display limits in pb unit", required=False)
parser.add_argument("-fb", "--fb", action="store_true", default=False,
                    help="Display limits in fb unit", required=False)
parser.add_argument("-y", "--Year", dest="year",
                    default="17", help="which year's datasetes")
parser.add_argument("-i", "--input", dest="input",
                    default="/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/HZGamma_SigModel_UL/Combine_results", help="limit files path")
parser.add_argument("-o", "--output", dest="output",
                    default="/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/HZGamma_SigModel_UL/Combine_results", help="plot out path")
parser.add_argument("--observe", action="store_true", default=False,
                    help="plot observed limits?", required=False)
args = parser.parse_args()

ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Preliminary"
#CMS_lumi.cmsText = ""
#CMS_lumi.extraText = ""
#CMS_lumi.extraText = "Private"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


# GET limits from root file
def getLimits(file_name):
    file = TFile(file_name)

    if file:
        print 'file opened'
    tree = file.Get("limit")
    if tree:
        print 'tree opened'

    limits = []
    print'tree = ', tree
    for quantile in tree:
        print'quantile = ', quantile
        
        limits.append(tree.limit)

    if args.observe:
        return limits[:6]
    else:
        return limits[:5]


# PLOT upper limits
def plotUpperLimits(cats, limit, year):
    N = len(cats)
    print str(N) + ' categories'

    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    observed = TGraph(N)      # observed line

    print'limit = ', limit

    x_mass = np.array([2., 4., 6., 8., 10.])
    x_mass_err = np.array([1., 1., 1., 1., 1.])

    limit_yellow_up = np.zeros(shape=(N))
    limit_green_up = np.zeros(shape=(N))
    limit_medium = np.zeros(shape=(N))
    limit_green_dn = np.zeros(shape=(N))
    limit_yellow_dn = np.zeros(shape=(N))

    for i in range(N):
        
        limit_yellow_up[i] = limit[cats[i]][4] - limit[cats[i]][2]
        limit_green_up[i] = limit[cats[i]][3] - limit[cats[i]][2]
        limit_medium[i] = limit[cats[i]][2]
        limit_green_dn[i] = limit[cats[i]][2] - limit[cats[i]][1] 
        limit_yellow_dn[i] = limit[cats[i]][2] - limit[cats[i]][0]

    print "DEBUG", limit_medium
    print "DEBUG", limit_green_dn, limit_green_up
    print "DEBUG", limit_yellow_dn, limit_yellow_up
    green = ROOT.TGraphAsymmErrors(N, x_mass, limit_medium, x_mass_err, x_mass_err, limit_green_dn, limit_green_up)
    yellow = ROOT.TGraphAsymmErrors(N, x_mass, limit_medium, x_mass_err, x_mass_err, limit_yellow_dn, limit_yellow_up)
    medium = ROOT.TGraphAsymmErrors(N, x_mass, limit_medium, x_mass_err, x_mass_err, np.zeros(shape=(N)), np.zeros(shape=(N)))

    W = 800
    H = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c", "c", 100, 100, W, H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(L/W)
    c.SetRightMargin(R/W)
    c.SetTopMargin(T/H)
    c.SetBottomMargin(B/H)
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()

    yellow.GetYaxis().CenterTitle()
    yellow.GetYaxis().SetTitleSize(0.04)
    yellow.GetXaxis().SetTitleSize(0.05)
    yellow.GetXaxis().SetLabelSize(0.04)
    yellow.GetYaxis().SetLabelSize(0.04)
    yellow.GetYaxis().SetTitleOffset(1.3)
    yellow.GetXaxis().SetNdivisions(6)
    yellow.GetYaxis().CenterTitle(True)
    #yellow.GetXaxis().SetTitle("m_{a} [GeV]")
    yellow.GetYaxis().SetTitle("95% CL limits on #sigma(pp#rightarrow H#rightarrow ll#gamma)/#sigma_{SM}")
    yellow.SetMinimum(1.)
    yellow.SetMaximum(200.)
    yellow.GetXaxis().SetLimits(0, 12)
    yellow.GetXaxis().ChangeLabel(1,-1,-1,-1,-1,-1," ")
    yellow.GetXaxis().ChangeLabel(2,-1,-1,-1,-1,-1,"cat0")
    yellow.GetXaxis().ChangeLabel(3,-1,-1,-1,-1,-1,"cat1")
    yellow.GetXaxis().ChangeLabel(4,-1,-1,-1,-1,-1,"cat2")
    yellow.GetXaxis().ChangeLabel(5,-1,-1,-1,-1,-1,"cat3")
    yellow.GetXaxis().ChangeLabel(6,-1,-1,-1,-1,-1,"Combined")
    yellow.GetXaxis().ChangeLabel(7,-1,-1,-1,-1,-1," ")

    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('a2')

    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('same2')

    medium.SetLineColor(1)
    medium.SetLineWidth(2)
    medium.SetLineStyle(2)
    medium.Draw('sameP')

    CMS_lumi.lumiText_posY = 0.0075
    CMS_lumi.lumiText_posX = 0.0
    CMS_lumi.CMSText_posX = 0.
    CMS_lumi.CMS_lumi(c, 4, 0, year)
    #ROOT.gPad.SetTicks(1, 1)
    #frame.Draw('sameaxis')

    # yboost = 0.075
    yboost = -0.1

    x1 = 0.65
    x2 = x1 + 0.24
    y2 = 0.95 + yboost
    y1 = 0.80 + yboost
    legend = TLegend(x1, y1, x2, y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextFont(42)
    #legend.AddEntry(median, "AsymptoticLimits CL_{s} expected", 'L')

    legend.AddEntry(median, "Median expected", 'L')
    legend.AddEntry(green, "68% CL_{s} expected", 'f')
    legend.AddEntry(yellow, "95% CL_{s} expected", 'f')
    legend.Draw()

    
    label = TLatex()
    label.SetNDC()
    label.SetTextAngle(0)
    label.SetTextColor(kBlack)
    label.SetTextFont(42)
    label.SetTextSize(0.045)
    label.SetLineWidth(2)
    label.DrawLatex(0.7,0.7 + yboost,"STAT only");
    

    c.SaveAs("{}/UpperLimit.pdf".format(args.output))
    c.SaveAs("{}/UpperLimit.png".format(args.output))
    c.Close()


# main
def main():

    verbose = 0
    year = args.year
    f = TFile()

    values = []
    limits = {}

    cats = ['cat0', 'cat1', 'cat2', 'cat3', 'allCats']

    for cat in cats:
        file_name = "{}/higgsCombine{}.AsymptoticLimits.mH125.root".format(args.input, cat)
        if verbose: print file_name

        limits[cat] = getLimits(file_name)
        if verbose: print limit[m]

    plotUpperLimits(cats, limits, year)


if __name__ == '__main__':
    main()
