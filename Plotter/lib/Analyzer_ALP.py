import os
import sys
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *
from array import array
#from ROOT import TColor
#####################################################################

def getVariableHistsEventsNumber(Tree,varName,sample,cut):

    #Make a canvas to do the work on
    canvas = TCanvas(varName,varName,1000,800)

    if sample in ['Electron', 'Muon']:
        if sample == 'Electron':
            Tree.Draw("{0}>>tree{0}".format(varName),"{0}>0.&&abs(l1_id)==11".format(varName))
        else:
            Tree.Draw("{0}>>tree{0}".format(varName),"{0}>0.&&abs(l1_id)==13".format(varName))
        Hist = gDirectory.Get("tree{0}".format(varName))

        canvas.Clear()
    else:
        #Extract the relevant variable from the trees
        Tree.Draw("{0}>>tree{0}".format(varName),"{0}".format(cut))
        Hist = gDirectory.Get("tree{0}".format(varName))

        canvas.Clear()

    return Hist.Integral()

def getVariableHistsEventsNumber_weight(Tree,varName,sample,cut):

    #Make a canvas to do the work on
    canvas = TCanvas(varName,varName,1000,800)

    #Extract the relevant variable from the trees
    Tree.Draw("{0}>>tree{0}".format(varName),"factor*pho1SFs*pho2SFs*({0})".format(cut))
    Hist = gDirectory.Get("tree{0}".format(varName))

    canvas.Clear()

    return Hist.Integral()

def CountCutFlow(ana_cfg, lumi, outPath, version):

    cutName = ["total Number of events", "HLT", "2 leptons", "isolation", "lep tight", "m_{ll} > 50GeV", "2 photons", "charge Iso", "neutral Iso","HOverE", "\\Delta R(l,\\gamma)>0.4","115 \\leq m_{H} \\leq 135","mva","trigger"]
    varName=["1", "passChaHadIso", "passNeuHadIso","passHOverE","passdR_gl","(H_m>115&&H_m<135)","passBDT","event_passedTrig"]

    outTable = open(outPath + 'cutFlowTable_'+ version + '.txt',"w")
    outTable.write("\\begin{landscape}\n")
    outTable.write("\\renewcommand{\\arraystretch}{1.5}\n")
    outTable.write("\\begin{sidewaystable}[tp]\n")
    outTable.write("\\centering\n")
    outTable.write("\\fontsize{6.5}{8}\\selectfont\n")

    Ncolum = "\\begin{tabular}{|c|"
    name_line = "lumi = " + lumi + " fb^{-1}"
    value = {}
    value_weight = {}

    line_weight = "event weight"
    line_xs = "cross section (pb)"
    weight = {}
    for sample in ana_cfg.samp_names:

        files = TFile(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
        filesTree = files.Get("passedEvents")

        ntuples = TChain("passedEvents","chain_" + sample)
        ntuples.Add(ana_cfg.sample_loc + '/ALP_%s.root' %sample)

        Ncolum = Ncolum + "c|"
        name_line = name_line + "&" + sample
    

        ntuples.GetEvent(1)
        if sample == 'DYJetsToLL':
            xs_value = 6435.0
        elif sample == 'data':
            xs_value = 1.0
        else:
            xs_value = 0.1
        
        if sample == 'data':
            weight_value = 1.0
        else:
            weight_value = xs_value*float(lumi)*1000.0/float(files.nEvents_ntuple.GetEntries())
        line_weight = line_weight + "&" + str(weight_value)
        weight[sample] = weight_value

        value[sample] = {}
        value_weight[sample] = {}
            
        
        line_xs = line_xs + "&" + str(xs_value)

        #value[sample][cutName[0]] = files.nEvents_ntuple.GetEntries()
        value[sample][cutName[0]] = files.nEvents_total.GetBinContent(1)
        value[sample][cutName[1]] = files.nEvents_trig.GetEntries()
        value[sample][cutName[2]] = files.Z_e_nocut.GetEntries() + files.Z_mu_nocut.GetEntries()
        value[sample][cutName[3]] = files.Z_e_lIso.GetEntries() + files.Z_mu_lIso.GetEntries()
        value[sample][cutName[4]] = files.Z_e_lIso_lTight.GetEntries() + files.Z_mu_lIso_lTight.GetEntries()
        value[sample][cutName[5]] = files.Z_50.GetEntries()


        #varName=["H_twopho>0.","Z_Ceta>-90","Z_pho_veto>-90", "passChaHadIso", "passNeuHadIso","passdR_gl","dR_pho<0.15", "passHOverE", "passIeIe", "passPhoIso","passH_m"]
        cut = ""
        for i in range(len(varName)):
            cut = cut + "&" + varName[i]
            cut = cut.lstrip("&")
            value[sample][cutName[i+6]] = getVariableHistsEventsNumber(filesTree, "H_m", sample, cut)
            value_weight[sample][cutName[i+6]] = getVariableHistsEventsNumber_weight(filesTree, "H_m", sample, cut)

    '''
    value_weight['totalMC'] = {}
    value['totalMC'] = {}
    for cut in cutName:
        value_weight['totalMC'][cut] = 0.
        value['totalMC'][cut] = 0.
        for sample in ana_cfg.samp_names:
            if sample in ana_cfg.sig_names + ['data']: continue
            value['totalMC'][cut] += value[sample][cut]
            if cutName.index(cut)<=5:
                value_weight['totalMC'][cut] += value[sample][cut] * weight[sample]
            else:
                value_weight['totalMC'][cut] += value_weight[sample][cut]
    '''

    #Ncolum = Ncolum + "c|" + "}\n"
    #name_line = name_line + "&totalMC" + "\\cr\\hline\n"

    #line_xs = line_xs + "&0" + "\\cr\\hline\n"
    #line_weight = line_weight + "&" + "\\cr\\hline\n"

    Ncolum = Ncolum + "}\n"
    name_line = name_line + "\\cr\\hline\n" #FIXMe

    line_xs = line_xs + "\\cr\\hline\n"
    line_weight = line_weight + "\\cr\\hline\n"

    outTable.write("\\resizebox{\\textwidth}{60mm}{\n")
    outTable.write("\\setlength{\\tabcolsep}{2mm}{\n")
    outTable.write(Ncolum)
    outTable.write("\\hline\n")
    outTable.write(name_line)

    outTable.write(line_xs)
    outTable.write(line_weight)
    outTable.write("\\hline\n")

    for i in range(len(cutName)):
        line_w = cutName[i]
        line_bw = "(befor weighted)"
        line_eff = "(cut efficiency)"
        #for sample in ana_cfg.samp_names + ['totalMC']:
        for sample in ana_cfg.samp_names: #FIXMe
            if sample != 'data':
                line_bw = line_bw + "&" + "(" + str(value[sample][cutName[i]]) + ")"
            else:
                line_bw = line_bw + "&"
            if i <= 5:
                if sample == 'totalMC':
                    line_w = line_w + "&" + str(value[sample][cutName[i]])
                else:
                    line_w = line_w + "&" + str(value[sample][cutName[i]] * weight[sample])
            else:
                line_w = line_w + "&" + str(value_weight[sample][cutName[i]])

            if i==0:
                line_eff = line_eff + "& 1.0"
            elif i <6:
                line_eff = line_eff + "&" + str((value[sample][cutName[i]] * weight[sample])/(value[sample][cutName[i-1]] * weight[sample]))
            elif i <7:
                line_eff = line_eff + "&" + str(value_weight[sample][cutName[i]]/(value[sample][cutName[i-1]] * weight[sample]))
            else:
                line_eff = line_eff + "&" + str(value_weight[sample][cutName[i]]/value_weight[sample][cutName[i-1]])

        outTable.write(line_w+"\\cr\n")
        outTable.write(line_bw+"\\cr\n")
        outTable.write(line_eff+"\\cr\\hline\n")
        if i == 4 or i == 5: outTable.write("\\hline\n")


    outTable.write("\\end{tabular}}}\n")
    outTable.write("\\end{sidewaystable}\n")
    outTable.write("\\end{landscape}\n")


def CountCutFlow_less(ana_cfg, lumi, outPath):

    cutName = ["total Number of events", "HLT", "all", "isolation", "lep tight", "m_{ll} > 50GeV", "find photon (at least 2 photons pt \\ge 10GeV)", "photon eta (\\vert\\eta|<1.4442 \\cup 1.566<|\\eta|<2.5)", "conversion veto", "charge Iso", "neutral Iso", "\\Delta R(l,\\gamma)>0.4", "\\Delta R(\\gamma\\gamma)<0.15", "sigmaIetaIeta", "HOverE", "photon Iso", "100\\leqslant m_{H}\\leqslant 180"]

    outTable = open(outPath + 'cutFlowTable_data_less.txt',"w")
    outTable.write("\\begin{landscape}\n")
    outTable.write("\\renewcommand{\\arraystretch}{1.5}\n")
    outTable.write("\\begin{table}[tp]\n")
    outTable.write("\\centering\n")
    outTable.write("\\fontsize{6.5}{8}\\selectfont\n")

    Ncolum = "\\begin{tabular}{|c|"
    name_line = "lumi = " + lumi + " fb^{-1}"
    value = {}

    line_weight = "event weight"
    line_xs = "cross section (pb)"
    weight = {}
    name = ''
    for sample in ana_cfg.samp_names:
        files = TFile(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
        filesTree = files.Get("passedEvents")

        ntuples = TChain("passedEvents","chain_" + sample)
        ntuples.Add(ana_cfg.sample_loc + '/ALP_%s.root' %sample)



        #weight_line = weight_line + "&" +

        ntuples.GetEvent(1)
        line_weight = line_weight + "&" + str(ntuples.event_weight)
        weight[sample] = ntuples.event_weight

        value[sample] = {}

        if sample in ['data', 'QCD']:
            line_xs = line_xs + "&" + "1.0"

            value[sample][cutName[0]] = 0
            value[sample][cutName[1]] = 0
            value[sample][cutName[2]] = 0
            value[sample][cutName[3]] = 0
            value[sample][cutName[4]] = 0
            value[sample][cutName[5]] = 0

        else:
            line_xs = line_xs + "&" + str(files.cross_section.GetBinContent(1)/(files.Events_weight.GetBinContent(1)/weight[sample]))

            value[sample][cutName[0]] = files.nEvents_ntuple.GetEntries()
            value[sample][cutName[1]] = files.nEvents_trig.GetEntries()
            value[sample][cutName[2]] = files.Z_e_nocut.GetEntries() + files.Z_mu_nocut.GetEntries()
            value[sample][cutName[3]] = files.Z_e_lIso.GetEntries() + files.Z_mu_lIso.GetEntries()
            value[sample][cutName[4]] = files.Z_e_lIso_lTight.GetEntries() + files.Z_mu_lIso_lTight.GetEntries()
            value[sample][cutName[5]] = files.Z_50.GetEntries()

        varName=["H_twopho","Z_Ceta","Z_pho_veto", "Z_CIso", "Z_NIso","Z_dR","Z_dR_pho", "Z_IeIe", "Z_HOE", "Z_PIso","Z_dR_pho_Cmh"]
        for i in range(len(varName)):
            value[sample][cutName[i+6]] = getVariableHistsEventsNumber(filesTree, varName[i], sample)

        if sample in ['WGGJets','WWToLNuQQ','WWW','WWZ','WZZ']:
            name = name + '/' + sample
            if sample !='WZZ': continue
        else:
            name = sample
        Ncolum = Ncolum + "c|"
        if sample =='WZZ':
            name_line = name_line + "&" + name.lstrip('TTTo2L2Nu/')
        else:
            name_line = name_line + "&" + name

    Ncolum = Ncolum + "}\n"
    name_line = name_line + "\\cr\\hline\n"

    outTable.write("\\resizebox{\\textwidth}{35mm}{\n")
    outTable.write("\\setlength{\\tabcolsep}{2mm}{\n")
    outTable.write(Ncolum)
    outTable.write("\\hline\n")
    outTable.write(name_line)

    outTable.write("\\hline\n")

    for i in range(len(cutName)):
        line_w = cutName[i]
        N = 0.0
        for sample in ana_cfg.samp_names:

            if sample in ['WGGJets','WWToLNuQQ','WWW','WWZ','WZZ']:
                N = N + value[sample][cutName[i]] * weight[sample]
                if sample !='WZZ': continue
                line_w = line_w + "&" + str(N)
            else:
                line_w = line_w + "&" + str(value[sample][cutName[i]] * weight[sample])
        outTable.write(line_w+"\\cr\\hline\n")
        if i == 4 or i == 5: outTable.write("\\hline\n")


    outTable.write("\\end{tabular}}}\n")
    outTable.write("\\end{table}\n")
    outTable.write("\\end{landscape}\n")


def CountCutFlow_mva(ana_cfg, lumi, outPath):

    if ana_cfg.mass == 'M1':
        cutName = ["total Number of events", "HLT", "all", "isolation", "lep tight", "m_{ll} > 50GeV", "find photon (at least 2 photons pt \\ge 10GeV)", "photon eta (\\vert\\eta|<1.4442 \\cup 1.566<|\\eta|<2.5)", "conversion veto", "charge Iso", "neutral Iso", "\\Delta R(l,\\gamma)>0.4", "\\Delta R(\\gamma\\gamma)<0.15", "HOverE", "mva"]
    elif ana_cfg.mass == 'M5':
        cutName = ["total Number of events", "HLT", "all", "isolation", "lep tight", "m_{ll} > 50GeV", "find photon (at least 2 photons pt \\ge 10GeV)", "photon eta (\\vert\\eta|<1.4442 \\cup 1.566<|\\eta|<2.5)", "conversion veto", "charge Iso", "neutral Iso", "\\Delta R(l,\\gamma)>0.4", "0.1<\\Delta R(\\gamma\\gamma)<0.5", "HOverE", "mva"]
    elif ana_cfg.mass == 'M15':
        cutName = ["total Number of events", "HLT", "all", "isolation", "lep tight", "m_{ll} > 50GeV", "find photon (at least 2 photons pt \\ge 10GeV)", "photon eta (\\vert\\eta|<1.4442 \\cup 1.566<|\\eta|<2.5)", "conversion veto", "charge Iso", "neutral Iso", "\\Delta R(l,\\gamma)>0.4", "0.2<\\Delta R(\\gamma\\gamma)<2", "HOverE", "mva"]
    else:
        cutName = ["total Number of events", "HLT", "all", "isolation", "lep tight", "m_{ll} > 50GeV", "find photon (at least 2 photons pt \\ge 10GeV)", "photon eta (\\vert\\eta|<1.4442 \\cup 1.566<|\\eta|<2.5)", "conversion veto", "charge Iso", "neutral Iso", "\\Delta R(l,\\gamma)>0.4", "0.6<\\Delta R(\\gamma\\gamma)<6", "HOverE", "mva"]
    outTable = open(outPath + 'cutFlowTable_mva_'+ana_cfg.mass+'.txt',"w")
    outTable.write("\\begin{landscape}\n")
    outTable.write("\\renewcommand{\\arraystretch}{1.5}\n")
    outTable.write("\\begin{table}[tp]\n")
    outTable.write("\\centering\n")
    outTable.write("\\fontsize{6.5}{8}\\selectfont\n")

    Ncolum = "\\begin{tabular}{|c|"
    name_line = "lumi = " + lumi + " fb^{-1}"
    value = {}
    value_weight = {}

    line_weight = "event weight"
    line_xs = "cross section (pb)"
    weight = {}
    for sample in ana_cfg.samp_names:
        files = TFile(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
        filesTree = files.Get("passedEvents")

        ntuples = TChain("passedEvents","chain_" + sample)
        ntuples.Add(ana_cfg.sample_loc + '/ALP_%s.root' %sample)

        Ncolum = Ncolum + "c|"
        name_line = name_line + "&" + sample
        #weight_line = weight_line + "&" +

        ntuples.GetEvent(1)
        line_weight = line_weight + "&" + str(ntuples.event_weight)
        weight[sample] = ntuples.event_weight

        value[sample] = {}
        value_weight[sample] = {}

        line_xs = line_xs + "&" + str(files.cross_section.GetBinContent(1)/(files.Events_weight.GetBinContent(1)/weight[sample]))

        value[sample][cutName[0]] = files.nEvents_ntuple.GetEntries()
        value[sample][cutName[1]] = files.nEvents_trig.GetEntries()
        value[sample][cutName[2]] = files.Z_e_nocut.GetEntries() + files.Z_mu_nocut.GetEntries()
        value[sample][cutName[3]] = files.Z_e_lIso.GetEntries() + files.Z_mu_lIso.GetEntries()
        value[sample][cutName[4]] = files.Z_e_lIso_lTight.GetEntries() + files.Z_mu_lIso_lTight.GetEntries()
        value[sample][cutName[5]] = files.Z_50.GetEntries()


        varName=["H_twopho","Z_Ceta","Z_pho_veto", "Z_CIso", "Z_NIso","Z_dR","Z_dR_pho", "Z_HOE", "Z_m"]
        for i in range(len(varName)):

            value[sample][cutName[i+6]] = getVariableHistsEventsNumber(filesTree, "H_twopho", sample, cut)
            value_weight[sample][cutName[i+6]] = getVariableHistsEventsNumber_weight(filesTree, "H_twopho", sample, cut)



    Ncolum = Ncolum + "}\n"
    name_line = name_line + "\\cr\\hline\n"

    line_xs = line_xs + "\\cr\\hline\n"
    line_weight = line_weight + "\\cr\\hline\n"

    outTable.write("\\resizebox{\\textwidth}{55mm}{\n")
    outTable.write("\\setlength{\\tabcolsep}{2mm}{\n")
    outTable.write(Ncolum)
    outTable.write("\\hline\n")
    outTable.write(name_line)

    outTable.write(line_xs)
    outTable.write(line_weight)
    outTable.write("\\hline\n")

    for i in range(len(cutName)):
        line_w = cutName[i]
        line_bw = "(befor weighted)"
        for sample in ana_cfg.samp_names:
            if sample != 'data':
                line_bw = line_bw + "&" + "(" + str(value[sample][cutName[i]]) + ")"
            else:
                line_bw = line_bw + "&"
            if i <= 5:
                line_w = line_w + "&" + str(value[sample][cutName[i]] * weight[sample])
            else:
                line_w = line_w + "&" + str(value_weight[sample][cutName[i]])
        outTable.write(line_w+"\\cr\n")
        outTable.write(line_bw+"\\cr\\hline\n")
        if i == 4 or i == 5: outTable.write("\\hline\n")


    outTable.write("\\end{tabular}}}\n")
    outTable.write("\\end{table}\n")
    outTable.write("\\end{landscape}\n")

def CountCutFlow_mva_less(ana_cfg, lumi, outPath):

    cutName = ["total Number of events", "HLT", "all", "isolation", "lep tight", "m_{ll} > 50GeV", "find photon (at least 2 photons pt \\ge 10GeV)", "photon eta (\\vert\\eta|<1.4442 \\cup 1.566<|\\eta|<2.5)", "conversion veto", "mva", "\\Delta R(l,\\gamma)>0.4", "\\Delta R(\\gamma\\gamma)<1.0", "100\\leqslant m_{H}\\leqslant 180"]

    outTable = open(outPath + 'cutFlowTable_mva_less.txt',"w")
    outTable.write("\\begin{landscape}\n")
    outTable.write("\\renewcommand{\\arraystretch}{1.5}\n")
    outTable.write("\\begin{table}[tp]\n")
    outTable.write("\\centering\n")
    outTable.write("\\fontsize{6.5}{8}\\selectfont\n")

    Ncolum = "\\begin{tabular}{|c|"
    name_line = "lumi = " + lumi + " fb^{-1}"
    value = {}

    line_weight = "event weight"
    line_xs = "cross section (pb)"
    weight = {}
    name = ''
    for sample in ana_cfg.samp_names:
        files = TFile(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
        filesTree = files.Get("passedEvents")

        ntuples = TChain("passedEvents","chain_" + sample)
        ntuples.Add(ana_cfg.sample_loc + '/ALP_%s.root' %sample)



        #weight_line = weight_line + "&" +

        ntuples.GetEvent(1)
        line_weight = line_weight + "&" + str(ntuples.event_weight)
        weight[sample] = ntuples.event_weight

        value[sample] = {}

        if sample in ['data', 'QCD']:
            line_xs = line_xs + "&" + "1.0"

            value[sample][cutName[0]] = 0
            value[sample][cutName[1]] = 0
            value[sample][cutName[2]] = 0
            value[sample][cutName[3]] = 0
            value[sample][cutName[4]] = 0
            value[sample][cutName[5]] = 0

        else:
            line_xs = line_xs + "&" + str(files.cross_section.GetBinContent(1))

            value[sample][cutName[0]] = files.nEvents_total.GetBinContent(1)
            value[sample][cutName[1]] = files.nEvents_trig.GetEntries()
            value[sample][cutName[2]] = files.Z_e_nocut.GetEntries() + files.Z_mu_nocut.GetEntries()
            value[sample][cutName[3]] = files.Z_e_lIso.GetEntries() + files.Z_mu_lIso.GetEntries()
            value[sample][cutName[4]] = files.Z_e_lIso_lTight.GetEntries() + files.Z_mu_lIso_lTight.GetEntries()
            value[sample][cutName[5]] = files.Z_50.GetEntries()

        varName=["H_twopho","Z_Ceta","Z_pho_veto","Z_pho_veto_mva","Z_dR","Z_dR_pho","Z_dR_pho_Cmh"]
        for i in range(7):
            value[sample][cutName[i+6]] = getVariableHistsEventsNumber(filesTree, varName[i], sample)

        if sample in ['WGGJets','WWToLNuQQ','WWW','WWZ','WZZ']:
            name = name + '/' + sample
            if sample !='WZZ': continue
        else:
            name = sample
        Ncolum = Ncolum + "c|"
        if sample =='WZZ':
            name_line = name_line + "&" + name.lstrip('TTTo2L2Nu/')
        else:
            name_line = name_line + "&" + name

    Ncolum = Ncolum + "}\n"
    name_line = name_line + "\\cr\\hline\n"

    outTable.write("\\resizebox{\\textwidth}{25mm}{\n")
    outTable.write("\\setlength{\\tabcolsep}{2mm}{\n")
    outTable.write(Ncolum)
    outTable.write("\\hline\n")
    outTable.write(name_line)

    outTable.write("\\hline\n")

    for i in range(len(cutName)):
        line_w = cutName[i]
        N = 0.0
        for sample in ana_cfg.samp_names:

            if sample in ['WGGJets','WWToLNuQQ','WWW','WWZ','WZZ']:
                N = N + value[sample][cutName[i]] * weight[sample]
                if sample !='WZZ': continue
                line_w = line_w + "&" + str(N)
            else:
                line_w = line_w + "&" + str(value[sample][cutName[i]] * weight[sample])
        outTable.write(line_w+"\\cr\\hline\n")
        if i == 4 or i == 5: outTable.write("\\hline\n")


    outTable.write("\\end{tabular}}}\n")
    outTable.write("\\end{table}\n")
    outTable.write("\\end{landscape}\n")


def PIso2D(canv, hitos2D_sig, hitos2D_bkg, isEB):

    canv.SetBottomMargin(0.1)
    canv.cd()

    upper_pad = TPad("upperpad", "upperpad", 0,0.505, 1,1)
    upper_pad.Draw()
    upper_pad.cd()
    hitos2D_sig.GetXaxis().SetTitle("Pflow photon isolation")
    hitos2D_sig.GetYaxis().SetTitle("photon Pt")
    hitos2D_sig.GetYaxis().SetTitleOffset(1.0)
    hitos2D_sig.Draw("COLZ")

    global l1
    if isEB:
        l1 = TLine(2.044,0,2.16451,30)
    else:
        l1 = TLine(3.032,0,3.143,30)
    l1.SetLineWidth(3)
    l1.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad", "lowerpad_", 0, 0.025, 1,0.495)
    lower_pad.Draw()
    lower_pad.cd()
    hitos2D_bkg.GetXaxis().SetTitle("Pflow photon isolation")
    hitos2D_bkg.GetYaxis().SetTitle("photon Pt")
    hitos2D_bkg.GetYaxis().SetTitleOffset(1.0)
    hitos2D_bkg.Draw("COLZ")

    global l2
    if isEB:
        l2 = TLine(2.044,0,2.16451,30)
    else:
        l2 = TLine(3.032,0,3.143,30)
    l2.SetLineWidth(3)
    l2.Draw("SAME")


def plot2D_mgg(canv, hitos2D_sig):

    canv.SetBottomMargin(0.1)
    canv.cd()

    upper_pad = TPad("upperpad", "upperpad", 0.05,0.05, 0.95,0.95)
    upper_pad.SetLeftMargin(0.13)
    upper_pad.SetRightMargin(0.2)
    upper_pad.SetTopMargin(0.085)
    upper_pad.Draw()
    upper_pad.cd()
    hitos2D_sig.GetXaxis().SetTitle("M_{#gamma#gamma} (GeV)")
    hitos2D_sig.GetXaxis().SetTitleSize(0.05)
    hitos2D_sig.GetYaxis().SetTitle("mvaVal")
    hitos2D_sig.GetYaxis().SetTitleSize(0.05)
    hitos2D_sig.GetYaxis().SetTitleOffset(1.0)
    gStyle.SetOptStat(0)
    '''
    stops = [ 0.00, 0.50, 1.00 ]
    red   = [ 1.00, 0.00, 0.00 ]
    green = [ 0.00, 1.00, 0.00 ]
    blue  = [ 0.00, 0.00, 1.00 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, 50)
    gStyle.SetNumberContours(50)
    '''
    
    
    hitos2D_sig.Draw("COLZ")
    lt = TLatex();
    lt.DrawLatexNDC(0.3,0.8,"Correlation: {:.3g}".format(hitos2D_sig.GetCorrelationFactor()))

    #canv.cd()
    #lower_pad = TPad("lowerpad", "lowerpad_", 0, 0.025, 1,0.495)
    #lower_pad.Draw()
    #lower_pad.cd()
    #hitos2D_bkg.GetXaxis().SetTitle("M_{llgg}")
    #hitos2D_bkg.GetXaxis().SetTitleSize(0.15)
    #hitos2D_bkg.GetYaxis().SetTitle("mvaVal")
    #hitos2D_bkg.GetYaxis().SetTitleSize(0.15)
    #hitos2D_bkg.GetYaxis().SetTitleOffset(1.0)
    #hitos2D_bkg.Draw("COLZ")

def plot2D_mllgg(canv, hitos2D_sig):

    canv.SetBottomMargin(0.1)
    canv.cd()

    upper_pad = TPad("upperpad", "upperpad", 0.05,0.05, 0.95,0.95)
    upper_pad.Draw()
    upper_pad.cd()
    hitos2D_sig.GetXaxis().SetTitle("M_{ll#gamma#gamma}")
    hitos2D_sig.GetXaxis().SetTitleSize(0.05)
    hitos2D_sig.GetYaxis().SetTitle("mvaVal")
    hitos2D_sig.GetYaxis().SetTitleSize(0.05)
    hitos2D_sig.GetYaxis().SetTitleOffset(1.0)
    gStyle.SetOptStat(0)
    '''
    stops = [ 0.00, 0.50, 1.00 ]
    red   = [ 1.00, 0.00, 0.00 ]
    green = [ 0.00, 1.00, 0.00 ]
    blue  = [ 0.00, 0.00, 1.00 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, 50)
    gStyle.SetNumberContours(50)
    '''
    
    
    hitos2D_sig.Draw("COLZ")
    lt = TLatex();
    lt.DrawLatexNDC(0.3,0.8,"Correlation: {:.3g}".format(hitos2D_sig.GetCorrelationFactor()))


def plot2D_CONT(canv, hitos2D_sig, hitos2D_bkg):

    canv.SetBottomMargin(0.1)
    canv.cd()

    upper_pad = TPad("upperpad", "upperpad", 0,0.505, 1,1)
    upper_pad.Draw()
    upper_pad.cd()
    hitos2D_sig.GetXaxis().SetTitle("Mh+MZ")
    #hitos2D_sig.GetXaxis().SetTitleSize(0.15)
    hitos2D_sig.GetYaxis().SetTitle("Mh")
    #hitos2D_sig.GetYaxis().SetTitleSize(0.15)
    hitos2D_sig.GetYaxis().SetTitleOffset(1.0)
    hitos2D_sig.Draw("CONT1")

    canv.cd()
    lower_pad = TPad("lowerpad", "lowerpad_", 0, 0.025, 1,0.495)
    lower_pad.Draw()
    lower_pad.cd()
    hitos2D_bkg.GetXaxis().SetTitle("Mh+MZ")
    #hitos2D_bkg.GetXaxis().SetTitleSize(0.15)
    hitos2D_bkg.GetYaxis().SetTitle("Mh")
    #hitos2D_bkg.GetYaxis().SetTitleSize(0.15)
    hitos2D_bkg.GetYaxis().SetTitleOffset(1.0)
    hitos2D_bkg.Draw("CONT1")



def CountYield(ana_cfg, histos, sys_name):

    file = open('./BDTSys/BDTSys.txt', 'a')
    print('Printing event yield: ')

    '''
    title = '{0:15}\t\t'.format('sample')
    for sys_name in sys_names:
        if sys_name == 'norm':
            title = title + '{0:15}\t\t'.format(sys_name)
        else:
            title = title + '{0:15}\t\t{0:15}\t\t'.format(sys_name + '_up', sys_name + '_dn')
    print title
    file.write(title + '\n')

    for sample in ana_cfg.samp_names:
        line = '{0:15}\t'.format(sample)
        for sys_name in sys_names:
            if sys_name == 'norm':
                line = line + '{0:15f}\t\t'.format(histos[sample]['norm'].Integral())
            else:
                line = line + '{0:15f}\t\t{0:15f}\t\t'.format(histos[sample][sys_name]['up'].Integral(), histos[sample][sys_name]['dn'].Integral())
        print line
        file.write(line + '\n')

    '''
    #title = 'sample \t\t normal'
    if sys_name in ['ShowerShape', 'pho_scale', 'pho_smear', 'lep_scale', 'lep_smear']:
        title = 'sample' + ' \t\t ' + sys_name + '_up' + ' \t\t ' + sys_name + '_dn'
        #title = '{0:15}  {1:15}  {2:15}'.format('sample', sys_name + '_up', sys_name + '_dn')
    else:
        title = 'sample \t\t normal'
        #title = '{0:15}  {1:15}'.format('sample', 'normal')
    print(title)
    file.write(title + '\n')
    for sample in ana_cfg.bkg_names:
        if sys_name in ['ShowerShape', 'pho_scale', 'pho_smear', 'lep_scale', 'lep_smear']:
            line = '%s' %(sample) + ' \t\t ' + str(histos[sample][sys_name]['up'].Integral()) + ' \t\t ' + str(histos[sample][sys_name]['dn'].Integral())
            #line = '{0:15}  {1:15d}  {2:15d}'.format(sample, histos[sample][sys_name]['up'].Integral(), histos[sample][sys_name]['dn'].Integral())
        else:
            line = '%s \t\t %f' %(sample, histos[sample]['norm'].Integral())
            #line = '{0:15}  {1:15d}'.format(sample, histos[sample]['norm'].Integral())
        print(line)
        file.write(line + '\n')



def CountNormalizationYield(ana_cfg, histos, sys_names, channel, year):

    if year == "-2016":
        file = open('./BDTSys/NormalizationSys_16APV_'+channel+'.txt', 'a')
    else:
        file = open('./BDTSys/NormalizationSys_'+year.lstrip('20')+'_'+channel+'.txt', 'a')
    print('Printing event yield: ')
    file.write('year:'+year+'_'+channel+ '\n')


    title = '{0:15}\t\t'.format('sample')
    for sys_name in sys_names:
            title = title + '{0:25}\t\t\t'.format(sys_name)
    print(title)
    file.write(title + '\n')

    for sample in ana_cfg.sig_names:
        line = '{0:10}\t'.format(sample)
        for sys_name in sys_names:
                line = line + '{0:15f} (unc:{1:5f})\t\t'.format(histos[sample][sys_name].Integral(),(histos[sample][sys_name].Integral()-histos[sample][sys_names[0]].Integral())/histos[sample][sys_names[0]].Integral())
        print(line)
        file.write(line + '\n')



def CountBDTSys(ana_cfg, histos, sys_names, channel, year):

    if year == "-2016":
        file = open('./BDTSys/BDTSys_16APV_'+channel+'.txt', 'a')
    else:
        file = open('./BDTSys/BDTSys_'+year.lstrip('20')+'_'+channel+'.txt', 'a')

    print('Printing event yield: ')

    file.write("year:"+ year +'_'+channel + '\n')

    title = '{0:15}\t\t'.format('sample')
    for sys_name in sys_names:
        if sys_name=='norm':
            title = title + '{0:25}\t\t\t'.format(sys_name)
        else:
            title = title + '{0:25}\t\t\t{1:25}\t\t\t'.format(sys_name+'_up',sys_name+'_dn')
    print(title)
    file.write(title + '\n')



    for sample in ana_cfg.sig_names:
        norm = histos[sample][sys_names[0]].Integral()
        line = '{0:10}\t'.format(sample)
        for sys_name in sys_names:
            if sys_name=='norm':
                line = line + '{0:15f} \t\t'.format(histos[sample][sys_name].Integral())
            else:
                line = line + '{0:15f} (unc:{1:5f})\t\t{2:15f} (unc:{3:5f})\t\t'.format(histos[sample][sys_name]['up'].Integral(),(histos[sample][sys_name]['up'].Integral()-norm)/norm,histos[sample][sys_name]['dn'].Integral(),(histos[sample][sys_name]['dn'].Integral()-norm)/norm)
        print(line)
        file.write(line + '\n')
