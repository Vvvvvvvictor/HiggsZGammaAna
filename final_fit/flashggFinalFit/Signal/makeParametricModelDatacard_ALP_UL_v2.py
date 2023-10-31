#!/usr/bin/env python

# Script adapted from original by Matt Kenzie.
# Used for Dry Run of Dec 2015 Hgg analysis.

###############################################################################
## IMPORTS ####################################################################
###############################################################################
from optparse import OptionParser
import os
import sys
import copy
import math
###############################################################################

###############################################################################
## PARSE ROOT MACROS  #########################################################
###############################################################################
import ROOT as r
# if options.quadInterpolate:
#  r.gROOT.ProcessLine(".L quadInterpolate.C+g")
#  from ROOT import quadInterpolate
r.gROOT.ProcessLine(
    ".L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
#r.gROOT.ProcessLine(".L ../libLoopAll.so")
###############################################################################

###############################################################################
## WSTFileWrapper  ############################################################
###############################################################################


class WSTFileWrapper:
    # self.fnList = [] # filename list
    # self.fileList = [] #file list
    # self.wsList = [] #workspace list

    def __init__(self, files, wsname):
        self.fnList = files.split(",")  # [1]
        self.fileList = []
        self.wsList = []  # now list of ws names...
        # print files
        for fn in self.fnList:  # [2]
            f = r.TFile.Open(fn)
            self.fileList.append(f)
            thing = f.Get(wsname)
            self.wsList.append(self.fileList[-1].Get(wsname))
            f.Close()

    def convertTemplatedName(self, dataName):
        theProcName = ""
        theDataName = ""
        tpMap = {"GG2H": "ggh", "VBF": "vbf", "TTH": "tth", "QQ2HLNU": "wh", "QQ2HLL": "zh",
                 "WH2HQQ": "wh", "ZH2HQQ": "zh", "testBBH": "bbh", "testTHQ": "th", "testTHW": "th"}
        for stxsProc in tpMap:
            if dataName.startswith(stxsProc):
                theProcName = stxsProc
                theDataName = dataName.replace(stxsProc, tpMap[stxsProc], 1)
        return [theDataName, theProcName]

    def data(self, dataName):
        thePair = self.convertTemplatedName(dataName)
        newDataName = thePair[0]
        newProcName = thePair[1]
        result = None
        complained_yet = 0
        for i in range(len(self.fnList)):
            if self.fnList[i] != "current file":
                if newProcName not in self.fnList[i] and newProcName != "":
                    continue
                this_result_obj = self.wsList[i].data(newDataName)
                if (result and this_result_obj and (not complained_yet)):
                    complained_yet = True
                if this_result_obj:  # [3]
                    result = this_result_obj
        return result

    def var(self, varName):
        result = None
        complained_yet = 0
        for i in range(len(self.fnList)):
            this_result_obj = self.wsList[i].var(varName)
            if this_result_obj:  # [3]
                result = this_result_obj

        return result

###############################################################################


###############################################################################
## OPTION PARSING  ############################################################
###############################################################################
parser = OptionParser()
parser.add_option("-i", "--infilename",
                  help="Input file (binned signal from flashgg)")
parser.add_option("-s", "--sysinfilename", help="Input sys files")
parser.add_option("--NormSys", help="Input normalization sys files")
parser.add_option("--BDTSys", help="Input BDT sys files")
parser.add_option("--channel", help="channel name")
parser.add_option("--interp", help="Input interpolation files")
parser.add_option("-o", "--outfilename", default="cms_hgg_datacard.txt",
                  help="Name of card to print (default: %default)")
parser.add_option("-p", "--procs", default="ggh,vbf,wh,zh,tth",
                  help="String list of procs (default: %default)")
parser.add_option("-c", "--cats", default="UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2",
                  help="Flashgg Categories (default: %default)")
parser.add_option("--uepsfilename", default="",
                  help="input files for calculating UEPS systematics; leave blank to use most recent set")
parser.add_option("--batch", default="LSF",
                  help="Batch system  (default: %default)")
parser.add_option("--photonCatScales", default="HighR9EE,LowR9EE,HighR9EB,LowR9EB",
                  help="String list of photon scale nuisance names - WILL NOT correlate across years (default: %default)")
parser.add_option("--photonCatScalesCorr", default="MaterialCentral,MaterialForward,FNUFEE,FNUFEB,ShowerShapeHighR9EE,ShowerShapeHighR9EB,ShowerShapeLowR9EE,ShowerShapeLowR9EB",
                  help="String list of photon scale nuisance names - WILL correlate across years (default: %default)")
parser.add_option("--photonCatSmears", default="HighR9EE,LowR9EE,HighR9EBRho,LowR9EBRho,HighR9EBPhi,LowR9EBPhi",
                  help="String list of photon smearing nuisance names - WILL NOT correlate across years (default: %default)")
parser.add_option("--photonCatSmearsCorr", default="",
                  help="String list of photon smearing nuisance names - WILL correlate across years (default: %default)")
#parser.add_option("--globalScales",default="NonLinearity:0.001,Geant4:0.0005,LightColl:0.0005,Absolute:0.0001",help="String list of global scale nuisances names with value separated by a \':\' - WILL NOT correlate across years (default: %default)")
parser.add_option("--globalScales", default="NonLinearity:0.001,Geant4:0.0005",
                  help="String list of global scale nuisances names with value separated by a \':\' - WILL NOT correlate across years (default: %default)")
parser.add_option("--globalScalesCorr", default="",
                  help="String list of global scale nuisances names with value separated by a \':\' - WILL correlate across years (default: %default)")
parser.add_option("--toSkip", default="",
                  help="proc:cat which are to skipped e.g ggH_hgg:11,qqH_hgg:12 etc. (default: %default)")
parser.add_option("--isMultiPdf", default=False, action="store_true")
parser.add_option("--submitSelf", default=False, action="store_true",
                  help="Tells script to submit itself to the batch")
parser.add_option("--justThisSyst", default="",
                  help="Only calculate the line corresponding to thsi systematic")
parser.add_option("--simplePdfWeights", default=False, action="store_true",
                  help="Condense pdfWeight systematics into 1 line instead of full shape systematic")
parser.add_option("--scaleFactors",
                  help="Scale factor for spin model pass as e.g. gg_grav:1.351,qq_grav:1.027")
parser.add_option("--quadInterpolate", type="int", default=0,
                  help="Do a quadratic interpolation of flashgg templates back to 1 sigma from this sigma. 0 means off (default: %default)")
parser.add_option("-m", "--mass", type="int", default=125,
                  help="Mass at which to calculate the systematic variations (default: %default)")
parser.add_option("--mA", type="int", default=5,
                  help="ALP Mass at which to calculate the systematic variations (default: %default)")
parser.add_option("--intLumi", type="float", default=3.71,
                  help="Integrated Lumi (default: %default)")
parser.add_option("--newGghScheme", default=False, action="store_true",
                  help="Use new WG1 scheme for ggH theory uncertainties")
parser.add_option("--doSTXS", default=False,
                  action="store_true", help="Use STXS Stage 0 processes")
(options, args) = parser.parse_args()
allSystList = []

###############################################################################


###############################################################################
## FILE I/O ###################################################################
###############################################################################
#inFile = r.TFile.Open(options.infilename)
outFile = open(options.outfilename, 'w')
###############################################################################

###############################################################################
## PROCS HANDLING & DICT ######################################################
###############################################################################
# convert flashgg style to combine style process
combProc = {'ggH_16_ele': 'ggH_16_ele_hza', 'ggH_16APV_ele': 'ggH_16APV_ele_hza', 'ggH_17_ele': 'ggH_17_ele_hza',
            'ggH_18_ele': 'ggH_18_ele_hza', 'ggH_16_mu': 'ggH_16_mu_hza', 'ggH_16APV_mu': 'ggH_16APV_mu_hza', 'ggH_17_mu': 'ggH_17_mu_hza',
            'ggH_18_mu': 'ggH_18_mu_hza', 'ggh': 'ggH_hZa', 'bkg_mass': 'bkg_mass'}
flashggProc = {'ggH_hZa': 'ggh', 'bkg_mass': 'bkg_mass'}
procId = {'ggH_16_ele_hza': -7, 'ggH_16APV_ele_hza': -6, 'ggH_17_ele_hza': -5, 'ggH_18_ele_hza': -4, 'ggH_16_mu_hza': -3, 'ggH_16APV_mu_hza': -2, 'ggH_17_mu_hza': -1, 'ggH_18_mu_hza': 0, 'bkg_mass': 1}
bkgProcs = ['bkg_mass']  # what to treat as background

# split procs vector
options.procs += ',bkg_mass'
options.procs = [combProc[p] for p in options.procs.split(',')]
options.toSkip = options.toSkip.split(',')
inFiles = options.infilename.split(',')
options.cats = options.cats.split(',')

print options.procs

###############################################################################

###############################################################################
## CATEGORISE TAGS FOR CONSIDERATION ##########################################
###############################################################################

###############################################################################

###############################################################################
## PHOTON SMEAR/SCALE SYSTEMATICS ## ##########################################
###############################################################################

###############################################################################

###############################################################################
## OPEN WORKSPACE AND EXTRACT INFO # ##########################################
sqrts = 13
print inFiles
norm = {}
sumEntries = {}
for iFile in range(len(inFiles)):

    inWS = WSTFileWrapper(inFiles[iFile], "wsig_13TeV")  # FIXED
    intL = inWS.var('IntLumi').getVal()
    SqrtS = inWS.var('IntLumi').getVal()
    print "[INFO] Get Intlumi from file, value : ", intL, " pb^{-1}", " SqrtS ", SqrtS

    fff = r.TFile.Open(inFiles[iFile])  # bing
    www = fff.Get("wsig_13TeV")  # bing
    dataNOMINAL = www.data("sig_%s_mass_m%d_%s" %
                           ("ggh", options.mass, "cat0"))  # bing
    NORM = www.function("hggpdfsmrel_13TeV_%s_%s_norm" %
                        ("ggh", "cat0"))  # bing

    # print "hggpdfsmrel_13TeV_%s_%s_norm"%("ggh", "cat0")
    # print "sig_%s_mass_m%d_%s"%("ggh",options.mass,"cat0")
    norm[options.procs[iFile]] = NORM.getVal()
    sumEntries[options.procs[iFile]] = dataNOMINAL.sumEntries()
    print "[[ INFO ]] norm: " + str(NORM.getVal()) + " ; sumEntries: " + str(dataNOMINAL.sumEntries())
###############################################################################


###############################################################################
## SHAPE SYSTEMATIC SETUP  ####################################################
###############################################################################
#file_ext = 'mva'
file_ext = 'data'

dataFile = 'CMS-HGG_mva_13TeV_multipdf.root'
bkgFile = dataFile
dataWS = 'multipdf'
bkgWS = 'multipdf'
#sigFile = 'CMS-HGG_%s_%dTeV_sigfit.root'%(file_ext,sqrts)
sigFile = 'CMS-HGG_sigfit_%s_$PROC_$CAT.root' % (file_ext)
# print "making sigfile " ,sigFile
sigWS = 'wsig_%dTeV' % (sqrts)
# file detaisl: for FLashgg always use unbinned signal and multipdf
fileDetails = {}
fileDetails['data_obs'] = [dataFile, dataWS, 'roohist_data_mass_$CHANNEL']
fileDetails['bkg_mass'] = [bkgFile, bkgWS,
                           'CMS_hgg_$CHANNEL_%dTeV_bkgshape' % sqrts]

years = ['16_ele', '16APV_ele', '17_ele', '18_ele', '16_mu', '16APV_mu', '17_mu', '18_mu']
for year in years:
    fileDetails['ggH_'+year+'_hza'] = [sigFile.replace('$PROC', "ggh").replace(
        ".root", "_"+year+".root"), sigWS, 'hggpdfsmrel_%dTeV_ggh_$CHANNEL' % sqrts]

###############################################################################

#############################################################################

###############################################################################
## GENERAL ANALYSIS SYSTEMATIC SETUP  #########################################
###############################################################################

#lumiSyst = 0.025  # Correct for Moriond17
lumiSyst = {'ggH_16_ele_hza':0.010, 'ggH_16APV_ele_hza':0.010, 'ggH_17_ele_hza':0.020, 'ggH_18_ele_hza':0.015, 'ggH_16_mu_hza':0.010, 'ggH_16APV_mu_hza':0.010, 'ggH_17_mu_hza':0.020, 'ggH_18_mu_hza':0.015}
# Printing Functions


def printLumiSyst():
    print '[INFO] Lumi...'
    for year in years[0:4]:
        year = year.rstrip("_"+year.split("_")[-1])
        #outFile.write('%-35s   lnN   ' % ('lumi_%dTeV' % sqrts))
        outFile.write('%-35s   lnN   ' % ('lumi_%dTeV_Uncorrelated_%s' % (sqrts,year)))
        for c in options.cats:
            for p in options.procs:
                if '%s:%s' % (p, c) in options.toSkip:
                    continue
                if p in bkgProcs:
                    outFile.write('- ')
                elif year+"_" in p:
                    outFile.write('%5.3f ' % (1.+lumiSyst[p]))
                else:
                    outFile.write('- ')
        outFile.write('\n')
    
    outFile.write('%-35s   lnN  1.006 1.006 1.009 1.020 1.006 1.006 1.009 1.020  -  ' % ('lumi_%dTeV_Correlated' % (sqrts)))
    outFile.write('\n')
    outFile.write('%-35s   lnN  - - 1.006 1.002 - - 1.006 1.002  -  ' % ('lumi_%dTeV_Correlated_1718' % (sqrts)))
    outFile.write('\n')
    outFile.write('\n')
###############################################################################

###############################################################################
##  DATACARD PREAMBLE TOOLS ###################################################
###############################################################################


def printPreamble():
    print '[INFO] Making Preamble...'
    outFile.write(
        'CMS-HZa datacard for parametric model - run2 %dTeV \n' % sqrts)
    outFile.write(
        'Auto-generated by flashggFinalFits/Datacard/makeParametricModelDatacardFLASHgg.py\n')
    outFile.write('Run with: combine\n')
    outFile.write('---------------------------------------------\n')
    outFile.write('imax *\n')
    outFile.write('jmax *\n')
    outFile.write('kmax *\n')
    outFile.write('---------------------------------------------\n')
    outFile.write('\n')
###############################################################################

###############################################################################
##  SHAPE SYSTEMATICS TOOLS ###################################################
###############################################################################


def printFileOptions():
    print '[INFO] File opts...'
    for typ, info in fileDetails.items():
        # print typ
        # print info
        for c in options.cats:
            file = info[0].replace('$CAT', '%s' % c)
            wsname = info[1]
            pdfname = info[2].replace('$CHANNEL', '%s' % c)
            if typ not in options.procs and typ != 'data_obs':
                continue
            #outFile.write('shapes %-10s %-15s %-30s %-30s\n'%(typ,'%s_%dTeV'%(c,sqrts),file.replace(".root","_%s_%s.root"%(typ,c)),wsname+':'+pdfname))
            outFile.write('shapes %-10s %-15s %-30s %-30s\n' %
                          (typ, '%s_%dTeV' % (c, sqrts), file, wsname+':'+pdfname))
    outFile.write('\n')
###############################################################################

###############################################################################
##  PROCESS/BIN LINES TOOLS ###################################################
###############################################################################


def printObsProcBinLines():
    print '[INFO] Rates...'

    outFile.write('%-15s ' % 'bin')
    for c in options.cats:
        outFile.write('%s_%dTeV ' % (c, sqrts))
    outFile.write('\n')

    outFile.write('%-15s ' % 'observation')
    for c in options.cats:
        outFile.write('-1 ')
    outFile.write('\n')

    outFile.write('%-15s ' % 'bin')
    for c in options.cats:
        for p in options.procs:
            if '%s:%s' % (p, c) in options.toSkip:
                continue
            outFile.write('%s_%dTeV ' % (c, sqrts))
    outFile.write('\n')

    outFile.write('%-15s ' % 'process')
    for c in options.cats:
        for p in options.procs:
            if '%s:%s' % (p, c) in options.toSkip:
                continue
            outFile.write('%s ' % p)
    outFile.write('\n')

    outFile.write('%-15s ' % 'process')
    for c in options.cats:
        for p in options.procs:
            if '%s:%s' % (p, c) in options.toSkip:
                continue
            outFile.write('%d ' % procId[p])
    outFile.write('\n')

    outFile.write('%-15s ' % 'rate')
    for c in options.cats:
        for p in options.procs:
            if '%s:%s' % (p, c) in options.toSkip:
                continue
            if p in bkgProcs:
                outFile.write('1.0 ')
            else:
                scale = 1.
                outFile.write('%8.7f ' % (sumEntries[p]/norm[p]))  # bing
    outFile.write('\n')
    outFile.write('\n')
###############################################################################

###############################################################################
##  DISCRETE SYST LINES TOOLS #################################################
###############################################################################


def printMultiPdf():
    if options.isMultiPdf:
        for c in options.cats:
            outFile.write('pdfindex_%s_%dTeV  discrete\n' % (c, sqrts))
###############################################################################


massList = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
if options.mA in massList:
    m = massList.index(options.mA)
else:
    for m_nom in massList:
        if m_nom<10:continue
        if abs(options.mA-m_nom)<=2:
            m = massList.index(m_nom)
            break
        else:
            continue


def printSyst():

    NormSysFiles = options.NormSys.split(',')
    if options.BDTSys is not None:
        BDTSysFiles = options.BDTSys.split(',')

    File_NormSys = open(NormSysFiles[0], 'r')

    NormSys_up = {}
    NormSys_dn = {}
    for line in File_NormSys.readlines()[1:]:
        if line.split()[0] == 'sample': 
            NormSys_names = line.split()[2:]
            NormSys_name = []
            for i in range(len(NormSys_names)/2):
                NormSys_up[NormSys_names[i*2][0:-3]] = {}
                NormSys_dn[NormSys_names[i*2][0:-3]] = {}
                NormSys_name.append(NormSys_names[i*2][0:-3])
                for year in years:
                    NormSys_up[NormSys_names[i*2][0:-3]]['ggH_'+year+'_hza'] = []
                    NormSys_dn[NormSys_names[i*2][0:-3]]['ggH_'+year+'_hza'] = []
                

    for f_norm in NormSysFiles:
        
        File_NormSys = open(f_norm, 'r')

        for line in File_NormSys.readlines():
            if line.split(':')[0]=="year":
                if line.split(':')[1].lstrip().rstrip('\n').split('_')[0] == "-2016":
                    year = line.split(':')[1].lstrip().rstrip('\n').replace("-2016","16APV")
                else:
                    year = line.split(':')[1].lstrip().rstrip('\n').lstrip('20')
            elif line.split()[0] == 'sample': 
                continue
            else:
                for j in range(len(NormSys_name)):
                    NormSysUp = abs(float(line.split()[3:][1::2][j*2].split(":")[1].rstrip(")")))
                    NormSysDn = abs(float(line.split()[3:][1::2][j*2+1].split(":")[1].rstrip(")")))

                    NormSys_up[NormSys_name[j]]['ggH_'+year+'_hza'].append(1+NormSysUp)
                    NormSys_dn[NormSys_name[j]]['ggH_'+year+'_hza'].append(1-NormSysDn)


    #for f_BDT in BDTSysFiles:
    #    File_BDTSys = open(f_BDT, 'r')

    BDTSys_up = {}
    BDTSys_dn = {}
    if options.BDTSys is not None:
        File_BDTSys = open(BDTSysFiles[0], 'r')
        for line in File_BDTSys.readlines()[1:]:
            if line.split()[0] == 'sample': 
                BDTSys_names = line.split()[2:]
                BDTSys_name = []
                for i in range(len(BDTSys_names)/2):
                    BDTSys_up[BDTSys_names[i*2][0:-3]] = {}
                    BDTSys_dn[BDTSys_names[i*2][0:-3]] = {}
                    BDTSys_name.append(BDTSys_names[i*2][0:-3])
                    for year in years:
                        BDTSys_up[BDTSys_names[i*2][0:-3]]['ggH_'+year+'_hza'] = []
                        BDTSys_dn[BDTSys_names[i*2][0:-3]]['ggH_'+year+'_hza'] = []
                    

        for f_BDT in BDTSysFiles:
            
            File_BDTSys = open(f_BDT, 'r')

            for line in File_BDTSys.readlines():
                if line.split(':')[0]=="year":
                    if line.split(':')[1].lstrip().rstrip('\n').split('_')[0] == "-2016":
                        year = line.split(':')[1].lstrip().rstrip('\n').replace("-2016","16APV")
                    else:
                        year = line.split(':')[1].lstrip().rstrip('\n').lstrip('20')
                elif line.split()[0] == 'sample': 
                    continue
                else:
                    for j in range(len(BDTSys_name)):
                        BDTSysUp = abs(float(line.split()[2:][1::2][j*2].split(":")[1].rstrip(")")))
                        BDTSysDn = abs(float(line.split()[2:][1::2][j*2+1].split(":")[1].rstrip(")")))

                        #print BDTSys#[BDTSys_name[j]]
                        BDTSys_up[BDTSys_name[j]]['ggH_'+year+'_hza'].append(1+BDTSysUp)
                        BDTSys_dn[BDTSys_name[j]]['ggH_'+year+'_hza'].append(1-BDTSysDn)



    print '[INFO] Normalization Sys...'
    
    for sys in NormSys_up:
        if 'lep' in sys:
            for chan in ['ele','mu']:
                outFile.write('%-35s   lnN   ' % (sys.replace("lep",chan)))
                for c in options.cats:
                    for p in options.procs:
                        if '%s:%s' % (p, c) in options.toSkip:
                            continue
                        if p in bkgProcs:
                            outFile.write('- ')
                        elif chan in p:
                            outFile.write('%5.4f/%5.4f ' % (NormSys_up[sys][p][m], NormSys_dn[sys][p][m]))
                        else:
                            outFile.write('- ')
                outFile.write('\n')      
            
            #outFile.write('\n')
        else:
            outFile.write('%-35s   lnN   ' % (sys))
            for c in options.cats:
                for p in options.procs:
                    if '%s:%s' % (p, c) in options.toSkip:
                        continue
                    if p in bkgProcs:
                        outFile.write('- ')
                    else:
                        outFile.write('%5.4f/%5.4f ' % (NormSys_up[sys][p][m], NormSys_dn[sys][p][m]))
                outFile.write('\n')
        
    outFile.write('\n')
    print '[INFO] BDT Sys...'
    for sys in BDTSys_up:
        if 'lep' in sys:
            for chan in ['ele','mu']:
                outFile.write('%-35s   lnN   ' % (sys.replace("lep",chan)))
                for c in options.cats:
                    for p in options.procs:
                        if '%s:%s' % (p, c) in options.toSkip:
                            continue
                        if p in bkgProcs:
                            outFile.write('- ')
                        elif chan in p:
                            outFile.write('%5.4f/%5.4f ' % (BDTSys_up[sys][p][m], BDTSys_dn[sys][p][m]))
                        else:
                            outFile.write('- ')
                outFile.write('\n')
        else:
            outFile.write('%-35s   lnN   ' % (sys))
            for c in options.cats:
                for p in options.procs:
                    if '%s:%s' % (p, c) in options.toSkip:
                        continue
                    if p in bkgProcs:
                        outFile.write('- ')
                    else:
                        outFile.write('%5.4f/%5.4f ' % (BDTSys_up[sys][p][m], BDTSys_dn[sys][p][m]))
            outFile.write('\n')

    outFile.write('\n')

def printShapeSys():
    print '[INFO] Shape Sys...'
    ShapeSysFiles = options.sysinfilename.split(',')
    #print ShapeSysFiles

    for f in ShapeSysFiles[0:8:4]:
        file_sys = open(f, 'r')
        lines = file_sys.readlines()
        for l in range(-4, 0):
            lines[l] = lines[l].rstrip('\n')
            sys_shapes = filter(None, lines[l].split("\t"))
            #print sys_shapes
            #print sys_shapes
            sys_val = [s for s in sys_shapes[1:] if float(s) != 0.0]
            #print sys_val[0]
            line_sys = 'CMS_hgg_nuisance_' + sys_shapes[0] + '\tparam\t0.0\t' + str(1.0)

            outFile.write(line_sys+'\n')

    line = 'pdfindex_0_13TeV \t\t\t discrete'
    outFile.write(line+'\n')


def printRate():
    interpFiles = open(options.interp)
    lines = interpFiles.readlines()
    outFile.write('\n')
    for line in lines:
        line_rate = 'yield_' + line.split()[0] + ' rateParam cat0_13TeV ggH_' + line.split()[0] + '_hza ' + line.split()[int(options.mA)].rstrip(',').rstrip(']').lstrip('[')
        outFile.write(line_rate+'\n')

    outFile.write('\n')
    for line in lines:
        line_rate = 'nuisance edit freeze yield_' + line.split()[0]
        outFile.write(line_rate+'\n')



###############################################################################
## MAIN #######################################################################
###############################################################################
# __main__ here
# preamble


printPreamble()
# shape systematic files
printFileOptions()
# obs proc/tag bins
printObsProcBinLines()

printLumiSyst()

printSyst()

printShapeSys()

if not options.interp is None:
    printRate()


