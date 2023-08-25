import os
import ROOT

#ROOT.PyConfig.IgnoreCommandLineOptions = True
#ROOT.gROOT.CloseFiles()
#ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2000;")
#ROOT.gROOT.SetMustClean(True)
ROOT.gStyle.SetOptStat("0")
ROOT.gStyle.SetTitleXSize(0.07)
ROOT.gStyle.SetTitleYSize(0.06)
ROOT.gStyle.SetPadColor(0)
ROOT.gStyle.SetCanvasColor(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.ROOT.EnableImplicitMT()

from optparse import OptionParser

print("\t===========================================================")
print("\tNOW!!! We will draw the most beautiful PICTURE in the world")
print("\tIt worth coutless $$$$")
print("\tLET'S GO!!!!!!")
print("\t===========================================================")

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

Legend_stack =["Fake photons1", "Fake photons2", "Fake photons3", "Fake photons4", "SM ZG", "VBS ZG", "TTG", "TG"]
Legend_line = ["data"]
Legend_pos = [0.75, 0.95]
ratio = 1.1

X_axis = [65, 0, 65]
x_name = "H_pt"
x_title = "Zero jets: P_[T]^[ll#gamma]"
input_dir = "/afs/cern.ch/user/j/jiehan/private/hmumuml"
input_stack_file = [["data_fake", "mc_true"], ["data_fake", "mc_true"], ["data_fake", "mc_true"], ["data_fake", "mc_true"], ["zg"], ["llajj"], ["ttg"], ["tg"]]
input_line_file = [["data"]]

print("\t===========================================================")
print("\tFinish basic setting")
print("\t===========================================================")

can = ROOT.Tcanvas("can", "", 800, 600)
pad1 = ROOT.TPad("c_1","",0,0,1,0.3)
pad2 = ROOT.TPad("c_2","", 0,0.25,1,0.95)





