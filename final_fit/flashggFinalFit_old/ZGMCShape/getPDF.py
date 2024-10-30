from ROOT import *

file_data = TFile("./InputData/two_jet_data/HZGamma_data_bkg_workspace_cat0.root")
file_bkg = TFile("./InputData/HZGamma_data_bkg_workspace_ZG.root")

ws_name = 'CMS_hzg_workspace'
inWS_data = file_data.Get(ws_name)
inWS_bkg = file_bkg.Get(ws_name)

nbins = 65
inWS_data.var("CMS_hzg_mass").setBins(nbins)
mass_ = inWS_bkg.var("CMS_hzg_mass")
mass_.SetTitle("m_{Zg}")
mass_.setUnit("GeV")
mass_.setRange(100.0, 180.0)
mass_.setBins(nbins)
intLumi_ = inWS_data.var("IntLumi")


datasets_data = inWS_data.data("data_mass_cat0")
datasets_bkg = inWS_bkg.data("data_mass_cat0")
mass_.setBins(nbins)
datasets_bkg.Print()
datasets_data.Print()


hist_bkg = datasets_bkg.binnedClone()
hist_bkg.Print()
mass_.Print()

ArgSet = RooArgSet("args")
ArgSet.add(mass_)
histpdf_bkg = RooHistPdf("histpdf_bkg", "histpdf_bkg", ArgSet, hist_bkg, 0)

canv = TCanvas("c1","c1",800,1000)
frame = mass_.frame()

datasets_data.plotOn(frame)
histpdf_bkg.plotOn(frame)

frame.Draw()
canv.SaveAs('DataMC_woweight_ZG.png')

##########################

nbins = 60000
mass_.setBins(nbins)

mass_.setRange(100.0, 180.0)
#hist_bkg_toy = histpdf_bkg.generate(ArgSet,50000)
hist_bkg_toy = histpdf_bkg.generate(ArgSet,1000000000)
hist_bkg_toy.Print()

kest1 = RooKeysPdf("ZGMCShape", "ZGMCShape", mass_, hist_bkg_toy, RooKeysPdf.MirrorBoth,1)

canv = TCanvas("c1","c1",800,1000)
mass_.setRange(105.0, 170.0)
frame = mass_.frame(100)

datasets_bkg.plotOn(frame)
kest1.plotOn(frame)

frame.Draw()
#canv.SaveAs('bkg_smooth.png')

w = RooWorkspace("w", "workspace")
getattr(w,'import')(kest1)
w.Print()
w.writeToFile("ZGMCShape.root")

f = TFile("ZGMCShape.root")
w = f.Get("w")
w.Print()
model = w.pdf("ZGMCShape")
x = w.var("CMS_hzg_mass")
model.Print("t")
canv = TCanvas("c1","c1",800,1000)
xframe = x.frame()
model.plotOn(xframe)
xframe.Draw()
canv.SaveAs('workspace.png')