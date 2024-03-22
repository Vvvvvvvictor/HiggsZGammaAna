from ROOT import *
file_bkg = TFile("./InputData/HZGamma_data_bkg_workspace_ZG.root")

ws_name = 'CMS_hzg_workspace'
inWS_bkg = file_bkg.Get(ws_name)

nbins = 65
mass_ = inWS_bkg.var("CMS_hzg_mass")
mass_.SetTitle("m_{ll\gamma}")
mass_.setUnit("GeV")
mass_.setRange(100.0, 180.0)
mass_.setBins(nbins)

datasets_bkg = inWS_bkg.data("data_mass_cat0")

f = TFile("ZGMCShape.root")
w = f.Get("w")
w.Print()
model = w.pdf("ZGMCShape")
x = w.var("CMS_hzg_mass")
x.SetTitle("m_{ll\gamma}")
x.setUnit("GeV")
x.setRange(100.0, 180.0)
x.setBins(nbins)
model.Print("t")

legend = TLegend(0.5,0.65,0.88,0.88)
legend.SetFillStyle(1001)
legend.SetBorderSize(1)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.SetShadowColor(0)
legend.SetFillColor(kWhite)

canv = TCanvas("c1","c1",800,1000)
xframe = x.frame()
datasets_bkg.plotOn(xframe)
model.plotOn(xframe)
xframe.Draw()
legend.AddEntry(datasets_bkg,"SM ZG MC")
legend.AddEntry(model,"Core Function","l")
legend.Draw("same")
canv.SaveAs('CoreFunction.png')
canv.SaveAs('CoreFunction.pdf')