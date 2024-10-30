from ROOT import *

file_data = TFile("./InputData/two_jet_data/HZGamma_data_bkg_workspace_cat0.root")

nbins = 65
ws_name = 'CMS_hzg_workspace'
inWS_data = file_data.Get(ws_name)
mass_ = inWS_data.var("CMS_hzg_mass")
mass_.SetTitle("m_{Zg}")
mass_.setUnit("GeV")
mass_.setBins(nbins)
datasets_data = inWS_data.data("data_mass_cat0")

f = TFile("ZGMCShape.root")
w = f.Get("w")
w.Print()
model = w.pdf("ZGMCShape")
x = w.var("CMS_hzg_mass")
model.Print("t")

'''
f1 = RooRealVar("f1","f1",0.9,0.0001,0.9999)
f2 = RooRealVar("f2","f2",0.4,0.0001,0.9999)
f3 = RooRealVar("f3","f3",0.1,0.0001,0.9999)
p1 = RooRealVar("p1","p1",-0.0005,-.002, -0.000001)
p2 = RooRealVar("p2","p2",-0.001,-.029,-0.000001)
#p2 = RooRealVar("p2","p2",0.001,.000000029,1.)
p3 = RooRealVar("p3","p3",-0.005,-.3,0.)

exp1 = RooExponential("exp1", "exp1", mass_, p1)
exp2 = RooExponential("exp2", "exp2", mass_, p2)
exp3 = RooExponential("exp3", "exp3", mass_, p3)
'''

#exp = exp1
#exp = RooAddPdf("exp", "exp", exp1, exp2, f1)
#exp = RooAddPdf("exp", "exp", mass_, RooArgList(exp1,exp2,exp3), RooArgList(f1,f2))


###### exp
'''
f1 = RooRealVar("f1","f1",0.5,0.,1.)
f2 = RooRealVar("f2","f2",0.4,0.,1.)
p1 = RooRealVar("p1","p1",-0.01,-.1,.1)
p2 = RooRealVar("p2","p2",-0.001,-.1,.1)
p3 = RooRealVar("p3","p3",-0.005,-.1,.1)

#pdf_exp = RooExponential("exp","exp",mass_,p1)
#exp = RooFormulaVar("exp", "exp","exp(@1*@0)",RooArgList(mass_,p1))
exp = RooFormulaVar("exp", "exp","@3*exp(@1*@0)+(1-@3)*exp(@2*@0)",RooArgList(mass_,p1,p2,f1))
#pdf_exp = RooFormulaVar("pdf_exp","@4*exp(@1*@0)+@5*exp(@2*@0)+(1.-@4-@5)*exp(@3*@0)",RooArgSet(mass_,p1,p2,p3,f1,f2))
'''

'''
###### pow

f1 = RooRealVar("f1","f1",1.,-10.,10.)
f2 = RooRealVar("f2","f2",1.,-10.,10.)
f3 = RooRealVar("f3","f3",1.,-10.,10.)
p1 = RooRealVar("p1","p1",-1.,-10.,5.)
p2 = RooRealVar("p2","p2",-1.,-10.,5.)
p3 = RooRealVar("p3","p3",-1.,-10.,5.)

#pows = RooGenericPdf("pow","pow", "@2(@0)^(@1)", RooArgList(mass_,p1,f1))
pows = RooFormulaVar("pow","pow", "@2*(@0)^(@1)+@4*(@0)^(@3)", RooArgList(mass_,p1,f1,p2,f2))
#pows = RooGenericPdf("pow","pow", "@2*(@0)^(@1)+@4*(@0)^(@3)+@6*(@0)^(@5)", RooArgList(mass_,p1,f1,p2,f2,p3,f3))
'''

###### pow


f1 = RooRealVar("f1","f1",0.1,-100.,100.)
f2 = RooRealVar("f2","f2",0.5,-100.,100.)
f3 = RooRealVar("f3","f3",0.01,-100.,100.)
f4 = RooRealVar("f4","f4",0.5,-100.,100.)

#lau = RooGenericPdf("lau","lau", "@1*(@0)^(-4)+@2*(@0)^(-5)", RooArgList(mass_,f1,f2))
lau = RooGenericPdf("lau","lau", "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)", RooArgList(mass_,f1,f2,f3))
#lau = RooGenericPdf("lau","lau", "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)+@4*(@0)^(-6)", RooArgList(mass_,f1,f2,f3,f4))


#model_final = RooEffProd("model_final", "model_final", model, exp)
model_final = RooEffProd("model_final", "model_final", model, lau)
#model_final = RooProdPdf("model_final", "model_final", RooArgList(model, exp), 1)

model_final.fitTo(datasets_data)

canv = TCanvas("c1","c1",800,1000)
frame = mass_.frame(65)


#histpdf_bkg.plotOn(frame)
datasets_data.plotOn(frame)
model_final.plotOn(frame,RooFit.LineColor(kRed))
#model.plotOn(frame,RooFit.LineColor(kRed),RooFit.Components("bkgxGau"),RooFit.LineStyle(kDashed))
#datasets_data.plotOn(frame)
#datasets_bkg.plotOn(frame)

frame.Draw()
canv.SaveAs('fit_results_lau.png')