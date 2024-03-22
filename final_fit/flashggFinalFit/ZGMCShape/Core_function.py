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
params_f1 = [0.01,0.1,0.5]
params_f2 = [0.01,0.1,0.5]
params_p1 = [0.02,0.01,-0.05,-0.01]
params_p2 = [0.02,0.01,-0.05,-0.01]
'''

params_f1 = [-100,1,10]
params_f2 = [-100,1,10]
params_p1 = [-100,1,10]
params_p2 = [-100,1,10]
model_final = {}

canv = TCanvas("c1","c1",800,1000)
frame = mass_.frame(65)

legend = TLegend(0.5,0.9,0.5,0.9)

model.plotOn(frame,RooFit.LineColor(kRed))

for p_f1 in params_f1:
    for p_f2 in params_f2:
        for p_p1 in params_p1:
            for p_p2 in params_p2:

                f1 = RooRealVar("f1","f1",p_f1,0.,1.)
                f2 = RooRealVar("f2","f2",p_f2,0.,1.)
                p1 = RooRealVar("p1","p1",p_p1,-.1,.1)
                p2 = RooRealVar("p2","p2",p_p2,-.1,.1)
                #p3 = RooRealVar("p3","p3",-0.005,-.1,.1)

                #exp = RooFormulaVar("exp", "exp","@3*exp(@1*@0)+@4*exp(@2*@0)",RooArgList(mass_,p1,p2,f1,f2))
                lau = RooGenericPdf("lau","lau", "@1*(@0)^(-4)+@2*(@0)^(-5)+@3*(@0)^(-3)+@4*(@0)^(-6)", RooArgList(mass_,p1,p2,f1,f2))
                #model_final["{}_{}_{}_{}".format(p_f1,p_f2,p_p1,p_p2)] = RooEffProd("model_final", "model_final", model, exp)
                model_final = RooEffProd("model_final", "model_final", model, lau)
                
                lau.plotOn(frame)
                legend.AddEntry(model_final, "f1:{}_f2:{}_p1:{}_p2:{}".format(p_f1,p_f2,p_p1,p_p2), "l")


'''
canv = TCanvas("c1","c1",800,1000)
frame = mass_.frame(65)

model.plotOn(frame,RooFit.LineColor(kRed))

for p_f1 in params_f1:
    for p_f2 in params_f2:
        for p_p1 in params_p1:
            for p_p2 in params_p2:
                model_final["{}_{}_{}_{}".format(p_f1,p_f2,p_p1,p_p2)].plotOn(frame)

'''
frame.Draw()
#legend.Draw()
canv.SaveAs('Core_functions_lau_lau.png')