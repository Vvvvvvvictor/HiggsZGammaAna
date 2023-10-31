from ROOT import *


lumi = {'16':16.81, '16APV':19.52, '17':41.48, '18':59.83}
years = ['16','16APV','17','18']
sysList = ['normal', 'phoScale_up', 'phoScale_dn', 'phoSmear_up', 'phoSmear_dn', 'lepScale_up', 'lepScale_dn', 'lepSmear_up', 'lepSmear_dn']
colors = [kRed, kBlue, kBlue, kGreen, kGreen, kYellow, kYellow, kMagenta, kMagenta]
a_masses = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
frac = {"coeff0":"frac_g0_ggh_cat0_rv_13TeV","coeff1":"hggpdfsmrel_13TeV_ggh_cat0_rv_13TeV_recursive_fraction_gaus_g1_ggh_cat0_rv_13TeV","coeff2":"hggpdfsmrel_13TeV_ggh_cat0_rv_13TeV_recursive_fraction_gaus_g2_ggh_cat0_rv_13TeV"}

param = {}
for a in a_masses:
    param[a] = {}
    for year in years:
        param[a][year] = {}

        file_param = open("../fit_parm_mu/paramDump_ggh_cat0_M"+str(a)+"_"+year+".txt",'r')
        file_param_root = TFile("../fit_parm_mu/paramDump_ggh_cat0_M"+str(a)+"_"+year+".root")
        file_param_root.wsig_13TeV.var("MH").setVal(125)
        #file_param_root.wsig_13TeV.func()

        for line in file_param.readlines()[1:-3]:
            param_name = line.split("ggh_cat0__rv_13TeV")[1].split()[0]
            param_value = line.split("ggh_cat0__rv_13TeV")[1].split()[4]
            
            param[a][year][param_name] = float(param_value)
    
            if param_name in frac:
                param[a][year][param_name] = file_param_root.wsig_13TeV.function(frac[param_name]).getVal()

param_float = {'Scale':'mean_float', 'Smear':'sigma_float'}

def sumGaussians(mass, param, scale):

    global MH, dm0, mean0, sigma0, gauss0, s0
    global dm1, mean1, sigma1, gauss1, coeff0, s1
    global dm2, mean2, sigma2, gauss2, coeff1, s2
    global dm3, mean3, sigma3, gauss3, coeff2, s3
    global mean_float, sigma_float

    MH = RooRealVar("MH","MH",125.0)
    mean_float = RooRealVar("mean_float","mean_float",0.0,-1.0,1.0)
    sigma_float = RooRealVar("sigma_float","sigma_float",0.0,-1.0,1.0)
    
    dm0 = RooRealVar("dm0","dm0",param['dm0'],-3.,3.)
    mean0 = RooFormulaVar("mean0","mean0","(@0+@1)*(1+@2)",RooArgList(MH,dm0,mean_float))
    s0 = RooRealVar("s0","s0",param['sigma0'],0.01,20.)
    sigma0 = RooFormulaVar("sigma0","sigma0","@0*(1+@1)",RooArgList(s0,sigma_float))
    gauss0 = RooGaussian("gaus0","gaus0",mass,mean0,sigma0)
    #coeff0 = RooRealVar("coeff0","coeff0",0.01,-0.005,0.005)

    dm1 = RooRealVar("dm1","dm1",param['dm1'],-3.,3.)
    #mean1 = RooFormulaVar("mean1","mean1","@0+@1",RooArgList(MH,dm1))
    #sigma1 = RooRealVar("sigma1","sigma1",param['sigma1'],0.01,20.)
    mean1 = RooFormulaVar("mean1","mean1","(@0+@1)*(1+@2)",RooArgList(MH,dm1,mean_float))
    s1 = RooRealVar("s1","s1",param['sigma1'],0.01,20.)
    sigma1 = RooFormulaVar("sigma1","sigma1","@0*(1+@1)",RooArgList(s1,sigma_float))
    gauss1 = RooGaussian("gaus1","gaus1",mass,mean1,sigma1)
    coeff0 = RooRealVar("coeff0","coeff0",param['coeff0'],0.,1.0)

    dm2 = RooRealVar("dm2","dm2",param['dm2'],-3.,3.)
    #mean2 = RooFormulaVar("mean2","mean2","@0+@1",RooArgList(MH,dm2))
    #sigma2 = RooRealVar("sigma2","sigma2",param['sigma2'],0.01,20.)
    mean2 = RooFormulaVar("mean2","mean2","(@0+@1)*(1+@2)",RooArgList(MH,dm2,mean_float))
    s2 = RooRealVar("s2","s2",param['sigma2'],0.01,20.)
    sigma2 = RooFormulaVar("sigma2","sigma2","@0*(1+@1)",RooArgList(s2,sigma_float))
    gauss2 = RooGaussian("gaus2","gaus2",mass,mean2,sigma2)
    coeff1 = RooRealVar("coeff1","coeff1",param['coeff1'],0.,1.0)

    dm3 = RooRealVar("dm3","dm3",param['dm3'],-3.,3.)
    #mean3 = RooFormulaVar("mean3","mean3","@0+@1",RooArgList(MH,dm3))
    #sigma3 = RooRealVar("sigma3","sigma3",param['sigma3'],0.01,20.)
    mean3 = RooFormulaVar("mean3","mean3","(@0+@1)*(1+@2)",RooArgList(MH,dm3,mean_float))
    s3 = RooRealVar("s3","s3",param['sigma3'],0.01,20.)
    sigma3 = RooFormulaVar("sigma3","sigma3","@0*(1+@1)",RooArgList(s3,sigma_float))
    gauss3 = RooGaussian("gaus3","gaus3",mass,mean3,sigma3)
    coeff2 = RooRealVar("coeff2","coeff2",param['coeff2'],0.,1.0);

    MH.setConstant(kTRUE)
    coeff0.setConstant(kTRUE)
    coeff1.setConstant(kTRUE)
    coeff2.setConstant(kTRUE)

    s0.setConstant(kTRUE)
    s1.setConstant(kTRUE)
    s2.setConstant(kTRUE)
    s3.setConstant(kTRUE)

    dm0.setConstant(kTRUE)
    dm1.setConstant(kTRUE)
    dm2.setConstant(kTRUE)
    dm3.setConstant(kTRUE)

    ### 1 for scale, 2 for smear, 0 for normal
    if scale == 1:
        sigma_float.setConstant(kTRUE)
    elif scale == 2:
        mean_float.setConstant(kTRUE)
    else:
        sigma_float.setConstant(kTRUE)
        mean_float.setConstant(kTRUE)

    func = RooAddPdf("nSum", "nSum", RooArgList(gauss0,gauss1,gauss2,gauss3), RooArgList(coeff0,coeff1,coeff2))
    

    #gauss0.Print()
    #func.Print()

    return func
    



def main():

    path_file = '/afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/'
    path_out = '/afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_sys_UL/shapeSys/fit_run2_UL_ScaleSmear_mu/ScaleSmear_plot/'
    ws_name = 'CMS_hza_workspace'

    for a_mass in a_masses:
        for year in years:
            file = TFile(path_file+'ALP_sig_Am'+str(a_mass)+'_Hm125_'+year+'_workspace.root')
            #file.ls()
            inWS = file.Get(ws_name)

            mass_ = inWS.var("CMS_hza_mass")
            mass_.SetTitle("m_{Za}")
            mass_.setUnit("GeV")
            mass_.setRange(110.0, 140.0)
            intLumi_ = inWS.var("IntLumi")

            data = {}
            model = {}
            fitResult = {}
            h = {}

            canv = TCanvas("c1","c1",800,1000)
            frame = mass_.frame()
            legend = TLegend(0.55,0.55,0.89,0.89)
            #legend.SetHeader("The Legend Title","C")

            print param[a_mass][year]
            

            for sys in sysList:
                if 'Scale' in sys:
                    s = 1
                elif 'Smear' in sys:
                    s = 2
                else:
                    s = 0
                data[sys] = inWS.data("ggh_125_13TeV_cat0_"+sys)
                model[sys] = sumGaussians(mass_, param[a_mass][year], s)
                model[sys].Print()
                data[sys].Print()
                fitResult[sys] = model[sys].fitTo(data[sys],RooFit.SumW2Error(kFALSE),RooFit.Save(kTRUE))
                #fitResult[sys] = gauss1.fitTo(data[sys],RooFit.SumW2Error(kFALSE),RooFit.Save(kTRUE))
                
                
                if sys == 'normal': data[sys].plotOn(frame)
                #data[sys].plotOn(frame, RooFit.LineColor(colors[sysList.index(sys)]))
                model[sys].plotOn(frame, RooFit.LineColor(colors[sysList.index(sys)]))
                #if sys == 'normal': data[sys].plotOn(frame)

                h[sys] = TH1F("h_"+sys,"h_"+sys,1,0,1)
                h[sys].SetLineColor(colors[sysList.index(sys)])
                legend.AddEntry(h[sys],sys,'L')
                

            #legend.Draw()
            frame.SetTitle("ScaleSmear_"+str(a_mass)+"_"+year)
            frame.Draw()
            legend.Draw()
            canv.SaveAs(path_out+"ScaleSmear_m"+str(a_mass)+"_"+year+".png")


            
            file = open(path_out+'ScaleSmear_m'+str(a_mass)+'_'+year+'.dat', 'a')
            title = '{0:15}\t'.format('#paramters')
            for sys in sysList:
                title = title + '{0:15}\t'.format(sys)

            file.write(title + '\n')

            for scsm in param_float:
                par = param_float[scsm]
                line = '#' + '{0:15}\t'.format(par)
                for sys in sysList:
                    if sys == 'normal':
                        if scsm == 'Scale':
                            par_val = fitResult['phoScale_up'].floatParsInit().find(par).getVal()
                        else:
                            par_val = fitResult['phoScale_up'].constPars().find(par).getVal()
                        par_val_norm = par_val
                    else:
                        if scsm in sys:
                            par_val = fitResult[sys].floatParsFinal().find(par).getVal()
                        else:
                            par_val = fitResult[sys].constPars().find(par).getVal()

                    line = line + '{0:.5f}\t\t\t'.format(par_val)
                file.write(line + '\n')

            

            
            sysType = ['pho', 'lep']
            
            for scsm in param_float:
                line = "photonCat" + scsm +"s="
                for p in sysType:
                    line = line + p + '_' + scsm + '_M' + str(a_mass)+'_ChanMu' + ','
                line = line.rstrip(',')
                file.write(line + '\n')

            L = ['diphotonCat=0', 'proc=ggh']
            for line in L:
                file.write(line + '\n')


            for scsm in param_float:
                for p in sysType:
                    line = p + '_' + scsm + '_M' + str(a_mass)+'_ChanMu' +'\t'

                    par_val_norm = -99.9
                    par_val_change = -99.9

                    par = param_float[scsm]
                    for sys in sysList:

                        if sys == 'normal':
                            if scsm == 'Scale':
                                par_val = fitResult['phoScale_up'].floatParsInit().find(par).getVal()
                            else:
                                par_val = fitResult['phoScale_up'].constPars().find(par).getVal()
                            par_val_norm = par_val
                        else:
                            if p not in sys: continue
                            if scsm in sys:
                                par_val = fitResult[sys].floatParsFinal().find(par).getVal()
                            else:
                                par_val = fitResult[sys].constPars().find(par).getVal()

                        if abs(par_val_norm - par_val) > par_val_change:
                            par_val_change = abs(par_val_norm - par_val)
                
                    
                    if scsm == 'Scale':
                        line = line + '{0:.5f}\t0.00000\t0.00000'.format(par_val_change)
                    else:
                        line = line + '\t0.00000\t{0:.5f}\t0.00000'.format(par_val_change)

                    file.write(line + '\n')


            

main()
