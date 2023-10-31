from ROOT import *

path = '/afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_runII/M'

masses = [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]#[7,8,9,10,11,12,13,14,15]

for m in masses:
    files = {}
    files[0] = TFile(path+str(m)+'/higgsCombine_M{0}_0.025.HybridNew.mH125.quant0.025.root'.format(m))
    files[1] = TFile(path+str(m)+'/higgsCombine_M{0}_0.16.HybridNew.mH125.quant0.160.root'.format(m))
    files[2] = TFile(path+str(m)+'/higgsCombine_M{0}_0.5.HybridNew.mH125.quant0.500.root'.format(m))
    files[3] = TFile(path+str(m)+'/higgsCombine_M{0}_0.84.HybridNew.mH125.quant0.840.root'.format(m))
    files[4] = TFile(path+str(m)+'/higgsCombine_M{0}_0.975.HybridNew.mH125.quant0.975.root'.format(m))
    files[5] = TFile(path+str(m)+'/higgsCombine_M{0}_Observed.HybridNew.mH125.root'.format(m))
    files[6] = TFile(path+str(m)+'/higgsCombine_observed.AsymptoticLimits.mH125.root')

    limits = []

    for f in files:
        mychain = files[f].Get('limit')
        entries = mychain.GetEntries()
        for jentry in range(entries):
            nb = mychain.GetEntry(jentry)
            limits.append(mychain.limit)
        files[f].Close()
    
    print '*******************'
    print '***** m(a)={0} ******'.format(m)
    print '*******************'
    print '     full CLs        Asymptotic'
    print '0    ',limits[0], '      ',limits[6]
    print '1    ',limits[1], '      ',limits[7]
    print '2    ',limits[2], '      ',limits[8]
    print '3    ',limits[3], '      ',limits[9]
    print '4    ',limits[4], '      ',limits[10]
    print '5    ',limits[5], '      ',limits[11]
    #print 'ALP mass: ', m, 'full CLs limits: ', limits[0:6], 'Asymptotic limits: ', limits[6:]