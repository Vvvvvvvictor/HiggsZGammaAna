#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.22.2018
#
#
#
#
#
#
import os
from ROOT import TFile, TH1F, TH1, gROOT, TH2F
import ctypes
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import uproot
# from root_pandas import *
import json
from categorizer import *
from tqdm import tqdm
from pdb import set_trace

pd.options.mode.chained_assignment = None

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'two_jet', help = 'Region to process')
    parser.add_argument('-va', '--variable', action = 'store', choices = ['bdt', 'NN'], default = 'bdt', help = 'MVA variable to use')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('--nscanvbf', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 4, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 4, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-es', '--estimate', action = 'store', choices = ['fullSim', 'fullSimrw', 'data_sid'], default = 'fullSimrw', help = 'Method to estimate significance')
    parser.add_argument('-f', '--nfold', type = int, default = 1, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = 30, help='early stopping rounds.')

    parser.add_argument('-s', '--shield', action='store', type=int, default=-1, help='Which variables needs to be shield')

    return  parser.parse_args()

def gettingsig(region, variable, boundaries_VBF, boundaries, transform, estimate):

    nbin = len(boundaries[0]) + len(boundaries_VBF[0])

    yields = pd.DataFrame({'sig': [0.]*nbin,
                          'sig_err': [0.]*nbin,
                          'data_sid': [0.]*nbin,
                          'data_sid_err': [0.]*nbin,
                          'bkgmc_sid': [0.]*nbin,
                          'bkgmc_sid_err': [0.]*nbin,
                          'bkgmc_cen': [0.]*nbin,
                          'bkgmc_cen_err': [0.]*nbin,
                          'VBF': [0.]*nbin,
                          'VBF_err': [0.]*nbin,
                          'ggH': [0.]*nbin,
                          'ggH_err': [0.]*nbin})

    for category in ['sig', 'VBF', 'ggH', "data_sid", "bkgmc_sid", "bkgmc_cen"]:

        # Construct the file path
        file_path = f'/eos/home-j/jiehan/root/outputs/test/{region}/{"bkgmc" if "bkgmc" in category else "Data" if "data" in category else "VBF_M125" if category == "VBF" else "ggH_M125" if category == "ggH" else category}.root'

        # Open the ROOT file
        with uproot.open(file_path) as f:
            # Specify the TTree name (assuming it's called 'test')
            tree = f[region]

            # Specify the columns you want to read
            columns = [
                f"{variable}_score{'_t' if transform else ''}",
                f"vbf_score{'_t' if transform else ''}",
                "H_mass",
                "weight",
                "event",
            ]

            # Read the data in chunks
            for data in tqdm(tree.iterate(columns, library='pd', step_size=500000), desc=f"Loading {category}", bar_format="{desc}: {percentage:3.0f}%|{bar:20}{r_bar}"):
    
                if 'sid' in category:
                    data = data[(data.H_mass >= 100) & (data.H_mass <= 180) & ((data.H_mass < 120) | (data.H_mass > 130))]
                    data['w'] = data.weight

                else:
                    data = data[(data.H_mass >= 120) & (data.H_mass <= 130)]
                    data['w'] = data.weight
        
                for i in range(len(boundaries)):
        
                    data_s = data[data.event % len(boundaries) == i]
        
                    data_ss = data_s[data_s[f"vbf_score{'_t' if transform else ''}"] < boundaries_VBF[i][0]]
        
                    for j in range(len(boundaries[i])):
        
                        data_sss = data_ss[data_ss[f"{variable}_score{'_t' if transform else ''}"] >= boundaries[i][j]]
                        if j != len(boundaries[i]) - 1: data_sss = data_sss[data_sss[f"{variable}_score{'_t' if transform else ''}"] < boundaries[i][j+1]]
        
                        # set_trace()
                        yields[category][j] += sum(data_sss.w)
                        yields[category+'_err'][j] = np.sqrt(yields[category+'_err'][j]**2 + sum(data_sss.w**2))
        
        
                    for j in range(len(boundaries_VBF[i])):
        
                        data_ss = data_s[data_s[f"vbf_score{'_t' if transform else ''}"] >= boundaries_VBF[i][j]]
                        if j != len(boundaries_VBF[i]) - 1: data_ss = data_ss[data_ss[f"vbf_score{'_t' if transform else ''}"] < boundaries_VBF[i][j+1]]
        
                        yields[category][j + len(boundaries[i])] += sum(data_ss.w)
                        yields[category+'_err'][j + len(boundaries[i])] = np.sqrt(yields[category+'_err'][j]**2 + sum(data_ss.w**2))


    if estimate == "data_sid":
        yields['bkg'] = yields['data_sid']*0.20
        yields['bkg_err'] = yields['data_sid_err']*0.20
    elif estimate == "fullSimrw":
        yields['bkg'] = yields['data_sid']*yields['bkgmc_cen']/yields['bkgmc_sid']
        yields['bkg_err'] = yields['data_sid_err']*yields['bkgmc_cen']/yields['bkgmc_sid']
    elif estimate == "fullSim":
        yields['bkg'] = yields['bkgmc_cen']
        yields['bkg_err'] = yields['bkgmc_cen_err']
    
    zs = calc_sig(yields.sig, yields.bkg, yields.sig_err, yields.bkg_err)
    yields['z'] = zs[0]
    yields['u'] = zs[1]
    yields['VBF purity [%]'] = yields['VBF']/yields['sig']*100
    yields['s/b'] = yields['sig']/yields['bkg']

    print(yields)

    z = np.sqrt((yields['z']**2).sum())
    u = np.sqrt((yields['z']**2 * yields['u']**2).sum())/z

    print(f'Significance:  {z:.4f} +/- {abs(u):.4f}')

    z_ggH = np.sqrt((yields['z']**2)[:4].sum())
    z_VBF = np.sqrt((yields['z']**2)[4:].sum())
    print(z_ggH, z_VBF)

    return z, abs(u), yields

def projectionY(hist, x1, x2, histname, ybin=100, y1=0, y2=1):

    hist_p = TH1F(histname, histname, ybin, y1, y2)
    for y in range(ybin):
        val, err = hist_integral(hist, x1, x2, y+1, y+1)
        hist_p.SetBinContent(y+1, val)
        hist_p.SetBinError(y+1, err)
    return hist_p

def categorizing(region, variable, sigs, bkgs, nscan, nscanvbf, minN, transform, nbin, nvbf, floatB, n_fold, fold, earlystop, estimate):

    # get inputs
    f_sig = TFile('/eos/home-j/jiehan/root/outputs/test/%s/sig.root' % (region))
    t_sig = f_sig.Get(region)

    if estimate in ["fullSim", "fullSimrw"]:
        f_bkgmc = TFile('/eos/home-j/jiehan/root/outputs/test/%s/bkgmc.root' % (region))
        t_bkgmc = f_bkgmc.Get(region)
    if estimate in ["fullSimrw", "data_sid"]:
        f_data_sid = TFile('/eos/home-j/jiehan/root/outputs/test/%s/Data.root' % (region))
        t_data_sid = f_data_sid.Get(region)

    # filling signal histograms
    h_sig_raw = TH2F('h_sig_raw', 'h_sig_raw', nscanvbf, 0., 1., nscan, 0., 1.) 
    h_sig=TH1F('h_sig','h_sig',nscanvbf,0,1)

    t_sig.Draw(f"{variable}_score{'_t' if transform else ''}:vbf_score{'_t' if transform else ''}>>h_sig_raw", "weight*%f*((H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_sig.Draw(f"vbf_score{'_t' if transform else ''}>>h_sig", "weight*%f*((H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,n_fold,fold if n_fold != 1 else 1))

    # filling bkg histograms
    if estimate in ["fullSim", "fullSimrw"]:
        h_bkgmc_cen_raw = TH2F('h_bkgmc_cen_raw', 'h_bkgmc_cen_raw', nscanvbf, 0., 1., nscan, 0., 1.)
        h_bkgmc_cen = TH1F('h_bkgmc_cen', 'h_bkgmc_cen', nscanvbf, 0., 1.)
        t_bkgmc.Draw(f"{variable}_score{'_t' if transform else ''}:vbf_score{'_t' if transform else ''}>>h_bkgmc_cen_raw", "weight*%f*((H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
        t_bkgmc.Draw(f"vbf_score{'_t' if transform else ''}>>h_bkgmc_cen", "weight*%f*((H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    if estimate in ["fullSimrw"]:
        h_bkgmc_sid_raw = TH2F('h_bkgmc_sid_raw', 'h_bkgmc_sid_raw', nscanvbf, 0., 1., nscan, 0., 1.)
        h_bkgmc_sid = TH1F('h_bkgmc_sid', 'h_bkgmc_sid', nscanvbf, 0., 1.)
        t_bkgmc.Draw(f"{variable}_score{'_t' if transform else ''}:vbf_score{'_t' if transform else ''}>>h_bkgmc_sid_raw", "weight*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
        t_bkgmc.Draw(f"vbf_score{'_t' if transform else ''}>>h_bkgmc_sid", "weight*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    if estimate in ["fullSimrw", "data_sid"]:
        h_data_sid_raw = TH2F('h_data_sid_raw', 'h_data_sid_raw', nscanvbf, 0., 1., nscan, 0., 1.)
        h_data_sid = TH1F('h_data_sid', 'h_data_sid', nscanvbf, 0., 1.)
        t_data_sid.Draw(f"{variable}_score{'_t' if transform else ''}:vbf_score{'_t' if transform else ''}>>h_data_sid_raw", "weight*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
        t_data_sid.Draw(f"vbf_score{'_t' if transform else ''}>>h_data_sid", "weight*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=122&&H_mass<=128)&&(event%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    fitboundary = 0.5
    fitboundary_g = 0.5

    if estimate == "data_sid":
        h_data_sid.Scale(0.20)
        cgz = categorizer(h_sig, h_data_sid)
    elif estimate == "fullSimrw":
        cgz = categorizer(h_sig, h_data_sid, h_bkg_rw_num=h_bkgmc_cen, h_bkg_rw_den=h_bkgmc_sid)
    elif estimate == "fullSim":
        cgz = categorizer(h_sig, h_bkgmc_cen)
    cgz.smooth(int(fitboundary * nscanvbf + 1), nscanvbf)

    zmax, boundaries_VBF, boundaries = 0, -1, -1
    for vb in tqdm(range(5, 95), ncols=90):

        boundaries_VBF, zv, uv = cgz.fit(vb, nscanvbf, nvbf, minN=minN, earlystop=earlystop)
        if boundaries_VBF == -1: break 

        h_sig_g = projectionY(h_sig_raw, 1, vb-1, f"h_sig_g_{vb}", nscan, 0, 1)
        if estimate in ["fullSim", "fullSimrw"]:
            h_bkgmc_cen_g = projectionY(h_bkgmc_cen_raw, 1, vb-1, f"h_bkgmc_cen_g_{vb}", nscan, 0, 1)
        if estimate in ["fullSimrw"]:
            h_bkgmc_sid_g = projectionY(h_bkgmc_sid_raw, 1, vb-1, f"h_bkgmc_sid_g_{vb}", nscan, 0, 1)
        if estimate in ["fullSimrw", "data_sid"]:
            h_data_sid_g = projectionY(h_data_sid_raw, 1, vb-1, f"h_data_sid_g_{vb}", nscan, 0, 1)

        if estimate == "data_sid":
            h_data_sid_g.Scale(0.20)
            cgz_g = categorizer(h_sig_g, h_data_sid_g)
        elif estimate == "fullSimrw":
            cgz_g = categorizer(h_sig_g, h_data_sid_g, h_bkg_rw_num=h_bkgmc_cen_g, h_bkg_rw_den=h_bkgmc_sid_g)
        elif estimate == "fullSim":
            cgz_g = categorizer(h_sig_g, h_bkgmc_cen_g)

        h_sig_g.Delete()
        if estimate in ["fullSimrw", "data_sid"]: h_data_sid_g.Delete()
        if estimate in ["fullSimrw"]: h_bkgmc_sid_g.Delete()
        if estimate in ["fullSim", "fullSimrw"]: h_bkgmc_cen_g.Delete()

        cgz_g.smooth(int(fitboundary_g * nscan + 1), nscan)
        boundaries, zg, ug = cgz_g.fit(1, nscan, nbin, minN=minN, floatB=floatB, earlystop=earlystop)
        if boundaries == -1: continue
 
        z = sqrt(zv**2 + zg**2)
        if z >= zmax:
            zmax, boundaries_VBF_max, boundaries_max = z, boundaries_VBF, boundaries

    boundaries_VBF_values = [(i-1.)/nscanvbf for i in boundaries_VBF_max]
    boundaries_values = [(i-1.)/nscan for i in boundaries_max]

    print('=========================================================================')
    print('Fold number ', fold)
    print('VBF boundaries: ', boundaries_VBF_values)
    print('ggF boundaries: ', boundaries_values)
    print('Total significane by fit: ', zmax)


    # get the significance estimated from the unfitted histogram
    nbkg = []
    for i in range(len(boundaries_VBF_max)):
        if estimate in ["fullSim", "fullSimrw"]:
            nbkgmc_cen, dbkgmc_cen = hist_integral(h_bkgmc_cen_raw, boundaries_VBF_max[i], boundaries_VBF_max[i+1]-1 if i != len(boundaries_VBF_max)-1 else nscanvbf+1, 1, nscan+1)
        if estimate in ["fullSimrw"]:
           nbkgmc_sid, dbkgmc_sid = hist_integral(h_bkgmc_sid_raw, boundaries_VBF_max[i], boundaries_VBF_max[i+1]-1 if i != len(boundaries_VBF_max)-1 else nscanvbf+1, 1, nscan+1)
        if estimate in ["fullSimrw", "data_sid"]:
            ndata_sid, ddata_sid = hist_integral(h_data_sid_raw, boundaries_VBF_max[i], boundaries_VBF_max[i+1]-1 if i != len(boundaries_VBF_max)-1 else nscanvbf+1, 1, nscan+1)
        if estimate == "data_sid":
            nbkg.append((ndata_sid*0.20, ddata_sid*0.20))
        elif estimate == "fullSimrw":
            nbkg.append((ndata_sid*nbkgmc_cen/nbkgmc_sid, ddata_sid*nbkgmc_cen/nbkgmc_sid))
        elif estimate == "fullSim":
            nbkg.append((nbkgmc_cen, dbkgmc_cen))
    for i in range(len(boundaries_max)):
        if estimate in ["fullSim", "fullSimrw"]:
            nbkgmc_cen, dbkgmc_cen = hist_integral(h_bkgmc_cen_raw, 1, boundaries_VBF_max[0]-1, boundaries_max[i], boundaries_max[i+1]-1 if i != len(boundaries_max)-1 else nscan+1)
        if estimate in ["fullSimrw"]:
           nbkgmc_sid, dbkgmc_sid = hist_integral(h_bkgmc_sid_raw, 1, boundaries_VBF_max[0]-1, boundaries_max[i], boundaries_max[i+1]-1 if i != len(boundaries_max)-1 else nscan+1)
        if estimate in ["fullSimrw", "data_sid"]:
            ndata_sid, ddata_sid = hist_integral(h_data_sid_raw, 1, boundaries_VBF_max[0]-1, boundaries_max[i], boundaries_max[i+1]-1 if i != len(boundaries_max)-1 else nscan+1)
        if estimate == "data_sid":
            nbkg.append((ndata_sid*0.20, ddata_sid*0.20))
        elif estimate == "fullSimrw":
            nbkg.append((ndata_sid*nbkgmc_cen/nbkgmc_sid, ddata_sid*nbkgmc_cen/nbkgmc_sid))
        elif estimate == "fullSim":
            nbkg.append((nbkgmc_cen, dbkgmc_cen))

    nsig = []
    for i in range(len(boundaries_VBF_max)):
        nsig.append(hist_integral(h_sig_raw, boundaries_VBF_max[i], boundaries_VBF_max[i+1]-1 if i != len(boundaries_VBF_max)-1 else nscanvbf+1, 1, nscan+1))
    for i in range(len(boundaries_max)):
        nsig.append(hist_integral(h_sig_raw, 1, boundaries_VBF_max[0]-1, boundaries_max[i], boundaries_max[i+1]-1 if i != len(boundaries_max)-1 else nscan+1))

    zs = [ calc_sig(nsig[i][0], nbkg[i][0], nsig[i][1], nbkg[i][1]) for i in range(len(nsig))]
    z, u = sum_z(zs)

    print('Significance from raw event yields in categorization set: ', z)
    print('=========================================================================')

    return boundaries_VBF_max, boundaries_VBF_values, boundaries_max, boundaries_values, zmax, z

def hist_integral(hist,i,j,k=-999,l=-999):
    err = ctypes.c_double() #Double()
    if i>j or k>l:
        n=0
        err=0
    elif k == -999 and l == -999: n = hist.IntegralAndError(i,j,err)
    else: n = hist.IntegralAndError(i,j,k,l,err)
    return n, err.value

def main():

    gROOT.SetBatch(True)
    TH1.SetDefaultSumw2(1)

    args=getArgs()
    shield = args.shield

    sigs = ['ggH_M125','VBF_M125'] #,'WminusH','WplusH','ZH','ttH']

    bkgs = ["ZGToLLG", "DYJetsToLL", "EWKZ2J"] #, "WGToLNuG", "ZG2JToG2L2J", "TT", "TTGJets", "TGJets", "ttWJets", "ttZJets", "WW", "WZ", "ZZ"]

    region = args.region

    variable = args.variable

    if not args.skip:
        siglist=''
        for sig in sigs:
            if os.path.isfile('/eos/home-j/jiehan/root/outputs/test/%s/%s.root'% (region,sig)): siglist+=' /eos/home-j/jiehan/root/outputs/test/%s/%s.root'% (region,sig)
        os.system("hadd -f /eos/home-j/jiehan/root/outputs/test/%s/sig.root"%(region)+siglist)

    if not args.skip:
        bkglist=''
        for bkg in bkgs:
            if os.path.isfile('/eos/home-j/jiehan/root/outputs/test/%s/%s.root'% (region,bkg)): bkglist+=' /eos/home-j/jiehan/root/outputs/test/%s/%s.root'% (region,bkg)
        os.system("hadd -f /eos/home-j/jiehan/root/outputs/test/%s/bkgmc.root"%(region)+bkglist)

    nscan=args.nscan
    nscanvbf = args.nscanvbf

    n_fold = args.nfold
    boundaries=[]
    boundaries_values=[]
    boundaries_VBF=[]
    boundaries_VBF_values=[]
    smaxs = []
    smaxs_raw = []
    for j in range(n_fold):
        out = categorizing(region, variable, sigs, bkgs, nscan, nscanvbf, args.minN, args.transform, args.nbin if not args.floatB else args.nbin + 1, args.vbf, args.floatB, n_fold, j, args.earlystop, estimate=args.estimate)
        boundaries_VBF.append(out[0])
        boundaries_VBF_values.append(out[1])
        boundaries.append(out[2])
        boundaries_values.append(out[3])
        smaxs.append(out[4])
        smaxs_raw.append(out[5])

    smax = sum(smaxs)/n_fold
    print('Averaged significance: ', smax)
    smax_raw = sum(smaxs_raw)/n_fold
    print('Averaged raw significance: ', smax_raw)

    s, u, yields = gettingsig(region, variable, boundaries_VBF_values, boundaries_values, args.transform, estimate=args.estimate)

    outs={}
    outs['boundaries'] = boundaries
    outs['boundaries_values'] = boundaries_values
    outs['boundaries_VBF'] = boundaries_VBF
    outs['boundaries_VBF_values'] = boundaries_VBF_values
    outs['smax'] = smax
    outs['smax_raw'] = smax_raw
    outs['significance'] = s
    outs['Delta_significance'] = u
    outs['nscan'] = nscan
    outs['nscanvbf'] = nscanvbf
    outs['transform'] = args.transform
    outs['floatB'] = args.floatB
    outs['nfold'] = n_fold
    outs['minN'] = args.minN
    outs['fine_tuned'] = False
    outs['variable'] = variable
    outs['estimate'] = args.estimate

    print(outs,'\n==================================================\n')

    '''
    if not os.path.isdir('significances/%s'%region):
        print(f'INFO: Creating output folder: "significances/{region}"')
        os.makedirs("significances/%s"%region)
    '''

    with open('/eos/home-j/jiehan/root/outputs/test/significances/%d_%s_2D_%d_%d.json' % (shield+1, region, args.vbf, args.nbin), 'w') as json_file:
        # json.dump(outs, json_file) 
        for i in list(yields['z']):
            json_file.write('\n{:.4f}'.format(i))
        json_file.write('\n{:.4f}'.format(s))

if __name__ == '__main__':
    main()
