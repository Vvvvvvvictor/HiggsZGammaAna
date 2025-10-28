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
import os
from ROOT import TFile, TH1F, TH1, gROOT, TH2, TH2F
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import uproot
# from root_pandas import *
import json
from categorizer import *
from pdb import set_trace

pd.options.mode.chained_assignment = None

hzgml_PATH="/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/hzgml"

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH'], default = 'zero_jet', help = 'Region to process')
    parser.add_argument('-i', '--input', action = 'store', default = '/eos/home-j/jiehan/root/outputs', help = 'Path of root file for categorization')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 10, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = int, default = 10, help = 'minimum number of events in mass window')
    parser.add_argument('--minSB', type = int, default = 10, help = 'minimum number of signal and background events')
    parser.add_argument('-v', '--variable', action = 'store', choices = ['bdt', 'NN'], default = 'bdt', help = 'MVA variable to use')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'se validation samples for categroization')
    parser.add_argument('-t', '--transform', type = bool, default = False, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-es', '--estimate', action = 'store', choices = ['fullSim', 'fullSimrw', 'data_sid'], default = 'fullSimrw', help = 'Method to estimate significance')

    parser.add_argument('-f', '--nfold', type = int, default = 1, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = -1, help='early stopping rounds.')

    parser.add_argument('-s', '--shield', action='store', type=int, default=-1, help='Which variables needs to be shielded')
    parser.add_argument('-a', '--add', action='store', type=int, default=-1, help='Which variables needs to be added')
    
    parser.add_argument('--score_min', type=float, default=0.0, help='Minimum boundary value for BDT/NN score range')
    parser.add_argument('--score_max', type=float, default=1.0, help='Maximum boundary value for BDT/NN score range')

    return  parser.parse_args()

def gettingsig(input_path, region, variable, boundaries, transform, estimate):

    nbin = len(boundaries[0])

    yields = pd.DataFrame({'sig': [0.]*nbin,
                          'sig_err': [0.]*nbin,
                          'sig_tot': [0.]*nbin,
                          'sig_tot_err': [0.]*nbin,
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

    for category in ['sig', 'VBF', 'ggH', "data_sid", "bkgmc_sid", "bkgmc_cen", "sig_tot"]:

        # for data in tqdm(read_root(f'{input_path}/{region}/{"bkgmc" if "bkgmc" in category else "data" if "data" in category else "sig" if "sig" in category else category}.root', key=region, columns=[f"{variable}_score{'_t' if transform else ''}", 'H_mass', 'weight_corr', 'event'], chunksize=500000), desc=f'Loading {category}', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
        # Define the file path
        file_path = f'{input_path}/{region}/{"bkgmc" if "bkgmc" in category else "Data" if "data" in category else "sig" if "sig" in category else "VBF_M125" if "VBF" in category else "ggH_M125" if "ggH" in category else category}.root'

        # Open the file
        file = uproot.open(file_path)

        # Define the columns to read
        columns = [f"{variable}_score{'_t' if transform else ''}", 'H_mass', 'weight_corr', 'event']

        # Iterate over the data in chunks
        for data in tqdm(file[region].iterate(columns, library='pd', step_size=500000), desc=f'Loading {category}', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
    
            if 'sid' in category:
                data = data[(data.H_mass >= 100) & (data.H_mass <= 180) & ((data.H_mass < 120) | (data.H_mass > 130))]
                data['w'] = data.weight_corr

            elif 'tot' in category:
                data = data[(data.H_mass >= 100) & (data.H_mass <= 180)]
                data['w'] = data.weight_corr

            else:
                data = data[(data.H_mass >= 120) & (data.H_mass <= 130)]
                data['w'] = data.weight_corr
           
            # print(data)
    
            for i in range(len(boundaries)):
    
                data_s = data[data.event % len(boundaries) == i]
    
                for j in range(len(boundaries[i])):
    
                    data_ss = data_s[data_s[f"{variable}_score{'_t' if transform else ''}"] >= boundaries[i][j]]
                    if j != len(boundaries[i]) - 1: data_ss = data_ss[data_ss[f"{variable}_score{'_t' if transform else ''}"] < boundaries[i][j+1]]
    
                    yields.loc[j, category] += data_ss.w.sum()
                    yields.loc[j, category+'_err'] = np.sqrt(yields[category+'_err'][j]**2 + (data_ss.w**2).sum())

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
    yields['u'] = zs[1] / zs[0]
    yields['VBF purity [%]'] = yields['VBF']/yields['sig']*100

    z = np.sqrt((yields['z']**2).sum())
    u = np.sqrt((zs[0]**2 * zs[1]**2).sum())/z

    print(f'Significance:  {z:.4f} +/- {abs(u):.4f}')

    return z, abs(u), yields

def categorizing(input_path, region, variable, sigs, bkgs, nscan, minN, transform, nbin, floatB, n_fold, fold, earlystop, minSB=10, estimate='fullSimrw', score_min=0.0, score_max=1.0):

    f_sig = TFile('%s/%s/sig.root' % (input_path, region))
    t_sig = f_sig.Get(region)
 
    if estimate in ["fullSim", "fullSimrw"]:
        f_bkgmc = TFile('%s/%s/bkgmc.root' % (input_path, region))
        t_bkgmc = f_bkgmc.Get(region)
    if estimate in ["fullSimrw", "data_sid"]:
        f_data_sid = TFile('%s/%s/Data.root' % (input_path, region))
        t_data_sid = f_data_sid.Get(region)


#################################################################################################

    # nmass = 10

    # h_sig = TH2F('h_sig','h_sig',nscan,score_min,score_max, nmass, 120, 130)
    # h_bkg = TH2F('h_bkg','h_bkg',nscan,score_min,score_max, nmass, 120, 130)

    # t_sig.Draw(f"H_mass:{variable}_score{'_t' if transform else ''}>>h_sig", "weight_corr*%f*((event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    # # filling bkg histograms
    # if estimate in ["fullSim", "fullSimrw"]:
    #     h_bkgmc_cen = TH2F('h_bkgmc_cen', 'h_bkgmc_cen', nscan, score_min, score_max, nmass, 120, 130)
    #     t_bkgmc.Draw(f"H_mass:{variable}_score{'_t' if transform else ''}>>h_bkgmc_cen", "weight_corr*%f*((event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    # if estimate in ["fullSimrw"]:
    #     h_bkgmc_sid = TH1F('h_bkgmc_sid', 'h_bkgmc_sid', nscan, score_min, score_max)
    #     t_bkgmc.Draw(f"{variable}_score{'_t' if transform else ''}>>h_bkgmc_sid", "weight_corr*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=120&&H_mass<=130)&&(event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    # if estimate in ["fullSimrw", "data_sid"]:
    #     h_data_sid = TH1F('h_data_sid', 'h_data_sid', nscan, score_min, score_max)
    #     t_data_sid.Draw(f"{variable}_score{'_t' if transform else ''}>>h_data_sid", "weight_corr*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=120&&H_mass<=130)&&(event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    # if estimate == "data_sid":
    #     h_data_sid.Scale(0.20)
    #     cgz = categorizer_shape(h_sig, h_data_sid)
    # elif estimate == "fullSimrw":
    #     cgz = categorizer_shape(h_sig, h_bkgmc_cen, h_bkg_rw_num=h_data_sid, h_bkg_rw_den=h_bkgmc_sid)
    # elif estimate == "fullSim":
    #     cgz = categorizer_shape(h_sig, h_bkgmc_cen)
    
    # cgz.smooth(50, nscan, SorB='B', function='Epoly2')
    # '''
    # uncomment upper line to fit a function to the BDT distribution. Usage: categorizer.smooth(left_bin_to_fit, right_bin_to_fit, SorB='S' (for signal) or 'B' (for bkg), function='Epoly2', printMessage=False (switch to "True" to print message))
    # TODO: first parameter must >= 1, 0 is banned.
    # '''

    # cgz.smooth_sim(1, SorB='B')

#################################################################################################  
 
    h_sig = TH1F('h_sig','h_sig',nscan,score_min,score_max)
    h_bkg = TH1F('h_bkg','h_bkg',nscan,score_min,score_max)

    t_sig.Draw(f"{variable}_score{'_t' if transform else ''}>>h_sig", "weight_corr*%f*((H_mass>=120&&H_mass<=130)&&(event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    # filling bkg histograms
    if estimate in ["fullSim", "fullSimrw"]:
        h_bkgmc_cen = TH1F('h_bkgmc_cen', 'h_bkgmc_cen', nscan, score_min, score_max)
        t_bkgmc.Draw(f"{variable}_score{'_t' if transform else ''}>>h_bkgmc_cen", "weight_corr*%f*((H_mass>=120&&H_mass<=130)&&(event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    if estimate in ["fullSimrw"]:
        h_bkgmc_sid = TH1F('h_bkgmc_sid', 'h_bkgmc_sid', nscan, score_min, score_max)
        t_bkgmc.Draw(f"{variable}_score{'_t' if transform else ''}>>h_bkgmc_sid", "weight_corr*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=120&&H_mass<=130)&&(event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    if estimate in ["fullSimrw", "data_sid"]:
        h_data_sid = TH1F('h_data_sid', 'h_data_sid', nscan, score_min, score_max)
        t_data_sid.Draw(f"{variable}_score{'_t' if transform else ''}>>h_data_sid", "weight_corr*%f*((H_mass>=100&&H_mass<=180)&&!(H_mass>=120&&H_mass<=130)&&(event%%%d!=%d))"%(n_fold/(n_fold-1) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    if estimate == "data_sid":
        h_data_sid.Scale(0.20)
        cgz = categorizer(h_sig, h_data_sid)
    elif estimate == "fullSimrw":
        cgz = categorizer(h_sig, h_bkgmc_cen, h_bkg_rw_num=h_data_sid, h_bkg_rw_den=h_bkgmc_sid)
    elif estimate == "fullSim":
        cgz = categorizer(h_sig, h_bkgmc_cen)
    
    # cgz.smooth(80, nscan, SorB='B', function='Epoly2')
    '''
    uncomment upper line to fit a function to the BDT distribution. Usage: categorizer.smooth(left_bin_to_fit, right_bin_to_fit, SorB='S' (for signal) or 'B' (for bkg), function='Epoly2', printMessage=False (switch to "True" to print message))
    TODO: first parameter must >= 1, 0 is banned.
    '''

    cgz.smooth_sim(1, SorB='B')

#################################################################################################
        
    bmax, zmax, umax = cgz.fit(1, nscan, nbin, minN=minN, minSB=minSB, floatB=floatB, earlystop=earlystop, pbar=True)
    umax = umax / zmax if zmax != 0 else 0
    print(bmax)
    boundaries = bmax
    boundaries_values = [np.round(score_min + (i-1)*(score_max-score_min)/nscan, 3) for i in boundaries]
    print('=========================================================================')
    print(f'Fold number {fold}')
    print(f'The maximal significance:  {zmax:.5f} +/- {umax:.5f}')
    print('Boundaries: ', boundaries_values)
    print('=========================================================================')

    return boundaries, boundaries_values, zmax
    
def main():

    gROOT.SetBatch(True)
    TH1.SetDefaultSumw2(1)

    args=getArgs()
    shield = args.shield
    add = args.add

    input_path = args.input + '/val'

    sigs = ['ggH_M125','VBF_M125'] #,'WminusH_M125','WplusH_M125','ZH_M125','ttH_M125', 'WH_M125']
    # sigs = ['ggH','VBF','WminusH','WplusH','ZH','ttH']

    # bkgs = ['data_fake', 'mc_med', 'mc_true']
    bkgs = ["ZGToLLG", "DYJetsToLL_ext", "EWKZ2J"] #, "ZG2JToG2L2J", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "WWG", "WZG", "ZZG", "ttZJets", "ttWJets"]
    if args.floatB and args.nbin == 16:
        print('ERROR: With floatB option, the maximun nbin is 15!!')
        quit()

    region = args.region

    nscan=args.nscan

    variable = args.variable

    if not args.skip:
        siglist=''
        for sig in sigs:
            if os.path.isfile('%s/%s/%s.root'% (input_path, region,sig)): siglist+=' %s/%s/%s.root'% (input_path, region,sig)
        os.system("hadd -f %s/%s/sig.root"%(input_path, region)+siglist)

    if not args.skip:
        bkglist=''
        for bkg in bkgs:
            if os.path.isfile('%s/%s/%s.root'% (input_path, region,bkg)): bkglist+=' %s/%s/%s.root'% (input_path,region,bkg)
        os.system("hadd -f %s/%s/bkgmc.root"%(input_path, region)+bkglist)


    n_fold = args.nfold
    boundaries=[]
    boundaries_values=[]
    smaxs = []
    for j in range(n_fold):
        bound, bound_value, smax = categorizing(input_path, region, variable, sigs, bkgs, nscan, args.minN, args.transform, args.nbin if not args.floatB else args.nbin + 1, args.floatB, n_fold, j, args.earlystop, minSB=args.minSB, estimate=args.estimate, score_min=args.score_min, score_max=args.score_max)
        boundaries.append(bound)
        boundaries_values.append(bound_value)
        smaxs.append(smax)

    # boundaries_values.append([0.00, 0.38, 0.71, 0.93])

    smax = sum(smaxs)/n_fold
    print('Averaged significance: ', smax)

    input_path = args.input + '/test'
    if not args.skip:
        siglist=''
        for sig in sigs:
            if os.path.isfile('%s/%s/%s.root'% (input_path, region,sig)): siglist+=' %s/%s/%s.root'% (input_path, region,sig)
        os.system("hadd -f %s/%s/sig.root"%(input_path, region)+siglist)
        bkglist=''
        for bkg in bkgs:
            if os.path.isfile('%s/%s/%s.root'% (input_path, region,bkg)): bkglist+=' %s/%s/%s.root'% (input_path,region,bkg)
        os.system("hadd -f %s/%s/bkgmc.root"%(input_path, region)+bkglist)

    s, u, yields = gettingsig(input_path, region, variable, boundaries_values, args.transform, estimate=args.estimate)

    outs={}
    outs['boundaries'] = boundaries
    outs['boundaries_values'] = boundaries_values
    outs['smax'] = smax
    outs['significance'] = s
    outs['Delta_significance'] = u
    outs['nscan'] = nscan
    outs['transform'] = args.transform
    outs['floatB'] = args.floatB
    outs['nfold'] = n_fold
    outs['minN'] = args.minN
    outs['minSB'] = args.minSB
    outs['fine_tuned'] = False
    outs['variable'] = variable
    outs['estimate'] = args.estimate
    outs['score_min'] = args.score_min
    outs['score_max'] = args.score_max
    outs.update(yields.to_dict())

    print(outs, '\n================================================\n')

    if not os.path.isdir('%s/significances/%s'%(input_path,region)):
        print(f'INFO: Creating output folder: "{input_path}/significances/{region}"')
        os.makedirs("%s/significances/%s"%(input_path,region))

    with open('%s/significances/bin_boundaries_1D_%s.txt'%(input_path,region), 'w') as json_file:
        for i in boundaries_values:
            # json_file.write('{:d} '.format(len(i)))
            for j in i:
                json_file.write('{:.2f} '.format(j))
        json_file.write('{:.2f}\n'.format(args.score_max))
        for i in list(yields['z']):
            json_file.write('{:.4f} '.format(i))
        json_file.write('%.4f %.4f' % (s, u))
    with open('%s/significances/%d_%d_%s_1D_%d.json' % (input_path, shield+1, add+1, region, args.nbin), 'w') as json_file:
        json.dump(outs, json_file, indent=4)
        # for i in list(yields['z']):
        #     json_file.write('{:.4f}\n'.format(i))
        # json_file.write('\n{:.4f}'.format(s))

    return



if __name__ == '__main__':
    main()
