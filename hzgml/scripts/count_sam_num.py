import uproot as up
import pandas as pd
import numpy as np
import uproot

def eta_sel(data, type):
        if type == 1:
            return (data.loc[:, 'gamma_eta']>-1.5) & (data.loc[:, 'gamma_eta']<1.5)
        if type == 0:
            return (data.loc[:, 'gamma_eta']<-1.5) | (data.loc[:, 'gamma_eta']>1.5)

def read_root(file_path, bar_id, lep_id, pt_range, prompt):
    if prompt:
        branches = uproot.open(file_path+':inclusive;1').arrays(['gamma_pt','gamma_eta', 'gamma_chiso',"gamma_mvaID",'weight','Z_lead_lepton_id', 'Z_sublead_lepton_pt', 'gamma_genPartFlav', 'n_iso_photons', 'weight'], library='pd')
    else:
        branches = uproot.open(file_path+':inclusive;1').arrays(['gamma_pt','gamma_eta', 'gamma_chiso',"gamma_mvaID",'weight','Z_lead_lepton_id', 'Z_sublead_lepton_pt', 'weight'], library='pd')

    com_cut =  branches.loc[:,'gamma_pt'] > 0

    bar_id_cut = eta_sel(branches, bar_id)
    lep_id_cut = branches.loc[:,'Z_lead_lepton_id'] == lep_id
    pt_cut = (branches.loc[:,'gamma_pt']<pt_range[1]) & (branches.loc[:,'gamma_pt']>=pt_range[0])

    com_cut = bar_id_cut & lep_id_cut & pt_cut & (branches.loc[:,'Z_sublead_lepton_pt'] > 15) & com_cut
    if prompt:
        com_cut = com_cut & (branches.loc[:,'gamma_genPartFlav'] == 1)

    if bar_id:
        low_mva = -0.4
        med_mva = 0.42
        med_chiso = 1.141
    else:
        low_mva = -0.59
        med_mva = 0.14
        med_chiso = 1.051

    low_cut = (branches.loc[:,'gamma_chiso']<med_chiso) & (branches.loc[:,'gamma_mvaID']<low_mva)
    med_cut = (branches.loc[:,'gamma_chiso']<med_chiso) & (branches.loc[:,'gamma_mvaID']>=low_mva) & (branches.loc[:,'gamma_mvaID']<med_mva)
    high_cut = (branches.loc[:,'gamma_mvaID']>=med_mva)

    low = branches.loc[com_cut & low_cut, 'weight'].sum()
    med = branches.loc[com_cut & med_cut, 'weight'].sum()
    high = branches.loc[com_cut & high_cut, 'weight'].sum()

    del branches
    # print(np.array([low, med, high]))
    return np.array([low[0], med[0], high[0]])

# if __name__ == '_main_':
path = '/afs/cern.ch/user/j/jiehan/private/hmumuml/skimmed_ntuples/'
signal = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
prompt = ["ZGToLLG", "TTGJets", "TGJets", "ZG2JToG2L2J", "DYJetsToLL"]
data = ["data"]
for bar_id in [0,1]:
    if bar_id:
        pt_ranges = [15,17,20,25,30,35,40,5000]
    else:
        pt_ranges = [15,17,20,25,30,40,5000]
    for lep_id in [11,13]:
        print(bar_id, " ", lep_id)
        for pt_i in range(len(pt_ranges)-1):
            pt_range = []
            pt_range.append(pt_ranges[pt_i])
            pt_range.append(pt_ranges[pt_i+1])
            # print(bar_id, " ", lep_id, " ", pt_range)
            
            sum_low = np.array([0., 0., 0.])
            sum_med = np.array([0., 0., 0.])
            sum_high = np.array([0., 0., 0.])
            for i in signal:
                sum_low += read_root(path+i+'/2017.root', bar_id, lep_id, pt_range, 0)
            # print("Signal samples! low:{:.2f}, medium:{:.2f}, high:{:.2f}".format(sum_low[0], sum_low[1], sum_low[2]))
            print("$[{:d},{:d})$ & {:.2f} & {:.2f} & {:.2f}".format(pt_range[0], pt_range[1], sum_low[0], sum_low[1], sum_low[2]), end='')
            
            for i in prompt:
                sum_med += read_root(path+i+'/2017.root', bar_id, lep_id, pt_range, 1)
            # print("Prompt MC samples! low:{:.2f}, medium:{:.2f}, high:{:.2f}".format(sum_med[0], sum_med[1], sum_med[2]))
            print(" & {:.2f} & {:.2f} & {:.2f}".format(sum_med[0], sum_med[1], sum_med[2]), end='')

            for i in data:
                sum_high += read_root(path+i+'/2017.root', bar_id, lep_id, pt_range, 0)
            # print("DATA samples! low:{:.2f}, medium:{:.2f}, high:{:.2f}".format(sum_high[0], sum_high[1], sum_high[2]))
            print(" & {:.0f} & {:.0f} & {:.0f} \\\\".format(sum_high[0], sum_high[1], sum_high[2]))