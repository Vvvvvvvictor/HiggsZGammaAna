from numpy import *

a_masses = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
years = ['16','16APV','17','18']
#a_masses = [1]
#years = ['16']

path = '/afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_sys_UL/shapeSys/'

obj = ['pho', 'lep']
val_pho_smear = {}
val_pho_scale = {}
val_lep_smear = {}
val_lep_scale = {}

for m in a_masses:
    val_pho_smear[m] = {}
    val_pho_scale[m] = {}
    val_lep_smear[m] = {}
    val_lep_scale[m] = {}
    for y in years:
        val_pho_smear[m][y] = []
        val_pho_scale[m][y] = []
        val_lep_smear[m][y] = []
        val_lep_scale[m][y] = []
        f = open(path+'fit_run2_UL_ScaleSmear_ele/ScaleSmear_plot/M{0}/ScaleSmear_m{0}_{1}.dat'.format(m,y))
        lines = []

        for l in f.readlines()[-4:]:
            lines.append(l)
        
        val_pho_smear[m][y].append(float(lines[0].split('\t')[-2]))
        val_pho_scale[m][y].append(float(lines[2].split('\t')[1]))
        val_lep_smear[m][y].append(float(lines[1].split('\t')[-2]))
        val_lep_scale[m][y].append(float(lines[3].split('\t')[1]))

        f.close()

for y in years:
    mean_pho_smear = []
    mean_pho_scale = []
    mean_lep_smear = []
    mean_lep_scale = []
    for m in a_masses:
        mean_pho_smear.append(val_pho_smear[m][y][0])
        mean_pho_scale.append(val_pho_scale[m][y][0])
        mean_lep_smear.append(val_lep_smear[m][y][0])
        mean_lep_scale.append(val_lep_scale[m][y][0])
    
    #print mean_pho_scale
    #print mean_lep_scale
    print mean_lep_smear

    #print "year:",y,"\tpho_smear:",mean(mean_pho_smear),"pho_scale:",mean(mean_pho_scale),"lep_smear:",mean(mean_lep_smear),"lep_scale:",mean(mean_lep_scale)
    print "year:",y,"\tpho_smear:",round(mean(mean_pho_smear),3),"pho_scale:",round(mean(mean_pho_scale),3),"lep_smear:",round(mean(mean_lep_smear),3),"lep_scale:",round(mean(mean_lep_scale),3)
