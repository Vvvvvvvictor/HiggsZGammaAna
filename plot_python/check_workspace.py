import ROOT
import os

file = ROOT.TFile.Open("/afs/cern.ch/work/z/zewang/public/jet01_run2_CoreFunc_v2/HZGamma_data_sig_Hm125_workspace_cat0.root")
workspace = file.Get("CMS_hzg_workspace")
var = workspace.var("CMS_hzg_mass")
dataset = workspace.data("ggh_125_13TeV_untag_cat0")
datahist = workspace.data("ggh_125_13TeV_untag_cat0_hist")
# for i in range(dataset.numEntries()):
#     print(f"dataset values = {dataset.get(i).getRealValue('CMS_hzg_weight')}")
for i in range(datahist.numEntries()):
    print(f"datahist values = {datahist.get(i).getRealValue('CMS_hzg_weight')}")