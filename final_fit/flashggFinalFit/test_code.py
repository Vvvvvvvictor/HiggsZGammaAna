import ROOT
import uproot
from pdb import set_trace

inputfile = ROOT.TFile.Open("/eos/home-j/jiehan/root/ws_data/Data_2017/Data_2017.root") #/eos/user/j/jiehan/root/ws_cor_syst/signal_2017/ggH_M125_2017.root
workspace = inputfile.Get("CMS_hzg_workspace")
dataset = workspace.data("Data_13TeV_ggH0") #ggH_125_13TeV_ttHl
for i in range(10):
    cms_hzg_mass = dataset.get(i).getRealValue("CMS_hzg_mass")
    weight = dataset.get(i).getRealValue("weight")
    print(cms_hzg_mass, weight)