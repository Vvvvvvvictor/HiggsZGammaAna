import ROOT

bkg_name = "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Background/HZGamma_BkgModel_UL_xingchen_chi2/fit_results_run2_cat3/CMS-HGG_mva_13TeV_multipdf_cat3.root"
file_bkg = ROOT.TFile(bkg_name)
#data_name = "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/InputData/xingchen/output_file/HZGamma_data_bkg_workspace_cat3.root"
#file_data = ROOT.TFile(data_name)

pdf_bkg = file_bkg.multipdf.pdf("CMS_hgg_untag_cat3_bkgshape")
hist_bkg = file_bkg.multipdf.data("roohist_data_mass_cat0")
file_bkg.multipdf.Print()
stat = ROOT.RooChi2Var("a", "a", pdf_bkg, hist_bkg)
#stat.getVal()