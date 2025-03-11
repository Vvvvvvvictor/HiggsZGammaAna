import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

# # Load the data
# data = uproot.open("/eos/user/j/jiehan/root/outputs/two_jet/data.root")['test'].arrays(library='pd')

# # Plot a 2D hotmap of two variables
# plt.figure()
# plt.hist2d(data['gamma_mvaID'], data['bdt_score_t'], bins=100)
# plt.xlabel('gamma_mvaID')
# plt.ylabel('bdt_score')
# plt.colorbar()
# plt.savefig('two_jet.png')

# # ==================================
# # Test yields in categories
# # ==================================
# cats = ['zero_to_one_jet', "two_jet", "VH", "ZH", "ttH_had", "ttH_lep"]
# # years = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]
# years = ["2016preVFP", "2016postVFP", "2017", "2018"]
# cats_yields = {'zero_to_one_jet': 0, "two_jet": 0, "VH": 0, "ZH": 0, "ttH_had": 0, "ttH_lep": 0}
# for year in years:
#     data = uproot.open(f"/eos/home-j/jiehan/root/skimmed_ntuples_run2/ggH_M125/{year}.root")
#     for cat in cats:
#         cat_data = data[cat].arrays(['weight', 'H_mass'], library='pd').query('H_mass > 120 & H_mass < 130')
#         cats_yields[cat] += cat_data['weight'].sum()
# print(cats_yields)

# for year in years:
#     data = uproot.open(f"/eos/home-j/jiehan/root/skimmed_ntuples_run2/ggH_M125/{year}.root")
#     data_inc = data["inclusive"].arrays(['z_mumu', 'z_ee', 'H_mass'], library='pd').query('H_mass > 100')
#     print(f'Year {year} inclusive z_mumu: {data_inc["z_mumu"].sum()}, z_ee: {data_inc["z_ee"].sum()}')
#     for cat in cats:
#         cat_data = data[cat].arrays(['z_mumu', 'z_ee', 'H_mass'], library='pd')
#         if "jet" in cat:
#             cat_data = cat_data.query('H_mass > 100')
#         else:
#             cat_data = cat_data.query('H_mass > 100 & (H_mass < 120 | H_mass > 130)')
#         print(f'{year} {cat} z_mumu: {cat_data["z_mumu"].sum()}, z_ee: {cat_data["z_ee"].sum()}')
#     print("================================================================================================")

# ===================================
# Print run, lumi, event and store in files(int number)
# ===================================
# cats = ["two_jet"]
# years = ['2022preEE']
# types = ["mu", "ele"]
# config = '356381 283 222634940'.split(' ')
# for year in years:
#     data = uproot.open(f"/eos/home-j/jiehan/root/skimmed_ntuples_run3/Data/{year}.root")
#     data = data['inclusive'].arrays(library='pd')
#     for i in data.columns:
#         temp = data.query(f'run == {config[0]} & luminosityBlock == {config[1]} & event == {config[2]}')
#         print(i, temp[i].values[0])
    # for cat in cats:
    #     cat_data = data[cat].arrays(['run', 'luminosityBlock', 'event', "z_mumu", "z_ee"], library='pd')
    #     with open(f'{year}_{cat}_{types[0]}.txt', 'w') as f:
    #         temp = cat_data.query('z_mumu > 0').reset_index(drop=True)
    #         for i in range(len(temp)):
    #             f.write(f'{int(temp["run"][i])} {int(temp["luminosityBlock"][i])} {int(temp["event"][i])}\n')
    #     with open(f'{year}_{cat}_{types[1]}.txt', 'w') as f:
    #         temp = cat_data.query('z_ee > 0').reset_index(drop=True)
    #         for i in range(len(temp)):
    #             f.write(f'{int(temp["run"][i])} {int(temp["luminosityBlock"][i])} {int(temp["event"][i])}\n')

# # ==================================
# # Show MET_pt and n_b_jets
# # ==================================
# years = ['2016preVFP'] #, '2016postVFP', '2017', '2018']
# for year in years:
#     data = uproot.open(f"/eos/home-j/jiehan/root/skimmed_ntuples_run2/Data/{year}.root")['inclusive'].arrays(['MET_pt', 'n_b_jets', 'n_jets'], library='pd')
#     data = data.query('n_jets < 2')
#     hist, _ = np.histogram(data['MET_pt'], bins=20, range=(0, 200))
#     for i, h in enumerate(hist):
#         print(f'{year} MET_pt {i*10}-{(i+1)*10}: {h}')

# =====================================
# Compare two log file
# =====================================

# with open('2022preEE_two_jet_ele.txt', 'r') as f:
#     data1 = f.readlines()
# with open('daje.txt', 'r') as f:
#     data2 = f.readlines()

# # check if each line in data1 is in data2
# print("In HiggsDNA not in DAJE:")
# for d1 in data1:
#     if d1 not in data2:
#         print(d1, end='')

# # check if each line in data2 is in data1
# print("In DAJE not in HiggsDNA:")
# for d2 in data2:
#     if d2 not in data1:
#         print(d2, end='')

# =====================================
# Get weight with syst
# =====================================
cats = ["inclusive"]
years = ['2017']
types = ["mumu", "ee"]
for year in years:
    # data = uproot.open(f"/eos/home-j/jiehan/root/cutflow/ggH_M125/{year}.root")
    data = pd.read_parquet(f"/eos/home-j/jiehan/parquet/nanov9/cutflow/ggH_M125_{year}/merged_nominal.parquet")
    print(data.columns)
    for cat in cats:
        # 加入所有满足这个格式的变量 'weight*central'
        variables = ['H_mass', 'z_mumu', 'z_ee']
        # for i in data[cat].keys():
        #     if i.startswith('weight') and i.endswith('central'):
        #         variables.append(i)
        # cat_data = data[cat].arrays(variables, library='pd')
        if cat == "inclusive":
            cat_data = data
        else:
            exit()
        for t in types:
            temp = cat_data.query(f'z_{t} > 0')
            weight = temp["weight_central"]
            print(f'{year} {cat} {t} lumi weight: {weight.sum()}')
            weight = weight * temp["weight_L1_prefiring_sf_central"]
            print(f'{year} {cat} {t} prefire weight: {weight.sum()}')
            weight = weight * temp["weight_pu_reweight_sf_central"]
            print(f'{year} {cat} {t} pileup weight: {weight.sum()}')
            weight = weight * temp["weight_btag_deepjet_wp_sf_SelectedJet_central"]
            print(f'{year} {cat} {t} btag weight: {weight.sum()}')
            weight = weight * temp["weight_electron_wplid_sf_SelectedElectron_central"]
            print(f'{year} {cat} {t} ele id weight: {weight.sum()}')
            weight = weight * temp["weight_muon_looseid_sf_SelectedMuon_central"]
            print(f'{year} {cat} {t} muon id weight: {weight.sum()}')

# =====================================
# Have a look at the hgg data
# =====================================
# data = uproot.open("/eos/user/c/chpan/hgg_ntuples/root_fiducial_final/root_fiducial/GluGluHtoGG_M-125_preEE/output_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_in.root")['DiphotonTree/ggh_125_13TeV_best_resolution'].arrays(library='pd')
# print(data.columns)


# =====================================
# Check ntuples
# =====================================
# for cat in ["VBF0", "VBF1", "VBF2", "VBF3"]:
#     cat_sum = 0
#     for proc in ["ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH"]:
#         for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
#             data = uproot.open(f"/eos/home-j/jiehan/root/fitting_signal/{proc}_M125_{year}/output_{proc}_M125.root")[f"DiphotonTree/{proc}_125_13TeV_{cat}"].arrays(['weight'], library='pd')
#             print(f'{year} {proc} {cat} weight: {data["weight"].sum()}')
#             cat_sum += data["weight"].sum()
#     print(f'{cat} total weight: {cat_sum}\n')

# for cat in ["VBF0", "VBF1", "VBF2", "VBF3"]:
#     cat_sum = 0
#     # data = uproot.open(f"/eos/home-j/jiehan/root/fitting_bkg/Data/output_Data_Run2.root")[f"DiphotonTree/Data_13TeV_{cat}"].arrays(['weight'], library='pd')
#     # print(f'data {cat} weight: {data["weight"].sum()}')
#     for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
#         data = uproot.open(f"/eos/home-j/jiehan/root/fitting_bkg/Data_{year}/output_Data_{year}.root")[f"DiphotonTree/Data_13TeV_{cat}"].arrays(['weight'], library='pd')
#         print(f'{year} data {cat} weight: {data["weight"].sum()}')
#         cat_sum += data["weight"].sum()
#     print(f'{cat} total weight: {cat_sum}\n')
