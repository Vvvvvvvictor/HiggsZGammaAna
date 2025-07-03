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

# # =====================================
# # Generate Event Count Table for ggH M125 2017
# # =====================================
# print("\nEvent Count Table for ggH M125 2017 (Electron Channel)")
# print("======================================================================")

# file_path = "/eos/home-j/jiehan/root/skimmed_ntuples_run2/ggH_M125/2017.root"
# tree_name = "inclusive"

# # Define all columns that will be used to speed up loading
# columns_to_load = [
#     'z_ee', 'n_iso_photons', 'Z_mass', 'gamma_pt', 'H_mass',
#     'Z_cos_theta', 'lep_cos_theta', 'lep_phi', 'gamma_mvaID',
#     'l1g_deltaR', 'l2g_deltaR', 'H_relpt', 'gamma_ptRelErr',
#     'Z_lead_lepton_eta', 'Z_sublead_lepton_eta', 'weight_central', # Add weight for potential weighted counts
#     'Z_lead_lepton_deltaphi', 'Z_sublead_lepton_deltaphi', 'gamma_eta'
# ]

# try:
#     df_full = uproot.open(file_path)[tree_name].arrays(columns_to_load, library='pd')
# except Exception as e:
#     print(f"Error loading file: {e}")
#     df_full = pd.DataFrame() # Empty dataframe to avoid further errors

# if not df_full.empty:
#     # Helper function to print rows in the desired format
#     def print_row(label, count_or_df, is_weighted=False):
#         count = 0
#         if isinstance(count_or_df, pd.DataFrame):
#             if is_weighted and 'weight_central' in count_or_df.columns:
#                 count = count_or_df['weight_central'].sum()
#             else:
#                 count = len(count_or_df)
#         else:
#             count = count_or_df
        
#         # Format count to match example (integer or float if weighted)
#         count_str = f"{count:.0f}" if not is_weighted or count == int(count) else f"{count:.2f}"
#         print(f"{label.ljust(60)} & {count_str}")

#     # Apply selections sequentially
#     df = df_full.copy()

#     df = df[df['z_ee'] > 0]  # Z_ee channel
#     print_row(r"$Z_ee channel$", df)

#     baseline_df = df.copy()
#     print_row("Event filters (baseline for subsequent)", baseline_df)

#     # # Cuts relative to baseline
#     # print_row(r"baseline + cosTheta>0", baseline_df[baseline_df['Z_cos_theta'] > 0])
#     # print_row(r"baseline + cosTheta>0.5", baseline_df[baseline_df['Z_cos_theta'] > 0.5])
    
#     # print_row(r"baseline + costheta>0", baseline_df[baseline_df['lep_cos_theta'] > 0])
#     # print_row(r"baseline + costheta>0.5", baseline_df[baseline_df['lep_cos_theta'] > 0.5])

#     # print_row(r"baseline + phi (psi)>0", baseline_df[baseline_df['lep_phi'] > 0])
#     # print_row(r"baseline + phi (psi)>0.5", baseline_df[baseline_df['lep_phi'] > 0.5])

#     # print_row(r"baseline + Photonmva>0.6", baseline_df[baseline_df['gamma_mvaID'] > 0.6])
#     # print_row(r"baseline + Photonmva>0.75", baseline_df[baseline_df['gamma_mvaID'] > 0.75])

#     # print_row(r"baseline + mindR>1", baseline_df[baseline_df['l2g_deltaR'] > 1]) # l2g_deltaR is min(dR(l,g))
#     # print_row(r"baseline + mindR>1.2", baseline_df[baseline_df['l2g_deltaR'] > 1.2])

#     # print_row(r"baseline + maxdR>1.5", baseline_df[baseline_df['l1g_deltaR'] > 1.5]) # l1g_deltaR is max(dR(l,g))
#     # print_row(r"baseline + maxdR>2", baseline_df[baseline_df['l1g_deltaR'] > 2])
    
#     # # Ensure H_mass is not zero for H_relpt calculation if not already filtered
#     # temp_df_ptmass = baseline_df[baseline_df['H_mass'] != 0]
#     # print_row(r"baseline + pT mass>0.25", temp_df_ptmass[temp_df_ptmass['H_relpt'] > 0.25])
#     # print_row(r"baseline + pT mass>0.5", temp_df_ptmass[temp_df_ptmass['H_relpt'] > 0.5])

#     # print_row(r"baseline + $\gamma$ resolution>0.02", baseline_df[baseline_df['gamma_ptRelErr'] > 0.02])
#     # print_row(r"baseline + $\gamma$ resolution>0.05", baseline_df[baseline_df['gamma_ptRelErr'] > 0.05])

#     # print_row(r"baseline + lead lep $\eta$>0", baseline_df[baseline_df['Z_lead_lepton_eta'] > 0])
#     # print_row(r"baseline + lead lep $\eta$>0.5", baseline_df[baseline_df['Z_lead_lepton_eta'] > 0.5])

#     # print_row(r"baseline + sublead lep $\eta$>0", baseline_df[baseline_df['Z_sublead_lepton_eta'] > 0])
#     # print_row(r"baseline + sublead lep $\eta$>0.5", baseline_df[baseline_df['Z_sublead_lepton_eta'] > 0.5])


#     print_row(r"baseline + lead lepton $\eta$>1", baseline_df[baseline_df['Z_lead_lepton_eta'] > 1])
#     print_row(r"baseline + lead lepton $\eta$>2", baseline_df[baseline_df['Z_lead_lepton_eta'] > 2])

#     print_row(r"baseline + sublead lepton $\eta$>1", baseline_df[baseline_df['Z_sublead_lepton_eta'] > 1])
#     print_row(r"baseline + sublead lepton $\eta$>2", baseline_df[baseline_df['Z_sublead_lepton_eta'] > 2])

#     print_row(r"baseline + photon $\eta$>1", baseline_df[baseline_df['gamma_eta'] > 1])
#     print_row(r"baseline + photon $\eta$>2", baseline_df[baseline_df['gamma_eta'] > 2])

#     print_row(r"baseline + lead lep photon dphi>1", baseline_df[baseline_df['Z_sublead_lepton_deltaphi'] > 1])
#     print_row(r"baseline + lead lep photon dphi>2", baseline_df[baseline_df['Z_sublead_lepton_deltaphi'] > 2])

#     print_row(r"baseline + sublead lep photon dphi>1", baseline_df[baseline_df['Z_lead_lepton_deltaphi'] > 1])
#     print_row(r"baseline + sublead lep photon dphi>2", baseline_df[baseline_df['Z_lead_lepton_deltaphi'] > 2])

# print("======================================================================")

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
years = ['2023postBPix']
types = ["mumu", "ee"]
for year in years:
    # data = uproot.open(f"/eos/home-j/jiehan/root/cutflow/ggH_M125/{year}.root")
    data = pd.read_parquet(f"/eos/home-j/jiehan/parquet/cutflow_ggf/ggH_M125_{year}/merged_nominal.parquet")
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

        temp = cat_data
        weight = temp["weight_central"]
        print(f'{year} {cat} total events: {len(temp)}')
        print(f'{year} {cat} lumi weight: {weight.sum()}')
        print(f'{year} {cat} prefire weight: {weight.sum()}')
        weight = weight * temp["weight_pu_reweight_sf_central"]
        print(f'{year} {cat} pileup weight: {weight.sum()}')
        weight = weight * temp["weight_btag_deepjet_wp_sf_SelectedJet_central"]
        print(f'{year} {cat} btag weight: {weight.sum()}')
        weight = weight * temp["weight_electron_iso_sf_SelectedElectron_central"] * temp["weight_electron_reco_sf_SelectedElectron_central"] * temp["weight_electron_wplid_sf_SelectedElectron_central"]
        print(f'{year} {cat} electron weight: {weight.sum()}')
        weight = weight * temp["weight_muon_looseid_sf_SelectedMuon_central"] * temp["weight_muon_iso_sf_SelectedMuon_central"] * temp["weight_muon_looseid_sf_SelectedMuon_central"]
        print(f'{year} {cat} muon weight: {weight.sum()}')
        weight = weight * temp["weight_photon_csev_sf_Photon_central"] * temp["weight_photon_id_sf_Photon_central"]
        print(f'{year} {cat} photon weight: {weight.sum()}')
        weight = weight * temp["weight_electron_hlt_sf_SelectedElectron_central"] * temp["weight_muon_hlt_sf_SelectedMuon_central"]
        print(f'{year} {cat} hlt weight: {weight.sum()}')
        print(weight)

        # for t in types:
        #     temp = cat_data.query(f'z_{t} > 0')
        #     weight = temp["weight_central"]
        #     print(f'{year} {cat} {t} total events: {len(temp)}')
        #     print(f'{year} {cat} {t} lumi weight: {weight.sum()}')
        #     weight = weight * temp["weight_L1_prefiring_sf_central"]
        #     print(f'{year} {cat} {t} prefire weight: {weight.sum()}')
        #     weight = weight * temp["weight_pu_reweight_sf_central"]
        #     print(f'{year} {cat} {t} pileup weight: {weight.sum()}')
        #     weight = weight * temp["weight_btag_deepjet_wp_sf_SelectedJet_central"]
        #     print(f'{year} {cat} {t} btag weight: {weight.sum()}')
        #     weight = weight * temp["weight_electron_iso_sf_SelectedElectron_central"] * temp["weight_electron_reco_sf_SelectedElectron_central"] * temp["weight_electron_wplid_sf_SelectedElectron_central"]
        #     print(f'{year} {cat} {t} electron weight: {weight.sum()}')
        #     weight = weight * temp["weight_muon_looseid_sf_SelectedMuon_central"] * temp["weight_muon_iso_sf_SelectedMuon_central"] * temp["weight_muon_looseid_sf_SelectedMuon_central"]
        #     print(f'{year} {cat} {t} muon weight: {weight.sum()}')
        #     weight = weight * temp["weight_photon_csev_sf_Photon_central"] * temp["weight_photon_id_sf_Photon_central"]
        #     print(f'{year} {cat} {t} photon weight: {weight.sum()}')
        #     weight = weight * temp["weight_electron_hlt_sf_SelectedElectron_central"] * temp["weight_muon_hlt_sf_SelectedMuon_central"]

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


# # =====================================
# # Check ntuples
# # =====================================

# path = "/eos/home-p/pelai/HZgamma/Root_Dataset/run3/NanoV12"
# data_name = ["Data"]
# bkg_name = ["DYGto2LG_10to50", "DYGto2LG_50to100", "DYJetsToLL", "ZG2JToG2L2J", "TT"]
# years = ["2022preEE", "2022postEE"]

# data, bkg = pd.DataFrame(), pd.DataFrame()
# for year in years:
#     for name in data_name:
#         temp = uproot.open(f"{path}/{name}/{year}.root")['two_jet'].arrays(library='pd')
#         data = pd.concat([data, temp])
#     for name in bkg_name:
#         temp = uproot.open(f"{path}/{name}/{year}.root")['two_jet'].arrays(library='pd')
#         bkg = pd.concat([bkg, temp])

# print(r"$n_{jets} \geq 2$ and $n_{b-jets} = 0$ and $n_{mu} = 2")
# selection = "n_muons==2"
# print(f"sf: {data.query(selection)['weight'].sum() / bkg.query(selection)['weight'].sum()}")
# selection = "(mass_jj>400) & (delta_eta_jj>2.5)"
# print(r"After applying $m_{jj} > 400$ GeV and $\Delta\eta_{jj} > 2.5$")
# print(f"sf: {data.query(selection)['weight'].sum() / bkg.query(selection)['weight'].sum()}")
# print(r"After applying $p_{T}^{j1} > 35$ GeV and $p_{T}^{j2} > 25$ GeV")
# selection += " & (jet_1_pt>35) & (jet_2_pt>25)"
# print(f"sf: {data.query(selection)['weight'].sum() / bkg.query(selection)['weight'].sum()}")

# # =================================
# # Check the data
# # =================================
# path = "/eos/home-j/jiehan/root/fitting_signal/"
# years = ["2016preVFP", "2016postVFP", "2017", "2018"]
# samples = ["ggH", "VBF", "ZH", "WH", "ttH"]
# cats = ["VBF0", "VBF1", "VBF2", "VBF3"]

# for cat in cats:
#     yields = 0
#     for sample in samples:
#         for year in years:
#             for flav in ('ele', 'mu'):
#                 # print(f"{path}/{sample}_M125_{year}/output_{sample}_M125.root")
#                 data = uproot.open(f"{path}/{sample}_M125_{year}/output_{sample}_M125.root")[f"DiphotonTree/{sample}_125_{flav}_13TeV_{cat}"].arrays(['weight'], library='pd')
#                 signal = data['weight'].sum()
#                 count = data['weight'].count()
#                 print(f'{year} {sample} {cat} {flav} yields: {signal:.3f}({count})')
#                 yields += data['weight'].sum()
#     print(f'{cat} yields: {yields}\n')