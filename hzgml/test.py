import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

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
# # Teset yields in categories
# # ==================================
# cats = ['zero_to_one_jet', "two_jet", "VH", "ZH", "ttH_had", "ttH_lep"]
# years = ['2016preVFP']
# for year in years:
#     data = uproot.open(f"/eos/home-j/jiehan/root/skimmed_ntuples_run2/Data/{year}.root")
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
cats = ['zero_to_one_jet', "two_jet"]
years = ['2016preVFP']
for year in years:
    data = uproot.open(f"/eos/home-j/jiehan/root/skimmed_ntuples_run2/Data/{year}.root")
    for cat in cats:
        cat_data = data[cat].arrays(['run', 'luminosityBlock', 'event'], library='pd')
        with open(f'{year}_{cat}_run_lumi_event.txt', 'w') as f:
            for i in range(len(cat_data['run'])):
                f.write(f'{int(cat_data["run"][i])} {int(cat_data["luminosityBlock"][i])} {int(cat_data["event"][i])}\n')

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