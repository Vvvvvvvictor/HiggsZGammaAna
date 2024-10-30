import numpy as np
import re
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

training_variables = ["ALL", "Z_cos_theta", "lep_phi", "lep_cos_theta", "Z_lead_lepton_eta", "Z_sublead_lepton_eta", "gamma_eta", "l2g_deltaR", "l1g_deltaR", "gamma_mvaID", "gamma_ptRelErr", "H_relpt", "delta_eta_jj", "pt_balance", "delta_phi_zgjj", "delta_phi_jj", "photon_zeppenfeld", "jet_1_pt", "jet_2_pt", "jet1G_deltaR" , "jet2G_deltaR"]
log_path = '/eos/user/j/jiehan/root/outputs/significances/'

significances, significance_errs = [], []
for i in range(21):
    with open(log_path + f'{i}_0_two_jet_1D_4.json') as f:
        data = f.readlines()
        for line in data:
            if '"significance"' in line:
                significances.append(float(re.search(r"[-+]?\d*\.\d+", line).group()))
            if '"Delta_significance"' in line:
                significance_errs.append(float(re.search(r"[-+]?\d*\.\d+", line).group()))
             
significances = np.array(significances)
mean_significance = np.mean(significances)
significance_errs = np.array(significance_errs)
plt.errorbar(significances, training_variables, xerr=significance_errs, fmt='o')
plt.axvline(x=mean_significance, color='r', linestyle='--', label=f'Mean Significance: {mean_significance:.2f}')
plt.xlabel('Significance')
plt.ylabel('Training Variables')
plt.tight_layout()
plt.savefig('/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/significance_training_variables.png')
