import uproot
import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = uproot.open("/eos/user/j/jiehan/root/outputs/two_jet/data.root")['test'].arrays(library='pd')

# Plot a 2D hotmap of two variables
plt.figure()
plt.hist2d(data['gamma_mvaID'], data['bdt_score_t'], bins=100)
plt.xlabel('gamma_mvaID')
plt.ylabel('bdt_score')
plt.colorbar()
plt.savefig('two_jet.png')
