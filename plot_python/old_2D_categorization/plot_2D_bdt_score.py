import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import uproot
import mplhep as hep
plt.style.use(hep.style.CMS)

filepath = "/eos/user/j/jiehan/root/outputs/two_jet/"

# Load data from ROOT files into pandas dataframes
data1 = uproot.open(filepath + "ggH_M125.root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight"], library="pd")
data2 = uproot.open(filepath + "VBF_M125.root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight"], library="pd")
data3 = uproot.open(filepath + "DYJetsToLL.root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight"], library="pd")
data4 = uproot.open(filepath + "ZGToLLG.root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight"], library="pd")
data5 = uproot.open(filepath + "ZG2JToG2L2J.root")["two_jet"].arrays(["bdt_score_t", "vbf_score_t", "weight"], library="pd")

# Create a figure with a grid layout
fig = plt.figure(constrained_layout=True, figsize=(8, 8), dpi=200)
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 2], width_ratios=[2, 1])

# Create the heatmap in the bottom-left position
ax_heatmap = fig.add_subplot(gs[1, 0])
sns.histplot(data=data1, x="vbf_score_t", y="bdt_score_t", weights="weight", bins=20, pmax=0.9, cmap='Reds', cbar=False, ax=ax_heatmap, alpha=0.5)
sns.histplot(data=data2, x="vbf_score_t", y="bdt_score_t", weights="weight", bins=20, pmax=0.9, cmap='Greens', cbar=False, ax=ax_heatmap, alpha=0.5)
sns.histplot(data=data3, x="vbf_score_t", y="bdt_score_t", weights="weight", bins=20, pmax=0.9, cmap='Oranges', cbar=False, ax=ax_heatmap, alpha=0.5)
sns.histplot(data=data4, x="vbf_score_t", y="bdt_score_t", weights="weight", bins=20, pmax=0.9, cmap='Blues', cbar=False, ax=ax_heatmap, alpha=0.5)
sns.histplot(data=data5, x="vbf_score_t", y="bdt_score_t", weights="weight", bins=20, pmax=0.9, cmap='Purples', cbar=False, ax=ax_heatmap, alpha=0.5)

# draw two lines
with open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/best_boundaries.txt", "r") as f:
    pos = list(map(float, f.readline().split(",")))
bin_x = [0,1]
bin_y = [0,1]
for p in pos[:-1]:
    ax_heatmap.axhline(p, color='black', linestyle='--')
    bin_y.insert(-1, p)
ax_heatmap.axvline(pos[-1], color='black', linestyle='--')
bin_x.insert(-1, pos[-1])
print(bin_x, bin_y)
for i in range((len(bin_y)-1) * (len(bin_x)-1)):
    ax_heatmap.annotate(f"{i}", 
                        ((bin_x[i//(len(bin_y)-1)]+bin_x[i//(len(bin_y)-1)+1])/2, 
                         (bin_y[i%(len(bin_y)-1)]+bin_y[i%(len(bin_y)-1)+1])/2), 
                        color='black', fontsize=16, ha='center', va='center')

ax_heatmap.set_xlabel('VBF score_t')
ax_heatmap.set_ylabel('BDT score_t')
ax_heatmap.set_xlim(0, 1)
ax_heatmap.set_ylim(0, 1)

# Create the histogram for x-axis projection in the top-left position
x_projection3, _ = np.histogram(data3["vbf_score_t"], bins=50, weights=data3["weight"], density=True)
x_projection4, bins = np.histogram(data4["vbf_score_t"], bins=50, weights=data4["weight"], density=True)
x_projection5, _ = np.histogram(data5["vbf_score_t"], bins=50, weights=data5["weight"], density=True)
x_projection1, _ = np.histogram(data1["vbf_score_t"], bins=50, weights=data1["weight"], density=True)
x_projection2, _ = np.histogram(data2["vbf_score_t"], bins=50, weights=data2["weight"], density=True)
ax_hist_x = fig.add_subplot(gs[0, 0], sharex=ax_heatmap)
hep.histplot(x_projection2, bins, color='green', histtype='fill', alpha=0.5, label='VBF', ax=ax_hist_x)
hep.histplot(x_projection5, bins, color='purple', histtype='fill', alpha=0.5, label='ZG2JToG2L2J', ax=ax_hist_x)
hep.histplot(x_projection1, bins, color='red', histtype='fill', alpha=0.5, label='ggH', ax=ax_hist_x)
hep.histplot(x_projection4, bins, color='blue', histtype='fill', alpha=0.5, label='ZGToLLG', ax=ax_hist_x)
hep.histplot(x_projection3, bins, color='orange', histtype='fill', alpha=0.5, label='DYJetsToLL', ax=ax_hist_x)
# ax_hist_x.set_title('vbf score_t')
ax_hist_x.set_ylabel('a.u.')
ax_hist_x.set_yscale('log')

# Create the histogram for y-axis projection in the bottom-right position
y_projection3, _ = np.histogram(data3["bdt_score_t"], bins=50, weights=data3["weight"], density=True)
y_projection4, bins = np.histogram(data4["bdt_score_t"], bins=50, weights=data4["weight"], density=True)
y_projection5, _ = np.histogram(data5["bdt_score_t"], bins=50, weights=data5["weight"], density=True)
y_projection1, _ = np.histogram(data1["bdt_score_t"], bins=50, weights=data1["weight"], density=True)
y_projection2, _ = np.histogram(data2["bdt_score_t"], bins=50, weights=data2["weight"], density=True)
ax_hist_y = fig.add_subplot(gs[1, 1], sharey=ax_heatmap)
hep.histplot(y_projection2, bins, color='green', histtype='fill', alpha=0.5, label='VBF', ax=ax_hist_y, orientation='horizontal')
hep.histplot(y_projection1, bins, color='red', histtype='fill', alpha=0.5, label='ggH', ax=ax_hist_y, orientation='horizontal')
hep.histplot(y_projection5, bins, color='purple', histtype='fill', alpha=0.5, label='ZG2JToG2L2J', ax=ax_hist_y, orientation='horizontal')
hep.histplot(y_projection4, bins, color='blue', histtype='fill', alpha=0.5, label='ZGToLLG', ax=ax_hist_y, orientation='horizontal')
hep.histplot(y_projection3, bins, color='orange', histtype='fill', alpha=0.5, label='DYJetsToLL', ax=ax_hist_y, orientation='horizontal')
# ax_hist_y.set_title('bdt score_t')
ax_hist_y.set_xlabel('a.u.')
ax_hist_y.set_xscale('log')

# add a legend in gs[0, 1]
ax_legend = fig.add_subplot(gs[0, 1])
ax_legend.axis('off')
# add objects to legend
hep.histplot([0], [0,1], color='red', histtype='fill', alpha=0.5, label='ggH', ax=ax_legend)
hep.histplot([0], [0,1], color='green', histtype='fill', alpha=0.5, label='VBF', ax=ax_legend)
hep.histplot([0], [0,1], color='orange', histtype='fill', alpha=0.5, label='DYJetsToLL', ax=ax_legend)
hep.histplot([0], [0,1], color='blue', histtype='fill', alpha=0.5, label='ZGToLLG', ax=ax_legend)
hep.histplot([0], [0,1], color='purple', histtype='fill', alpha=0.5, label='ZG2JToG2L2J', ax=ax_legend)
ax_legend.legend(loc='center')

# Save the combined plot as an image file
plt.savefig('/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/vbf_2j_heatmap.png')