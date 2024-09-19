import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

def get_yield(file_name):
    input_file = uproot.open(file_name)
    yields = []
    for tree in ["zero_jet", "one_jet", "two_jet", "ZH", "ttH_had", "ttH_lep"]:
        # only store the sum of weights for each channel
        sam_yield = np.sum(input_file[tree].arrays(["weight"], library="pd")["weight"])
        yields.append(sam_yield)
    print(file_name, ":", yields)
    return np.array(yields)

PATH = "/eos/home-j/jiehan/root/skimmed_ntuples/"
prods = ["ggH_M125", "VBF_M125", "WminusH_M125", "WplusH_M125", "ZH_M125", "ttH_M125"]
prod_labels = ["ggH", "VBF", r"$W^-$H", r"$W^+$H", "ZH", "ttH"]
chans = ["zero_jet", "one_jet", "two_jet", "ZH", "ttH_had", "ttH_lep"]
years = ["2016preVFP", "2016postVFP", "2017", "2018"] # ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE"]
prod_yield = []
for prod in prods:
    # add up all yields for each production mode
    yield_sum = []
    for year in years:
        yield_sum.append(get_yield(PATH+f"{prod}/{year}.root"))
    yield_sum = np.sum(yield_sum, axis=0)
    prod_yield.append(yield_sum)
    
prod_ratio = list(np.array(prod_yield)/np.sum(prod_yield, axis=0))

# plot prod_ratio in each prods
fig, ax = plt.subplots()
hep.histplot(
    prod_ratio, 
    bins=np.linspace(-0.5, 5.5, 7), 
    histtype='fill', 
    stack=True,
    color=['b', 'g', 'r', 'c', 'm', 'y'],
    ax=ax,
    label=prod_labels
)
ax.set_xticks(np.arange(6))
ax.set_xticklabels(chans)
ax.set_ylabel("Ratio")
ax.set_ylim(0, 1.2)
ax.legend(ncol=3)
plt.savefig("prod_ratio_label_chans.png")
plt.clf()
    
# transpose prod_yield 
prod_ratio = list(np.array(prod_yield).T/np.sum(prod_yield, axis=1))

# plot prod_ratio in each chans
fig, ax = plt.subplots()
hep.histplot(
    prod_ratio, 
    bins=np.linspace(-0.5, 5.5, 7), 
    histtype='fill', 
    stack=True,
    color=['b', 'g', 'r', 'c', 'm', 'y'],
    ax=ax,
    label=chans
)
ax.set_xticks(np.arange(6))
ax.set_xticklabels(prod_labels)
ax.set_ylabel("Ratio")
ax.set_ylim(0, 1.2)
ax.legend(ncol=3)
plt.savefig("prod_ratio_label_prods.png")