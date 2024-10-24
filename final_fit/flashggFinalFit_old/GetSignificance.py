from array import array
from ROOT import *
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import math
import numpy as np

significances = []
significances_global = []
N_indep = 4.

for m in range(1,31):
	path = "/afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_runII/M"+str(m)
	f = TFile(path+'/higgsCombine_observed.Significance.mH125.root')

	tree = f.Get("limit")
	for quantile in tree:
		significances.append(tree.limit)
		N_trails = 1.0+np.sqrt(math.pi/2.)*N_indep*tree.limit
		p_local = 1.0 - norm.cdf(tree.limit,loc=0,scale=1)
		#p_global = p_local*N_trails
		p_global = 1.0-pow((1.-p_local),N_trails)
		print "sig_local:", tree.limit, "p_local:", p_local, "N_trails:", N_trails, "p_global:", p_global, "Z_global:", norm.ppf(1.0 - p_global,loc=0,scale=1)
		significances_global.append(norm.ppf(1.0 - p_global,loc=0,scale=1))
		#print tree.limit

	f.Close()

#print "p_count: ", p_count, "nToys: ", nToys, "p_global: ", float(p_count)/float(nToys), "Z_global", norm.ppf(1.0 - float(p_count)/float(nToys),loc=0,scale=1)

x = range(1,31)	
plt.xlim(xmin=0,xmax=31)
plt.xlabel('m(gg) GeV')
plt.ylabel('Significance')
plt.plot(x, significances, 'o-', c='blue')

plt.grid()
plt.savefig('local_p.png')
plt.close('all')
######

plt.xlim(xmin=0,xmax=31)
plt.xlabel('m(gg) GeV')
plt.ylabel('Global significance')
plt.plot(x, significances_global, 'o-', c='blue')

plt.grid()
plt.savefig('global_p.png')
plt.close('all')
