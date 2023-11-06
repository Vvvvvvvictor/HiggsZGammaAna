from array import array
from ROOT import *
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


nToys = 1000
path = "/afs/cern.ch/work/z/zewang/private/flashggfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/ALP_SigModel_param_UL/fit_results_runII/M3/"
p_local = 2.6

p_toys = {}
p_toys_max = []
p_count = 0
p_count_iToys = []
N = []
c0 = 0.5

for toy_i in range(1,nToys+1):
	f = open(path+"toy_"+str(toy_i)+".txt","r")
	p_toys[toy_i] = []
	N_temp = 0.0

	for l in f.readlines():
		if 'Significance:' in l:
			p_toys[toy_i].append(float(l.split(' ')[1]))
	
	#for p_t in p_toys[toy_i]:


	p_toys_max.append(max(p_toys[toy_i]))
	if max(p_toys[toy_i]) >= p_local:
		p_count = p_count + 1
		p_count_iToys.append(toy_i)

print "p_count: ", p_count, "nToys: ", nToys, "p_global: ", float(p_count)/float(nToys), "Z_global", norm.ppf(1.0 - float(p_count)/float(nToys),loc=0,scale=1)

x = range(95,181)	
plt.xlim(xmin=95,xmax=180)
plt.xlabel('m(llgg) GeV')
plt.ylabel('significance')
plt.plot([0., 180.], [p_local, p_local], c='black', linestyle='--')
plt.plot(x, p_toys[p_count_iToys[0]], 'o-', c='blue')

plt.grid()
plt.savefig('LookElesWhere.png')
plt.close('all')
######

plt.hist(p_toys_max, bins=40, normed=0, facecolor="aqua", edgecolor="black", alpha=0.3)
plt.xlabel('p max')
plt.ylabel('nToys')
plt.axvline(p_local)
plt.savefig('p_Max.png')
plt.close('all')
