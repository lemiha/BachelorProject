from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

# Initialize variables
hmc_names = np.sort(glob.glob("hmc/corrh*.dat"))
hmcn_names = np.sort(glob.glob("hmcn/corrh*.dat"))
smd_names = np.sort(glob.glob("smd/corrs*.dat"))
smdn_names = np.sort(glob.glob("smdn/corrs*.dat"))

rvalues = np.array([2, 4, 6, 8, 10])
xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps = (1.0/xvalues)*(1.0/xvalues)

avesfile = "data_singlefile.csv"
hmcaves = np.genfromtxt(avesfile, delimiter=',')[:,0]
hmcnaves = np.genfromtxt(avesfile, delimiter=',')[:,1]
smdaves = np.genfromtxt(avesfile, delimiter=',')[:,2]
smdnaves = np.genfromtxt(avesfile, delimiter=',')[:,3]

dataerrors = 'data_errors2.csv'
hmc_error = np.genfromtxt(dataerrors, delimiter=',')[:,0]
hmcn_error = np.genfromtxt(dataerrors, delimiter=',')[:,1]
smd_error = np.genfromtxt(dataerrors, delimiter=',')[:,2]
smdn_error = np.genfromtxt(dataerrors, delimiter=',')[:,3]

# Plot data
plt.errorbar(rvalues, hmcaves, yerr=hmc_error, marker = 'o', capsize=6, markersize = 5, label = 'HMC w. Metro')
plt.errorbar(rvalues, hmcnaves, yerr=hmcn_error, marker = 'o', capsize=6, markersize = 5, label = 'HMCN w/o. Metro')
plt.errorbar(rvalues, smdaves, yerr=smd_error, marker = 'o', capsize=6, markersize = 5, label = 'SMD w. Metro')
plt.errorbar(rvalues, smdnaves, yerr=smdn_error, marker = 'o', capsize=6, markersize = 5, label = 'SMDN w/o. Metro')# # plt.ylim((1.5,2.7))

plt.xlabel(r'$r$', fontsize=13)
plt.ylabel(r'$\left\langle C(r) \right\rangle$', fontsize=13)
plt.legend(loc = 'upper right')

plt.show()