from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

# Initialize variables
hmc_names = np.sort(glob.glob("hmc/h*.dat"))
hmcn_names = np.sort(glob.glob("hmcn/h*.dat"))
smd_names = np.sort(glob.glob("smd/s*.dat"))
smdn_names = np.sort(glob.glob("smdn/s*.dat"))

dataaves = 'data_aves.csv'
dataerrors = 'data_errors.csv'

hmc_ave = np.genfromtxt(dataaves, delimiter=',')[:,0]
hmcn_ave = np.genfromtxt(dataaves, delimiter=',')[:,1]
smd_ave = np.genfromtxt(dataaves, delimiter=',')[:,2]
smdn_ave = np.genfromtxt(dataaves, delimiter=',')[:,3]

hmc_error = np.genfromtxt(dataerrors, delimiter=',')[:,0]
hmcn_error = np.genfromtxt(dataerrors, delimiter=',')[:,1]
smd_error = np.genfromtxt(dataerrors, delimiter=',')[:,2]
smdn_error = np.genfromtxt(dataerrors, delimiter=',')[:,3]

xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps = (1.0/xvalues)*(1.0/xvalues)

plt.xlabel(r'$\epsilon^2$', fontsize=13)
plt.ylabel(r'$\left\langle E \right\rangle$', fontsize=13)
plt.errorbar(eps, hmc_ave, yerr=hmc_error, marker = 'o', capsize=6, markersize = 5, label = 'HMC w. Metro')
plt.errorbar(eps, hmcn_ave, yerr=hmcn_error, marker = 'o', capsize=6, markersize = 5, label = 'HMCN w/o. Metro')
plt.errorbar(eps, smd_ave, yerr=smd_error, marker = 'o', capsize=6, markersize = 5, label = 'SMD w. Metro')
plt.errorbar(eps, smdn_ave, yerr=smdn_error, marker = 'o', capsize=6, markersize = 5, label = 'SMDN w/o. Metro')
plt.legend(loc = 'upper left')

plt.show()