from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

# Initialize variables
xy_names = np.sort(glob.glob("hmc/h*.dat"))
xyn_names = np.sort(glob.glob("hmcn/h*.dat"))
xysmd_names = np.sort(glob.glob("smd/s*.dat"))
xysmdn_names = np.sort(glob.glob("smdn/s*.dat"))

dataaves = 'data_aves.csv'
dataerrors = 'data_errors.csv'

xy_ave = np.genfromtxt(dataaves, delimiter=',')[:,4]
xyn_ave = np.genfromtxt(dataaves, delimiter=',')[:,5]
xysmd_ave = np.genfromtxt(dataaves, delimiter=',')[:,6]
xysmdn_ave = np.genfromtxt(dataaves, delimiter=',')[:,7]

xy_error = np.genfromtxt(dataerrors, delimiter=',')[:,4]
xyn_error = np.genfromtxt(dataerrors, delimiter=',')[:,5]
xysmd_error = np.genfromtxt(dataerrors, delimiter=',')[:,6]
xysmdn_error = np.genfromtxt(dataerrors, delimiter=',')[:,7]

xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps = (1.0/xvalues)*(1.0/xvalues)

plt.xlabel(r'$\epsilon^2$', fontsize=13)
plt.ylabel(r'$\left\langle E \right\rangle$', fontsize=13)
plt.errorbar(eps, xy_ave, yerr=xy_error, marker = 'o', capsize=6, markersize = 5, label = 'XYHMC w. Metro')
plt.errorbar(eps, xyn_ave, yerr=xyn_error, marker = 'o', capsize=6, markersize = 5, label = 'XYHMCN w/o. Metro')
plt.errorbar(eps, xysmd_ave, yerr=xysmd_error, marker = 'o', capsize=6, markersize = 5, label = 'XYSMD w. Metro')
plt.errorbar(eps, xysmdn_ave, yerr=xysmdn_error, marker = 'o', capsize=6, markersize = 5, label = 'XYSMDN w/o. Metro')
plt.legend(loc = 'upper left')

plt.show()