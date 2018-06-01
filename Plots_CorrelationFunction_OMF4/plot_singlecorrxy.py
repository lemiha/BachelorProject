from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

# Initialize variables
xy_names = np.sort(glob.glob("xy/corrx*.dat"))
xyn_names = np.sort(glob.glob("xyn/corrx*.dat"))
xysmd_names = np.sort(glob.glob("xysmd/corrx*.dat"))
xysmdn_names = np.sort(glob.glob("xysmdn/corrx*.dat"))

rvalues = np.array([2, 4, 6, 8, 10])
xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps = (1.0/xvalues)*(1.0/xvalues)

avesfile = "data_singlefile.csv"
xyaves = np.genfromtxt(avesfile, delimiter=',')[:,4]
xynaves = np.genfromtxt(avesfile, delimiter=',')[:,5]
xysmdaves = np.genfromtxt(avesfile, delimiter=',')[:,6]
xysmdnaves = np.genfromtxt(avesfile, delimiter=',')[:,7]

dataerrors = 'data_errors2.csv'
xy_error = np.genfromtxt(dataerrors, delimiter=',')[:,4]
xyn_error = np.genfromtxt(dataerrors, delimiter=',')[:,5]
xysmd_error = np.genfromtxt(dataerrors, delimiter=',')[:,6]
xysmdn_error = np.genfromtxt(dataerrors, delimiter=',')[:,7]

# Plot data
plt.errorbar(rvalues, xyaves, yerr=xy_error, marker = 'o', capsize=6, markersize = 5, label = 'XYHMC w. Metro')
plt.errorbar(rvalues, xynaves, yerr=xyn_error, marker = 'o', capsize=6, markersize = 5, label = 'XYHMCN w/o. Metro')
plt.errorbar(rvalues, xysmdaves, yerr=xysmd_error, marker = 'o', capsize=6, markersize = 5, label = 'XYSMD w. Metro')
plt.errorbar(rvalues, xysmdnaves, yerr=xysmdn_error, marker = 'o', capsize=6, markersize = 5, label = 'XYSMDN w/o. Metro')

plt.xlabel(r'$r$', fontsize=13)
plt.ylabel(r'$\left\langle C(r) \right\rangle$', fontsize=13)
plt.legend(loc = 'upper right')

plt.show()