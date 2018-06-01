from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

# names = np.sort(glob.glob("hmc/h*.dat"))
# names = np.sort(glob.glob("hmcn/h*.dat"))
# names = np.sort(glob.glob("smd/s*.dat"))
# names = np.sort(glob.glob("smdn/s*.dat"))
names = np.sort(glob.glob("xy/x*.dat"))
# names = np.sort(glob.glob("xyn/x*.dat"))
# names = np.sort(glob.glob("xysmd/x*.dat"))
# names = np.sort(glob.glob("xysmdn/x*.dat"))

xy_error = np.genfromtxt('data_errorxy.csv', delimiter='\n')

kappa = np.genfromtxt('data_xy.csv', delimiter=',')[:,0]
susp = np.genfromtxt('data_xy.csv', delimiter=',')[:,1]

plt.errorbar(kappa, susp, yerr=xy_error, marker = 'o', capsize=5, markersize = 5, label = 'HMC w. Metro')
plt.xlabel(r'$1/\kappa$', fontsize=13)
plt.ylabel(r'$C_V$', fontsize=13)
plt.legend(loc = 'upper left')

plt.show()