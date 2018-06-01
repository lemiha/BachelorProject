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
xy_names = np.sort(glob.glob("xy/x*.dat"))
xyn_names = np.sort(glob.glob("xyn/x*.dat"))
xysmd_names = np.sort(glob.glob("xysmd/x*.dat"))
xysmdn_names = np.sort(glob.glob("xysmdn/x*.dat"))

hmc_ave = []
hmcn_ave = []
smd_ave = []
smdn_ave = []
xy_ave = []
xyn_ave = []
xysmd_ave = []
xysmdn_ave = []

# xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps = (1.0/xvalues)*(1.0/xvalues)

# Calculate the averages
for i in range(len(xvalues)):
    hmc_ave.append(np.average(np.genfromtxt(hmc_names[i], delimiter='  ')[:,1]))
    hmcn_ave.append(np.average(np.genfromtxt(hmcn_names[i], delimiter='  ')[:,1]))
    smd_ave.append(np.average(np.genfromtxt(smd_names[i], delimiter='  ')[:,1]))
    smdn_ave.append(np.average(np.genfromtxt(smdn_names[i], delimiter='  ')[:,1]))
    xy_ave.append(np.average(np.genfromtxt(xy_names[i], delimiter='  ')[:,1]))
    xyn_ave.append(np.average(np.genfromtxt(xyn_names[i], delimiter='  ')[:,1]))
    print(xysmd_names[i])
    xysmd_ave.append(np.average(np.genfromtxt(xysmd_names[i], delimiter='  ')[:,1]))
    xysmdn_ave.append(np.average(np.genfromtxt(xysmdn_names[i], delimiter='  ')[:,1]))

np.savetxt("data_aves.csv",  np.c_[hmc_ave, hmcn_ave, smd_ave, smdn_ave, xy_ave, xyn_ave, xysmd_ave, xysmdn_ave],
        delimiter=',', fmt='%1.6f', header="hmc_ave, hmcn_ave, smd_ave, smdn_ave, xy_ave, xyn_ave, xysmd_ave, xysmdn_ave")