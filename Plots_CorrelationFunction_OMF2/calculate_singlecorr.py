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
xy_names = np.sort(glob.glob("xy/corrx*.dat"))
xyn_names = np.sort(glob.glob("xyn/corrx*.dat"))
xysmd_names = np.sort(glob.glob("xysmd/corrx*.dat"))
xysmdn_names = np.sort(glob.glob("xysmdn/corrx*.dat"))

rvalues = np.array([2, 4, 6, 8, 10])
xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps = (1.0/xvalues)*(1.0/xvalues)
hmcaves = []
hmcnaves = []
smdaves = []
smdnaves = []
xyaves = []
xynaves = []
xysmdaves = []
xysmdnaves = []

xstart = 1
for i in range(len(rvalues)):
    hmcaves.append(np.average(np.genfromtxt(hmc_names[xstart], delimiter='  ')[:,i]))
    hmcnaves.append(np.average(np.genfromtxt(hmcn_names[xstart], delimiter='  ')[:,i]))
    smdaves.append(np.average(np.genfromtxt(smd_names[xstart], delimiter='  ')[:,i]))
    smdnaves.append(np.average(np.genfromtxt(smdn_names[xstart], delimiter='  ')[:,i]))
    xyaves.append(np.average(np.genfromtxt(xy_names[xstart], delimiter='  ')[:,i]))
    xynaves.append(np.average(np.genfromtxt(xyn_names[xstart], delimiter='  ')[:,i]))
    xysmdaves.append(np.average(np.genfromtxt(xysmd_names[xstart], delimiter='  ')[:,i]))
    xysmdnaves.append(np.average(np.genfromtxt(xysmdn_names[xstart], delimiter='  ')[:,i]))

np.savetxt("data_singlefile.csv",  np.c_[hmcaves, hmcnaves, smdaves, smdnaves, xyaves, xynaves, xysmdaves, xysmdnaves],
        delimiter=',', fmt='%1.6f', header="hmcaves, hmcnaves, smdaves, smdnaves, xyaves, xynaves, xysmdaves, xysmdnaves")
    