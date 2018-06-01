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

all_names = np.array([hmc_names, hmcn_names, smd_names, smdn_names, xy_names, xyn_names, xysmd_names, xysmdn_names])
all_time = []

# Gather the run times in an array
for names in all_names:
    tmp = []
    for name in names:
        fileHandle = open(name, "r")
        lineList = fileHandle.readlines()
        tmp.append(float(lineList[-1].split("=", 2)[2]))
    all_time.append(tmp)

all_time = np.asarray(all_time)
fileHandle.close()

np.savetxt("data_time.csv",  np.c_[all_time[0], all_time[1], all_time[2], all_time[3], all_time[4], all_time[5], all_time[6], all_time[7]],
        delimiter=',', fmt='%1.6f', header="hmc_time, hmcn_time, smd_time, smdn_time, xy_time, xyn_time, xysmd_time, xysmdn_time")





