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

all_names = np.array([hmc_names, smd_names, xy_names, xysmd_names])
all_acc = []

# Gather the acc. ratios in an array
for names in all_names:
    tmp = []
    for name in names:
        fileHandle = open(name, "r")
        lineList = fileHandle.readlines()
        tmp.append(float(lineList[-1].split("=")[1].split(" ")[0]))
    all_acc.append(tmp)

all_acc = np.asarray(all_acc)

np.savetxt("data_acceptance.csv",  np.c_[all_acc[0], all_acc[1], all_acc[2], all_acc[3]],
        delimiter=',', fmt='%1.6f', header="hmc_ave, smd_ave, xy_ave, xysmd_ave")
