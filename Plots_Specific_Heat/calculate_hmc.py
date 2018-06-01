from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

names = np.sort(glob.glob("hmc/h*.dat"))

kappa = []
susp = []

for name in names:
    fileHandle = open(name, "r")
    lineList = fileHandle.readlines()
    tmp = lineList[1].split(" ")[2]
    kappa.append(float(tmp.split("=",1)[1]))

kappa = np.asarray(kappa)
kappa = 1/kappa

for i in range(len(names)):
    tmp = np.genfromtxt(names[i], delimiter='  ')[:,1]
    tmpSquared = tmp**2
    susp.append(np.mean(tmpSquared) - np.mean(tmp)**2)

np.savetxt("data_hmc.csv",  np.c_[kappa, susp],
        delimiter=',', fmt='%1.6f', header="kappa, susp")