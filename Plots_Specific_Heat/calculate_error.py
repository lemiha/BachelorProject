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

hmc_error = []
hmcn_error = []
smd_error = []
smdn_error = []
xy_error = []
xyn_error = []
xysmd_error = []
xysmdn_error = []

# Calculate the averages
nspace = 500 # Bin size
for i in range(len(hmc_names)):
    tmphmc = np.genfromtxt(hmc_names[i], delimiter='  ')[:,1]
    hmcinterval = np.linspace(0, len(tmphmc)-1, num = nspace)
    hmcstd = []
    for j in range(1, len(hmcinterval)):
        start = int(hmcinterval[j-1])
        end = int(hmcinterval[j])
        hmcstd.append(np.mean(tmphmc[start:end]))
    hmc_error.append(np.std(hmcstd)/np.sqrt(nspace))

for i in range(len(hmcn_names)):
    tmphmcn = np.genfromtxt(hmcn_names[i], delimiter='  ')[:,1]
    hmcninterval = np.linspace(0, len(tmphmcn)-1, num = nspace)
    hmcnstd = []
    for j in range(1, len(hmcninterval)):
        start = int(hmcninterval[j-1])
        end = int(hmcninterval[j])
        hmcnstd.append(np.mean(tmphmcn[start:end]))
    hmcn_error.append(np.std(hmcnstd)/np.sqrt(nspace))

for i in range(len(smd_names)):
    tmpsmd = np.genfromtxt(smd_names[i], delimiter='  ')[:,1]
    smdinterval = np.linspace(0, len(tmpsmd)-1, num = nspace)
    smdstd = []
    for j in range(1, len(smdinterval)):
        start = int(smdinterval[j-1])
        end = int(smdinterval[j])
        smdstd.append(np.mean(tmpsmd[start:end]))
    smd_error.append(np.std(smdstd)/np.sqrt(nspace))

for i in range(len(smdn_names)):
    tmpsmdn = np.genfromtxt(smdn_names[i], delimiter='  ')[:,1]
    smdninterval = np.linspace(0, len(tmpsmdn)-1, num = nspace)
    smdnstd = []
    for j in range(1, len(smdninterval)):
        start = int(smdninterval[j-1])
        end = int(smdninterval[j])
        smdnstd.append(np.mean(tmpsmdn[start:end]))
    smdn_error.append(np.std(smdnstd)/np.sqrt(nspace))

for i in range(len(xy_names)):
    tmpxy = np.genfromtxt(xy_names[i], delimiter='  ')[:,1]
    xyinterval = np.linspace(0, len(tmpxy)-1, num = nspace)
    xystd = []
    for j in range(1, len(xyinterval)):
        start = int(xyinterval[j-1])
        end = int(xyinterval[j])
        xystd.append(np.mean(tmpxy[start:end]))
    xy_error.append(np.std(xystd)/np.sqrt(nspace))

for i in range(len(xyn_names)):
    tmpxyn = np.genfromtxt(xyn_names[i], delimiter='  ')[:,1]
    xyninterval = np.linspace(0, len(tmpxyn)-1, num = nspace)
    xynstd = []
    for j in range(1, len(xyninterval)):
        start = int(xyninterval[j-1])
        end = int(xyninterval[j])
        xynstd.append(np.mean(tmpxyn[start:end]))
    xyn_error.append(np.std(xynstd)/np.sqrt(nspace))

for i in range(len(xysmd_names)):
    tmpxysmd = np.genfromtxt(xysmd_names[i], delimiter='  ')[:,1]
    xysmdinterval = np.linspace(0, len(tmpxysmd)-1, num = nspace)
    xysmdstd = []
    for j in range(1, len(xysmdinterval)):
        start = int(xysmdinterval[j-1])
        end = int(xysmdinterval[j])
        xysmdstd.append(np.mean(tmpxysmd[start:end]))
    xysmd_error.append(np.std(xysmdstd)/np.sqrt(nspace))

for i in range(len(xysmdn_names)):
    tmpxysmdn = np.genfromtxt(xysmdn_names[i], delimiter='  ')[:,1]
    xysmdninterval = np.linspace(0, len(tmpxysmdn)-1, num = nspace)
    xysmdnstd = []
    for j in range(1, len(xysmdninterval)):
        start = int(xysmdninterval[j-1])
        end = int(xysmdninterval[j])
        xysmdnstd.append(np.mean(tmpxysmdn[start:end]))
    xysmdn_error.append(np.std(xysmdnstd)/np.sqrt(nspace))

np.savetxt("data_errorhmc.csv",  np.c_[hmc_error],
        delimiter=',', fmt='%1.6f', header="hmc_error")

np.savetxt("data_errorxy.csv",  np.c_[xy_error],
        delimiter=',', fmt='%1.6f', header="xy_error")