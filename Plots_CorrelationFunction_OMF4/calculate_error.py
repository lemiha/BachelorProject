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

hmc_error = []
hmcn_error = []
smd_error = []
smdn_error = []
xy_error = []
xyn_error = []
xysmd_error = []
xysmdn_error = []

# Calculate the averages
bins = 500
for i in range(len(hmc_names)):
    print(i)
    tmphmc = np.genfromtxt(hmc_names[i], delimiter='  ')[:,3]
    hmcinterval = np.linspace(0, len(tmphmc)-1, num = bins)
    hmcstd = []
    for j in range(1, len(hmcinterval)):
        start = int(hmcinterval[j-1])
        end = int(hmcinterval[j])
        hmcstd.append(np.mean(tmphmc[start:end]))
    hmc_error.append(np.std(hmcstd)/np.sqrt(bins))

    tmphmcn = np.genfromtxt(hmcn_names[i], delimiter='  ')[:,3]
    hmcninterval = np.linspace(0, len(tmphmcn)-1, num = bins)
    hmcnstd = []
    for j in range(1, len(hmcninterval)):
        start = int(hmcninterval[j-1])
        end = int(hmcninterval[j])
        hmcnstd.append(np.mean(tmphmcn[start:end]))
    hmcn_error.append(np.std(hmcnstd)/np.sqrt(bins))

    tmpsmd = np.genfromtxt(smd_names[i], delimiter='  ')[:,3]
    smdinterval = np.linspace(0, len(tmpsmd)-1, num = bins)
    smdstd = []
    for j in range(1, len(smdinterval)):
        start = int(smdinterval[j-1])
        end = int(smdinterval[j])
        smdstd.append(np.mean(tmpsmd[start:end]))
    smd_error.append(np.std(smdstd)/np.sqrt(bins))

    tmpsmdn = np.genfromtxt(smdn_names[i], delimiter='  ')[:,3]
    smdninterval = np.linspace(0, len(tmpsmdn)-1, num = bins)
    smdnstd = []
    for j in range(1, len(smdninterval)):
        start = int(smdninterval[j-1])
        end = int(smdninterval[j])
        smdnstd.append(np.mean(tmpsmdn[start:end]))
    smdn_error.append(np.std(smdnstd)/np.sqrt(bins))

    tmpxy = np.genfromtxt(xy_names[i], delimiter='  ')[:,3]
    xyinterval = np.linspace(0, len(tmpxy)-1, num = bins)
    xystd = []
    for j in range(1, len(xyinterval)):
        start = int(xyinterval[j-1])
        end = int(xyinterval[j])
        xystd.append(np.mean(tmpxy[start:end]))
    xy_error.append(np.std(xystd)/np.sqrt(bins))

    tmpxyn = np.genfromtxt(xyn_names[i], delimiter='  ')[:,3]
    xyninterval = np.linspace(0, len(tmpxyn)-1, num = bins)
    xynstd = []
    for j in range(1, len(xyninterval)):
        start = int(xyninterval[j-1])
        end = int(xyninterval[j])
        xynstd.append(np.mean(tmpxyn[start:end]))
    xyn_error.append(np.std(xynstd)/np.sqrt(bins))

    tmpxysmd = np.genfromtxt(xysmd_names[i], delimiter='  ')[:,3]
    xysmdinterval = np.linspace(0, len(tmpxysmd)-1, num = bins)
    xysmdstd = []
    for j in range(1, len(xysmdinterval)):
        start = int(xysmdinterval[j-1])
        end = int(xysmdinterval[j])
        xysmdstd.append(np.mean(tmpxysmd[start:end]))
    xysmd_error.append(np.std(xysmdstd)/np.sqrt(bins))

    tmpxysmdn = np.genfromtxt(xysmdn_names[i], delimiter='  ')[:,3]
    xysmdninterval = np.linspace(0, len(tmpxysmdn)-1, num = bins)
    xysmdnstd = []
    for j in range(1, len(xysmdninterval)):
        start = int(xysmdninterval[j-1])
        end = int(xysmdninterval[j])
        xysmdnstd.append(np.mean(tmpxysmdn[start:end]))
    xysmdn_error.append(np.std(xysmdnstd)/np.sqrt(bins))

# Save data to csv
np.savetxt("data_errors.csv",  np.c_[hmc_error, hmcn_error, smd_error, smdn_error, xy_error, xyn_error, xysmd_error, xysmdn_error],
        delimiter=',', fmt='%1.6f', header="hmc_error, hmcn_error, smd_error, smdn_error, xy_error, xyn_error, xysmd_error, xysmdn_error")

# Plot data
# plt.plot(eps, hmc_ave, marker = 'o', label = 'HMC w. Metro')
# plt.plot(eps, hmcn_ave, marker = 'o', label = 'HMCN w/o. Metro')
# plt.plot(eps, smd_ave, marker = 'o', label = 'SMD w. Metro')
# plt.plot(eps, smdn_ave, marker = 'o', label = 'SMDN w/o. Metro')

# # plt.ylim((1.5,2.7))
# plt.legend(loc = 'upper left')

# plt.show()