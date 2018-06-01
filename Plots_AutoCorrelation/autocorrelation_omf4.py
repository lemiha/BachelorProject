import numpy as np
import matplotlib.pyplot as plt
import glob

# Initialize variables
hmc_names = np.sort(glob.glob("hmc_omf4/autoh*.dat"))
hmcn_names = np.sort(glob.glob("hmcn_omf4/autoh*.dat"))
smd_names = np.sort(glob.glob("smd_omf4/autos*.dat"))
smdn_names = np.sort(glob.glob("smdn_omf4/autos*.dat"))
xy_names = np.sort(glob.glob("xy_omf4/autox*.dat"))
xyn_names = np.sort(glob.glob("xyn_omf4/autox*.dat"))
xysmd_names = np.sort(glob.glob("xysmd_omf4/autox*.dat"))
xysmdn_names = np.sort(glob.glob("xysmdn_omf4/autox*.dat"))

all_names = np.array([hmc_names, hmcn_names, smd_names, smdn_names, xy_names, xyn_names, xysmd_names, xysmdn_names])
all_tauint = []

for names in all_names:
    tmp = []
    for name in names:
        autoc = np.genfromtxt(name, delimiter='  ')[: ,0]
        tauint = np.genfromtxt(name, delimiter='  ')[: ,1]
        for i in range(len(autoc)):
            if autoc[i] < 0:
                tmp.append(tauint[i-1])
                break
    all_tauint.append(tmp)

np.savetxt("data_tauint_omf4.csv",  np.c_[all_tauint[0], all_tauint[1], all_tauint[2], all_tauint[3], all_tauint[4], all_tauint[5], all_tauint[6], all_tauint[7]],
        delimiter=',', fmt='%1.6f', header="hmc_tauint, hmcn_tauint, smd_tauint, smdn_tauint, xy_tauint, xyn_tauint, xysmd_tauint, xysmdn_tauint")




# plt.plot(range(len(tauint)), tauint)

# plt.show()