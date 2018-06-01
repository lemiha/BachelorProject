import numpy as np
import matplotlib.pyplot as plt
import glob

# Initialize variables
hmc_names_omf2 = np.sort(glob.glob("hmc_omf2/autoh*.dat"))
hmc_names_omf4 = np.sort(glob.glob("hmc_omf4/autoh*.dat"))

autoc = np.genfromtxt(hmc_names_omf4[4], delimiter='  ')[: ,0]
tauint = np.genfromtxt(hmc_names_omf4[4], delimiter='  ')[: ,1]


# hmcn_names = np.sort(glob.glob("hmcn_omf2/autoh*.dat"))
# smd_names = np.sort(glob.glob("smd_omf2/autos*.dat"))
# smdn_names = np.sort(glob.glob("smdn_omf2/autos*.dat"))
# xy_names = np.sort(glob.glob("xy_omf2/autox*.dat"))
# xyn_names = np.sort(glob.glob("xyn_omf2/autox*.dat"))
# xysmd_names = np.sort(glob.glob("xysmd_omf2/autox*.dat"))
# xysmdn_names = np.sort(glob.glob("xysmdn_omf2/autox*.dat"))
cutoff = 3000

# Autocorrelation
plt.figure(0)
plt.xlabel('Iterations ' r'$t$', fontsize=13)
plt.ylabel(r'$\tau_{\rm{int}}(t)$', fontsize=13)
plt.plot(range(cutoff), tauint[0:cutoff])

# Autocorrelation time
plt.figure(1)
plt.xlabel('Iterations ' r'$t$', fontsize=13)
plt.ylabel(r'$R(t)$', fontsize=13)
plt.plot(range(cutoff), autoc[0:cutoff])

plt.show()