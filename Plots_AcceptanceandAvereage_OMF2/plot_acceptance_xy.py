from timeit import default_timer as timer
from subprocess import call
import matplotlib.pyplot as plt
import math
import numpy as np
import glob

# Initialize variables
dataacc = 'data_acceptance.csv'
xy_acc = np.genfromtxt(dataacc, delimiter=',')[:,2]
xysmd_acc = np.genfromtxt(dataacc, delimiter=',')[:,3]

xvalues = np.array([2, 3, 4, 5, 6, 7, 10, 15, 20, 30, 40, 50])
eps2 = (1.0/xvalues)*(1.0/xvalues)
eps = (1.0/xvalues)

all_acc = np.asarray([xy_acc, xysmd_acc])
all_acceps = 1/(eps*(all_acc))

# Plot data
plt.figure(0)
plt.xlabel(r'$\epsilon^2$', fontsize=13)
plt.ylabel(r'$Acc. \mathrm{in} \ \%$', fontsize=13)
plt.plot(eps2, all_acc[0], marker = 'o', label = 'XYHMC w. Metro')
plt.plot(eps2, all_acc[1], marker = 'o', label = 'XYSMD w. Metro')
plt.legend(loc = 'lower left')

plt.figure(1)
plt.plot(eps2, all_acceps[0], marker = 'o', label = 'XYHMC w. Metro')
plt.plot(eps2, all_acceps[1], marker = 'o', label = 'XYSMD w. Metro')
plt.xlabel(r'$\epsilon^2$', fontsize=13)
plt.ylabel(r'$1/(\epsilon(acc.))$', fontsize=13)
plt.legend(loc = 'upper right')

plt.show()





