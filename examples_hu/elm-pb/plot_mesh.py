from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
import matplotlib.pyplot as plt
from boututils.file_import import file_import
from matplotlib.ticker import FixedFormatter, FormatStrFormatter, AutoLocator, AutoMinorLocator

g = file_import("data/cbm18_dens8.grid_nx68ny64.nc")

x_temp=g.get("Rxy")
y_temp=g.get("Zxy")

#print x_temp.shape
nR = x_temp.shape[0]
nZ = x_temp.shape[1]

for r in np.arange(nR):
	plt.plot(x_temp[r,:], y_temp[r,:], color='black', linewidth=0.2)

for z in np.arange(nZ):
	plt.plot(x_temp[:,z], y_temp[:,z], color='black', linewidth=0.2)

plt.savefig('mesh.pdf')
plt.show()
