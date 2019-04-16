from __future__ import division
from builtins import range
from past.utils import old_div
from boutdata.collect import collect
import numpy as np
import matplotlib.pyplot as plt
from boututils.file_import import file_import
from matplotlib.ticker import FixedFormatter, FormatStrFormatter, AutoLocator, AutoMinorLocator


#read grid data
g = file_import('data/slab.grd.nc')
x_temp=g.get("Rxy")
y_temp=g.get("Zxy")

#read simulation data
z_temp=collect("Ni",path="data")


nx = z_temp.shape[1]
ny = z_temp.shape[2]



time = 10
plot = "1d"

fig = plt.figure(figsize=(9,6), dpi=100)


if plot == "2d":
	ax2d = fig.add_subplot(111, aspect='equal')
	x=np.zeros((nx,ny))
	y=np.zeros((nx,ny))

	z=z_temp[10,:,:]

	for i in np.arange(nx):
		for j in np.arange(ny):
			x[i,j]=x_temp[i,j]
			y[i,j]=y_temp[i,j]


	levels = MaxNLocator(nbins=60).tick_values(z.min(), z.max())
	cmap = plt.get_cmap('jet')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)


	cf = ax2d.contourf(x,y,z_temp[time,:,:],levels=levels,cmap=cmap,antialiased=True)
	ax2d.axis([x.min(),x.max(),y.min(),y.max()])
	cbar = fig.colorbar(cf)
	fig.savefig('contour2d.png')

elif plot == "1d":
	ax1d = fig.add_subplot(111)
	x=np.zeros(nx)
	for i in np.arange(nx):
		x[i] = i
	ntime = (z_temp.shape)[0] - 1
	for time in np.linspace(0, ntime, 5):
		ax1d.plot(x, z_temp[time,:,0], label = "time " + str(int(time)))

	#plt.xlim((0,20))
	plt.legend()
	fig.savefig('line1d.png')


plt.show()
