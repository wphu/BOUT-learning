#!/usr/bin/env python
"""
An animated image
"""

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#python 34
from netCDF4 import Dataset
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from boutdata import *

#read grid data
f=Dataset('data/slab.grd.nc')
x_temp=f.variables['Rxy']
y_temp=f.variables['Zxy']

#read simulation data
z_temp=collect("Ni",path="data")
print "222",type(z_temp),z_temp.shape

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
		x[i]=x_temp[i,0]

	ax1d.plot(x, z_temp[time,:,0])
	fig.savefig('line1d.png')

plt.show()
