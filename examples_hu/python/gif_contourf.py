#!/usr/bin/env python
"""
An animated image
"""
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

f=Dataset('data/bout.grd.east50.nc')
x_temp=f.variables['Rxy']
y_temp=f.variables['Zxy']
print "111",x_temp[:,:]

x=np.zeros((50,64))
y=np.zeros((50,64))
z=np.zeros((50,64))
z_temp=collect("Vi",path="data")
print "222",type(z_temp)
#print z_temp.shape
z=z_temp[100,:,:,0]
for i in np.arange(50):
	for j in np.arange(64):
		x[i,j]=x_temp[i,j]
		y[i,j]=y_temp[i,j]


levels = MaxNLocator(nbins=60).tick_values(z.min(), z.max())
#cmap = plt.get_cmap('PiYG')
#cmap = plt.get_cmap('spectral')
cmap = plt.get_cmap('jet')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#y,x=np.mgrid[slice(1,ny+dy,dy),slice(1,nx+dx,dx)]

fig = plt.figure(figsize=(6,12))

#x = np.linspace(0, 2 * np.pi, 120)
#y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []
#for i in range(200):
#	im = plt.pcolormesh(x,y,z_temp[i,:,:,0],cmap=cmap,norm=norm)
	#im = plt.contourf(x,y,z_temp[i,:,:,0],levels=levels,cmap=cmap,antialiased=True)
	#im = plt.contourf(x,y,z_temp[i,:,:,0],cmap=cmap,norm=norm)
#	ims.append([im])
#plt.colorbar()
#plt.title('Te/eV')
def init():
	plt.contourf(x,y,z_temp[0,:,:,0],levels=levels,cmap=cmap,antialiased=True)
	#if(nframe==0):
	plt.colorbar()
	plt.title('Te/eV')

def animate(nframe):
	plt.contourf(x,y,z_temp[nframe,:,:,0],levels=levels,cmap=cmap,antialiased=True)
	#if(nframe==0):
	#plt.colorbar()
	#plt.title('Te/eV')
#ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,repeat_delay=1000)
ani = animation.FuncAnimation(fig,animate,frames=200, interval=50, blit=True,init_func=init)

#ani.save('dynamic_images.mp4')
#anim = animation.FuncAnimation(fig, ims,
#                            frames=100, interval=20, blit=True)
 
# this is how you save your animation to file:
ani.save('bout_movie.gif', writer='imagemagick_file', fps=30)
#plt.show()
