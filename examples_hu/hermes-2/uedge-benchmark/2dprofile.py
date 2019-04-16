from __future__ import division
from builtins import range
from past.utils import old_div
from boutdata.collect import collect
import numpy as np
import matplotlib.pyplot as plt
from boututils.file_import import file_import
from matplotlib.ticker import FixedFormatter, FormatStrFormatter, AutoLocator, AutoMinorLocator

g = file_import("uedge.grd_Up_Ni_Tei_2d.nc")
var="Pe"
unit="eV"
t=900



x_temp=g.get("Rxy")
y_temp=g.get("Zxy")


nx=g.get("nx")
ny=g.get("ny")
Ni_x=g.get("Ni_x")
Ti_x=g.get("Ti_x")
Te_x=g.get("Te_x")
Vi_x=g.get("Vi_x")


x=np.zeros((nx,ny))
y=np.zeros((nx,ny))
z=np.zeros((nx,ny))
z_temp=collect(var,path="./")


z=z_temp[t,:,:,0]*Te_x
print(z_temp.shape)

for i in np.arange(nx):
	for j in np.arange(ny):
		x[i,j]=x_temp[i,j]
		y[i,j]=y_temp[i,j]



cmap = plt.get_cmap('jet')


fig = plt.figure(figsize=(6,12))


plt.contourf(x,y,z,cmap=cmap,antialiased=True)
plt.colorbar()
plt.axis([x.min(),x.max(),y.min(),y.max()])
plt.title(var+'/'+unit)


plt.savefig('contour2d.png')
plt.show()
