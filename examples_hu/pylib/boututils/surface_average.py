from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from boututils import *
from .idl_tabulate import idl_tabulate
from bunch import bunchify


#from gen_surface import gen_surface


# Perform surface average
#
# var - [x,y,z] or [t,x,y,z]
#
#
# area=area : Average by flux-surface area = (B/Bp)*dl * R*dz
#
# By default, averages over poloidal angle such that
# surface_average(nu) = q
#

def surface_average ( var, g, area=None):

    s = np.ndim(var)



    if s == 4 :
        nx = np.shape(var)[1]
        ny = np.shape(var)[2]
        nt = np.shape(var)[0]

        result = np.zeros((nx,nt))
        for t in range (nt):

            result[:,t] = surface_average(var[t,:,:,:], g, area=area)

        return result
    elif s != 3 :
        print("ERROR: surface_average var must be 3 or 4D")
        return 0


  # 3D [x,y,z]
    nx = np.shape(var)[0]
    ny = np.shape(var)[1]
#    nz = np.shape(var)[2]

# Use bunch to create grid structure
    grid=bunchify(g)


  # Calculate poloidal angle from grid
    theta = np.zeros((nx,ny))

  #status = gen_surface(mesh=grid) ; Start generator
    xi = -1
    yi = np.arange(0,ny,dtype=int)
    last = 0
    while True:
    #yi = gen_surface(last=last, xi=xi, period=periodic)
        xi = xi + 1
        if xi == nx-1 :
            last = 1

        dtheta = 2.*np.pi / np.float(ny)
        r = grid.Rxy[xi,yi]
        z = grid.Zxy[xi,yi]
        n = np.size(r)

        dl = old_div(np.sqrt( deriv(r)**2 + deriv(z)**2 ), dtheta)
        if area:
            dA = (old_div(grid.Bxy[xi,yi],grid.Bpxy[xi,yi]))*r*dl
            A = int_func(np.arange(n),dA)
            theta[xi,yi] = 2.*np.pi*A/A[n-1]
        else:
            nu = dl * (grid.Btxy[xi,yi]) / ((grid.Bpxy[xi,yi]) * r )
            theta[xi,yi] = int_func(np.arange(n)*dtheta,nu)
            theta[xi,yi] = 2.*np.pi*theta[xi,yi] / theta[xi,yi[n-1]]

        if last==1 : break

    vy = np.zeros(ny)
    result = np.zeros(nx)
    for x in range(nx) :
        for y in range(ny) :
            vy[y] = np.mean(var[x,y,:])

        result[x] = old_div(idl_tabulate(theta[x,:], vy), (2.*np.pi))

    return result
