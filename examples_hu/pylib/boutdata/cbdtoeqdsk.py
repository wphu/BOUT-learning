from __future__ import print_function
from boututils.file_import import file_import
from bunch import Bunch
import numpy as np


def cbmtogeqdsk(g):
    gg=Bunch()
    gg.r=g['Rxy']
    gg.z=g['Zxy']
    gg.psi=g['psi']
    gg.pres=g['mu0p']
    gg.qpsi=g['qsafe']
    gg.fpol=g['f']
    gg.nx=g['nx']
    gg.ny=g['ny']
    i=np.argwhere(g['mu0p']==0)
    
    gg.simagx=gg.psi.min()
    gg.sibdry=gg.psi[i[0]]
    gg.xlim=0
    gg.ylim=0
    gg.nlim=0

    return gg

if __name__ == '__main__':
    gfile='../cbm18_dens8.dskgato.cdl'
    g=file_import(gfile)
    print(cbmtogeqdsk(g))
