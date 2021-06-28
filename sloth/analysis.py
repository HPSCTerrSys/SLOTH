"""analysis - submodule of SLOTH

author: Nikals WAGNER
e-mail: n.wagner@fz-juelich.de
version: 2021-06-24

Description:
tooBox.py is aimed to hold standalone functions used / developed / found
within the realm of SLOTH development.
[...]

Check examples/ to see how this works. 
"""
import sys
import os
import numpy as np
import netCDF4 as nc
import datetime as dt
from . import pars_ParFlowTCL as ppfl
from . import ParFlow_IO as pio
import matplotlib.pyplot as plt

def calc_wtd(press, cellDepths):
    totalColumnDepth = ht.sum(cellDepths)
    dim = press.ndim
    if dim==3:
        wtd = totalColumnDepth - (press[0] + cellDepths[0] * 0.5)
    elif dim==4:
        wtd = totalColumnDepth - (press[:,0,...] + cellDepths[0] * 0.5)
    else:
        print(f'ERROR: dim of press input is not supported --> EXIT!')
        sys.exit()
    # Above could be also solved by below oneliner:
    # >> wtd = totalColumnDepth - (press[...,0,:,:] + cellDepth[0] * 0.5)
    # but maybe this is not as intuitive as the used if-clause?
    wtd = ht.where((wtd < 0), 0, wtd)
    return wtd

def calc_gwr(spw, wtd, cellCenterDepth3D):
    nz, ny, nx = cellCenterDepth3D.shape
    gwr = ht.zeros((ny, nx))
    # wtd_z_index should contain the z-index of that level which does
    # contribute to the ground water recharge. To my (NWR) unserstanding
    # this is the index of the next cell center above the WTD
    # 1. Calculate diff between WTD and cellCenterDepth
    tmp_1 = wtd-cellCenterDepth3D
    # 2. 'mask' those cells where WTD is above (masking = set to large value)
    tmp_1 = ht.where(tmp_1 <= 0, 1000, tmp_1)
    # 3. find next / closest cellCenter above WTD
    wtd_z_index = ht.argmin(tmp_1,axis=0)
    # 4. get spw at wtd_z_index which is gwr
    for z in range(nz):
        gwr = ht.where(wtd_z_index==z, spw[z], gwr)
    return gwr

