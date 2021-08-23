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
import heat as ht
import numpy as np
import netCDF4 as nc
import datetime as dt
from . import pars_ParFlowTCL as ppfl
from . import ParFlow_IO as pio
import matplotlib.pyplot as plt

def calc_wtd(press, cellDepths):
    """ calculate the water table depth

    assuming dimensions:
    dim=3D: (z, y, x)
    dim=4D: (t, z, y, x)
    """
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

def calc_gwr_v1(spw, wtd, cellCenterDepth3D):
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

def calc_gwr_v2(spw, wtd_z_index):
    """
    see also:
    https://www.mi.fu-berlin.de/en/math/groups/ag-numerik/download/dateien/VL_Engelhardt_3.pdf
    """
    spwShape = spw.shape
    nx = spwShape[-1]
    ny = spwShape[-2]
    nz = spwShape[-3]
    gwr = ht.zeros((ny, nx))
    # wtd_z_index is the index corrospond to first cell outsid the 
    # groundwaterbody which is to to my (NWR) unserstanding that lvl
    # contributing to ground-water-recharge
    for z in range(nz):
        gwr = ht.where(wtd_z_index==z, spw[z], gwr)
    return gwr

def get_3Dgroundwaterbody_mask(satur):
    ''' calculating a 3D mask of the groundwater-body

    The groundwater-body is that volume of a soil column which is continuously 
    fully saturated seen from the very bottom of the model domain (bedrock).

    Return:
      gwb_mask    = ht-ndarray
        A HeAT-array of the same shape as 'satur' input, holding 1 at pixels
        belonging to the groundwater-body, and 0 for other pixel
      wtd_z_index = ht-ndarray
        A HeAT-array of the same 2D shape as 'satur' (no z dimension),
        holding the index of the first (counted from bottom/bedrock to 
        top/surface) unsaturated layer for each soil column. If there are N 
        layers with the model and the index has a value of N+1, than the entire 
        column is saturated!
    '''
    nz, ny, nx = satur.shape
    gwb        = ht.zeros((nz+1, ny, nx))
    gwb_mask   = ht.zeros((nz, ny, nx))

    # set unsaturated pixel to 0 and saturated pixel to 1
    # add 'zero' level on top ( ht.zeros((nz+1, ny, nx)) )
    # to guarantee finding a min with ht.argmin()
    gwb[:-1] = ht.where(satur>=1, 1, 0)
    # get index of first unsaturated cell in columne, countet from model bottom
    # this index corrospond to first cell outsid the groundwaterbody
    wtd_z_index = ht.argmin(gwb,axis=0)
    for z in range(nz):
        # with `z<wtd_z_index` we do get fully saturated cells only, NOT 
        # including the first unsaturated layer (as we want to!)
        gwb_mask[z] = ht.where(z<wtd_z_index, 1, 0)

    return gwb_mask, wtd_z_index

