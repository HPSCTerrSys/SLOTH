"""analysis - submodule of SLOTH

Description:
analysis.py is aimed to hold standalone functions used for analysis purpose.
Basically this are functions calculating some variabels.

[...]
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

def vanGenuchten(refP, sSat, sRes, nVanG, aVanG):
    """  Calculates the degree of saturation as a function of the pressure head.

    The degree of saturation is calculated as a function of the pressure head
    according to M. Th. van Genuchten.
    Name: A Closed‐form Equation for Predicting the Hydraulic Conductivity of
            Unsaturated Soils
    DOI: https://doi.org/10.2136/sssaj1980.03615995004400050002x

    Parameters:
    refP:       Pressure head [L]
    sSat:       Relative saturated water content [-]
    sRes:       Relative residual saturation [-]
    nVanG:      Non linearity coefficient of the soil [-]
    aVanG:      Air entry values of the soil [L^-1]

    Returns:
    vanG:

    """
    mVanG = 1 - (1 / nVanG)

    vanG = ( (sSat - sRes) / ( 1 + (aVanG * np.absolute(refP))**(nVanG) )**mVanG ) + sRes

    vanG = np.where(refP>=0., 1., vanG) # avoid unsaturated values where water is ponding

    return vanG

