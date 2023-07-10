"""analysis - submodule of SLOTH

Description:
analysis.py is aimed to hold standalone functions used for analysis purpose.
Basically this are functions calculating some variabels.

[TBE]
"""
import sys
import os
import numpy as np
try:
    import heat as ht
except ImportError:
    print(f'ERROR SLOTH.analysis --> Could not import heat, so heat versions of analysis functions are not available.')
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt


def calc_wtd_HT(press, cellDepths):
    """
    Calculate the water table depth.

    Parameters
    ----------
    press : htdarray
        Array of pressure values. The dimensions depend on the scenario:
        - For 3D: (z, y, x)
        - For 4D: (t, z, y, x)

    cellDepths : ndarray
        Array of cell depths.

    Returns
    -------
    wtd : ndarray
        Array of water table depths.

    Notes
    -----
    The function calculates the water table depth based on the given pressure values
    and cell depths. The water table depth is determined by subtracting the lower most 
    pressure values and cell depths from the total column depth.

    - For 3D input, the water table depth is calculated as:
        wtd = totalColumnDepth - (press[0] + cellDepths[0] * 0.5)

    - For 4D input, assuming a time series of pressure values, the water table depth is calculated as:
        wtd = totalColumnDepth - (press[:, 0, ...] + cellDepths[0] * 0.5)

    Negative water table depths are set to 0 using `np.where()`.

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

def calc_wtd(press, cellDepths):
    """
    Calculate the water table depth.

    Parameters
    ----------
    press : ndarray
        Array of pressure values. The dimensions depend on the scenario:
        - For 3D: (z, y, x)
        - For 4D: (t, z, y, x)

    cellDepths : ndarray
        Array of cell depths.

    Returns
    -------
    wtd : ndarray
        Array of water table depths.

    Notes
    -----
    The function calculates the water table depth based on the given pressure values
    and cell depths. The water table depth is determined by subtracting the lower most 
    pressure values and cell depths from the total column depth.

    - For 3D input, the water table depth is calculated as:
        wtd = totalColumnDepth - (press[0] + cellDepths[0] * 0.5)

    - For 4D input, assuming a time series of pressure values, the water table depth is calculated as:
        wtd = totalColumnDepth - (press[:, 0, ...] + cellDepths[0] * 0.5)

    Negative water table depths are set to 0 using `np.where()`.

    """
    totalColumnDepth = np.sum(cellDepths)
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
    wtd = np.where((wtd < 0), 0, wtd)
    return wtd

def get_3Dgroundwaterbody_mask_HT(satur):
    """
    Calculate a 3D mask of the groundwater body.

    The groundwater body refers to the volume of a soil column that is continuously
    fully saturated, starting from the very bottom of the model domain (bedrock).

    Parameters
    ----------
    satur : htdarray
        Input array representing the saturation level of the soil column.
        Shape: (nz, ny, nx)

    Returns
    -------
    gwb_mask : ndarray
        Numpy array of the same shape as the 'satur' input, where each pixel
        belonging to the groundwater body is assigned a value of 1, and 0 otherwise.
        Shape: (nz, ny, nx)

    wtd_z_index : ndarray
        Numpy array of the same 2D shape as the 'satur' input (without the z dimension).
        It holds the index of the first unsaturated layer, counted from the bottom/bedrock
        to the top/surface, for each soil column. If the index has a value of N+1, where N
        is the number of layers in the model, the entire column is saturated.

    """
    nz, ny, nx = satur.shape
    gwb        = ht.zeros((nz+1, ny, nx))
    gwb_mask   = ht.zeros((nz, ny, nx))

    # set unsaturated pixel to 0 and saturated pixel to 1
    # add 'zero' level on top ( np.zeros((nz+1, ny, nx)) )
    # to guarantee finding a min with np.argmin()
    gwb[:-1] = ht.where(satur>=1, 1, 0)
    # get index of first unsaturated cell in columne, countet from model bottom
    # this index corrospond to first cell outsid the groundwaterbody
    wtd_z_index = ht.argmin(gwb,axis=0)
    for z in range(nz):
        # with `z<wtd_z_index` we do get fully saturated cells only, NOT 
        # including the first unsaturated layer (as we want to!)
        gwb_mask[z] = ht.where(z<wtd_z_index, 1, 0)

    return gwb_mask, wtd_z_index
def get_3Dgroundwaterbody_mask(satur):
    """
    Calculate a 3D mask of the groundwater body.

    The groundwater body refers to the volume of a soil column that is continuously
    fully saturated, starting from the very bottom of the model domain (bedrock).

    Parameters
    ----------
    satur : ndarray
        Input array representing the saturation level of the soil column.
        Shape: (nz, ny, nx)

    Returns
    -------
    gwb_mask : ndarray
        Numpy array of the same shape as the 'satur' input, where each pixel
        belonging to the groundwater body is assigned a value of 1, and 0 otherwise.
        Shape: (nz, ny, nx)

    wtd_z_index : ndarray
        Numpy array of the same 2D shape as the 'satur' input (without the z dimension).
        It holds the index of the first unsaturated layer, counted from the bottom/bedrock
        to the top/surface, for each soil column. If the index has a value of N+1, where N
        is the number of layers in the model, the entire column is saturated.

    """
    nz, ny, nx = satur.shape
    gwb        = np.zeros((nz+1, ny, nx))
    gwb_mask   = np.zeros((nz, ny, nx))

    # set unsaturated pixel to 0 and saturated pixel to 1
    # add 'zero' level on top ( np.zeros((nz+1, ny, nx)) )
    # to guarantee finding a min with np.argmin()
    gwb[:-1] = np.where(satur>=1, 1, 0)
    # get index of first unsaturated cell in columne, countet from model bottom
    # this index corrospond to first cell outsid the groundwaterbody
    wtd_z_index = np.argmin(gwb,axis=0)
    for z in range(nz):
        # with `z<wtd_z_index` we do get fully saturated cells only, NOT 
        # including the first unsaturated layer (as we want to!)
        gwb_mask[z] = np.where(z<wtd_z_index, 1, 0)

    return gwb_mask, wtd_z_index

def vanGenuchten(refP, sSat, sRes, nVanG, aVanG):
    """
    Calculate the degree of saturation as a function of the pressure head.

    The degree of saturation is calculated as a function of the pressure head
    according to the M. Th. van Genuchten equation.

    Reference:
    M. Th. van Genuchten. "A Closed‐form Equation for Predicting the Hydraulic Conductivity
    of Unsaturated Soils." Soil Science Society of America Journal, 1980.
    DOI: https://doi.org/10.2136/sssaj1980.03615995004400050002x

    Parameters
    ----------
    refP : ndarray
        Pressure head. [L]

    sSat : float
        Relative saturated water content. [-]

    sRes : float
        Relative residual saturation. [-]

    nVanG : float
        Non-linearity coefficient of the soil. [-]

    aVanG : float
        Air entry value of the soil. [L^-1]

    Returns
    -------
    vanG : ndarray
        Degree of saturation.

    Notes
    -----
    This function calculates the degree of saturation as a function of the pressure head
    using the van Genuchten equation. The van Genuchten equation relates the degree of
    saturation to the pressure head and soil-specific parameters.

    The equation is given by:

        vanG = ((sSat - sRes) / (1 + (aVanG * np.absolute(refP)) ** nVanG) ** mVanG) + sRes

    Where:
    - sSat is the relative saturated water content,
    - sRes is the relative residual saturation,
    - nVanG is the non-linearity coefficient of the soil,
    - aVanG is the air entry value of the soil.

    The coefficient mVanG is calculated as mVanG = 1 - (1 / nVanG).

    To avoid unsaturated values where water is ponding, the output vanG is adjusted to 1
    wherever refP is greater than or equal to 0.

    """
    mVanG = 1 - (1 / nVanG)

    vanG = ( (sSat - sRes) / ( 1 + (aVanG * np.absolute(refP))**(nVanG) )**mVanG ) + sRes

    vanG = np.where(refP>=0., 1., vanG) # avoid unsaturated values where water is ponding

    return vanG

