""" Example script to apply SanityCheck().

In many cases it is needed to quickly inspect some data-sets. 
[ADD SOME DESCRIPTION HERE]
"""
import numpy as np
import netCDF4 as nc
import sys
import os

sloth_path='../'
sys.path.append(sloth_path)
import sloth.SanityCheck


###############################################################################
### Define some paths, filenames, etc
###############################################################################

dataRootDir = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType    = 'postpro'
dataYear    = '1979'
varName     = 'T_S'
fileNames   = [
        f'{dataRootDir}/{datasetName}/{procType}/1979_12/{varName}_ts.nc',
        f'{dataRootDir}/{datasetName}/{procType}/1980_01/{varName}_ts.nc',
        f'{dataRootDir}/{datasetName}/{procType}/1980_02/{varName}_ts.nc'
        ]

###############################################################################
### read in data and save as ndarray
###############################################################################
var      = []
var_mask = []
for fileName in fileNames:
    with nc.Dataset(f'{fileName}', 'r') as nc_file:
        nc_var   = nc_file.variables[f'{varName}']
        tmp_var  = nc_var[...]
        tmp_var_mask  = tmp_var.mask
        var.append(tmp_var)
        var_mask.append(tmp_var_mask)
var      = np.concatenate(var, axis=0)
var_mask = np.concatenate(var_mask, axis=0)
# Merging array in above style is removing the masked attached to the array.
# So readd the maske:
var      = np.ma.masked_where(var_mask, var)


###############################################################################
### Start SanityCheck
###############################################################################
# define some title for plot, which can be passed via functionarguments
# see funciton definition for full potential
fig_title    = f'Sanity-Check for {varName} --- {dataYear} (DJF)'
figname      = f'./examples_SanityCheck_Season.pdf'
minax_title  = f'{varName} min'
maxax_title  = f'{varName} max'
kinax_title  = f'{varName} mean'
hisax_title  = f'{varName} mean - distribution'

# For mor detailed information about how plot_SanityCheck() does work, see
# sloth/SanityCheck.py --> plot_SanityCheck()
sloth.SanityCheck.plot_SanityCheck(data=var, 
        # below is optional
        kind='mean', figname=figname,
        lowerP=2, upperP=98, interactive=False,
        fig_title=fig_title, minax_title=minax_title, maxax_title=maxax_title, 
        kinax_title=kinax_title, hisax_title=hisax_title, cmapName='Spectral_r')
