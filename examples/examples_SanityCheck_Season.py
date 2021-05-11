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
import sloth


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
        # netCDF4 returns masked-array, which is nice, but one should stick to one
        # way handling missing values and I do prefere np.nan, so filling
        # var to get pure numpy ndarray
        # However we can keep the mask, to keep the information where masked values
        # are / were - can be helpfull
        tmp_var_mask  = tmp_var.mask
        # handling if no values are masked (one aspect why I prefere handling of
        # missing values via np.nan...
        if not tmp_var_mask.any():
            tmp_var_mask  = np.zeros(tmp_var.shape, dtype=bool)
        tmp_var = tmp_var.filled(fill_value=np.nan)
        var.append(tmp_var)
        var_mask.append(tmp_var_mask)
var      = np.concatenate(var, axis=0)
var_mask = np.concatenate(var_mask, axis=0)


###############################################################################
### Start SanityCheck 3D
###############################################################################
# define some title for plot, which can be passed via functionarguments
# see funciton definition for full potential
fig_title    = f'Sanity-Check for {varName} --- {dataYear} (DJF)'
figname      = f'./examples_SanityCheck_Season.pdf'
minax_title  = f'{varName} min'
maxax_title  = f'{varName} max'
kinax_title  = f'{varName} mean'
hisax_title  = f'{varName} mean - distribution'

# use 'plot_SanityCheck_3D' from script sloth//SanityCheck.py imported above
sloth.SanityCheck.plot_SanityCheck_3D(data=var, 
        # below is optional
        data_mask=var_mask, kind='mean', figname=figname,
        lowerP=2, upperP=98, interactive=False,
        fig_title=fig_title, minax_title=minax_title, maxax_title=maxax_title, 
        kinax_title=kinax_title, hisax_title=hisax_title, cmapName='Spectral_r')
