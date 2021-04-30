""" Example script to apply SanityCheck().

In many cases it is needed to quickly inspect some data-sets. 
[ADD SOME DESCRIPTION HERE]
"""
import numpy as np
import netCDF4 as nc
import sys
import os

src_path='../src/'
sys.path.append(src_path)
import SanityCheck 

###############################################################################
### Define some paths, filenames, etc
###############################################################################

dataRootDir = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType    = 'postpro'
dataYear    = '1979_01'
varName     = 'T_S'
fileName    = f'{dataRootDir}/{datasetName}/{procType}/{dataYear}/{varName}_ts.nc'

###############################################################################
### read in data and save as ndarray
###############################################################################

with nc.Dataset(f'{fileName}', 'r') as nc_file:
    nc_var   = nc_file.variables[f'{varName}']
    var      = nc_var[...]
    # netCDF4 returns masked-array, which is nice, but one should stick to one
    # way handling missing values and I do prefere np.nan, so filling
    # var to get pure numpy ndarray
    # However we can keep the mask, to keep the information where masked values
    # are / were - can be helpfull
    var_mask  = var.mask
    # handling if no values are masked (one aspect why I prefere handling of
    # missing values via np.nan...
    if not var_mask.any():
        var_mask  = np.zeros(var.shape, dtype=bool)
    var = var.filled(fill_value=np.nan)


###############################################################################
### Start SanityCheck 3D
###############################################################################
# define some title for plot, which can be passed via functionarguments
# see funciton definition for full potential
fig_title    = f'Sanity-Check for {varName} --- {dataYear}'
figname      = f'./examples_SanityCheck_Month.pdf'
minax_title  = f'{varName} min'
maxax_title  = f'{varName} max'
kinax_title  = f'{varName} mean'
hisax_title  = f'{varName} mean - distribution'

# use 'plot_SanityCheck_3D' from script ../src/SanityCheck.py imported above
SanityCheck.plot_SanityCheck_3D(data=var, 
        # below is optional
        data_mask=var_mask, kind='mean', figname=figname,
        lowerP=2, upperP=98, interactive=False,
        # below is even more optional (**kwargs)
        fig_title=fig_title, minax_title=minax_title, maxax_title=maxax_title, 
        kinax_title=kinax_title, hisax_title=hisax_title, cmapName='Spectral_r')
