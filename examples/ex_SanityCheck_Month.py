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
dataYear    = '1979_01'
varName     = 'T_S'
fileName    = f'{dataRootDir}/{datasetName}/{procType}/{dataYear}/{varName}_ts.nc'

###############################################################################
### read in data and save as ndarray
###############################################################################

with nc.Dataset(f'{fileName}', 'r') as nc_file:
    nc_var   = nc_file.variables[f'{varName}']
    var      = nc_var[...]

###############################################################################
### Start SanityCheck
###############################################################################
# define some title for plot, which can be passed via functionarguments
# see funciton definition for full potential
fig_title    = f'Sanity-Check for {varName} --- {dataYear}'
figname      = f'./ex_SanityCheck_Month.pdf'
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
        # below is even more optional (**kwargs)
        fig_title=fig_title, minax_title=minax_title, maxax_title=maxax_title, 
        kinax_title=kinax_title, hisax_title=hisax_title, cmapName='Spectral_r')
