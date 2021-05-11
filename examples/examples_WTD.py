import numpy as np
import netCDF4 as nc
import argparse
import glob
import os
import sys
import heat as ht
import matplotlib.pyplot as plt
import datetime as dt

catchyNAME_path='../'
sys.path.append(catchyNAME_path)
import catchyNAME


def calc_WTD(press, dz, dz_mult):
    """ calculate the water table depth

    assuming/supporting following dimensions of press:
    dim=3D: (z, y, x)
    dim=4D: (t, z, y, x)
    exit if not 3D or 4D!
    """
    cellDepth = dz*dz_mult
    totalColumnDepth = ht.sum(cellDepth)
    dim = press.ndim
    if dim==3:
        wtd = totalColumnDepth - (press[0] + cellDepth[0] * 0.5)
    elif dim==4:
        wtd = totalColumnDepth - (press[:,0,...] + cellDepth[0] * 0.5)
    else:
        print(f'ERROR: dim of press input is not supported --> EXIT!')
        sys.exit()
    # Above could be also solved by below oneliner:
    # >> wtd = totalColumnDepth - (press[...,0,:,:] + cellDepth[0] * 0.5)
    # but maybe this is not as intuitive as the used if-clause?
    return wtd

###############################################################################
### Define some paths, filenames, etc
###############################################################################
# Need for HeAT
split = None

dataRootDir   = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName   = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType      = 'postpro'
pressVarName  = 'press'
pressFiles    = [
          f'{dataRootDir}/{datasetName}/{procType}/1980_01/{pressVarName}.nc',
          f'{dataRootDir}/{datasetName}/{procType}/1980_02/{pressVarName}.nc'
          ]
maskFile      = f'{dataRootDir}/{datasetName}/{procType}/1980_01/mask.nc'

# read in also ParFlow namelist and indicator file.
# I want to use / show a function defined in ANalysisTool (need to rethink the name)
# called 'mappIndicator()' which does read the ParFlow namelist and map
# with indicator file to generate vanGenuchten parameter fields.
# I do not need those fields for WTD, but it is also possible to extract other
# variables out of the ParFlow namelist, as e.g. dz and dzmult.
# This way one can avoid to type those values by hand what is very error-prone.
pflname       = f'{dataRootDir}/{datasetName}/namelists/coup_oas.tcl'
indicatorfile = f'{dataRootDir}/{datasetName}/geo/parflow/ParFlow_SOIL_INDICATOR3_from_EUR03_x1592y1544z15_on_EUR11_x436y424z15.pfb'

###############################################################################
### Read in all files and data
###############################################################################
press = []
for pressFile in pressFiles:
    print(f'ht.load_netcdf --> {pressVarName}')
    tmp_press = ht.load_netcdf(pressFile, split=split, variable=pressVarName)
    press.append(tmp_press)
press = ht.concatenate(press, axis=0)
print(f'press.shape: {press.shape}')

mask = ht.load_netcdf(maskFile, split=split, variable='mask')
print(f'mask.shape: {mask.shape}')
# strange mask... but below is correct for this particular case
# 99999 is land and other is lake / sea
# also reduce form 4D (1,z,y,x) to 3D (z,y,x)
mask = ht.where(mask[0]==99999,0,1)  # 0 land, 1 sea
print(f'mask.shape: {mask.shape}')

# reading pfl-namelist and map indicator file against
print('reading pfl-namelist and map Indicator against')
indicatorMap = catchyNAME.toolBox.mappIndicator(ParFlowNamelist=pflname, IndicatorFile=indicatorfile)
dz_mult = ht.array(indicatorMap['dz_mult'], is_split=0, comm=ht.MPI_WORLD)
dz      = indicatorMap['dz']

###############################################################################
### Calculate WTD
###############################################################################
print('calculating WTD')
wtd = calc_WTD(press, dz=dz, dz_mult=dz_mult)


###############################################################################
#### Plot
###############################################################################
### Below is ugly but needed as SanityCheck is based on numpy instead of heat
print('convert to numpy...')
wtd_numpy = wtd.numpy()
mask_numpy = mask.numpy()
wtd_numpy = np.where(mask_numpy[0]==0,wtd_numpy,np.nan)

print('plot')
catchyNAME.SanityCheck.plot_SanityCheck_3D(data=wtd_numpy,
    kind='mean', figname='./examples_WTD.pdf',
    fig_title='WTD [m] (t,y,x)', minax_title='min', maxax_title='max', 
    kinax_title='mean', cmapName='seismic')

###############################################################################
#### Create netCDF file and fill with basic attributes
###############################################################################
netCDFFileName = catchyNAME.toolBox.createNetCDF('./WTD_testFile.nc', domain='EU11', 
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description='I want tot est and more interactivly check WTD calculation!',
    source='add source here',NBOUNDCUT=4)

###############################################################################
#### Create the actual variable we want to store the data at.
###############################################################################
with nc.Dataset(netCDFFileName, 'a') as nc_file:
    # Name of the variable: 'TestData'
    ncVar = nc_file.createVariable('wtd', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999,
                                    zlib=True)
    ncVar.standard_name = 'wtd'
    ncVar.long_name = 'water table depth'
    ncVar.units ='m'
    ncVar.grid_mapping = 'rotated_pole'

    ncTime = nc_file.createVariable('time', 'i2', ('time',))
    ncTime.standard_name = 'time'
    ncTime.units = 'days since 1984-01-01 00:00:00'
    ncTime.calendar = '365_day'

    ncVar[...] = wtd_numpy[...]
    ncTime[...] = np.arange(wtd_numpy.shape[0])