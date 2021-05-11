import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import sys
import os

sloth_path='../'
sys.path.append(sloth_path)
import sloth

dataRootDir = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType    = 'postpro'
dataDate    = '1979_01'
varName     = 'T_2M'
fileNameA = f'{dataRootDir}/{datasetName}/{procType}/{dataDate}/{varName}_ts.nc'
# procType    = 'simres'
# dataDate    = '1987_01'
# fileNameA = f'{dataRootDir}/{datasetName}/{procType}/{dataDate}/cosmo/lffd1987020100.nc'

fileNameB = '../data/example_WriteNetCDF/example_netCDF_EU11'

fileNameB = sloth.toolBox.createNetCDF(fileNameB, domain='EU11', fillValue=-123456, NBOUNDCUT=4)

with nc.Dataset(fileNameA, 'r') as nc_fileA:
	rlonA = nc_fileA.variables['rlon'][...]
	rlatA = nc_fileA.variables['rlat'][...]

with nc.Dataset(fileNameB, 'r') as nc_fileB:
	rlonB = nc_fileB.variables['rlon'][...]
	rlatB = nc_fileB.variables['rlat'][...]

print(f'np.all(rlonA==rlonB): {np.all(rlonA==rlonB)}')
print(f'np.all(rlatA==rlatB): {np.all(rlatA==rlatB)}')