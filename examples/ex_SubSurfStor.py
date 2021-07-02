import numpy as np
import netCDF4 as nc
import argparse
import glob
import os
import sys
import heat as ht
import matplotlib.pyplot as plt
import datetime as dt

sloth_path='../'
sys.path.append(sloth_path)
import sloth

def get_3Dgroundwaterbody_mask(satur):
    ''' calculating a 3D mask of the groundwater-body

    Return:
      gwb_mask    = ht-ndarray
        A HeAT-array of the same shape as 'satur' input, holding 1 at pixels
        belonging to the groundwater-body, and 0 for other pixel
      wtd_z_index = ht-ndarray
        A HeAT-array of the same spatial shape as 'satur' (no z dimension),
        holding the index of the first (first from bottom to top/surface) 
        unsaturated layer for each pixel. If there are N layers with the model
        and the index has a value of N+1, than the entire column is saturated!
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
          f'{dataRootDir}/{datasetName}/{procType}/1983_06/{pressVarName}.nc',
          f'{dataRootDir}/{datasetName}/{procType}/1983_07/{pressVarName}.nc'
          ]
maskFile      = f'{dataRootDir}/{datasetName}/{procType}/1980_06/mask.nc'
sstoragefile  = f'{dataRootDir}/{datasetName}/tmp_static/specific_storage.pfb'
porofile      = f'{dataRootDir}/{datasetName}/tmp_static/porosity.pfb'
permFileZ     = f'{dataRootDir}/{datasetName}/tmp_static/perm_z.pfb'
permFileY     = f'{dataRootDir}/{datasetName}/tmp_static/perm_y.pfb'
permFileX     = f'{dataRootDir}/{datasetName}/tmp_static/perm_x.pfb'
pflname       = f'{dataRootDir}/{datasetName}/namelists/coup_oas.tcl'
indicatorfile = f'{dataRootDir}/{datasetName}/geo/parflow/ParFlow_SOIL_INDICATOR3_from_EUR03_x1592y1544z15_on_EUR11_x436y424z15.pfb'
slopeFileX    = f'{dataRootDir}/{datasetName}/geo/parflow/slopex.pfb'
slopeFileY    = f'{dataRootDir}/{datasetName}/geo/parflow/slopey.pfb'

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

# read time-info, which are needed to properly set time-dim in netCDF out
with nc.Dataset(pressFiles[0], 'r') as nc_file:
    nc_time    = nc_file.variables['time']
    calendar   = nc_time.calendar
    time_units = nc_time.units

# At the end 'mask' is requiered as 3D ht-ndarray (z,y,x)
# !A value 0 has to indicate sean pixel!
# !A value 1 has to indicate land pixel!
mask = ht.load_netcdf(maskFile, split=split, variable='mask')
print(f'mask.shape (raw): {mask.shape}')
# --> reduce form 4D (1,z,y,x) to 3D (z,y,x)
# --> ensure 'mask' is a binary field (sometime there are nan values inside)
mask = ht.where(mask[0]==0,0,1)  # 0 sea, 1 land
print(f'mask.shape: {mask.shape}')

sstorage = sloth.extern.diagIO.read_pfb(sstoragefile, split=split)
print(f'sstorage.shape: {sstorage.shape}')

poro = sloth.extern.diagIO.read_pfb(porofile, split=split)
print(f'poro.shape: {poro.shape}')

permz = sloth.extern.diagIO.read_pfb(permFileZ, split=split)
permy = sloth.extern.diagIO.read_pfb(permFileY, split=split)
permx = sloth.extern.diagIO.read_pfb(permFileX, split=split)
print(f'permz.shape: {permz.shape}')
print(f'permy.shape: {permy.shape}')
print(f'permx.shape: {permx.shape}')

slopex = sloth.extern.diagIO.read_pfb(slopeFileX, split=split)
slopey = sloth.extern.diagIO.read_pfb(slopeFileY, split=split)
print(f'slopex.shape: {slopex.shape}')
print(f'slopey.shape: {slopey.shape}')

# reading pfl-namelist and map indicator file against
print('reading pfl-namelist and map Indicator against')
# For mor detailed information about how mappIndicator() does work, see
# sloth/toolBox.py --> mappIndicator()
indicatorMap = sloth.toolBox.mappIndicator(ParFlowNamelist=pflname, IndicatorFile=indicatorfile)
alpha   = ht.array(indicatorMap['alpha'], is_split=0, comm=ht.MPI_WORLD)
nvg     = ht.array(indicatorMap['nvg'], is_split=0, comm=ht.MPI_WORLD)
sres    = ht.array(indicatorMap['sres'], is_split=0, comm=ht.MPI_WORLD)
dzmult  = ht.array(indicatorMap['dz_mult'], is_split=0, comm=ht.MPI_WORLD)
dz      = indicatorMap['dz']
dy      = indicatorMap['dy']
dx      = indicatorMap['dx']
nz      = indicatorMap['nz']
ny      = indicatorMap['ny']
nx      = indicatorMap['nx']

perm = ht.zeros((3,nz,ny,nx),split=split)
perm[0]=permz
perm[1]=permy
perm[2]=permx

mannings = ht.full(alpha.shape, 5.5e-5, split=split)
###############################################################################
### Initialize Diagnostics
###############################################################################
diag = sloth.extern.Diagnostics.Diagnostics(Mask=mask, Perm=perm,
    Poro=poro, Sstorage=sstorage,
    Ssat=1., Sres=sres, Nvg=nvg, Alpha=alpha,
    Mannings=mannings, Slopex=slopex, Slopey=slopey,
    Dx=dx, Dy=dy, Dz=dz,
    Dzmult=dzmult,
    Nx=nx, Ny=ny, Nz=nz,
    Terrainfollowing=True, Split=split)

###############################################################################
### Calculate stuff
###############################################################################
print('calculating SubSurfStor')
sss         = []  # SubSurfaceStorage             (t,z,y,x) [m^3]
ss          = []  # SurfaceStorage                (t,y,x)   [m^3]
satur       = []  # SATURation                    (t,z,y,x) [-]
gwb_mask    = []  # GroundWaterBody MASK          (t,z,y,x) [bool] 
wtd_z_index = []  # WaterTabelDepth Z INDEX       (t,y,x)   [-]
satu_sss    = []  # SATUrated SubSurfaceStorage   (t,y,x)   [m^3]
unsa_sss    = []  # UNSAturated SubSurfaceStorage (t,y,x)   [m^3]
# Diagnostics.py can handle 3D (z,y,x) data only. So looping over time steps:
for t in range(press.shape[0]):
    # calculate saturation
    tmp_satur, krel = diag.VanGenuchten(press[t])
    # find groundwater body based on saturation
    tmp_gwb_mask, tmp_wtd_z_index = get_3Dgroundwaterbody_mask(tmp_satur)
    # calculate subsurface storage
    tmp_sss = diag.SubsurfaceStorage(press[t], tmp_satur)

    # calculate saturated subsurface storage
    tmp_satu_sss = ht.where(tmp_gwb_mask==1, tmp_sss, 0)
    tmp_satu_sss = ht.sum(tmp_satu_sss, axis=0)

    # calculate unsaturated subsurface storage
    tmp_unsa_sss = ht.where(tmp_gwb_mask==0, tmp_sss, 0)
    tmp_unsa_sss = ht.sum(tmp_unsa_sss, axis=0)

    # reduce z-axis of subsurface storage 
    tmp_sss = ht.sum(tmp_sss, axis=0)

    # calculate surface storage
    tmp_ss       = diag.SurfaceStorage(press[t,-1])

    # add variables calculated for the current time-step to list
    sss.append(tmp_sss)
    ss.append(tmp_ss)
    satur.append(tmp_satur)
    gwb_mask.append(tmp_gwb_mask)
    wtd_z_index.append(tmp_wtd_z_index)
    satu_sss.append(tmp_satu_sss)
    unsa_sss.append(tmp_unsa_sss)

# stack list entries to ht-ndarray
sss      = ht.stack(sss, axis=0)
print(f'sss.shape: {sss.shape}')
ss       = ht.stack(ss, axis=0)
print(f'ss.shape: {ss.shape}')
satur    = ht.stack(satur, axis=0)
print(f'satur.shape: {satur.shape}')
gwb_mask = ht.stack(gwb_mask, axis=0)
print(f'gwb_mask.shape: {gwb_mask.shape}')
wtd_z_index = ht.stack(wtd_z_index, axis=0)
print(f'wtd_z_index.shape: {wtd_z_index.shape}')
satu_sss = ht.stack(satu_sss, axis=0)
print(f'satu_sss.shape: {satu_sss.shape}')
unsa_sss = ht.stack(unsa_sss, axis=0)
print(f'unsa_sss.shape: {unsa_sss.shape}')

###############################################################################
#### Convert from HeAT to Numpy array, 
#### as those are more easy to mask and dump to proper netCDF files
###############################################################################
print('convert to numpy...')
sss_numpy      = sss.numpy()
ss_numpy       = ss.numpy()
satur_numpy    = satur.numpy()
gwb_mask_numpy = gwb_mask.numpy()
wtd_z_index_numpy = wtd_z_index.numpy()
satu_sss_numpy = satu_sss.numpy()
unsa_sss_numpy = unsa_sss.numpy()
mask_numpy     = mask.numpy()

print('apply land/sea mask...')
# There is a need to broadcast the mask fromn (z,y,x) to (t,yz,y,x) before, 
# as np.ma seems to have some problmes doing this automatically
mask4D    = np.broadcast_to(mask_numpy==0, satur_numpy.shape) # (t,z,y,x)
mask3D    = np.broadcast_to(mask_numpy[-1]==0, ss_numpy.shape) # (t,y,x)
sss_numpy = np.ma.masked_where(mask3D, sss_numpy)
ss_numpy  = np.ma.masked_where(mask3D, ss_numpy)
satur_numpy    = np.ma.masked_where(mask4D, satur_numpy)
gwb_mask_numpy = np.ma.masked_where(mask4D, gwb_mask_numpy)
wtd_z_index_numpy = np.ma.masked_where(mask3D, wtd_z_index_numpy)
satu_sss_numpy = np.ma.masked_where(mask3D, satu_sss_numpy)
unsa_sss_numpy = np.ma.masked_where(mask3D, unsa_sss_numpy)

###############################################################################
#### Create netCDF file 
#### with some basic attributes and store variables to
###############################################################################
print('storing variables in netCDF...')
# There is a need to broadcast the mask fromn (z,y,x) to (t,yz,y,x) before, 
# For mor detailed information about how createNetCDF() does work, see
# sloth/toolBox.py --> createNetCDF()
netCDFFileName = sloth.toolBox.createNetCDF('./SubSurfaceStorage.nc', domain='EU11', 
    nz=nz, timeCalendar=calendar, timeUnit=time_units,
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description='I want to test the sss calculation!',
    source='add source here',NBOUNDCUT=4)

with nc.Dataset(netCDFFileName, 'a') as nc_file:
    ncVar = nc_file.createVariable('sss', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'sss'
    ncVar.long_name     = 'subsurface storage'
    ncVar.units         = 'm^3'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = sss_numpy[...]


    ncVar = nc_file.createVariable('ss', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'ss'
    ncVar.long_name     = 'surface storage'
    ncVar.units         = 'm^3'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = ss_numpy[...]


    ncVar = nc_file.createVariable('gwb', 'f8', ('time', 'lvl', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'gwb'
    ncVar.long_name     = 'groundwater body mask'
    ncVar.units         = '--'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = gwb_mask_numpy[...]


    ncVar = nc_file.createVariable('wtd_zidx', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'wtd_zidx'
    ncVar.long_name     = 'water table depth z index'
    ncVar.units         = '--'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = wtd_z_index_numpy[...]


    ncVar = nc_file.createVariable('satu_sss', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'satu_sss'
    ncVar.long_name     = 'saturated subsurface storage'
    ncVar.units         = 'm^3'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = satu_sss_numpy[...]


    ncVar = nc_file.createVariable('unsa_sss', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'unsa_sss'
    ncVar.long_name     = 'unsaturated subsurface storage'
    ncVar.units         = 'm^3'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = unsa_sss_numpy[...]


    ncVar = nc_file.createVariable('satur', 'f8', ('time', 'lvl', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'satur'
    ncVar.long_name     = 'saturation'
    ncVar.units         = '--'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = satur_numpy[...]


    ncTime      = nc_file.variables['time']
    ncTime[...] = np.arange(satur_numpy.shape[0])
