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
          f'{dataRootDir}/{datasetName}/{procType}/1983_01/{pressVarName}.nc',
          f'{dataRootDir}/{datasetName}/{procType}/1983_07/{pressVarName}.nc'
          ]
maskFile      = f'{dataRootDir}/{datasetName}/{procType}/1980_01/mask.nc'
# specific_storage
sstoragefile  = f'{dataRootDir}/{datasetName}/tmp_static/specific_storage.pfb'
# porosity
# should be possible to get this from indicator files and namelist...!
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

mask = ht.load_netcdf(maskFile, split=split, variable='mask')
print(f'mask.shape: {mask.shape}')
# strange mask... but below is correct for this particular case
# 99999 is land and other is lake / sea
# also reduce form 4D (1,z,y,x) to 3D (z,y,x)
mask = ht.where(mask[0]==99999,0,1)  # 0 land, 1 sea
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
    Poro=poro,
    Sstorage=sstorage,
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
satur       = []  # SATURation                    (t,z,y,x) [-]
gwb_mask    = []  # GroundWaterBody MASK          (t,z,y,x) [bool] 
wtd_z_index = []  # WaterTabelDepth Z INDEX       (t,y,x)   [-]
satu_sss    = []  # SATUrated SubSurfaceStorage   (t,y,x)   [m^3]
unsa_sss    = []  # UNSAturated SubSurfaceStorage (t,y,x)   [m^3]
for t in range(press.shape[0]):
    tmp_satur, krel = diag.VanGenuchten(press[t])
    tmp_gwb_mask, tmp_wtd_z_index = get_3Dgroundwaterbody_mask(tmp_satur)
    tmp_sss = diag.SubsurfaceStorage(press[t], tmp_satur)

    tmp_satu_sss = ht.where(tmp_gwb_mask==1, tmp_sss, 0)
    tmp_satu_sss = ht.sum(tmp_satu_sss, axis=0)

    tmp_unsa_sss = ht.where(tmp_gwb_mask==0, tmp_sss, 0)
    tmp_unsa_sss = ht.sum(tmp_unsa_sss, axis=0)

    sss.append(tmp_sss)
    satur.append(tmp_satur)
    gwb_mask.append(tmp_gwb_mask)
    wtd_z_index.append(tmp_wtd_z_index)
    satu_sss.append(tmp_satu_sss)
    unsa_sss.append(tmp_unsa_sss)
sss      = ht.stack(sss, axis=0)
print(f'sss.shape: {sss.shape}')
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

# SWITCH to Numpy as np.ma is more powerfull
print('convert to numpy...')
sss_numpy      = sss.numpy()
satur_numpy    = satur.numpy()
gwb_mask_numpy = gwb_mask.numpy()
wtd_z_index_numpy = wtd_z_index.numpy()
satu_sss_numpy = satu_sss.numpy()
unsa_sss_numpy = unsa_sss.numpy()
mask_numpy     = mask.numpy()

###############################################################################
#### Plot
###############################################################################
### Below is ugly but needed as SanityCheck is based on numpy instead of heat
sss_numpy = np.where(mask_numpy[0]==0,sss_numpy,np.nan)
gwb_mask_numpy = np.where(mask_numpy[0]==0,gwb_mask_numpy,np.nan)
wtd_z_index_numpy = np.where(mask_numpy[0]==0,wtd_z_index_numpy,np.nan)
satu_sss_numpy = np.where(mask_numpy[0]==0,satu_sss_numpy,np.nan)
unsa_sss_numpy = np.where(mask_numpy[0]==0,unsa_sss_numpy,np.nan)

print('plot')
# For mor detailed information about how plot_SanityCheck_3D() does work, see
# sloth/SanityCheck.py --> plot_SanityCheck_3D()
sloth.SanityCheck.plot_SanityCheck_3D(data=wtd_z_index_numpy,
    kind='mean', figname='./examples_SubSurfStor.pdf',
    fig_title='wtd_z_index [-] (t,y,x)', minax_title='min', maxax_title='max', 
    kinax_title='mean', cmapName='tab20')

###############################################################################
#### Create netCDF file and fill with basic attributes
###############################################################################
"""
# For mor detailed information about how createNetCDF() does work, see
# sloth/toolBox.py --> createNetCDF()
netCDFFileName = sloth.toolBox.createNetCDF('./sss_testFile.nc', domain='EU11', 
    nz=nz, timeCalendar=calendar, timeUnit=time_units,
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description='I want to test the sss calculation!',
    source='add source here',NBOUNDCUT=4)

with nc.Dataset(netCDFFileName, 'a') as nc_file:
    # Name of the variable: 'TestData'
    ncVar = nc_file.createVariable('sss', 'f8', ('time', 'lvl', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'sss'
    ncVar.long_name     = 'subsurfacestorage'
    ncVar.units         ='m^3'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = sss_numpy[...]

    ncTime      = nc_file.variables['time']
    ncTime[...] = np.arange(sss_numpy.shape[0])
"""

# For mor detailed information about how createNetCDF() does work, see
# sloth/toolBox.py --> createNetCDF()
netCDFFileName = sloth.toolBox.createNetCDF('./gwb_testFile.nc', domain='EU11', 
    nz=nz, timeCalendar=calendar, timeUnit=time_units,
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description='I want to inspect the GWB!',
    source='add source here',NBOUNDCUT=4)

with nc.Dataset(netCDFFileName, 'a') as nc_file:
    # Name of the variable: 'TestData'
    ncVar = nc_file.createVariable('gwb', 'f8', ('time', 'lvl', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'gwb'
    ncVar.long_name     = 'groundwaterbody'
    ncVar.units         ='--'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = gwb_mask_numpy[...]

    ncTime      = nc_file.variables['time']
    ncTime[...] = np.arange(gwb_mask_numpy.shape[0])

# For mor detailed information about how createNetCDF() does work, see
# sloth/toolBox.py --> createNetCDF()
netCDFFileName = sloth.toolBox.createNetCDF('./satur_testFile.nc', domain='EU11', 
    nz=nz, timeCalendar=calendar, timeUnit=time_units,
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description='I want to inspect the satur!',
    source='add source here',NBOUNDCUT=4)

with nc.Dataset(netCDFFileName, 'a') as nc_file:
    # Name of the variable: 'TestData'
    ncVar = nc_file.createVariable('satur', 'f8', ('time', 'lvl', 'rlat', 'rlon',),
                                    fill_value=-9999, zlib=True)
    ncVar.standard_name = 'satur'
    ncVar.long_name     = 'saturation'
    ncVar.units         ='--'
    ncVar.grid_mapping  = 'rotated_pole'
    ncVar[...]          = satur_numpy[...]

    ncTime      = nc_file.variables['time']
    ncTime[...] = np.arange(satur_numpy.shape[0])
