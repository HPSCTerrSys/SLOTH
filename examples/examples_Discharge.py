import numpy as np
import netCDF4 as nc
import argparse
import glob
import os
import sys
import heat as ht
import matplotlib.pyplot as plt
import datetime as dt

src_path='../src/'
sys.path.append(src_path)
# ANalysisTool does contain the function 'get_intervalSlice()'
import ANalysisTool as ANT
import SanityCheck 
src_path='../extern/ana_parflow-diagnostics_pythonheat'
sys.path.append(src_path)
# ANalysisTool does contain the function 'get_intervalSlice()'
import Diagnostics
import IO as htio


###############################################################################
### Define some paths, filenames, etc
###############################################################################
# Need for HeAT
split = None

dataRootDir   = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName   = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType      = 'postpro'
pressVarName  = 'press'
pressFiles    = sorted(glob.glob(f'{dataRootDir}/{datasetName}/{procType}/1980_01/{pressVarName}.nc'))
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
### Read in data
###############################################################################
press = ht.load_netcdf(pressFiles[0], split=split, variable=pressVarName)
# press = []
# for pressFile in pressFiles:
#     print(f'ht.load_netcdf')
#     tmp_press = ht.load_netcdf(pressFile, split=split, variable=pressVarName)
#     press.append(tmp_press)
# press = ht.concatenate(press, axis=0)
print(f'press.shape: {press.shape}')

mask = ht.load_netcdf(maskFile, split=split, variable='mask')[0]
# strange mask, but should be correct for sample data-set
# 99999 is land and other is lake / sea
# needed by Diagnostics:
# 0 = pixel to mask (sea)
# 1 = pixel to keep (land)
mask = ht.where(mask==99999,1,0)
print(f'mask.shape: {mask.shape}')
SanityCheck.plot_SanityCheck_3D(data=mask.numpy(),
    kind='mean', figname='./Mask.pdf',
    fig_title='Q [L^3/T] (t,y,x)', minax_title='min', maxax_title='max', 
    kinax_title='mean', cmapName='Blues')
# sys.exit()

sstorage = htio.read_pfb(sstoragefile, split=split)
print(f'sstorage.shape: {sstorage.shape}')

poro = htio.read_pfb(porofile, split=split)
print(f'poro.shape: {poro.shape}')

permz = htio.read_pfb(permFileZ, split=split)
permy = htio.read_pfb(permFileY, split=split)
permx = htio.read_pfb(permFileX, split=split)
print(f'permz.shape: {permz.shape}')
print(f'permy.shape: {permy.shape}')
print(f'permx.shape: {permx.shape}')

slopex = htio.read_pfb(slopeFileX, split=split)
slopey = htio.read_pfb(slopeFileY, split=split)
print(f'slopex.shape: {slopex.shape}')
print(f'slopey.shape: {slopey.shape}')


###############################################################################
### Initialize Diagnostics and calculate discharge (Q)
###############################################################################
indicatorMap = ANT.toolBox.mappIndicator(ParFlowNamelist=pflname, IndicatorFile=indicatorfile)
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

diag = Diagnostics.Diagnostics(Mask=mask, Perm=perm,
    Poro=poro,
    Sstorage=sstorage,
    Ssat=1., Sres=sres, Nvg=nvg, Alpha=alpha,
    Mannings=mannings, Slopex=slopex, Slopey=slopey,
    Dx=dx, Dy=dy, Dz=dz, 
    Dzmult=dzmult, 
    Nx=nx, Ny=ny, Nz=nz,
    Terrainfollowing=True, Split=split)

net_lateral_overlandflow = ht.full(press[:,0,...].shape,0,split=split)
Toplayerpress = ht.full(press[:,0,...].shape,0,split=split)
discharge     = ht.full(press[:,0,...].shape,0,split=split)
for t in range(press.shape[0]):
    print(f'handle t={t}')
    Toplayerpress[t] = diag.TopLayerPressure(Press=press[t])
    flowx, flowy = diag.OverlandFlow(Toplayerpress=Toplayerpress[t])
    flowx = ht.where(mask[0]!=0,flowx,0)
    flowy = ht.where(mask[0]!=0,flowy,0)
    net_lateral_overlandflow[t] = diag._NetLateralOverlandFlow(overland_flow_x=flowx, overland_flow_y=flowy)
    discharge[t] = ht.abs(flowx)*dx + ht.abs(flowy)*dy

###############################################################################
#### Plot
###############################################################################
### Below is ugly but needed as SanityCheck is based on numpy instead of heat
print('convert to numpy...')
discharge_numpy     = discharge.numpy()
net_lateral_overlandflow_numpy     = net_lateral_overlandflow.numpy()
Toplayerpress_numpy = Toplayerpress.numpy()
mask_numpy = mask.numpy()
discharge_numpy = np.where(mask_numpy[0]!=0,discharge_numpy,np.nan)
Toplayerpress_numpy = np.where(mask_numpy[0]!=0,Toplayerpress_numpy,np.nan)
net_lateral_overlandflow_numpy = np.where(mask_numpy[0]!=0,net_lateral_overlandflow_numpy,np.nan)

print('plot')
SanityCheck.plot_SanityCheck_3D(data=discharge_numpy,
    kind='mean', figname='./examples_Discharge.pdf',
    fig_title='Q [L^3/T] (t,y,x)', minax_title='min', maxax_title='max', 
    kinax_title='mean', cmapName='Blues')
SanityCheck.plot_SanityCheck_3D(data=net_lateral_overlandflow_numpy,
    kind='mean', figname='./NetLatOverlandFlow.pdf',
    fig_title='Q [L^3/T] (t,y,x)', minax_title='min', maxax_title='max', 
    kinax_title='mean', cmapName='Blues')
SanityCheck.plot_SanityCheck_3D(data=Toplayerpress_numpy,
    kind='mean', figname='./TopLayerPRess.pdf',
    fig_title='pressure [m] (t,y,x)', minax_title='min', maxax_title='max', 
    kinax_title='mean', cmapName='Blues')

###############################################################################
#### Create netCDF file and fill with basic attributes
###############################################################################
netCDFFileName = ANT.toolBox.createNetCDF('./Q_testFile.nc', domain='EU11', 
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description='Discharge - to - store - with - testdataset!',
    source='add source here',NBOUNDCUT=4)

###############################################################################
#### Create the actual variable we want to store the data at.
###############################################################################
with nc.Dataset(netCDFFileName, 'a') as nc_file:
    # Name of the variable: 'TestData'
    ncVar = nc_file.createVariable('Q', 'f8', ('time', 'rlat', 'rlon',),
                                    fill_value=-9999,
                                    zlib=True)
    ncVar.standard_name = 'Q'
    ncVar.long_name = 'surface discharge'
    ncVar.units ='L^3/T'
    ncVar.grid_mapping = 'rotated_pole'

    ncTime = nc_file.createVariable('time', 'i2', ('time',))
    ncTime.standard_name = 'time'
    ncTime.units = 'days since 1980-01-01 00:00:00'
    ncTime.calendar = '365_day'

    ncVar[...] = discharge_numpy[...]
    ncTime[...] = np.arange(discharge_numpy.shape[0])
