import numpy as np
import sys
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import heat as ht

catchyNAME_path='../'
sys.path.append(catchyNAME_path)
import catchyNAME

ext_src='../extern/ana_parflow-diagnostics_pythonheat'
sys.path.append(ext_src)
import Diagnostics
import IO as htio


###############################################################################
### Define some paths, filenames, etc
###############################################################################
split = None
dataRootDir  = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName  = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType     = 'postpro'
pressVarName = 'press'
pressFile    = f'{dataRootDir}/{datasetName}/{procType}/1980_01/{pressVarName}.nc'
GRDCdataset  = 'GRDC'

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
press   = ht.load_netcdf(pressFile, split=split, variable=pressVarName)
SimLons = ht.load_netcdf(pressFile, split=split, variable='lon')
SimLats = ht.load_netcdf(pressFile, split=split, variable='lat')
print(f'press.shape: {press.shape}')

mask = ht.load_netcdf(maskFile, split=split, variable='mask')[0]
# 99999 is land and other is lake / sea
# for Diagnostics convert this to:
# 0 = pixel to mask (sea)
# 1 = pixel to keep (land)
mask = ht.where(mask==99999,1,0)
print(f'mask.shape: {mask.shape}')

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
indicatorMap = catchyNAME.toolBox.mappIndicator(ParFlowNamelist=pflname, IndicatorFile=indicatorfile)
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

Q  = ht.full((press[:,0,...].shape),0,split=split)
for t in range(press.shape[0]):
    print(f'handle t={t}')
    Toplayerpress = diag.TopLayerPressure(Press=press[t])
    flowx, flowy = diag.OverlandFlow(Toplayerpress=Toplayerpress)
    flowx = ht.where(mask[0]!=0,flowx,0)
    flowy = ht.where(mask[0]!=0,flowy,0)
    Q[t] = ht.abs(flowx)*dx + ht.abs(flowy)*dy
# convert from HeAT to numpy, as Diagnostics are based on HeAT
Q = Q.numpy()
SimMeanQ = np.mean(Q, axis=0)

###############################################################################
### Initialize GRDC dataset
###############################################################################
# provide list of individual GRDC files
GRDCfiles    = sorted(glob.glob(f'{dataRootDir}/{GRDCdataset}/*.mon'))
GRDC_example = catchyNAME.GRDCdataset.GRDCdataset(GRDCfiles=GRDCfiles)
# GRDC_example.filter_index(key='Country', value='DE')
GRDC_example.filter_index(key='GRDC-No', value='6122110')
# GRDC_example.filter_index_date(start='1979-01', end='1980-12', form='%Y-%m')
GRDC_example.read_files(start='1980-01-01', end='1980-01-31', form='%Y-%m-%d')

print(f'############################################################################')
print(GRDC_example.data)
print(f'Found {GRDC_example.data.shape[0]} stations with {GRDC_example.data.shape[1]} data-points for applied filter.')
print(f'Possible accessible data:')
print(f'GRDC_example.id.shape: {GRDC_example.id.shape} ({type(GRDC_example.id)})')
# print(f'GRDC_example.id: {GRDC_example.id}')
print(f'GRDC_example.data.shape: {GRDC_example.data.shape} ({type(GRDC_example.data)})')
# print(f'GRDC_example.data: {GRDC_example.data}')
print(f'GRDC_example.lats.shape: {GRDC_example.lats.shape} ({type(GRDC_example.lats)})')
# print(f'GRDC_example.lats: {GRDC_example.lats}')
print(f'GRDC_example.lons.shape: {GRDC_example.lons.shape} ({type(GRDC_example.lons)})')
# print(f'GRDC_example.lons: {GRDC_example.lons}')
print(f'GRDC_example.time.shape: {GRDC_example.time.shape} ({type(GRDC_example.time)})')
# print(f'GRDC_example.time: {GRDC_example.time}')
print(f'GRDC_example.meanArea.shape: {GRDC_example.meanArea.shape} ({type(GRDC_example.meanArea)})')
# print(f'GRDC_example.meanArea: {GRDC_example.meanArea}')
print(f'############################################################################')


###############################################################################
### Initialize Mapper-Object and map data on SimGrid
###############################################################################
# convert from HeAT to numpy, as Diagnostics are based on HeAT but mapper is 
# based on numpy
SimLons = SimLons.numpy()
SimLats = SimLats.numpy()
slopex  = slopex.numpy()
slopey  = slopey.numpy()
Mapper  = catchyNAME.mapper.mapper(SimLons=SimLons, SimLats=SimLats,
	                ObsLons=GRDC_example.lons, ObsLats=GRDC_example.lats,
	                ObsIDs=GRDC_example.id, 
	                SimMeanQ=SimMeanQ, ObsMeanQ=np.nanmean(GRDC_example.data, axis=1))

print(f'Map Catchment')
Mapper.ObsMeanArea = GRDC_example.meanArea
Mapper.MapBestCatchment(search_rad=3, slopex=slopex[0], slopey=slopey[0])
MapBestCatchment_X = Mapper.MapXIdx_fit
MapBestCatchment_Y = Mapper.MapYIdx_fit

print(f'MeanArea:')
print(f'Mapper.ObsMeanArea: {Mapper.ObsMeanArea}')
print(f'Mapper.SimMeanArea: {Mapper.SimMeanArea}')