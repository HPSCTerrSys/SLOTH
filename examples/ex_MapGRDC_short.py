import numpy as np
import sys
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import heat as ht

sloth_path='../'
sys.path.append(sloth_path)
import sloth.GRDCdataset
import sloth.mapper
import sloth.IO 
import sloth.PlotLib


###############################################################################
### Define some paths, filenames, etc
###############################################################################
dataRootDir   = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName   = 'ParFlow_MB3km'
procType      = 'postpro'
dischargeFile = f'{dataRootDir}/{datasetName}/{procType}/pf_discharge1997_2006.nc'
GRDCdataset   = 'GRDC'

CoordFile     = f'{dataRootDir}/{datasetName}/geo/geo_em.d03_MODIS3km1kmlakes_geo_cutdowntocosmo.nc' 
slopeFileX    = f'{dataRootDir}/{datasetName}/geo/ParFlow_MB3km_SLPX_x1592y1544.pfb'
slopeFileY    = f'{dataRootDir}/{datasetName}/geo/ParFlow_MB3km_SLPY_x1592y1544.pfb'

###############################################################################
### Read in data
###############################################################################
with nc.Dataset(dischargeFile,'r') as ncFile:
    Q = ncFile.variables['flow'][...]

print(f'Q.shape: {Q.shape}')
SimMeanQ = np.ma.mean(Q, axis=0)
print(f'SimMeanQ.shape: {SimMeanQ.shape}')

with nc.Dataset(CoordFile,'r') as ncFile:
    SimLons = ncFile.variables['XLONG_M'][0]
    SimLats = ncFile.variables['XLAT_M'][0]
print(f'SimLons.shape: {SimLons.shape}')
print(f'SimLats.shape: {SimLats.shape}')

slopex = sloth.IO.read_pfb(slopeFileX)[0]
slopey = sloth.IO.read_pfb(slopeFileY)[0]
print(f'slopex.shape: {slopex.shape}')
print(f'slopey.shape: {slopey.shape}')

###############################################################################
### Initialize GRDC dataset
###############################################################################
# provide list of individual GRDC files
GRDCfiles    = sorted(glob.glob(f'{dataRootDir}/{GRDCdataset}/*.mon'))
# For mor detailed information about how GRDCdataset() does work, see
# sloth/GRDCdataset.py --> GRDCdataset()
GRDC_example = sloth.GRDCdataset.GRDCdataset(GRDCfiles=GRDCfiles)
# GRDC_example.filter_index(key='Country', value='DE')
GRDC_example.filter_index(key='Catchment area', value=1000, operant='>')
# GRDC_example.filter_index(key='GRDC-No', value=[6233201], operant='in')
# GRDC_example.filter_index(key='GRDC-No', value=[6122110, 6119200, 6142520, 6335050, 6142660], operant='in')
# GRDC_example.dump_index(keys2dump=['GRDC-No', 'River', 'Country',
#                                    'Date start', 'Date end', 'Catchment area'])
# Filtering for time-period of interest is important, as read_files() is 'stupid' 
# and will crash if time-period is incomplete 
GRDC_example.filter_index_date(start='1997-01', end='2007-01', form='%Y-%m')
GRDC_example.read_files(start='1997-01-01', end='2006-12-31', form='%Y-%m-%d')

print(f'############################################################################')
print(GRDC_example.data)
print(f'Found {GRDC_example.data.shape[0]} stations with {GRDC_example.data.shape[1]} data-points for applied filter.')
print(f'Possible accessible data:')
print(f'GRDC_example.id.shape: {GRDC_example.id.shape} ({type(GRDC_example.id)})')
print(f'GRDC_example.data.shape: {GRDC_example.data.shape} ({type(GRDC_example.data)})')
print(f'GRDC_example.lats.shape: {GRDC_example.lats.shape} ({type(GRDC_example.lats)})')
print(f'GRDC_example.lons.shape: {GRDC_example.lons.shape} ({type(GRDC_example.lons)})')
print(f'GRDC_example.time.shape: {GRDC_example.time.shape} ({type(GRDC_example.time)})')
print(f'GRDC_example.meanArea.shape: {GRDC_example.meanArea.shape} ({type(GRDC_example.meanArea)})')
print(f'############################################################################')


###############################################################################
### Initialize Mapper-Object and map data on SimGrid
###############################################################################
# For mor detailed information about how mapper() does work, see
# sloth/mapper.py --> mapper()
Mapper  = sloth.mapper.mapper(SimLons=SimLons, SimLats=SimLats,
	                ObsLons=GRDC_example.lons, ObsLats=GRDC_example.lats,
	                ObsIDs=GRDC_example.id, 
	                SimMeanQ=SimMeanQ, ObsMeanQ=np.nanmean(GRDC_example.data, axis=1))

# For mor detailed information about how MapBestQ() does work, see
# sloth/mapper.py --> MapBestQ()
Mapper.MapBestQ(search_rad=2)
# For mor detailed information about how plot_MappedSubAreas() does work, see
# sloth/PlotLib.py --> plot_MappedSubAreas()
sloth.PlotLib.plot_MappedSubAreas(mapper=Mapper, fit_name='BestQ', search_rad=10)
print(f'Map Catchment')
Mapper.ObsMeanArea = GRDC_example.meanArea
# For mor detailed information about how MapBestCatchment() does work, see
# sloth/mapper.py --> MapBestCatchment()
Mapper.MapBestCatchment(search_rad=3, dy=3000, dx=3000, slopex=slopex, slopey=slopey)
# For mor detailed information about how plot_MappedSubAreas() does work, see
# sloth/PlotLib.py --> plot_MappedSubAreas()
sloth.PlotLib.plot_MappedSubAreas(mapper=Mapper, fit_name='BestCatchment', search_rad=10)
MapBestCatchment_X = Mapper.MapXIdx_fit
MapBestCatchment_Y = Mapper.MapYIdx_fit

print(f'MeanArea:')
print(f'Mapper.ObsMeanArea: {Mapper.ObsMeanArea}')
print(f'Mapper.SimMeanArea: {Mapper.SimMeanArea}')
