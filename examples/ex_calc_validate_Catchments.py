import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime

import sys
import glob

import sloth
import sloth.ParFlow_IO as pio

def fileInRange(f, startDate, endDate, dateForm='%Y_%m'):
    f = f.split('/')[-2]
    tmp_date = datetime.datetime.strptime(f, dateForm)
    if tmp_date < startDate or tmp_date > endDate:
        return False
    else:
        return True

def read_data(fileNames, varName, convCoef):
    # default 'not rotated' --> overwriten below
    pole_latitude  = 90
    pole_longitude = 180
    rlon1D = None
    rlat1D = None
    data = []
    dates = []
    calendar = []
    for fileName in fileNames:
        print(f'fileName: {fileName}')
        with nc.Dataset(fileName, 'r') as nc_file:
            data.append(nc_file.variables[varName][...]*convCoef)
            pole_longitude = float(nc_file.variables['rotated_pole'].grid_north_pole_longitude)
            pole_latitude  = float(nc_file.variables['rotated_pole'].grid_north_pole_latitude)
            rlat1D = nc_file.variables['rlat'][...]
            rlon1D = nc_file.variables['rlon'][...]

            nc_time   = nc_file.variables['time']
            tmp_dates     = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
            calendar.append(nc_time.calendar)
            dates.append(tmp_dates)


    if len(data) == 1:
        data = data[0]
        dates    = dates[0]
        calendar = calendar[0]
    else:
        data = np.concatenate(data, axis=0)
        dates = np.concatenate(dates, axis=0)
        # check if all calendar strings are equal, which is IMPORTANT!
        if calendar[:-1] == calendar[1:]:
            calendar = calendar[0]
        else:
            print(f'ERROR: red in different calendar strings --> EXIT')
            print(f'in function "read_dates()"')
            sys.exit(1)

    rlon2D, rlat2D = np.meshgrid(rlon1D, rlat1D)

    return data, rlon2D, rlat2D, pole_longitude, pole_latitude, dates, calendar

def calc_NSE(Qm, Qo, norm=False):
    Qom = Qo.mean()
    NSE = 1 - (np.abs(Qm-Qo)).sum() / (np.abs(Qo-Qom)).sum()
    #NSE = 1 - ((Qm-Qo)**2).sum() / ((Qo-Qom)**2).sum()
    if norm:
        NSE = 1. / (2.-NSE)

    return NSE

###############################################################################
### Define some paths, filenames, etc
###############################################################################
startDate        = datetime.datetime(1981, 1, 1) 
endDate          = datetime.datetime(2010, 12, 31)

dataRootDir      = '/p/scratch/cesmtst/tsmpforecast/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow/'
grdcRootDir      = f'{dataRootDir}/analysis/GRDC/data'
GRDCstationsFile = f'{grdcRootDir}/grdc_stations/GRDC_Stations.csv'

srRootDir     = f'{dataRootDir}/ctrl/postpro/tmp_postpro_out'
tmpFileNames  = sorted(glob.glob(f'{srRootDir}/????_??/sr.nc'))
srFileNames   = [f for f in tmpFileNames if fileInRange(f, startDate, endDate)]

CoordFile     = f'{dataRootDir}/geo/DE-0055_LAND-LAKE-SEA-MASK_regridded_DE-011.nc'
slopeFileX    = f'{dataRootDir}/geo/DE-0055_XSLOPE_TPS_MERIT_sea_streams_corr_611.pfb'
slopeFileY    = f'{dataRootDir}/geo/DE-0055_YSLOPE_TPS_MERIT_sea_streams_corr_611.pfb'

###############################################################################
### Read in data
###############################################################################
with nc.Dataset(CoordFile,'r') as ncFile:
    SimLons = ncFile.variables['lon'][...]
    SimLats = ncFile.variables['lat'][...]
print(f'SimLons.shape: {SimLons.shape}')
print(f'SimLats.shape: {SimLats.shape}')

slopex = pio.read_pfb(slopeFileX)[0]
slopey = pio.read_pfb(slopeFileY)[0]
print(f'slopex.shape: {slopex.shape}')
print(f'slopey.shape: {slopey.shape}')

sr, rlon2D, rlat2D, pole_longitude, pole_latitude, dates, calendar = read_data(fileNames=srFileNames, varName='sr', convCoef=1.)

###############################################################################
### Initialize GRDC dataset
###############################################################################
GRDCfiles    = sorted(glob.glob(f'{grdcRootDir}/2021-08-20_13-35/*_Q_Month.txt'))
# For mor detailed information about how GRDCdataset() does work, see
# sloth/GRDCdataset.py --> GRDCdataset()
GRDC_example = sloth.GRDCdataset.GRDCdataset(GRDCfiles=GRDCfiles, GRDCstationsFile=GRDCstationsFile)
# manually inspect which station I want --> dump keys and get IDs
#GRDC_example.filter_index(key='Country', value='DE', operant='==')
#GRDC_example.filter_index(key='Catchment area', value=10, operant='>')
#GRDC_example.filter_index(key='River', value=['RHINE RIVER', 'ELBE RIVER', 'DANUBE RIVER', 'WESER'], operant='in')
#GRDC_example.filter_index(key='Date end', value=2001, operant='>')
#GRDC_example.filter_index(key='GRDC-No', value=[6335020, 6340150, 6342900, 6337200], operant='in')
#GRDC_example.dump_index(keys2dump=['GRDC-No', 'River', 'Station', 'Country',
#                                    'Date start', 'Date end', 'Catchment area'])
#sys.exit()
GRDC_example.filter_index_date(start='1980-01', end='2011-01', form='%Y-%m')
GRDC_example.read_files(start=startDate.strftime('%Y-%m-%d'), 
        end=endDate.strftime('%Y-%m-%d'), 
        form='%Y-%m-%d')
#GRDC_example.filter_index_date(start='1961-01', end='2009-01', form='%Y-%m')
#GRDC_example.read_files(start='1961-01-01', end='2008-12-31', form='%Y-%m-%d')

print(f'############################################################################')
#print(GRDC_example.data)
print(f'Found {GRDC_example.data.shape[0]} stations with {GRDC_example.data.shape[1]} data-points for applied filter.')
print(f'Possible accessible data:')
print(f'GRDC_example.id.shape: {GRDC_example.id.shape} ({type(GRDC_example.id)})')
print(f'GRDC_example.data.shape: {GRDC_example.data.shape} ({type(GRDC_example.data)})')
print(f'GRDC_example.lats.shape: {GRDC_example.lats.shape} ({type(GRDC_example.lats)})')
print(f'GRDC_example.lons.shape: {GRDC_example.lons.shape} ({type(GRDC_example.lons)})')
print(f'GRDC_example.time.shape: {GRDC_example.time.shape} ({type(GRDC_example.time)})')
#print(f'GRDC_example.meanArea.shape: {GRDC_example.meanArea.shape} ({type(GRDC_example.meanArea)})')
print(f'############################################################################')
#sys.exit()

print(GRDC_example.GRDCindexObj)
##############################################################################
### Initialize Mapper-Object and map data on SimGrid
###############################################################################
# SimMeanQ is needed, so set dummy
dummy_SimMeanQ = np.zeros_like(SimLons)
Mapper  = sloth.mapper.mapper(SimLons=SimLons, SimLats=SimLats,
                        ObsLons=GRDC_example.lons, ObsLats=GRDC_example.lats,
                        ObsIDs=GRDC_example.id,
                        SimMeanQ=dummy_SimMeanQ, ObsMeanQ=np.nanmean(GRDC_example.data, axis=1))

print(f'Map Catchment')
#Mapper.ObsMeanArea = np.zeros_like(GRDC_example.time) # dummy!
Mapper.ObsMeanArea = GRDC_example.meanArea
# For mor detailed information about how MapBestCatchment() does work, see
# sloth/mapper.py --> MapBestCatchment()
Mapper.MapRaw()
#MapBestCatchment_X = Mapper.MapXIdx_raw
#MapBestCatchment_Y = Mapper.MapYIdx_raw
Mapper.MapBestCatchment(search_rad=2, dy=611, dx=611, slopex=slopex, slopey=slopey)
MapBestCatchment_X = Mapper.MapXIdx_fit
MapBestCatchment_Y = Mapper.MapYIdx_fit
print(f'MapBestCatchment_X: {MapBestCatchment_X}')
print(f'MapBestCatchment_Y: {MapBestCatchment_Y}')
print(f'MeanArea:')
print(f'Mapper.ObsMeanArea: {Mapper.ObsMeanArea}')
print(f'Mapper.SimMeanArea: {Mapper.SimMeanArea}')

#MapBestCatchment_X = [663, 996, 1285, 1427]
#MapBestCatchment_Y = [1265, 1440, 1415, 601]
stationCoords = zip(MapBestCatchment_Y, MapBestCatchment_X)
#stationCoords = [(1265, 663),(1440,996),(1415, 1285),(601, 1427)]

PFL_data  = []
GRDC_data = []
goodStations = 0
badStations  = 0
for idx, stationCoord in enumerate(stationCoords):
    y,x = stationCoord
    print(f'check Obs and Sim catchment')
    MeanAreaRatio = Mapper.ObsMeanArea[idx] / Mapper.SimMeanArea[idx]
    if not (0.9 < MeanAreaRatio < 1.1):
        print('    Catchment area does not fit in Obs and Sim')
        print(f'    Mapper.ObsMeanArea[idx]: {Mapper.ObsMeanArea[idx]}')
        print(f'    Mapper.SimMeanArea[idx]: {Mapper.SimMeanArea[idx]}')
        badStations += 1
        continue
    goodStations += 1
    print(f'start calc Catchment from -- x:{x} y:{y}')
    stationID   = GRDC_example.id[idx]
    tmpHeader   = GRDC_example.GRDCindexObj[0]
    tmpEntries  = GRDC_example.GRDCindexObj[1][idx]
    tmpMetadata = {tmpHeader[i]:tmpEntries[i] for i, _ in enumerate(tmpHeader) }
    
    catchmentMask = sloth.toolBox.calc_catchment(slopex=slopex, slopey=slopey, x=x, y=y)
    # a mask in sense of np.ma is: 1=to mask, 0=not to mask
    catchmentMask = np.where(catchmentMask==0, 1, 0)
    print('end calc Catchment')

    # get surface runoff (in [m/month]) for catchment
    catchmentMask_b = np.broadcast_to(catchmentMask, sr.shape )
    tmp_PFL_data        = np.ma.masked_where(catchmentMask_b, sr)
    tmp_PFL_data        = tmp_PFL_data.sum(axis=(1,2))
    tmp_PFL_data        = tmp_PFL_data.mean()
    # convert [m/month] --> [m^3/month]
    tmp_PFL_data *= 611*611
    # convert [m^3/month] --> [m^3/s]
    tmp_PFL_data /= 60*60*24*30.5 # 30.5 Tage average per month
    PFL_data.append(tmp_PFL_data)

    tmp_GRDC_data = GRDC_example.data[idx]
    tmp_GRDC_data = tmp_GRDC_data.mean()
    GRDC_data.append(tmp_GRDC_data)

    print(f'tmp_PFL_data: {tmp_PFL_data} -- tmp_GRDC_data: {tmp_GRDC_data}')

print(f'badStations: {badStations}')
print(f'goodStations: {goodStations}')
print(f'DEBUG: start writign to netCDF')

with open('./PFL_scatterValidation_1981-2010.npy', 'wb') as f:
    np.save(f, np.array(PFL_data))
with open('./GRDC_scatterValidation_1981-2010.npy', 'wb') as f:
    np.save(f, np.array(GRDC_data))

