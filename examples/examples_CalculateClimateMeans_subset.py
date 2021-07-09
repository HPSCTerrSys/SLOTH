""" Example script to calculate climatology means.

Calculating means with python is not a difficult task and even the choice  
of the interval we are handling (e.g. daily or monthly means) does not really 
matter. At the end we do simply have to apply e.g. np.mean() to our data and
we are done!
However, extracting the correct slices or intervals out of a given data set to 
apply the mean-function to, could be very nasty. 
The correct slicing depends on the actual data we are facing as e.g. the 
starting time step and the time resolution (hours, daily, monthly, or even 
yearly data?). Individual datasets need individual treatment.

This example make use of the function 'get_intervalSlice()' which aims to 
generalize the progress to get the correct slices or intervals.
Mandatory for this purpose is, the handled data does contain the correct 
time information, as e.g. the correct time-variables with netCDF files.
If this is the case the function does loop over all your dates and returns
a list of slices to extract the correct intervals out of your data. If you 
data-set does consist out of individual files or is represented by one big 
file does not matter.
"""
import numpy as np
import datetime as dt
import netCDF4 as nc
import matplotlib as mpl
import sys
import os
import glob
import cftime

sloth_path='../'
sys.path.append(sloth_path)
import sloth


###############################################################################
### Define some paths, filenames, options, etc
###############################################################################
dataRootDir = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType    = 'postpro'
# CLM style
varName = 'TSA'
fileName = f'{varName}.nc'
# COSMO style
#varName = 'T_2M'
#fileName = f'{varName}_ts.nc'
files = sorted(glob.glob(f'{dataRootDir}/{datasetName}/{procType}/*/{fileName}'))

# Climatology calculations are always based on comparisons of individual 
# intervals between different years. 
# If the interval we are interested in is 'month', we do compare the same 
# month between different years.
# If the interval we are interested in is 'day', we do compare the same 
# day between different years.
# etc.
# In any case we do have to know the Number of Intervals (NoI) a year does 
# contain. For daily based calculation this usually is NoI=365 (365 days a 
# year), for monthly based calculations this usually is NoI=12 (12 month a 
# year).
# However we could also handle summer months (JJA) only and do a monthly based
# calculation, in which case the Number of Intervals is NoI=3, as we do 
# investigate 3 months only.
meanInterval = 'day'
NoI = 92
# Define which dates to handle
validYears  = np.arange(1979,1982)
validMonths = [6,7,8]


###############################################################################
#### Read in data and calculate interval mean and calculate center time-step
###############################################################################
# Appending data to a list and convert this list into a ndarray once at the end 
# is much faster than appending to ndarray (np.append()) at each iteration 
# step. I guess this is because a list contains pointers and appending a 
# pointer is quiet fast, while np.append() open and restructure the
# entire array which is getting bigger and bigger while reading in more
# and more data / files.
tmp_IntervalMean = []
tmp_IntervalTime = []
tmp_timeCalendar = []
tmp_timeUnits    = []
# Loop over all files of our data-set.
for file in files:
    # Open the file and read in data and time
    with nc.Dataset(file, 'r') as nc_file:
        data  = nc_file.variables[varName]
        print(f'data.shape: {data.shape}')
        # add some chekc here if ndims is correct
        # Extract units, time, etc for later usage.
        dataUnits    = data.units
        nc_time      = nc_file.variables['time']
        calendar     = nc_time.calendar
        dates        = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
        timeValues   = nc_time[:]
        timeCalendar = nc_time.calendar
        timeUnits    = nc_time.units
        # Store timeUnits and timeCalendar to proper store output netCDF
        tmp_timeCalendar.append(timeCalendar)
        tmp_timeUnits.append(timeUnits)

        # Calculate the slices for the current file based on the choose meanInterval
        # For mor detailed information about how get_intervalSlice() does work, see
        # sloth/toolBox.py --> get_intervalSlice()
        dailySlices = sloth.toolBox.get_intervalSlice(dates=dates, sliceInterval=meanInterval)
        # Loop over all slices, mask 'missing' values with np.nan, and calculate 
        # related mean. The averaged data gets appended for each file and slice.
        for Slice in dailySlices:
            # Skip dates not in valid range
            tmp_dates    = dates[Slice]
            yearsValid   = [ date.year in validYears for date in tmp_dates]
            if not all(yearsValid):
                continue
            monthsValid   = [ date.month in validMonths for date in tmp_dates]
            if not all(monthsValid):
                continue

            tmp_time     = timeValues[Slice]
            tmp_time     = tmp_time.filled(fill_value=np.nan)
            tmp_var      = data[Slice]

            tmp_timeMean  = tmp_time.mean(axis=0, keepdims=True)
            tmp_timeMean  = nc.num2date(tmp_timeMean, units=timeUnits,calendar=timeCalendar)
            tmp_monthMean = tmp_var.mean(axis=0, keepdims=True, dtype=float)

            tmp_IntervalMean.append(tmp_monthMean)
            tmp_IntervalTime.append(tmp_timeMean)
            print(f'len(tmp_IntervalMean): {len(tmp_IntervalMean)}')


    print('##################')
    print('##################')
intervalMean = np.concatenate(tmp_IntervalMean, axis=0)
print(f'intervalMean.shape: {intervalMean.shape}')
intervalTime = np.concatenate(tmp_IntervalTime, axis=0)

# First create an empty array of same shape as intervalMean but with 
# t-dim (0-axis) = NoI
climaDim = [NoI]    
climaDim = climaDim + [dim for dim in intervalMean[0].shape]
clima = np.empty(climaDim)
for curr_interval in range(NoI):
    # The climat mean is simple the mean of all entries with NoI in distance.
    # NoI=12 --> for all entries with a distance of 12 --> month
    # NoI=365 --> for all entries with a distance of 365 --> days
    clima[curr_interval] = intervalMean[curr_interval::NoI].mean(axis=0, dtype=float)

###############################################################################
#### Prepare output netCDF file to appand different variables
###############################################################################
# For mor detailed information about how createNetCDF() does work, see
# sloth/toolBox.py --> createNetCDF()
saveDir = f'../data/example_ClimateMeans/'
# save everything for later use         
if not os.path.exists(f'{saveDir}'):
            os.makedirs(f'{saveDir}')

#### Store the ClimatMean values
###############################################################################
YearStart=f'{intervalTime[0].year}'
YearEnd=f'{intervalTime[-1].year}'
descriptionStr = [
        f'This file does contain climate mean values on a {meanInterval} base ',
        f'calculated for the years {YearStart} to {YearEnd}.',
        ]
descriptionStr = ''.join(descriptionStr)
saveFile=f'{saveDir}/ClimateMeans_{YearStart}-{YearEnd}_NoI-{NoI}_MeanInterval-{meanInterval}.nc'
netCDFFileName = sloth.toolBox.createNetCDF(saveFile, domain='EU11_TSMP',
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description=descriptionStr),
    source='---', NBOUNDCUT=4)

with nc.Dataset(netCDFFileName, 'a') as nc_file:
    ncVar = nc_file.createVariable(f'ClimatMean_{varName}', 'f8', ('time', 'rlat', 'rlon',),
                    fill_value=9999, zlib=True)
    ncVar.standard_name = f'climat_{varName}'
    ncVar.long_name = f'climat mean of {varName}'
    ncVar.units =f'{dataUnits}'
    ncVar.description = f'{meanInterval} based climat mean values'
    ncVar[...] = clima[...]
    
    ncTime = nc_file.createVariable('time', 'i4', ('time',),
                    fill_value=9999, zlib=True)
    ncTime.description = f'{meanInterval} number in year'
    # NOTE:
    # np.arange(clima.shape[0]) is not correct if calculating subset of data
    # only
    ncTime[...] = np.arange(clima.shape[0])

#### Store the Mean values
###############################################################################
descriptionStr = [
        f'This file does contain {meanInterval} mean values ',
        f'calculated for the years {YearStart} to {YearEnd}.',
        ]
descriptionStr = ''.join(descriptionStr)
# take the first read in timeUnits and timeCalendar for output netCDF
timeUnitsOut    = tmp_timeUnits[0]
timeCalendarOut = tmp_timeCalendar[0]
saveFile=f'{saveDir}/Means_{YearStart}-{YearEnd}_MeanInterval-{meanInterval}.nc'
netCDFFileName = sloth.toolBox.createNetCDF(saveFile, domain='EU11_TSMP',
    timeCalendar=timeCalendarOut, timeUnit=timeUnitsOut,
    author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
    institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
    description=descriptionStr),
    source='---', NBOUNDCUT=4)

with nc.Dataset(netCDFFileName, 'a') as nc_file:
    ncVar = nc_file.createVariable(f'mean_{varName}', 'f8', ('time', 'rlat', 'rlon',),
                    fill_value=9999, zlib=True)
    ncVar.standard_name = f'mean {varName}'
    ncVar.long_name = f'mean of {varName}'
    ncVar.units =f'{dataUnits}'
    ncVar.description = f'{meanInterval} based mean values'
    ncVar[...] = intervalMean[...]

    ncTime = nc_file.variables['time']
    ncTime[...] = nc.date2num(intervalTime, units=timeUnitsOut,calendar=timeCalendarOut)[...]

# plot if meanInterval='month'
if meanInterval == 'month':
    kwargs = {
            'title': 'Test climate plot',
            #'title': '\n'.join(tmp_titlesubstr),
            'infostr': True,
            #'var_vmax': 1,
            #'var_vmin': 0,
            'var_cmap': mpl.cm.get_cmap('jet'),
            'saveFile': f'./examples_CalculateClimateMeans.pdf',
            #'dpi': 100,
            'figsize': (10, 4),
            }
    sloth.PlotLib.plot_ClimateYearMonth(clima, **kwargs)
