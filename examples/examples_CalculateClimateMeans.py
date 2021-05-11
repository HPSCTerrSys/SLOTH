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
# In any case we do have to know how the Number of Intervals (NoI) a year does 
# contain. For daily based calculation this usually is NoI=365 (365 days a 
# year), for monthly based calculations this usually is NoI=12 (12 month a 
# year).
# However we could also handle summer months (JJA) only and do a monthly based
# calculation, in which case the Number of Intervals is NoI=3, as we do 
# investigate 3 months only.
meanInterval = 'day'
NoI = 365


###############################################################################
#### Read in data and calculate interval mean and calculate center time-step
###############################################################################
# Create a tmp lists to append data at. Append data to a list and convert this 
# list into a ndarray once at the end is much faster than appending to ndarray 
# (np.append()) at each iteration step.
# I guess this is because list contains pointers and appending a 
# pointer is quiet fast, while np.append() open and restructure the
# entire array which is getting bigger and bigger while reading in more
# and more data / files.
tmp_IntervalMean = []
tmp_IntervalTime = []
# Loop over all files of our data-set.
for file in files:
    # Open the file and read in data and time
    with nc.Dataset(file, 'r') as nc_file:
        data  = nc_file.variables[varName]
        print(f'data.shape: {data.shape}')
        # add some chekc here if ndims is correct
        # Extract data, units, time, etc for later usage.
        data         = data[...]
        nc_time      = nc_file.variables['time']
        calendar     = nc_time.calendar
        dates        = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
        timeValues   = nc_time[:]
        timeCalendar = nc_time.calendar
        timeUnits    = nc_time.units

    # Calculate the slices for the current file based on the choose meanInterval
    dailySlices = sloth.toolBox.get_intervalSlice(dates=dates, sliceInterval=meanInterval)
    # Loop over all slices, mask 'missing' values with np.nan, and calculate 
    # related mean. The averaged data gets appended for each file and slice.
    for Slice in dailySlices:
        tmp_time     = timeValues[Slice]
        tmp_time     = tmp_time.filled(fill_value=np.nan)
        tmp_var      = data[Slice]
        tmp_var_mask = tmp_var.mask
        if not tmp_var_mask.any():
            tmp_var_mask  = np.zeros(tmp_var.shape, dtype=bool)
        tmp_var       = tmp_var.filled(fill_value=np.nan)

        tmp_timeMean  = np.nanmean(tmp_time, axis=0, keepdims=True)
        # tmp_timeMean  = np.ma.mean(tmp_time, axis=0, keepdims=True)
        tmp_timeMean  = nc.num2date(tmp_timeMean, units=timeUnits,calendar=timeCalendar)
        tmp_monthMean = np.nanmean(tmp_var, axis=0, keepdims=True, dtype=float)
        # tmp_monthMean = np.ma.mean(tmp_var, axis=0, keepdims=True, dtype=float)

        tmp_IntervalMean.append(tmp_monthMean)
        tmp_IntervalTime.append(tmp_timeMean)
        print(f'len(tmp_IntervalMean): {len(tmp_IntervalMean)}')


    print('##################')
    print('##################')
intervalMean = np.concatenate(tmp_IntervalMean, axis=0)
print(f'intervalMean.shape: {intervalMean.shape}')
intervalTime = np.concatenate(tmp_IntervalTime, axis=0)
# save everything for later use         
if not os.path.exists(f'../data/example_ClimateMeans/'):
            os.makedirs(f'../data/example_ClimateMeans/')
with open(f'../data/example_ClimateMeans/intervalMean_{meanInterval}.npy', 'wb') as f:
    np.save(f, intervalMean)
with open(f'../data/example_ClimateMeans/intervalTime_{meanInterval}.npy', 'wb') as f:
    np.save(f, intervalTime)


# First create an empty array of same shape as intervalMean but with 
# t-dim (0-axis) = NoI
climaDim = [NoI]    
climaDim = climaDim + [dim for dim in intervalMean[0].shape]
clima = np.empty(climaDim)
    
for curr_day in range(NoI):
    # The climat mean is simple the mean of all entries with NoI in distance.
    clima[curr_day] = np.nanmean(intervalMean[curr_day::NoI], axis=0, dtype=float)        
# dump climate data
with open(f'../data/example_ClimateMeans/climate_{meanInterval}.npy', 'wb') as f:
    np.save(f, clima)

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