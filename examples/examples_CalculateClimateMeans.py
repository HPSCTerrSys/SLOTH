""" Example script to calculate climatology means.

Calculating means with python is not a difficult task and even the choice  
of the interval we are handling (e.g. daily- or monthly means) does not really 
matter. At the end we do simply have to apply e.g. np.mean() to our data and
we are done!
However, extracting the correct slices or intervals out of a given data set to 
apply the mean-function to, could be very nasty. 
The correct slicing depends on the actual data we are facing as e.g. the 
starting time steps and the time resolution (hours, daily, monthly, or even 
yearly data?). Individual datasets need individual treatment.

This example make use of the function 'get_intervalSlice()' which aims to 
generalize the progress to get the correct slices or intervals.
Mandatory for this purpose is, the handled data does contain the correct 
time information, as e.g. the correct time-variables with netCDF files.
If this is the case the function does loop over all your dates and returns
a list of slices to extract the correct intervals out of your data. If you 
data-set does consist out of individual files or is one big files does not
matter.
"""
import numpy as np
import datetime as dt
import netCDF4 as nc
import matplotlib as mpl
import sys
import os
import glob
import cftime
src_path='../src/'
sys.path.append(src_path)
# ANalysisTool does contain the named function 'get_intervalSlice()'
import ANalysisTool as ANT
import PlotLib 


###############################################################################
#### Some settings
###############################################################################
# You have to define the dumpinterval of your mode-output. 
dumpIntervall = dt.timedelta(hours=3)

# You have to define the interval you are interested in. Currently supported:
# 'day', 'month' 
meanInterval = 'month'

###############################################################################
#### CALCULATE CLIMATOLOGY
###############################################################################
# Climatology calculations are always based on comparisons of individual 
# intervals between different years. 
# If the interval we are interested in is 'month', we do compare the same 
# month between different years.
# If the interval we are interested in is 'day', we do compare the same 
# day between different years.
# etc.
# IN any case we do have to know how many intervals a year does cover. 
# For daily based calculation this usually is NoI=365, for monthly based 
# calculations this usually is NoI=12.
# However we could also handle summer months (JJA) only and do a monthly based
# calculation, in which case NoI=3 etc.!
NoI = 12

# set the root dir of your data
dataRootDir = f'/p/scratch/cjibg35/poshyvailo1/HiCam-CORDEX_EUR-11_MPI-ESM-LR_histo_r1i1p1_FZJ-IBG3-TSMP120EC_v00aJuwelsCpuProdTt-1949_2005/postpro/postpro_before_18032021'
# set the variable name of your data
varName = 'TSA' # T_2M
# set you file name
# for our postpor you have to distinguish between COSMO and CLM/PFL only. COSMO does contain a '_ts' in it filename
fileName = f'{varName}.nc' # CLM/PFL style
#fileName = f'{varName}_ts.nc' # COSMO style

# use wildcats here as for Linux terminal
# e.g. 1997_?? does contain every month of 1997
# e.g. 199*    does contain every month of every year starting with 199
# etc.
#files = sorted(glob.glob(f'{dataRootDir}/199?_??/{fileName}'))
date_start_clim = cftime.datetime(1965, 6, 1, calendar='noleap')
date_final_clim = cftime.datetime(1986, 1, 1, calendar='noleap')
files = []
for year in range(date_start_clim.year,date_final_clim.year+1):
     tmp_files = sorted(glob.glob(f'{dataRootDir}/{year}_??/{fileName}')) #all month within the year
     files += tmp_files
     print(tmp_files)


###############################################################################
#### READ IN ALL DATA AND CALCULATE DAILY MEAN and relate middle time-step
###############################################################################
# Already create variables we want to store our data at
intervalMean = None
tmp_IntervalMean = []
intervalTime = None
tmp_IntervalTime = []
# Loop over all files of our data-set.
for file in files:
    # Open the file
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

    # Calculate the slices for the current file based on the passed meanInterval
    dailySlices = ANT.toolBox.get_intervalSlice(dates=dates, dumpIntervall=dumpIntervall, sliceInterval=meanInterval)
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
        tmp_timeMean  = nc.num2date(tmp_timeMean, units=timeUnits,calendar=timeCalendar)
        tmp_monthMean = np.nanmean(tmp_var, axis=0, keepdims=True, dtype=float)

        #######################################################################
        # START -- IMPORTANT NOTE FROM AND FOR NWR:
        # appending individual intervals to list is MUCH faster than appending 
        # to ndarray (np.append) directly!
        # I guess this is because list contains pointers and appending a 
        # pointer is quiet fast, while np.append() open and restructure the
        # entire array which is getting slower and slower while reading in more
        # and more data / files.
        tmp_IntervalMean.append(tmp_monthMean)
        tmp_IntervalTime.append(tmp_timeMean)
        print(f'len(tmp_IntervalMean): {len(tmp_IntervalMean)}')

        # if intervalMean is None:
        #     intervalMean = tmp_monthMean
        #     intervalTime = tmp_timeMean
        # else:
        #     intervalMean = np.append(intervalMean, tmp_monthMean, axis=0)
        #     intervalTime = np.append(intervalTime, tmp_timeMean, axis=0)
        # print(f'intervalMean.shape: {intervalMean.shape}')
        # END -- IMPORTANT NOTE FROM AND FOR NWR:
        #######################################################################

    print('##################')
    print('##################')
intervalMean = np.concatenate(tmp_IntervalMean, axis=0)
print(f'intervalMean.shape: {intervalMean.shape}')
intervalTime = np.concatenate(tmp_IntervalTime, axis=0)
# save everything for later use 
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
with open(f'../data/example_ClimateMeans/climate_{meanInterval}_{date_start_clim.strftime("%Y_%m")}_{date_final_clim.strftime("%Y_%m")}.npy', 'wb') as f:
    np.save(f, clima)

# plot if meanInterval='month'
if meanInterval == 'month':
    kwargs = {
            'title': 'Test Clima plot',
            #'title': '\n'.join(tmp_titlesubstr),
            'infostr': True,
            #'var_vmax': 1,
            #'var_vmin': 0,
            'var_cmap': mpl.cm.get_cmap('jet'),
            'saveFile': f'./examples_CalculateClimateMeans.pdf',
            #'dpi': 100,
            'figsize': (10, 4),
            }
    PlotLib.plot_ClimateYearMonth(clima, **kwargs)