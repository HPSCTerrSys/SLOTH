import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import datetime as dt
import netCDF4 as nc
import argparse
import sys
import os
import glob
import cftime

'''
assuming time axis is first axis! (t,z,y,x) or (t,y,x)!
'''

def slice_months(data, dates, dumpIntervall, calendar):
    ''' 
    assuming at least monthly data!
    assuming at leas hourly steps / dumpIntervalls! 
    assuming no seconds in data!
    '''

    print(f'########################################################################')
    print(f'#### checking dates-dim matched data-first-dim')
    print(f'########################################################################')
    if dates.shape[0] != data.shape[0]:
        print(f'ERROR: first dim of data {data.shape[0]} does not match dim of dates {dates.shape[0]}')
        return False
    else:
        print(f'DONE')

    print(f'########################################################################')
    print(f'#### checking dt == dumpIntervall for all data-points')
    print(f'########################################################################')
    # checking if dt is equal to dumpInterval
    tmp_dates1   = dates.copy()
    tmp_dates2   = np.roll(dates, 1, axis=0)
    allEqualDump = (tmp_dates1 - tmp_dates2) == dumpIntervall
    # set first entry to True, as first step is always False due to roll()
    allEqualDump[0] = True
    if not allEqualDump.all():
        print('ERROR: calculated dt does not match dumpIntervall for all steps')
        return False
    else:
        print(f'DONE')
    
    print(f'########################################################################')
    print(f'#### Finding first of month')
    print(f'########################################################################')
    Offset = 0
    while True:
        # verify first date is first of month.
        # idea: first step - dumpIntervall (or 0.5*dumpIntervall) is first of month 0UTC
        print(f'Offset: {Offset}')
        tmp_first = dates[Offset]
        # get first of current month at midnight
        tmp_first_ref = tmp_first.replace(day=1, hour=0, minute=0, second=0)
        if (tmp_first - dumpIntervall) == tmp_first_ref:
            print(f'check step {Offset} is first of a month')
            break
        elif (tmp_first - (0.5*dumpIntervall)) == tmp_first_ref:
            print(f'check step {Offset} is first of a month')
            break
        else:
            print(f'ERROR: step {Offset} is not first step of month!')
            print(f'step: {tmp_first}')
            Offset += 1
            # hard break of loop
            if Offset > 50:
                break
            continue

    print(f'########################################################################')
    print(f'#### getting month series / slices')
    print(f'########################################################################')
    # there should be some clever and short solution for below loop ...!
    # now that we know dt=dumpIntervall and first step is first of month, loop over all data
    t_lower = Offset
    Slices = []
    for t in range(Offset, dates.shape[0]):
        currMonth = dates[t].month
        # I decided to go for the solution checking current month and next month
        # to catch the case if the dateset contains one month only!
        #print(f'BE LOW IS AH ACK') # the hack for LPOs COSMO output
        nextMonth = (dates[t] + 2*dumpIntervall).month
        #nextMonth = (dates[t] + dumpIntervall).month
        if nextMonth != currMonth:
            Slices.append(slice(t_lower,t+1,None))
            t_lower = t+1

    return Slices

def slice_days(data, dates, dumpIntervall, calendar):
    ''' 
    assuming at least monthly data!
    assuming at leas hourly steps / dumpIntervalls! 
    assuming no seconds in data!
    '''

    print(f'########################################################################')
    print(f'#### checking dates-dim matched data-first-dim')
    print(f'########################################################################')
    if dates.shape[0] != data.shape[0]:
        print(f'ERROR: first dim of data {data.shape[0]} does not match dim of dates {dates.shape[0]}')
        return False
    else:
        print(f'DONE')

    print(f'########################################################################')
    print(f'#### checking dt == dumpIntervall for all data-points')
    print(f'########################################################################')
    # checking if dt is equal to dumpInterval
    tmp_dates1   = dates.copy()
    tmp_dates2   = np.roll(dates, 1, axis=0)
    allEqualDump = (tmp_dates1 - tmp_dates2) == dumpIntervall
    # set first entry to True, as first step is always False due to roll()
    allEqualDump[0] = True
    if not allEqualDump.all():
        print('ERROR: calculated dt does not match dumpIntervall for all steps')
        return False
    else:
        print(f'DONE')
    
    print(f'########################################################################')
    print(f'#### Finding first of month at midnight')
    print(f'########################################################################')
    Offset = 0
    while True:
        # verify first date is first of month.
        # idea: first step - dumpIntervall (or 0.5*dumpIntervall) is first of month 0UTC
        print(f'Offset: {Offset}')
        tmp_first = dates[Offset]
        # get first of current month at midnight
        tmp_first_ref = tmp_first.replace(day=1, hour=0, minute=0, second=0)
        if (tmp_first - dumpIntervall) == tmp_first_ref:
            print(f'check step {Offset} is first of a month at midnight')
            break
        elif (tmp_first - (0.5*dumpIntervall)) == tmp_first_ref:
            print(f'check step {Offset} is first of a month at midnight')
            break
        else:
            print(f'ERROR: step {Offset} is not first step of month at midnight!')
            print(f'step: {tmp_first}')
            Offset += 1
            # hard break of loop
            if Offset > 50:
                break
            continue

    print(f'########################################################################')
    print(f'#### getting month series / slices')
    print(f'########################################################################')
    # there should be some clever and short solution for below loop ...!
    # now that we know dt=dumpIntervall and first step is first of month, loop over all data
    t_lower = Offset
    Slices = []
    for t in range(Offset, dates.shape[0]):
        currDay = dates[t].day
        # I decided to go for the solution checking current month and next month
        # to catch the case if the dateset contains one month only!
        #print(f'BE LOW IS AH ACK') # the hack for LPOs COSMO output
        #nextDay = (dates[t] + 2*dumpIntervall).day
        nextDay = (dates[t] + dumpIntervall).day
        if nextDay != currDay:
            Slices.append(slice(t_lower,t+1,None))
            t_lower = t+1

    return Slices

###############################################################################
###############################################################################
###############################################################################
#### you have to adjust things below only

# set your dumpintervall
dumpIntervall = dt.timedelta(hours=3)
# set the root dir of your data
dataRootDir = f'/p/scratch/cjibg35/poshyvailo1/HiCam-CORDEX_EUR-11_MPI-ESM-LR_histo_r1i1p1_FZJ-IBG3-TSMP120EC_v00aJuwelsCpuProdTt-1949_2005/postpro/'
# set the variable name of your data
varName = 'T_2M' # T_2M
# set you file name
# for our postpor you have to distinguish between COSMO and CLM/PFL only. COSMO does contain a '_ts' in it filename
fileName = f'{varName}_ts.nc' # CLM/PFL style
#fileName = f'{varName}_ts.nc' # COSMO style

# use wildcats here as for Linux terminal
# e.g. 1997_?? does contain every month of 1997
# e.g. 199*    does contain every month of every year starting with 199
# etc.
files = sorted(glob.glob(f'{dataRootDir}/199?_??/{fileName}'))

#### you have to adjust things above only
###############################################################################
###############################################################################
###############################################################################
###############################################################################


###############################################################################
#### READ IN ALL DATA AND CALCULATE MONTHLY MEAN and relate middle time-step
###############################################################################
try:
    dailyMean = np.load(f'DailyMeans_TestDump.npy')
    print(f'dailyMean.shape: {dailyMean.shape}')
    print(f'dailyMean.dtype: {dailyMean.dtype}')
    dailyTime = np.load(f'DailyTime_TestDump.npy', allow_pickle=True)
    print(f'dailyTime.shape: {dailyTime.shape}')
except FileNotFoundError:
    dailyMean = None
    dailyTime = None
    for file in files:
        with nc.Dataset(file, 'r') as nc_file:
            data  = nc_file.variables[varName]
            print(f'data.shape: {data.shape}')
            # add some chekc here if ndims is correct
            data         = data[...]
            nc_time      = nc_file.variables['time']
            calendar     = nc_time.calendar
            dates        = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
            timeValues   = nc_time[:]
            timeCalendar = nc_time.calendar
            timeUnits    = nc_time.units

        dailySlices = slice_days(data=data, dates=dates, calendar=calendar, dumpIntervall=dumpIntervall)
        print(dailySlices)
        # this loop is needed if there are more than one day in data and 
        # therefore dailySlices contains more than one entry...
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
            if dailyMean is None:
                dailyMean = tmp_monthMean
                dailyTime = tmp_timeMean
            else:
                dailyMean = np.append(dailyMean, tmp_monthMean, axis=0)
                dailyTime = np.append(dailyTime, tmp_timeMean, axis=0)
            print(f'dailyMean.shape: {dailyMean.shape}')
        print('##################')
        print('##################')
    # dump daily mean
    with open('DailyMeans_TestDump.npy', 'wb') as f:
        np.save(f, dailyMean)
    with open('DailyTime_TestDump.npy', 'wb') as f:
        np.save(f, dailyTime)

###############################################################################
#### CALCULATE CLIMATOLOGY
###############################################################################
# number of days (NoD) should be 365, but I want this flexible, e.g. if one 
# does pass data for one day (NoD=1) or similar.
NoD = 365
try:
    clima = np.load(f'dailyClima_TestDump.npy')
    print(f'dailyClima.shape: {clima.shape}')
    print(f'dailyClima.dtype: {clima.dtype}')
except FileNotFoundError:
    climaDim = [NoD]
    climaDim = climaDim + [dim for dim in dailyMean[0].shape]
    clima = np.empty(climaDim)
    for month in range(NoD):
        clima[month] = np.nanmean(dailyMean[month::NoD], axis=0, dtype=float)
    # dump clima
    with open('DailyClima_TestDump.npy', 'wb') as f:
        np.save(f, clima)

    print(f'dailyClima.shape: {clima.shape}')
    print(f'np.nanmax(clima): {np.nanmax(clima)}')
    print(f'np.nanmin(clima): {np.nanmin(clima)}')

###############################################################################
##### CALCULATE ANOMALIES
###############################################################################
try:
    dailyAnomalyDomain = np.load(f'DailyAnomalies_TestDump.npy')
    print(f'dailyAnomalyDomain.shape: {dailyAnomalyDomain.shape}')
    print(f'dailyAnomalyDomain.dtype: {dailyAnomalyDomain.dtype}')
except FileNotFoundError:
    dailyAnomalyDomain = np.empty(dailyMean.shape[0])
    for n, day in enumerate(dailyMean):
        idx_day = n%NoD
        dailyAnomalyDomain[n] = np.nanmean(day, dtype=float) - np.nanmean(clima[idx_day], dtype=float)
    # dump anomalies
    with open('DailyAnomalies_TestDump.npy', 'wb') as f:
        np.save(f, dailyAnomalyDomain)

###############################################################################
#### Plot stuff
###############################################################################
fig, ax = plt.subplots(figsize=(8,4))
# makes a grid on the background
ax.grid()
ax.plot(dailyAnomalyDomain, color='black')

# BELOW IS UGLY STYLE AND HACKY...!
# We need to draw the canvas, otherwise the labels won't be positioned and 
# won't have values yet.
# fig.canvas.draw()
# labels_idx = [item.get_text() for item in ax.get_xticklabels()]
# labels_idx = [ int(item) for item in labels_idx[1:-2]]
# new_labels = [ cftime.datetime.strftime(item, '%Y-%m') for item in monthlyTime[labels_idx] ]
# print(new_labels)
# ax.set_xticks(labels_idx)
# ax.set_xticklabels(new_labels)
labels_idx = np.arange(0,dailyTime.shape[0], NoD)
ax.set_xticks(labels_idx)
labels = [ cftime.datetime.strftime(item, '%Y') for item in dailyTime[labels_idx] ]
ax.set_xticklabels(labels)

x = np.arange(dailyAnomalyDomain.shape[0])
ax.fill_between(x, 0, dailyAnomalyDomain, where=dailyAnomalyDomain>0, facecolor='red', interpolate=True, alpha=0.75) 
ax.fill_between(x, 0, dailyAnomalyDomain, where=dailyAnomalyDomain<0, facecolor='blue', interpolate=True, alpha=0.75)
plt.savefig(f'DailyAnomalieTestPlot_{varName}.png', dpi=380) 
