import numpy as np
import matplotlib as mpl
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

def prep_monthMean(data, dates, dumpIntervall, calendar):
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
    # set first entry to True, as first step is always False due to toll()
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
        # idea: first step - dumpIntervall (or 0.5*dumpIntervall) is first of month
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
    # there should be some cleaver and short solution for below loop ...!
    # now that we know dt=dumpIntervall and first step is first of month, loop over all data
    t_lower = Offset
    Slices = []
    for t in range(Offset, dates.shape[0]):
        currMonth = dates[t].month
        # I decided to go for the solution checking current month and next month
        # to catch the case if the dateset contains one month only!
        nextMonth = (dates[t] + dumpIntervall).month
        if nextMonth != currMonth:
            Slices.append(slice(t_lower,t+1,None))
            t_lower = t+1

    return Slices

dumpIntervall = dt.timedelta(hours=3)
dataRootDir = f'/work/workflowuser01/tmp_postpro/'
#varName = 'WIND'
varName = 'TMAX_2M'
#fileName = f'{varName}.nc'
fileName = f'{varName}_ts.nc'

files = sorted(glob.glob(f'{dataRootDir}/*/{fileName}'))
#files = sorted(glob.glob(f'{dataRootDir}/1996_12/{fileName}'))

# READ IN ALL DATA AND CALCULATE MONTHLY MEAN
climaRaw = None
for file in files:
    with nc.Dataset(file, 'r') as nc_file:
        data = nc_file.variables[varName]
        print(f'data.shape: {data.shape}')
        # add some chekc here if ndims is correct
        data     = data[...]
        nc_time  = nc_file.variables['time']
        calendar = nc_time.calendar
        dates    = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)

    monthlySlices = prep_monthMean(data=data, dates=dates, calendar=calendar, dumpIntervall=dumpIntervall)
    print(monthlySlices)
    # this loop is needed if there are more than one month in data and 
    # therefore monthlySlices contains more than one entry...
    for Slice in monthlySlices:
        tmp_monthMean = np.nanmean(data[Slice], axis=0, keepdims=True)
        if climaRaw is None:
            climaRaw = tmp_monthMean
        else:
            climaRaw = np.append(climaRaw, tmp_monthMean, axis=0)
        print(f'climaRaw.shape: {climaRaw.shape}')
    print('##################')
    print('##################')

# CALCULATE CLIMATOLOGY
clima = np.empty()
for month in range(12):
    if clima is None:
        clima = 
