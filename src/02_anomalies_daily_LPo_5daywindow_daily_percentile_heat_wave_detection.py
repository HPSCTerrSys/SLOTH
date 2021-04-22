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
from scipy import signal
#from scipy.interpolate import interp1d
#from scipy.signal import savgol_filter
#from scipy.signal import filtfilt
#from scipy.signal import ellip
import scipy.ndimage as ndimage
from datetime import date

'''
assuming time axis is first axis! (t,z,y,x) or (t,y,x)!
original code is from N.Wagner, modified L.Poshyvailo-Strube 2021.
LPo: Some things are takes also from here https://github.com/ecjoliver/marineHeatWaves/blob/master/marineHeatWaves.py
'''

varName, vmin, vmax =('TSA', 10, 24) #Celsius
path_TSMP='/p/scratch/cjibg35/poshyvailo1/HiCam-CORDEX_EUR-11_MPI-ESM-LR_histo_r1i1p1_FZJ-IBG3-TSMP120EC_v00aJuwelsCpuProdTt-1949_2005/postpro/postpro_before_18032021'
year_start_clim=1965
year_final_clim=1985

year_start_hw=1972 #hw=heatwave year
year_final_hw=1972 #hw=heatwave
month_start_hw=6 #hw=heatwave year
month_final_hw=8 #hw=heatwave

#currently it works only for 5-d window
window_percentile = 5
minDuration=6 #of heat wave
#joinAcrossGaps=True
#maxGap=1

pixel_x=250
pixel_y=300

hw_investigation='yes'

mhw = {}
mhw['time_start'] = [] # datetime format; Start time of MHW [datetime format]
mhw['time_end'] = [] # datetime format;  End time of MHW [datetime format]
mhw['date_start'] = [] # datetime format
mhw['date_end'] = [] # datetime format
mhw['index_start'] = []
mhw['index_end'] = []
mhw['index_peak'] = []
mhw['date_peak'] = [] # datetime format
mhw['duration'] = [] # [days]
mhw['max_temp'] = [] # [days]
mhw['intensity_max'] = [] # [deg C]
mhw['intensity_mean'] = [] # [deg C]
#mhw['intensity_var'] = [] # [deg C]
mhw['intensity_cumulative'] = [] # [deg C]
mhw['intensity_max_relThresh'] = [] # [deg C]
mhw['intensity_mean_relThresh'] = [] # [deg C]
#mhw['intensity_var_relThresh'] = [] # [deg C]
mhw['intensity_cumulative_relThresh'] = [] # [deg C]
mhw['intensity_max_abs'] = [] # [deg C]
mhw['intensity_mean_abs'] = [] # [deg C]
#mhw['intensity_var_abs'] = [] # [deg C]
mhw['intensity_cumulative_abs'] = [] # [deg C]
#mhw['category'] = []
#mhw['rate_onset'] = [] # [deg C / day]
#mhw['rate_decline'] = [] # [deg C / day]

        
############################################
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
    print(f'#### getting daily series / slices')
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
dataRootDir = f'{path_TSMP}'
#dataRootDir = f'/p/scratch/cjibg35/poshyvailo1/HiCam-CORDEX_EUR-11_MPI-ESM-LR_histo_r1i1p1_FZJ-IBG3-TSMP120EC_v00aJuwelsCpuProdTt-1949_2005/postpro/'
# set the variable name of your data
#varName = 'T_2M' # T_2M
# set you file name
# for our postpor you have to distinguish between COSMO and CLM/PFL only. COSMO does contain a '_ts' in it filename
fileName = f'{varName}.nc' # CLM/PFL style
#fileName = f'{varName}_ts.nc' # COSMO style

# use wildcats here as for Linux terminal
# e.g. 1997_?? does contain every month of 1997
# e.g. 199*    does contain every month of every year starting with 199
# etc.
#files = sorted(glob.glob(f'{dataRootDir}/199?_??/{fileName}'))

files = []
for year in range(year_start_clim,year_final_clim+1):
    for month in range(month_start_hw, month_final_hw+1):
         month_file=str(month).zfill(2)
         tmp_files = sorted(glob.glob(f'{dataRootDir}/{year}_{month_file}/{fileName}')) #all month within the year
         files += tmp_files
         print(tmp_files)


#### you have to adjust things above only
###############################################################################
###############################################################################
###############################################################################
###############################################################################


###############################################################################
#### READ IN ALL DATA AND CALCULATE DAILY MEAN and relate middle time-step
###############################################################################
try:
    dailyMean = np.load(f'DailyMeans_TestDump.npy')
    print(f'dailyMean.shape: {dailyMean.shape}')
    print(f'dailyMean.dtype: {dailyMean.dtype}')
    dailyTime = np.load(f'DailyTime_TestDump.npy', allow_pickle=True)
    print(f'dailyTime.shape: {dailyTime.shape}')
    print('Daily mean:', dailyMean[:,pixel_y,pixel_x])       
    print('MIN of daily mean:', np.nanmin(dailyMean[:,pixel_y,pixel_x]))
    print('MAX of daily mean:', np.nanmax(dailyMean[:,pixel_y,pixel_x]))
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
    print('Daily mean:', dailyMean[2,pixel_y,pixel_x])
    print('MIN of daily mean:', np.nanmin(dailyMean[:,pixel_y,pixel_x]))
    print('MAX of daily mean:', np.nanmax(dailyMean[:,pixel_y,pixel_x]))
       
###############################################################################
#### CALCULATE CLIMATOLOGY
###############################################################################
# number of days (NoD) should be 365, but I want this flexible, e.g. if one 
# does pass data for one day (NoD=1) or similar.
sum_days=0
for month in range(month_start_hw, month_final_hw+1):
     if( month==1 or month==3 or month==5 or month==7 or month==8 or month==10 or month==12):
         nt = 31
     elif(month==2):
         nt = 28
     else:
         nt =30
     sum_days=nt+sum_days
     print('sum processed days = ', sum_days)

###############################################
NoD = sum_days #214 #(len(dailyTime)/(year_start_clim-year_final_clim+1), dtype=int)

month_start_hw_print=str(month_start_hw).zfill(2)
month_final_hw_print=str(month_final_hw).zfill(2)
try:
    clima = np.load(f'DailyClima_TestDump_'+str(year_start_clim)+'_'+str(month_start_hw).zfill(2) + '_' +str(year_final_clim)+'_'+  str(month_final_hw).zfill(2) + '.npy')        
    print(f'dailyClima.shape: {clima.shape}')
    print(f'dailyClima.dtype: {clima.dtype}')
    print('Clima:', clima[2,pixel_y,pixel_x])
except FileNotFoundError:
    climaDim = [NoD]    
    climaDim = climaDim + [dim for dim in dailyMean[0].shape]
    clima = np.empty(climaDim)
        
    for curr_day in range(NoD):
        clima[curr_day] = np.nanmean(dailyMean[curr_day::NoD], axis=0, dtype=float)        
    # dump clima
    with open('DailyClima_TestDump_'+str(year_start_clim)+'_'+str(month_start_hw).zfill(2) + '_' +str(year_final_clim)+'_'+  str(month_final_hw).zfill(2) + '.npy', 'wb') as f:
        np.save(f, clima)

    print(f'dailyClima.shape: {clima.shape}')
    print(f'np.nanmax(clima): {np.nanmax(clima)}')
    print(f'np.nanmin(clima): {np.nanmin(clima)}')
    print('Clima:', clima[:,pixel_y,pixel_x])

############### 90th percentile from all days of the selcted month (see at the beginning of the code), and the period covering climatologies
#clima_90p=np.percentile(dailyMean, 90, axis=0)    
#climaDim_90p=[NoD]
#climaDim_90p = climaDim_90p + [dim for dim in dailyMean[0].shape]
#clima_90p_full = np.empty(climaDim_90p)
#for curr_day in range(NoD):
#       clima_90p_full[curr_day]=clima_90p
#print('!!!!!!!!!!!!!!!!!!!!!!!!! 90 percentile', clima_90p_full[2,pixel_y,pixel_x])

############### 97.5th percentile from all days of the selcted month (see at the beginning of the code), and the period covering climatologies
#clima_975p=np.percentile(dailyMean, 97.5, axis=0)    
#climaDim_975p=[NoD]
#climaDim_975p = climaDim_975p + [dim for dim in dailyMean[0].shape]
#clima_975p_full = np.empty(climaDim_975p)
#for curr_day in range(NoD):
#       clima_975p_full[curr_day]=clima_975p
       
############### 99.5th percentile from all days of the selcted month (see at the beginning of the code), and the period covering climatologies
#clima_995p=np.percentile(dailyMean, 99.5, axis=0)    
#climaDim_995p=[NoD]
#climaDim_995p = climaDim_995p + [dim for dim in dailyMean[0].shape]
#clima_995p_full = np.empty(climaDim_995p)
#for curr_day in range(NoD):
#       clima_995p_full[curr_day]=clima_995p

       
############### 90th percentile for each day -- does not work so well
climaDim_90p = [NoD]    #dimension of 90th percentile variable
climaDim_90p = climaDim_90p + [dim for dim in dailyMean[0].shape]
clima_90p = np.empty(climaDim_90p)
print(f'dailyClima_90p.shape: {clima_90p.shape}')
print(dailyMean.shape)

############## calcualting percentile with 5-days window
for curr_day in range(NoD):
    if (curr_day == 0): #for the first day of the period -- take the first day +2 days after for calculating the percentile
        window_temp=3
        dim_temp=window_temp*(year_final_clim-year_start_clim+1)
        daily_mean_temp_dim=[dim_temp]
        daily_mean_temp_dim=daily_mean_temp_dim + [dim for dim in dailyMean[0].shape]
        daily_mean_temp = np.empty(daily_mean_temp_dim)
        print(f'daily_mean_temp: {daily_mean_temp.shape}')        
        for index_temp in range(window_temp):
             dim_temp=(index_temp+1)*(year_final_clim-year_start_clim+1)
             a_temp=index_temp*(year_final_clim-year_start_clim+1)            
             daily_mean_temp[a_temp:dim_temp]=dailyMean[index_temp::NoD]
             #print('daily_temp', daily_mean_temp[a_temp:dim_temp,pixel_y,pixel_x])
             #print('daily_MEAN!!!', dailyMean[0:3,pixel_y,pixel_x])
        clima_90p[curr_day]=np.percentile(daily_mean_temp, 90, axis=0)
        #print(clima_90p[0,250,300])
    
    if (curr_day == 1): # take the second day +2 days after +1 day before  for calculating the percentile
        window_temp=4
        dim_temp=window_temp*(year_final_clim-year_start_clim+1)
        daily_mean_temp_dim=[dim_temp]
        daily_mean_temp_dim=daily_mean_temp_dim + [dim for dim in dailyMean[0].shape]
        daily_mean_temp = np.empty(daily_mean_temp_dim)
        print(f'daily_mean_temp: {daily_mean_temp.shape}')        
        for index_temp in range(window_temp):
             dim_temp=(index_temp+1)*(year_final_clim-year_start_clim+1)
             a_temp=index_temp*(year_final_clim-year_start_clim+1)
             print(a_temp,dim_temp)
             daily_mean_temp[a_temp:dim_temp]=dailyMean[index_temp::NoD]
             #print('daily_temp', daily_mean_temp[a_temp:dim_temp,pixel_y,pixel_x])
             #print('daily_MEAN!!!', dailyMean[0:4,pixel_y,pixel_x])
        clima_90p[curr_day]=np.percentile(daily_mean_temp, 90, axis=0)
        print('90th percentile',clima_90p[curr_day,pixel_y,pixel_x])        

    if (curr_day > 1) and (curr_day < NoD-2): #5d window: take a day +-2 days for percentile calcualtions
        dim_temp=window_percentile*(year_final_clim-year_start_clim+1)
        daily_mean_temp_dim=[dim_temp]
        daily_mean_temp_dim=daily_mean_temp_dim + [dim for dim in dailyMean[0].shape]
        daily_mean_temp = np.empty(daily_mean_temp_dim)
        print(f'daily_mean_temp: {daily_mean_temp.shape}')        
        for index_temp in range(window_percentile):
             dim_temp=(index_temp+1)*(year_final_clim-year_start_clim+1)
             a_temp=index_temp*(year_final_clim-year_start_clim+1)
             b_temp=curr_day+index_temp-2
             print(a_temp,dim_temp,b_temp)
             daily_mean_temp[a_temp:dim_temp]=dailyMean[b_temp::NoD]
             #print('daily_temp', daily_mean_temp[a_temp:dim_temp,pixel_y,pixel_x])
             #print('daily_MEAN!!!', dailyMean[0:4,pixel_y,pixel_x])
        clima_90p[curr_day]=np.percentile(daily_mean_temp, 90, axis=0)
        print('90th percentile',clima_90p[curr_day,pixel_y,pixel_x])

    if (curr_day == NoD-2): # take the second day from the end +2 days before  +1 day after for calculating the percentile
        window_temp=4
        dim_temp=window_temp*(year_final_clim-year_start_clim+1)
        daily_mean_temp_dim=[dim_temp]
        daily_mean_temp_dim=daily_mean_temp_dim + [dim for dim in dailyMean[0].shape]
        daily_mean_temp = np.empty(daily_mean_temp_dim)
        print(f'daily_mean_temp: {daily_mean_temp.shape}')
        for index_temp in range(window_temp):
             dim_temp=(index_temp+1)*(year_final_clim-year_start_clim+1)
             a_temp=index_temp*(year_final_clim-year_start_clim+1)
             print(a_temp,dim_temp)
             b_temp=curr_day+index_temp-2
             daily_mean_temp[a_temp:dim_temp]=dailyMean[b_temp::NoD]
             #print('daily_temp', daily_mean_temp[a_temp:dim_temp,pixel_y,pixel_x])
             #print('daily_MEAN!!!', dailyMean[0:4,pixel_y,pixel_x])
        clima_90p[curr_day]=np.percentile(daily_mean_temp, 90, axis=0)
        print('90th percentile',clima_90p[curr_day,pixel_y,pixel_x])

    if (curr_day == NoD-1): # take the last day from he end +2 days before  +1 day after for calculating the percentile
        window_temp=3
        dim_temp=window_temp*(year_final_clim-year_start_clim+1)
        daily_mean_temp_dim=[dim_temp]
        daily_mean_temp_dim=daily_mean_temp_dim + [dim for dim in dailyMean[0].shape]
        daily_mean_temp = np.empty(daily_mean_temp_dim)
        print(f'daily_mean_temp: {daily_mean_temp.shape}')
        for index_temp in range(window_temp):
             dim_temp=(index_temp+1)*(year_final_clim-year_start_clim+1)
             a_temp=index_temp*(year_final_clim-year_start_clim+1)
             print(a_temp,dim_temp)
             b_temp=curr_day+index_temp-2
             daily_mean_temp[a_temp:dim_temp]=dailyMean[b_temp::NoD]
             #print('daily_temp', daily_mean_temp[a_temp:dim_temp,pixel_y,pixel_x])
             #print('daily_MEAN!!!', dailyMean[0:4,pixel_y,pixel_x])
        clima_90p[curr_day]=np.percentile(daily_mean_temp, 90, axis=0)
        print('90th percentile',clima_90p[curr_day,pixel_y,pixel_x])
        
#    clima_90p[curr_day]=np.percentile(dailyMean[curr_day::NoD], 90, axis=0)
#    clima_90p[curr_day]=np.percentile(dailyMean[curr_day::NoD], 90, axis=0)
    #print('!!!!!!!!!!!!!!!!!!!!!!!!!', clima_90p)
    #print('Percentile is calculated from this daily values corresponding to each year: ', dailyMean[curr_day::NoD,pixel_y,pixel_x])
    #print('-----------------')
    
#print('clima_90p: ', clima_90p[:,pixel_y,pixel_x])
clima_90p_full = clima_90p

###############################################################################
#HEATWAVES ANALYSIS
###############################################################################
if (hw_investigation=='yes'):
    dumpIntervall_hw = dt.timedelta(hours=3)
    dataRootDir = f'{path_TSMP}'
    fileName = f'{varName}.nc' # CLM/PFL style
    #fileName = f'{varName}_ts.nc' # COSMO style

    files_hw = []
    for year in range(year_start_hw, year_final_hw+1):
        for month in range(month_start_hw, month_final_hw+1):
             month_file=str(month).zfill(2)
             tmp_files_hw = sorted(glob.glob(f'{dataRootDir}/{year}_{month_file}/{fileName}')) #all month within the year
             files_hw += tmp_files_hw

    ###############################################################################
    #### READ IN ALL DATA AND CALCULATE DAILY MEAN and relate middle time-step
    ###############################################################################
    try:
        month_start_hw_print=str(month_start_hw).zfill(2)
        month_final_hw_print=str(month_final_hw).zfill(2)
        dailyMean_hw = np.load(f'DailyMeans_TestDump_{year_start_hw}_{month_start_hw_print}_{year_final_hw}_{month_final_hw_print}.npy')
        print(f'dailyMean_hw.shape: {dailyMean_hw.shape}')
        print(f'dailyMean_hw.dtype: {dailyMean_hw.dtype}')
        dailyTime_hw = np.load(f'DailyTime_TestDump_{year_start_hw}_{month_start_hw_print}_{year_final_hw}_{month_final_hw_print}.npy', allow_pickle=True)
        print(f'dailyTime_hw.shape: {dailyTime_hw.shape}')
    except FileNotFoundError:
        dailyMean_hw = None
        dailyTime_hw = None
        for file_hw in files_hw:
            with nc.Dataset(file_hw, 'r') as nc_file_hw:
                data_hw  = nc_file_hw.variables[varName]
                print(f'data_hw.shape: {data_hw.shape}')
                # add some chekc here if ndims is correct
                data_hw         = data_hw[...]
                nc_time_hw      = nc_file_hw.variables['time']
                calendar_hw     = nc_time_hw.calendar
                dates_hw        = nc.num2date(nc_time_hw[:],units=nc_time_hw.units,calendar=nc_time_hw.calendar)
                timeValues_hw   = nc_time_hw[:]
                timeCalendar_hw = nc_time_hw.calendar
                timeUnits_hw    = nc_time_hw.units

            dailySlices_hw = slice_days(data=data_hw, dates=dates_hw, calendar=calendar_hw, dumpIntervall=dumpIntervall_hw)
            print(dailySlices_hw)
            # this loop is needed if there are more than one day in data and 
            # therefore dailySlices contains more than one entry...
            for Slice_hw in dailySlices_hw:
                tmp_time_hw      = timeValues_hw [Slice_hw]
                tmp_time_hw      = tmp_time_hw .filled(fill_value=np.nan)
                tmp_var_hw       = data_hw[Slice_hw]
                tmp_var_mask_hw  = tmp_var_hw.mask
                if not tmp_var_mask_hw .any():
                    tmp_var_mask_hw  = np.zeros(tmp_var_hw .shape, dtype=bool)
                tmp_var_hw        = tmp_var_hw .filled(fill_value=np.nan)
                tmp_timeMean_hw   = np.nanmean(tmp_time_hw , axis=0, keepdims=True)
                tmp_timeMean_hw  = nc.num2date(tmp_timeMean_hw , units=timeUnits_hw,calendar=timeCalendar_hw)
                tmp_monthMean_hw  = np.nanmean(tmp_var_hw , axis=0, keepdims=True, dtype=float)
                if dailyMean_hw is None:
                    dailyMean_hw = tmp_monthMean_hw 
                    dailyTime_hw = tmp_timeMean_hw
                else:
                    dailyMean_hw = np.append(dailyMean_hw, tmp_monthMean_hw, axis=0) 
                    dailyTime_hw = np.append(dailyTime_hw, tmp_timeMean_hw, axis=0)
                print(f'dailyMean_hw.shape: {dailyMean_hw.shape}')
            print('##################')
            print('##################')
        # dump daily mean
        with open('DailyMeans_TestDump_'+str(year_start_hw)+'_'+str(month_start_hw).zfill(2) + '_' + str(year_final_hw)+ '_' +str(month_final_hw).zfill(2) + '.npy', 'wb') as f:
            np.save(f, dailyMean_hw)
        with open('DailyTime_TestDump_' + str(year_start_hw)+'_'+str(month_start_hw).zfill(2) + '_' + str(year_final_hw)+ '_' +str(month_final_hw).zfill(2) + '.npy', 'wb') as f:
            np.save(f, dailyTime_hw)


###################### HEAT WAVE indices calculation: https://github.com/ecjoliver/marineHeatWaves/blob/master/marineHeatWaves.py
#
# Find MHWs as exceedances above the threshold
#
# Time series of "True" when threshold is exceeded, "False" otherwise


################### filtering the noise in daily climatological data ##########
b, a = signal.butter(3, 0.1, btype='lowpass', analog=False) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
#low_passed = signal.filtfilt(b, a, noisy_signal)
#b, a = ellip(4, 0.01, 120, 0.125)  # Filter to be applied.
clima_plot=signal.filtfilt(b, a, clima[:,pixel_y,pixel_x], padlen=15)

################### 90th percentile data ##########
################### filtering the noise in daily 90th percentile data ##########
b, a = signal.butter(3, 0.1, btype='lowpass', analog=False) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
clima_90p_plot_temp=signal.filtfilt(b, a, clima_90p_full[:,pixel_y,pixel_x], padlen=15)

#I do here calculations with NOT SMOOTHED percentile
exceed_bool = dailyMean_hw[:,pixel_y,pixel_x] - clima_90p_full[:,pixel_y,pixel_x]
#exceed_bool = dailyMean_hw[:,pixel_y,pixel_x] - clima_90p_plot_temp[:,pixel_y,pixel_x]
exceed_bool[exceed_bool<=0] = False
exceed_bool[exceed_bool>0] = True

# Fix issue where missing temp vaues (nan) are counted as True
exceed_bool[np.isnan(exceed_bool)] = False
# Find contiguous regions of exceed_bool = True
events, n_events = ndimage.label(exceed_bool)
print('Number of events (groups of days) exceeding 90th percentile:', n_events)

# Find all MHW events of duration >= minDuration
for ev in range(1,n_events+1):
        print("yes")
        event_duration = (events == ev).sum()
        if event_duration < minDuration:
            continue
        #if event_duration >= minDuration:
            #print('Event duration', event_duration)            
        mhw['time_start'].append(dailyTime_hw[np.where(events == ev)[0][0]])
        mhw['time_end'].append(dailyTime_hw[np.where(events == ev)[0][-1]])
        print('Time start',mhw['time_start'])
        print('Time end',mhw['time_end'])      
#print('N events: ', n_events)
#print('----------------------')
#print( np.array(mhw['time_end'][0:])-np.array(mhw['time_start'][0:]))
#print(np.array(mhw['time_end'][0:-1]))
#print(dailyTime_hw)


#################
# Calculate heat wave properties
mhw['n_events'] = len(mhw['time_start'])
#categories = np.array(['Moderate', 'Strong', 'Severe', 'Extreme'])
print('Numbers of the heat-events exceeding min. length of ', str(minDuration), 'days: ', mhw['n_events'], ' event')
    
abs_intensity_array= np.zeros(NoD) #initialize empty array
rel_intensity_array= np.zeros(NoD) #initialize empty array

for ev in range(mhw['n_events']):
        #print(mhw['time_start'][ev]) #datetime.strftime
        print('################ Start of the heat event #', ev) 
        mhw['date_start'].append(mhw['time_start'][ev])
        mhw['date_end'].append(mhw['time_end'][ev])
        # Get SST series during MHW event, relative to both threshold and to seasonal climatology
        tt_start = np.where(dailyTime_hw==mhw['time_start'][ev])[0][0]
        tt_end = np.where(dailyTime_hw==mhw['time_end'][ev])[0][0]
        #print('tt_start', tt_start)
        #print('tt_end:', tt_end)
        print('date_start: ', cftime.datetime.strftime(mhw['time_start'][ev], '%Y-%m-%d'))
        print('date_end: ', cftime.datetime.strftime(mhw['time_end'][ev], '%Y-%m-%d'))      
        #print('dailyTime_hw: ', dailyTime_hw[tt_start:tt_end+1])
        #print('----------------------')
        mhw['index_start'].append(tt_start)
        mhw['index_end'].append(tt_end)
        
        temp_mhw = dailyMean_hw[tt_start:tt_end+1,pixel_y,pixel_x]
        #print('temp_mhv: ',temp_mhw)
        thresh_mhw = clima_90p_full[tt_start:tt_end+1,pixel_y,pixel_x]
        seas_mhw = clima_plot[tt_start:tt_end+1]        
        
        mhw_relSeas = temp_mhw - seas_mhw
        mhw_relThresh = temp_mhw - thresh_mhw #like in Vautard2013, I will call it 'absolute intensity'
        mhw_relThreshNorm = (temp_mhw - thresh_mhw) / (thresh_mhw - seas_mhw)
        #mhw_abs = temp_mhw
        #print('mhw_relThresh:' ,mhw_relThresh)
        #print('mhw_relThreshNorm: ', mhw_relThreshNorm)
        #print('----------------------')
        abs_intensity_array[tt_start:tt_end+1]=mhw_relThresh
        rel_intensity_array[tt_start:tt_end+1]=mhw_relThreshNorm
                
        # Find peak
        #tt_peak = np.argmax(mhw_relSeas)
        tt_peak = np.argmax(mhw_relThresh)
        #print('tt_peak: ', tt_peak)
        #mhw['time_peak'].append(mhw['time_start'][ev] + tt_peak)
        #mhw['date_peak'].append(mhw['time_start'][ev] + tt_peak)
        mhw['max_temp'].append(temp_mhw[tt_peak])
        mhw['date_peak'].append(dailyTime_hw[tt_start+tt_peak])           
        mhw['index_peak'].append(tt_start + tt_peak)
        print('max_temp: ', temp_mhw[tt_peak])
        print('date_peak: ', cftime.datetime.strftime(dailyTime_hw[tt_start+tt_peak], '%Y-%m-%d'))
        #print('index_peak: ', tt_start + tt_peak)
        #print('----------------------')
                        
        # MHW Duration
        mhw['duration'].append(len(mhw_relSeas))
        print('duration: ', len(mhw_relSeas), ' days')
        
        # MHW Intensity metrics
        mhw['intensity_max'].append(mhw_relSeas[tt_peak])
        mhw['intensity_mean'].append(mhw_relSeas.mean())
        #mhw['intensity_var'].append(np.sqrt(mhw_relSeas.var()))
        mhw['intensity_cumulative'].append(mhw_relSeas.sum())
        
        mhw['intensity_max_relThresh'].append(mhw_relThresh[tt_peak])
        mhw['intensity_mean_relThresh'].append(mhw_relThresh.mean())
        print('Maximum of abs. intensity:', mhw_relThresh[tt_peak], ' K')
        #mhw['intensity_var_relThresh'].append(np.sqrt(mhw_relThresh.var()))
        mhw['intensity_cumulative_relThresh'].append(mhw_relThresh.sum())
        
        #mhw['intensity_max_abs'].append(mhw_abs[tt_peak])
        #mhw['intensity_mean_abs'].append(mhw_abs.mean())
        #mhw['intensity_var_abs'].append(np.sqrt(mhw_abs.var()))
        #mhw['intensity_cumulative_abs'].append(mhw_abs.sum())
        print('################ End of the heat event #', ev)              
print('Abs. intensity:', abs_intensity_array)
print('Rel. intensity:', rel_intensity_array)

###############################################################################
#### Plot stuff
###############################################################################
fig, ax = plt.subplots(figsize=(8,6))
# makes a grid on the background
ax.grid()
dailyMean_hw_celsius=dailyMean_hw[:,pixel_y,pixel_x]-273.15        
ax.plot(dailyMean_hw_celsius, color='black', label='daily '+str(varName), linewidth=2)

################### filtering the noise in daily climatological data ##########
#b, a = signal.butter(3, 0.1, btype='lowpass', analog=False) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
#low_passed = signal.filtfilt(b, a, noisy_signal)
##b, a = ellip(4, 0.01, 120, 0.125)  # Filter to be applied.
#clima_plot=signal.filtfilt(b, a, clima[:,pixel_y,pixel_x], padlen=30)
clima_plot_celsius=clima_plot-273.15
#clima_plot=savgol_filter(clima[:,pixel_y,pixel_x], 11, 3) #savitzky_golay(clima[:,pixel_y,pixel_x], 51, 3) # window size 51, polynomial order 3
ax.plot(clima_plot_celsius, color='green', label='daily clim. '+ str(varName) + ': ' + str(year_start_clim) +'-'+str(year_final_clim) +' (smoothed)',  linewidth=2)
ax.plot(clima[:,pixel_y,pixel_x]-273.15, color='green', linestyle='dashed', label='daily clim. '+ str(varName) + ': ' + str(year_start_clim) +'-'+str(year_final_clim) +' (original)',  linewidth=2)

################### 90th percentile data ##########
################### filtering the noise in daily 90th percentile data ##########
#b, a = signal.butter(3, 0.1, btype='lowpass', analog=False) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
#clima_90p_plot_temp=signal.filtfilt(b, a, clima_90p_full[:,pixel_y,pixel_x], padlen=30)
#clima_90p_plot_celsius_filt=clima_90p_plot_temp-273.15
#ax.plot(clima_90p_plot_celsius_filt, color='darkviolet', linewidth=2, label='S_int: daily 90th percentile: ' + str(year_start_clim)+'-'+str(year_final_clim)+ ' (smoothed)')

clima_90p_plot_celsius=clima_90p_full[:,pixel_y,pixel_x]-273.15
#clima_plot=savgol_filter(clima[:,pixel_y,pixel_x], 11, 3) #savitzky_golay(clima[:,pixel_y,pixel_x], 51, 3) # window size 51, polynomial order 3
ax.plot(clima_90p_plot_celsius, color='darkviolet', linewidth=2, label='daily 5-day window 90th percentile from ' + str(year_start_clim)+'-'+str(year_final_clim))

################### 97.5th percentile data ##########
#clima_975p_plot_celsius=clima_975p_full[:,pixel_y,pixel_x]-273.15
#clima_plot=savgol_filter(clima[:,pixel_y,pixel_x], 11, 3) #savitzky_golay(clima[:,pixel_y,pixel_x], 51, 3) # window size 51, polynomial order 3
#ax.plot(clima_975p_plot_celsius, color='black', linestyle='dotted', label='S_deb: 97.5th percentile: ' + str(year_start_clim)+'-'+str(year_final_clim))

################### 99.5th percentile data ##########
#clima_995p_plot_celsius=clima_995p_full[:,pixel_y,pixel_x]-273.15
#clima_plot=savgol_filter(clima[:,pixel_y,pixel_x], 11, 3) #savitzky_golay(clima[:,pixel_y,pixel_x], 51, 3) # window size 51, polynomial order 3
#ax.plot(clima_995p_plot_celsius, color='black', linestyle='dashdot', label='S_pic: 99.5th percentile: ' + str(year_start_clim)+'-'+str(year_final_clim))

##############################
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

labels_idx = np.arange(0,dailyTime_hw.shape[0], 15) #NoD) #15
ax.set_xticks(labels_idx)
labels = [cftime.datetime.strftime(item, '%Y-%m-%d') for item in dailyTime_hw[labels_idx] ]
ax.set_xticklabels(labels)
x = np.arange(dailyTime_hw.shape[0])
  
#ax.legend()  
#ax.legend(loc='lower left',fancybox=True, shadow=True) # ncol=4)
ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

plt.title('TSMP '+str(varName) +' [X=' + str(pixel_x) +', Y=' + str(pixel_y) +']' +' heat waves investigation: '+str(year_start_hw),fontsize=12,fontweight='bold')                
plt.ylabel('Daily mean temperatures (' + varName +'), K', fontsize=11)
plt.ylim(vmin, vmax)
#plt.xlabel('Days',fontsize=11)  

ax.fill_between(x, clima_plot_celsius, dailyMean_hw_celsius, where=dailyMean_hw_celsius>clima_plot_celsius, facecolor='lightcoral', interpolate=True, alpha=1.0) #lightsalmon
ax.fill_between(x, clima_plot_celsius, dailyMean_hw_celsius, where=dailyMean_hw_celsius<clima_plot_celsius, facecolor='deepskyblue', interpolate=True, alpha=1.0) 
ax.fill_between(x, clima_90p_plot_celsius, dailyMean_hw_celsius, where=dailyMean_hw_celsius>clima_90p_plot_celsius, facecolor='darkred', interpolate=True, alpha=1.0)

plt.savefig(f'Daily_{varName}_diff_daily.png', dpi=380) 
plt.show()



###############################################################################
#### Plot stuff +1
###############################################################################
fig, ax = plt.subplots(figsize=(8,3))
# makes a grid on the background
ax.grid()
ax.plot(abs_intensity_array, color='blue', label='Absolute intensity', linewidth=2)

labels_idx = np.arange(0,dailyTime_hw.shape[0], 15) #NoD) #15
ax.set_xticks(labels_idx)
labels = [cftime.datetime.strftime(item, '%Y-%m-%d') for item in dailyTime_hw[labels_idx] ]
ax.set_xticklabels(labels)
x = np.arange(dailyTime_hw.shape[0])
  
#ax.legend()  
#ax.legend(loc='lower left',fancybox=True, shadow=True) # ncol=4)
ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

plt.title('TSMP '+str(varName) +' [X=' + str(pixel_x) +', Y=' + str(pixel_y) +']' +' heat waves investigation: '+str(year_start_hw),fontsize=12,fontweight='bold')                
plt.ylabel('Abs. intensity of heat wave, K', fontsize=11)
plt.savefig(f'Absolute intensity_{varName}_diff_daily.png', dpi=380) 
#plt.show()


###############################################################################
#### Plot stuff +1
###############################################################################
fig, ax = plt.subplots(figsize=(8,3))
# makes a grid on the background
ax.grid()
ax.plot(rel_intensity_array, color='blue', linestyle = 'dashed', label='Relative intensity', linewidth=2)

labels_idx = np.arange(0,dailyTime_hw.shape[0], 15) #NoD) #15
ax.set_xticks(labels_idx)
labels = [cftime.datetime.strftime(item, '%Y-%m-%d') for item in dailyTime_hw[labels_idx] ]
ax.set_xticklabels(labels)
x = np.arange(dailyTime_hw.shape[0])
  
#ax.legend()  
#ax.legend(loc='lower left',fancybox=True, shadow=True) # ncol=4)
ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

plt.title('TSMP '+str(varName) +' [X=' + str(pixel_x) +', Y=' + str(pixel_y) +']' +' heat waves investigation: '+str(year_start_hw),fontsize=12,fontweight='bold')                
plt.ylabel('Rel. intensity of heat wave', fontsize=11)
plt.savefig(f'Relative intensity_{varName}_diff_daily.png', dpi=380) 
#plt.show()
