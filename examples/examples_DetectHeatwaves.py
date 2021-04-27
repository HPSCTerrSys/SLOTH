import numpy as np
import sys
import os
import cftime
from scipy import signal
import scipy.ndimage as ndimage

src_path='../src/'
sys.path.append(src_path)
import PlotLib 

###############################################################################
#### General settings
###############################################################################
# set start and end date of climatology (already calculated)
date_start_clim = cftime.datetime(1965, 6, 1, calendar='noleap')
date_final_clim = cftime.datetime(1986, 1, 1, calendar='noleap')

# Year to investigate for HeatWaves (HW)
# currently only one year possible!
hwYears = [1973]#, 1973]

# HW investigation performed with this script is based on percentile 
# calculation. To 'smooth' the calculation the percentiles are calculated not
# for the time-series of a single day, but as kind of running-percentile for
# a 'window' of X days
window_percentile = 5

# A single hot day crossing the 90p threshold is not a HW, minDuration does
# define how many 'too hot days' in a row building a HW
minDuration = 6

# A singel pixel for which to perform the HW analysis 
# --> There should be an easy way to extend this to entire domain...
pixel_x=250
pixel_y=300

# why needed?
varName = 'TSA'

# Prepare a dict holding some stats of the HW 
mhw_def = {}
mhw_def['time_start'] = [] # cftime format; Start time of MHW [cftime format]
mhw_def['time_end'] = [] # cftime format;  End time of MHW [cftime format]
mhw_def['date_start'] = [] # cftime format
mhw_def['date_end'] = [] # cftime format
mhw_def['index_start'] = []
mhw_def['index_end'] = []
mhw_def['index_peak'] = []
mhw_def['date_peak'] = [] # cftime format
mhw_def['duration'] = [] # [days]
mhw_def['max_temp'] = [] # [days]
mhw_def['intensity_max'] = [] # [deg C]
mhw_def['intensity_mean'] = [] # [deg C]
#mhw_def['intensity_var'] = [] # [deg C]
mhw_def['intensity_cumulative'] = [] # [deg C]
mhw_def['intensity_max_relThresh'] = [] # [deg C]
mhw_def['intensity_mean_relThresh'] = [] # [deg C]
#mhw_def['intensity_var_relThresh'] = [] # [deg C]
mhw_def['intensity_cumulative_relThresh'] = [] # [deg C]
mhw_def['intensity_max_abs'] = [] # [deg C]
mhw_def['intensity_mean_abs'] = [] # [deg C]
#mhw_def['intensity_var_abs'] = [] # [deg C]
mhw_def['intensity_cumulative_abs'] = [] # [deg C]
#mhw_def['category'] = []
#mhw_def['rate_onset'] = [] # [deg C / day]
#mhw_def['rate_decline'] = [] # [deg C / day]

###############################################################################
##### Read in climatology calculated with 'examples_CalculateClimateMeans.py
###############################################################################
# The climatology calculate in 'examples_CalculateClimateMeans.py' is able to 
# generate different climatologies based on different mean-intervals (daily, 
# monthly). For HW investigation we always need daily mean values, wherefore 
# 'meanInterval' is not optional but set for consistency.
meanInterval = 'day'
# Similar for the Number Of Intervals (NoI) contained with the climatology 
# loaded.
# For full years and 'meanInteral' = month --> NoI 12 (12 month per year)
# For full years and 'meanInteral' = day   --> NoI 365 (365 days per year)
# However, one could also calculate the climatology for summer-months only 
# (JJA) in which case the number of intervals would reduce to NoI=3 or NoI=92
# In this case NoI is not optional and have to be set NoI=365 for a full year
# of daily means
NoI = 365
try:
    dailyMean = np.load(f'../data/example_ClimateMeans/intervalMean_{meanInterval}.npy')
    print(f'dailyMean.shape: {dailyMean.shape}')
    dailyTime = np.load(f'../data/example_ClimateMeans/intervalTime_{meanInterval}.npy', allow_pickle=True)
    print(f'dailyTime.dtype: {dailyTime.dtype}')
    clima = np.load(f'../data/example_ClimateMeans/climate_{meanInterval}_{date_start_clim.strftime("%Y_%m")}_{date_final_clim.strftime("%Y_%m")}.npy')
except FileNotFoundError:
    print(f'needed climateMeans files were not found: EXIT')
    sys.exit()


###############################################################################
#### preparing HeatWave investigation
###############################################################################
# HeatWave investigation performed in this script is based on the 90p 
# percentile of the climatology. 
# Create array holding the 90p percentiles, which shape is [365, Y, X] for 
# our case (NoI=365)
climaDim_90p = [NoI]
climaDim_90p = climaDim_90p + [dim for dim in dailyMean[0].shape]
clima_90p = np.empty(climaDim_90p)
print(f'intervalClima_90p.shape: {clima_90p.shape}')

# To 'smooth' the 90p percentile, a window around each day is taken into 
# account, that percentiles are not calculated out of the time-series of a 
# single day, but as kind of running-percentile for a 'window' of X days.
# A 'window' of e.g. 5 is actual a day increment by +- 2 days where incWindow 
# is calculated as incWindow = window_percentile//2 (integer devision)
incWindow = window_percentile // 2
for curr_day in range(NoI):
    print(f'handling day: {curr_day}')
    # tmp variable holding all data within the 'window'
    tmpWindow = None
    # Gathering the window.
    # +1 because python does exclude the last step
    for i in range(curr_day-incWindow, curr_day+incWindow+1):
        # skip if i is out of range for dailyMean
        if i < 0 or i > dailyMean.shape[0]:
            continue
        if tmpWindow is None:
            tmpWindow = dailyMean[i::NoI]
        else:
            tmpWindow = np.append(tmpWindow, dailyMean[i::NoI], axis=0)
    # Calculating the 'windowed' 90p percentile for each day a year. 
    clima_90p[curr_day]=np.percentile(tmpWindow, 90, axis=0)
    # print('highest 90p percentile over entire domain: ',np.nanmax(clima_90p[curr_day]))

###############################################################################
#### HEAT WAVE indices calculation
#### see also: https://github.com/ecjoliver/marineHeatWaves/blob/master/marineHeatWaves.py
###############################################################################
# A HW is detected if X days in a row exceeding the 90p threshold, wherefore X 
# is minDuration.
# Each year passed to investigate, is treating individually
for hwYear in hwYears:
    mhw = mhw_def
    ###############################################################################
    #### Extract year of interests for HW analysis
    ###############################################################################
    # Extract the year / period of interest to analysis for HWs
    boolSlice = [ hwYear == item.year for item in dailyTime]
    dailyMean_hw = dailyMean[boolSlice]
    dailyTime_hw = dailyTime[boolSlice]

    # Calculate where daily mean temperature exceeded the 90p threshold
    exceed_bool = dailyMean_hw - clima_90p
    exceed_bool[exceed_bool<=0] = False
    exceed_bool[exceed_bool>0] = True
    # Fix issue where missing temp values (nan) are counted as True
    exceed_bool[np.isnan(exceed_bool)] = False
    # Investigate HW for given pixel only
    # --> There HAST TO BE a way to calculate this for the entire domain easily... 
    exceed_bool = exceed_bool[:,pixel_y,pixel_x]

    # Find contiguous regions of exceed_bool = True
    events, n_events = ndimage.label(exceed_bool)
    print('Number of events (groups of days) exceeding 90th percentile:', n_events)

    # Find all MHW events of duration >= minDuration
    for ev in range(1,n_events+1):
            print("yes")
            event_duration = (events == ev).sum()
            if event_duration < minDuration:
                continue          
            mhw['time_start'].append(dailyTime_hw[np.where(events == ev)[0][0]])
            mhw['time_end'].append(dailyTime_hw[np.where(events == ev)[0][-1]])
    mhw['n_events'] = len(mhw['time_start'])
    #categories = np.array(['Moderate', 'Strong', 'Severe', 'Extreme'])
    print('Numbers of the heat-events exceeding min. length of ', str(minDuration), 'days: ', mhw['n_events'], ' event')
     
    ###############################################################################
    #### Calculate some HW properties as intensity etc.
    ###############################################################################
    # First filtering / smooth the noise in daily climatological data
    # --> see also: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
    b, a = signal.butter(3, 0.1, btype='lowpass', analog=False) 
    clima_smooth=signal.filtfilt(b, a, clima[:,pixel_y,pixel_x], padlen=15)


    abs_intensity_array= np.zeros(NoI) #initialize empty array
    rel_intensity_array= np.zeros(NoI) #initialize empty array
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
            thresh_mhw = clima_90p[tt_start:tt_end+1,pixel_y,pixel_x]
            seas_mhw = clima_smooth[tt_start:tt_end+1]        
            
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

    ###############################################################################
    #### Plot
    ###############################################################################
    vmin = 10
    vmax = 24
    PlotLib.plot_HeatWaveInvest(abs_intensity_array=abs_intensity_array, dailyTime_hw=dailyTime_hw, 
        varName=varName, pixel_x=pixel_x, pixel_y=pixel_y, hwYear=hwYear,
        rel_intensity_array=rel_intensity_array, dailyMean_hw=dailyMean_hw, clima_smooth=clima_smooth, 
        date_start_clim=date_start_clim, date_final_clim=date_final_clim, clima_90p=clima_90p,
        clima=clima, vmin=vmin, vmax=vmax)