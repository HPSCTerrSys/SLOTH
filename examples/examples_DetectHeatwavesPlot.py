""" Example script showing how to detect HeatWaves

[ADD SOME DESCRIPTIN HERE]
"""
import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import glob
import os
import cftime

sloth_path='../'
sys.path.append(sloth_path)
import sloth

###############################################################################
#### Define some paths, filenames, options, etc
###############################################################################
#rootDateDir       = '/p/project/cesmtst/poshyvailo1/SLOTH/data'
rootDateDir       = '../data'
#projectName       = 'IPSL-IPSL-CM5A-LR_REMO2015'
projectName       = 'TSMP_ERA5'
EvFileNamePattern = f'Events_????_{projectName}.nc'
EvFiles = sorted(glob.glob(f'{rootDateDir}/example_HWevents/{projectName}/{EvFileNamePattern}'))

###############################################################################
##### Read in all events -- 'examples_DetectHeatwaves_Domain.py'
###############################################################################
accumulated_events = []
days = np.arange(1,26)
eventsA = []
eventsB = []
fig6_amplitude = np.arange(10)
fig6_frequency = []
for EvFile in EvFiles:
    print(f'EvFile: {EvFile}')
    with nc.Dataset(EvFile, 'r') as nc_file:
        events = nc_file.variables['events'][...]
        # change dtype to int
        events = events.astype(int)
        Nevents = np.max(events)
        # MASK 'EVENTS' HERE WITH LAND SEA MASK OR WHATEVER YOU WANT
        # EXAMPLE BELOW for a LandSeaMask (LSM) holding zeros for sea
        # and ones for land pixel:
        # events = np.ma.masked_where(LSM==0, events)

        # 'events' does contain the 3D information if pixel does exceed the 90p
        # threshold. Further those pixels are labeled, that connected pixel 
        # share the same value. This way 'counting' 
        # (np.unique(EvStats, return_counts=True)) does return the 
        # length / duration of each individual event
        eventsLabel, eventsDuration = np.unique(events, return_counts=True)

        # Read in 1D variables for detailed information on HeatWaves which are
        # heat-events with a duration >= minDuration
        intensity_max_relThresh = nc_file.variables['intensity_max_relThresh'][...]
        print(f'intensity_max_relThresh.shape: {intensity_max_relThresh.shape:}')
        print(f'Nevents: {Nevents}')
        # Mask same values for 1D variables which were masked for 'events' (e.g.
        # according to a LandSeaMask).
        mask1D = np.full(Nevents, True)
        mask1D[eventsLabel] = False
        intensity_max_relThresh = np.ma.masked_where(mask1D, intensity_max_relThresh)
        # You can do this for any other 1D variable as well, as the mask1D
        # stays the same. For example:
        # Further1DVariable = np.ma.masked_where(mask1D, Further1DVariable)

        # Calculate / collect the two variables (eventsA and eventsB) needed
        # for 'Vautard Fig 5.'
        eventsA.append([ np.sum(eventsDuration==tmpThreshold) for tmpThreshold in days] )
        eventsB.append([ np.sum(eventsDuration>=tmpThreshold) for tmpThreshold in days] )

        # Count nonzero values / labels along axis=0, the time axis.
        # Those pixel does exceed the 90p threshold and can be used to plot
        # a 'number of hot days' map.
        accumulated_events.append(np.count_nonzero(events, axis=0))
    
        # Calculate / collect variables needed for 'Vautard Fig 6.'
        fig6_frequency.append(np.array([ np.sum(intensity_max_relThresh>=threshold) for threshold in fig6_amplitude ]))

eventsA = np.stack(eventsA)
eventsA = np.sum(eventsA, axis=0)
eventsB = np.stack(eventsB)
eventsB = np.sum(eventsB, axis=0)
accumulated_events = np.stack(accumulated_events)
accumulated_events = np.sum(accumulated_events, axis=0)
fig6_frequency = np.stack(fig6_frequency)
fig6_frequency = np.sum(fig6_frequency, axis=0)
# frequency: NumberOfEventsExceedingAmplitude / TotalNumberOfEvents 
fig6_frequency = fig6_frequency / fig6_frequency[0]
###############################################################################
#### Plot stuff Fig.5 Vautard
###############################################################################
# hardcoded value for landPixels, until this is available via events.nc
landPixels = 90000
# normalize by 'land-pixel'
eventsA_mean = eventsA / landPixels
# normalize by 'land-pixel'
eventsB_mean = eventsB / landPixels

fig, ax = plt.subplots(figsize=(9,5))
ax.grid()
ax.plot(days,eventsA_mean, ls='--', label='TSMP heatwaves (= duration from X-axis)', linewidth=2)
ax.plot(days,eventsB_mean, ls='--', label='TSMP heatwaves (>=duration from X-axis)', linewidth=2)
ax.legend(loc='upper center', fancybox=True, shadow=True)
plt.title(f'{projectName} heat waves investigation',fontsize=12,fontweight='bold')
plt.xlabel('Duration [days]', fontsize=11)
plt.ylabel('Mean numbers of events', fontsize=11)
plt.savefig(f'Test_{projectName}.pdf')

###############################################################################
#### Plot stuff on p.5 LPO
###############################################################################
fig, ax = plt.subplots(figsize=(9,5))
ax.set_title(f'{projectName} Accumulated hot days > 90p TSA prec. for JJA 1979')
img = ax.imshow(accumulated_events, origin='lower', cmap='Reds', interpolation=None)
fig.colorbar(img)
fig.savefig(f'p5_{projectName}.pdf')


###############################################################################
#### Plot Fig.6 Vautard
###############################################################################
fig, ax = plt.subplots(figsize=(9,5))
ax.grid()
ax.plot(fig6_amplitude,fig6_frequency, ls='--', label='TSMP heatwaves with duration > 6 days & >= amplitude from X-axis', linewidth=2)
ax.legend(loc='upper center', fancybox=True, shadow=True)
plt.title(f'{projectName} heat waves investigation',fontsize=12,fontweight='bold')
plt.xlabel('Amplitude (T exceeds 90p) [K]', fontsize=11)
plt.ylabel('Frequency of events', fontsize=11)
plt.savefig(f'p6_{projectName}.pdf')

"""    
###############################################################################
#### Plot stuff Fig.5 Vautard
###############################################################################
# set array of individual durations in days
days = np.arange(1,16)
# Find and sum up all events equal to some ref-value, what is equal of finding
# and counting events equal to a given duration.
# (--1--) >> np.sum(ALL_event_label_counts==tmpThreshold)
# Apply this to all values / durations / days and save as ndarray
# (--2--) >> np.array([ (--1--) for tmpThreshold in days] )
HeventsA = np.array([ np.sum(ALL_event_label_counts==tmpThreshold) for tmpThreshold in days] )
# normalize by 'land-pixel'
HeventsA_mean = HeventsA / landPixels
HeventsB = np.array([ np.sum(ALL_event_label_counts>=tmpThreshold) for tmpThreshold in days] )
HeventsB_mean = HeventsB / landPixels

fig, ax = plt.subplots(figsize=(9,5))
ax.grid()
ax.plot(days,HeventsA_mean, color='blue', label='TSMP heatwaves '+ varName + '(= duration from X-axis)', linewidth=2)  #varName
ax.plot(days,HeventsB_mean, color='red', label='TSMP heatwaves '+ varName + '(>=duration from X-axis)', linewidth=2)  #varName
#ax.plot(days, n_event_modified, color='red', label='TSMP heatwaves (>90th percentile of '+ varName +')', linewidth=2)  #varName, L504='>90p'
#ax.set_xticks(labels_idx)
ax.legend(loc='upper center', fancybox=True, shadow=True)
plt.title('TSMP '+str(varName) + ' heat waves investigation',fontsize=12,fontweight='bold')                
plt.xlabel('Duration [days]', fontsize=11)
plt.ylabel('Mean numbers of events', fontsize=11)
plt.savefig(f'Mean_numbers_heatwaves_{varName}.pdf') 
# plt.show()

###############################################################################
#### Plot stuff Fig.6a Vautard, 2013

###############################################################################
# set array of individual amplitudes 
amplitude = np.arange(0,9)
# get labels of all events exceeding a given minDuration. As described above
# the length of this array is equal to the number of events exceeding the 
# minDuration.
maxEvent  = ALL_event_labels[ALL_event_label_counts>=minDuration].shape[0]
frequencyGTamplitude = np.array([ len({k:v for k,v in HWevents.items() if v['duration'] >= minDuration and v['intensity_max_relThresh'] >= tmpThreshold}) for tmpThreshold in amplitude], dtype=float)
frequencyGTamplitude *= 1./maxEvent
fig, ax = plt.subplots(figsize=(9,5))
ax.grid()
ax.plot(amplitude,frequencyGTamplitude, color='red', label='TSMP heatwaves with duration >=6 days & >= amplitude from X-axis', linewidth=2)  #varName
ax.legend(loc='upper right', fancybox=True, shadow=True,  ncol=4)
plt.title('TSMP '+str(varName) + ' heat waves investigation',fontsize=12,fontweight='bold')                
plt.xlabel('Amplitude (T excess of 90th percentile) [K]', fontsize=11)
plt.ylabel('Frequency of events', fontsize=11)
plt.savefig(f'Mean_frequency_amplitude_{varName}.pdf') 
# plt.show()

# ###############################################################################
# #### Plot stuff Fig.7 Vautard, 2013
# ###############################################################################
# fig, ax = plt.subplots(figsize=(9,5))
# ax.grid()
# ax.plot(tot_year,freq_tot_hw_days, color='red', label='TSMP mean hot days (exceeding 90th percentile)', linewidth=2)  #varName
# ax.legend(loc='upper left', fancybox=True, shadow=True)
# plt.title('TSMP '+str(varName) + ' heat waves investigation: JJA '+str(year_start_hw) + '-' +str(year_final_hw),fontsize=12,fontweight='bold')                
# plt.xlabel('Year', fontsize=11)
# plt.ylabel('Hot days frquency', fontsize=11)
# plt.savefig(f'Mean_frequency_hot_days_{varName}.png', dpi=380) 
# plt.show()
# exit()
"""
