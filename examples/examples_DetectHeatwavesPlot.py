""" Example script showing how to detect HeatWaves

[ADD SOME DESCRIPTIN HERE]
"""
import numpy as np
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
rootDateDir       = '/p/project/cesmtst/poshyvailo1/SLOTH/data'
#projectName       = 'IPSL-IPSL-CM5A-LR_REMO2015'
projectName       = 'TSMP_MPI-ESM-LR'
refClimatology    = f'{rootDateDir}/example_ClimateMeans/climate_day_{projectName}_1976_2005.npy'
HwFileNamePattern = f'HWevents_????_{projectName}.npy'
EvFileNamePattern = f'Events_????_{projectName}.npy'

HwStatsFiles = sorted(glob.glob(f'{rootDateDir}/example_HWevents/{projectName}/{HwFileNamePattern}'))
EvStatsFiles = sorted(glob.glob(f'{rootDateDir}/example_HWevents/{projectName}/{EvFileNamePattern}'))
#EvStatsFiles = sorted(glob.glob(f'{rootDateDir}/example_HWevents/{projectName}/{EvFileNamePattern}'))
#print(f'HWstatsFiles: {HWstatsFiles}')

###############################################################################
##### Read in climatology calculated with 'examples_CalculateClimateMeans.py
##### those are only needed for 'landPixeld' --> how to better solve this?
###############################################################################
try:
    clima = np.load(refClimatology)

    # mask holding 1 where land
    mask = np.zeros_like(clima[0])
    mask[clima[0]!=np.nan] = 1
    landPixels = np.sum(mask)
except FileNotFoundError:
    print(f'ERROR: not all needed files were not found: EXIT')
    sys.exit()

###############################################################################
##### Read in HwStats created with 'examples_DetectHeatwaves_Domain.py'
###############################################################################
HWeventDurations = []
for HwStatsFile in HwStatsFiles:
    print(f'HwStatsFile: {HwStatsFile}')
    # HW were stored as python-dicts. Therefore it is needed to read those in
    # with the flag 'allow_pickle=True' 
    HwStats = np.load(HwStatsFile, allow_pickle=True)
    # Further 'np.load()' does return a numpy ndarray, that we need to extract 
    # the item to get a dict again.
    HwStats = HwStats.item()
    
    tmp_HWeventDurations = [HwStats[key]['duration'] for key in HwStats.keys()]
    HWeventDurations += tmp_HWeventDurations


###############################################################################
##### Read in all events created with 'examples_DetectHeatwaves_Domain.py'
###############################################################################
events_list = []
eventsDuration_list = []
for EvStatsFile in EvStatsFiles:
    print(f'EvStatsFile: {EvStatsFile}')
    EvStats = np.load(EvStatsFile)
    events_list.append(EvStats)
    
    # EvStats does contain the 3D information if pixel does exceed the 90p
    # TSA value. Further those are labeld, that connected pixel share the same
    # label. This way 'counting' (np.unique(EvStats, return_counts=True)) does
    # return the length / duration of each individual event
    tmp_eventsLabel, tmp_eventsDuration = np.unique(EvStats, return_counts=True)
    # remove (slice out) label 0, as is is 'no HW'
    tmp_eventsDuration = tmp_eventsDuration[1:]
    print(f'type(tmp_eventsDuration): {type(tmp_eventsDuration)}')
    eventsDuration_list.append(tmp_eventsDuration)
    print(f'EvStats.shape: {EvStats.shape}')
events         = np.concatenate(events_list, axis=0)
print(f'events.shape: {events.shape}')
eventsDuration = np.concatenate(eventsDuration_list, axis=0)
print(f'eventsDuration: \n{eventsDuration}')

###############################################################################
#### Plot stuff Fig.5 Vautard
###############################################################################
# set array of individual durations in days
days = np.arange(1,26)
eventsA = np.array([ np.sum(eventsDuration==tmpThreshold) for tmpThreshold in days] )
# normalize by 'land-pixel'
eventsA_mean = eventsA / landPixels
eventsB = np.array([ np.sum(eventsDuration>=tmpThreshold) for tmpThreshold in days] )
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
# Count non-zero values in 3D events. This way one gets the accumulated
# hot-days for each pixel over time.
accumulated_events = np.count_nonzero(events[:92], axis=0)
fig, ax = plt.subplots(figsize=(9,5))
ax.set_title(f'{projectName} Accumulated hot days > 90p TSA prec. for JJA 1979')
img = ax.imshow(accumulated_events, origin='lower', cmap='Reds', interpolation=None)
fig.colorbar(img)
fig.savefig(f'p5_{projectName}.pdf')


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
