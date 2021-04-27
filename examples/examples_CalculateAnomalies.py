import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
import glob
import cftime
src_path='../src/'
sys.path.append(src_path)
import PlotLib 

###############################################################################
##### Read in climate means
###############################################################################
meanInterval = 'month'
varName = 'TSA'
NoI = 12
date_start_clim = cftime.datetime(1965, 6, 1, calendar='noleap')
date_final_clim = cftime.datetime(1986, 1, 1, calendar='noleap')
try:
    intervalMean = np.load(f'../data/example_ClimateMeans/intervalMean_{meanInterval}.npy')
    print(f'intervalMean.shape: {intervalMean.shape}')
    intervalTime = np.load(f'../data/example_ClimateMeans/intervalTime_{meanInterval}.npy', allow_pickle=True)
    print(f'intervalTime.dtype: {intervalTime.dtype}')
    clima = np.load(f'../data/example_ClimateMeans/climate_{meanInterval}_{date_start_clim.strftime("%Y_%m")}_{date_final_clim.strftime("%Y_%m")}.npy')
except FileNotFoundError:
    print(f'needed climateMeans files were not found: EXIT')
    sys.exit()

###############################################################################
##### CALCULATE ANOMALIES
###############################################################################
try:
    dailyAnomalyDomain = np.load(f'IntervalAnomalies_TestDump.npy')
    print(f'dailyAnomalyDomain.shape: {dailyAnomalyDomain.shape}')
    print(f'dailyAnomalyDomain.dtype: {dailyAnomalyDomain.dtype}')
except FileNotFoundError:
    intervalAnomalyDomain = np.empty(intervalMean.shape[0])
    for n, interval in enumerate(intervalMean):
        idx_interval = n%NoI
        intervalAnomalyDomain[n] = np.nanmean(interval, dtype=float) - np.nanmean(clima[idx_interval], dtype=float)
    # dump anomalies
    with open(f'IntervalyAnomalies_{meanInterval}.npy', 'wb') as f:
        np.save(f, intervalAnomalyDomain)

###############################################################################
#### Plot stuff
###############################################################################
fig, ax = plt.subplots(figsize=(8,4))
# makes a grid on the background
ax.grid()
ax.plot(intervalAnomalyDomain, color='black')

labels_idx = np.arange(0,intervalTime.shape[0], NoI)
ax.set_xticks(labels_idx)
labels = [ cftime.datetime.strftime(item, '%Y') for item in intervalTime[labels_idx] ]
ax.set_xticklabels(labels)

x = np.arange(intervalAnomalyDomain.shape[0])
ax.fill_between(x, 0, intervalAnomalyDomain, where=intervalAnomalyDomain>0, facecolor='red', interpolate=True, alpha=0.75) 
ax.fill_between(x, 0, intervalAnomalyDomain, where=intervalAnomalyDomain<0, facecolor='blue', interpolate=True, alpha=0.75)
plt.savefig(f'./examples_CalculateAnomalies.pdf') 