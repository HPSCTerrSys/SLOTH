import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
import glob
import cftime


###############################################################################
##### Read in climate means
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
meanInterval = 'month'
NoI = 12

try:
    intervalMean = np.load(f'../data/example_ClimateMeans/intervalMean_{meanInterval}.npy')
    print(f'intervalMean.shape: {intervalMean.shape}')
    intervalTime = np.load(f'../data/example_ClimateMeans/intervalTime_{meanInterval}.npy', allow_pickle=True)
    print(f'intervalTime.dtype: {intervalTime.dtype}')
    clima = np.load(f'../data/example_ClimateMeans/climate_{meanInterval}.npy')
except FileNotFoundError:
    print(f'needed climateMeans files were not found: EXIT')
    sys.exit()

###############################################################################
##### CALCULATE ANOMALIES
###############################################################################
try:
    dailyAnomalyDomain = np.load(f'../data/example_ClimateMeans/IntervalyAnomalies_{meanInterval}.npy')
    print(f'dailyAnomalyDomain.shape: {dailyAnomalyDomain.shape}')
    print(f'dailyAnomalyDomain.dtype: {dailyAnomalyDomain.dtype}')
except FileNotFoundError:
    intervalAnomalyDomain = np.empty(intervalMean.shape[0])
    for n, interval in enumerate(intervalMean):
        idx_interval = n%NoI
        intervalAnomalyDomain[n] = np.nanmean(interval, dtype=float) - np.nanmean(clima[idx_interval], dtype=float)
    # dump anomalies
    with open(f'../data/example_ClimateMeans/IntervalyAnomalies_{meanInterval}.npy', 'wb') as f:
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