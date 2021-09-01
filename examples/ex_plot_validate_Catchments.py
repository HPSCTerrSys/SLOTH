import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import datetime

import sys
import glob

import sloth
import sloth.ParFlow_IO as pio

def calc_NSE(Qm, Qo, norm=False):
    Qom = Qo.mean()
    NSE = 1 - (np.abs(Qm-Qo)).sum() / (np.abs(Qo-Qom)).sum()
    #NSE = 1 - ((Qm-Qo)**2).sum() / ((Qo-Qom)**2).sum()
    if norm:
        NSE = 1. / (2.-NSE)

    return NSE

PFL_data  = np.load('./PFL_scatterValidation_1981-2010.npy')
GRDC_data = np.load('./GRDC_scatterValidation_1981-2010.npy')
NSE = calc_NSE(Qm=PFL_data, Qo=GRDC_data)
axmaxlim = max(np.max(PFL_data), np.max(GRDC_data))
xValues = np.arange(axmaxlim)
mainDiagValues  = 1.  * xValues + 0
lowerDiagValues = 0.8 * xValues + 0
upperDiagValues = 1.2 * xValues + 0

fig   = plt.figure(facecolor='w', edgecolor='k')
fig.set_size_inches(8.27,11.69/2.) # A4 / 2

ax = fig.add_axes([0.1,0.1,0.85,0.85])
ax.set_aspect(1)

ax.set_title(f'Mittlerer monatlicher Abfluss (1981 - 2010)')

ax.scatter(GRDC_data, PFL_data, s=2, color='blue')
ax.plot(mainDiagValues, color='black', linewidth=0.9)
ax.plot(lowerDiagValues, color='gray', linewidth=0.9, linestyle='--')
ax.plot(upperDiagValues, color='gray', linewidth=0.9, linestyle='--')

ax.set_xlabel(r'GRDC [$m^{3} s^{-1}$]')
ax.set_xlim([1,axmaxlim])
ax.set_xscale('log')
ax.set_ylabel(r'PFL [$m^{3} s^{-1}$]')
ax.set_ylim([1,axmaxlim])
ax.set_yscale('log')

infoStr = [f'number of Stations: {GRDC_data.shape[0]}',
        f'NSE: {NSE:0.2f}']
infoStr = '\n'.join(infoStr)
ax.text(0.99, 0.01, infoStr, transform=ax.transAxes, 
        verticalalignment='bottom', horizontalalignment='right')

idealPerformance = Line2D([0], [0], label='ideal performance', color='black')
error20p         = Line2D([0], [0], label='20% error', linestyle='--', color='gray')
discharge        = Line2D(range(1), range(1), label='discharge values',
        color="white", marker='o', markerfacecolor="blue")
ax.legend(handles=[idealPerformance, error20p, discharge])

fig.savefig(f'./ScatterValidation.png')
