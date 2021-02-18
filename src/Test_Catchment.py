import numpy as np
import sys
import os
import glob

import VAlidationTool as vat
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc

import time

import sys
sys.path.append("../extern/transform_rotated_pole")
import rotate
sys.path.append("../extern/parflow")
from ParFlow_IO import read_pfb


pfb_file_path = '/home/n.wagner/development/VAlidationTool_data/LPOCatchment'

slopex = read_pfb(f'{pfb_file_path}/slopex.pfb')
slopey = read_pfb(f'{pfb_file_path}/slopey.pfb')

# 329, 166 are taken from BIBIS GRDC dataset for ID 6842900
# 163	326
# 241	245 for ID 6457010

x = 326
y = 163
rad_x = 10
rad_y = 15

"""
start = time.time()
cmask_new = vat.toolBox.calc_catchment(slopex[0], slopey[0], x, y)
end = time.time()
print(f'vat.toolBox.calc_catchment: {end - start}')
start = time.time()
sys.exit()
cmask_old_majorslope = vat.toolBox.calc_catchment_old_majorslope(slopex[0], slopey[0], x, y)
end = time.time()
print(f'vat.toolBox.calc_catchment_old_majorslope: {end - start}')
start = time.time()
cmask_old_nozeroslope = vat.toolBox.calc_catchment_old_nonzeroslope(slopex[0], slopey[0], x, y)
end = time.time()
print(f'vat.toolBox.calc_catchment_old_nonzeroslope: {end - start}')
#cmask_old = np.load(f'{pfb_file_path}/catchment_mask_6842900.npy')
sys.exit()
"""
cmask = np.load('catchment_mask.npy')
cmask_old_majorslope = np.load('catchment_mask_old_majorslope.npy')
cmask_old_nonslope = np.load('catchment_mask_old_nonzeroslope.npy')

print(f'cmask==cmask_old_majorslope: {(cmask==cmask_old_majorslope).all()}')
print(f'cmask==cmask_old_nonslope: {(cmask==cmask_old_nonslope).all()}')
print(f'cmask==cmask: {(cmask==cmask).all()}')
sys.exit()
fig = plt.figure()
gs1 = fig.add_gridspec(nrows=2, ncols=2,
                   top=1.0, bottom=0.1)
axes = []
# index 0 ist am oberrand und index -1 ist am unterrand
axes.append(fig.add_subplot(gs1[0,0]))
axes.append(fig.add_subplot(gs1[0,1]))
axes.append(fig.add_subplot(gs1[1,0]))
img_cmaks = axes[0].imshow(cmask)
img_cmask_old_majorslope = axes[1].imshow(cmask_old_majorslope)
img_cmask_old_nonslope = axes[2].imshow(cmask_old_nonslope)
plt.show()

sys.exit()
fig = plt.figure()
gs1 = fig.add_gridspec(nrows=3, ncols=2,
                   top=1.0, bottom=0.1)
axes = []
# index 0 ist am oberrand und index -1 ist am unterrand
axes.append(fig.add_subplot(gs1[0,0]))
axes.append(fig.add_subplot(gs1[0,1]))
axes.append(fig.add_subplot(gs1[1,0]))
#axes.append(fig.add_subplot(gs1[1,1]))
axes.append('dummy')
axes.append(fig.add_subplot(gs1[2,0]))
axes.append(fig.add_subplot(gs1[2,1]))

slopex_bin = np.zeros_like(slopex[-1])
slopex_bin = np.where(slopex<0, -1, slopex_bin)
slopex_bin = np.where(slopex>0, 1, slopex_bin)
slopey_bin = np.zeros_like(slopey[-1])
slopey_bin = np.where(slopey<0, -1, slopey_bin)
slopey_bin = np.where(slopey>0, 1, slopey_bin)

cnew = axes[0].imshow(cmask_new[y-rad_y:y+rad_y, x-rad_x:x+rad_x], aspect='auto')
axes[0].scatter(rad_x, rad_y, c='green', marker='x', label='start')
plt.colorbar(cnew, ax=axes[0])
cold = axes[1].imshow(cmask_old[y-rad_y:y+rad_y, x-rad_x:x+rad_x], aspect='auto')
axes[1].scatter(rad_x, rad_y, c='green', marker='x', label='start')
plt.colorbar(cold, ax=axes[1])
cdif = axes[2].imshow((cmask_new-cmask_old)[y-rad_y:y+rad_y, x-rad_x:x+rad_x], aspect='auto')
axes[2].scatter(rad_x, rad_y, c='green', marker='x', label='start')
plt.colorbar(cdif, ax=axes[2])
slox = axes[4].imshow(slopex_bin[0, y-rad_y:y+rad_y, x-rad_x:x+rad_x], aspect='auto', vmin=-1, vmax=1, cmap=plt.cm.bwr)
axes[4].scatter(rad_x, rad_y, c='green', marker='x', label='start')
plt.colorbar(slox, ax=axes[4])
sloy = axes[5].imshow(slopey_bin[0, y-rad_y:y+rad_y, x-rad_x:x+rad_x], aspect='auto', vmin=-1, vmax=1, cmap=plt.cm.bwr)
axes[5].scatter(rad_x, rad_y, c='green', marker='x', label='start')
plt.colorbar(sloy, ax=axes[5])
plt.show()
# print( (cmask_old == cmask_new).all() )
# vat.toolBox.calc_catchment_old(slopex[0], slopey[0], x, y)
sys.exit()



