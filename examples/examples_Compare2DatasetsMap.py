""" Example script to show how to use different methods of this repo.

Test-Case: compare tow datasets, as e.g. modifed run vs referenz run.
"""
import numpy as np
import matplotlib as mpl
import netCDF4 as nc
import cftime
import sys
import glob
import copy

src_path='../src/'
sys.path.append(src_path)
import PlotLib 


###############################################################################
### Define some paths, filenames, etc
###############################################################################

dataRootDir = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName_D1 = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
datasetName_D2 = 'HiCam-CORDEX_EUR-11_MPI-ESM-LR_histo_r1i1p1_FZJ-IBG3'
procType    = 'postpro'
varName     = 'satur'
fileName    = f'{varName}.nc'

files_D1 = [
        f'{dataRootDir}/{datasetName_D1}/{procType}/1979_09/{fileName}',
        f'{dataRootDir}/{datasetName_D1}/{procType}/1979_10/{fileName}',
        f'{dataRootDir}/{datasetName_D1}/{procType}/1979_11/{fileName}',
        ]
files_D2 = [
        f'{dataRootDir}/{datasetName_D2}/{procType}/1979_09/{fileName}',
        f'{dataRootDir}/{datasetName_D2}/{procType}/1979_10/{fileName}',
        f'{dataRootDir}/{datasetName_D2}/{procType}/1979_11/{fileName}',
        ]


###############################################################################
### read in data and save as ndarray
###############################################################################
print(f'read in data and save as ndarray')
for idx in range(len(files_D1)):
    print(f'idx: {idx}')
    with nc.Dataset(f'{files_D1[idx]}', 'r') as nc_file_D1, nc.Dataset(f'{files_D2[idx]}', 'r') as nc_file_D2:
        tmp_var1  = nc_file_D1.variables[varName][:,-1,...]
        tmp_mask1 = tmp_var1.mask
        if not tmp_mask1.any():
            tmp_mask1  = np.zeros(tmp_var1.shape, dtype=bool)
            tmp_var1   = tmp_var1.filled(fill_value=np.nan)
            tmp_var1[tmp_var1 <= -1e+38] = np.nan
        tmp_var2  = nc_file_D2.variables[varName][:,-1,...]
        tmp_mask2 = tmp_var2.mask
        if not tmp_mask2.any():
            tmp_mask2  = np.zeros(tmp_var2.shape, dtype=bool)
            tmp_var2   = tmp_var2.filled(fill_value=np.nan)
            tmp_var2[tmp_var2 <= -1e+38] = np.nan

        var_units = nc_file_D2.variables[varName].units
        nc_time1 = nc_file_D1.variables['time']
        nc_time2 = nc_file_D2.variables['time']
        tmp_dates1 = nc.num2date(nc_time1[:],units=nc_time1.units,calendar=nc_time1.calendar)
        tmp_dates2 = nc.num2date(nc_time2[:],units=nc_time2.units,calendar=nc_time2.calendar)
        if idx == 0:
            var_D1   = tmp_var1
            dates_D1 = tmp_dates1
            mask_D1  = tmp_mask1
            var_D2   = tmp_var2
            dates_D2 = tmp_dates2
            mask_D2  = tmp_mask2
        else:
            var_D1    = np.append(var_D1, tmp_var1, axis=0)
            dates_D1  = np.append(dates_D1, tmp_dates1, axis=0)
            mask_D1   = np.append(mask_D1, tmp_mask1, axis=0)
            var_D2    = np.append(var_D2, tmp_var2, axis=0)
            dates_D2  = np.append(dates_D2, tmp_dates2, axis=0)
            mask_D2   = np.append(mask_D2, tmp_mask2, axis=0)

### check read in data
# think about to exit if D1.shape and D2.shape does not match
print(f'var_D1.shape: {var_D1.shape}')
print(f'mask_D1.shape: {mask_D1.shape}')
print(f'var_D1 max: {np.nanmax(var_D1)}')
print(f'var_D2 min: {np.nanmin(var_D2)}')
print(f'var_D1 has inf: {np.isinf(var_D1).any()}')
print(f'var_D2 has inf: {np.isinf(var_D2).any()}')

###############################################################################
### Do some calculation on read in data
###############################################################################
print(f'Do some calculation on read in data')
var_D1_ave = np.nanmean(var_D1, axis=0)
var_D2_ave = np.nanmean(var_D2, axis=0)

###############################################################################
### Define a special colormap
### seel also: https://matplotlib.org/stable/gallery/color/colormap_reference.html
###############################################################################
print(f'Define a special colormap')
cmap = copy.copy(mpl.cm.get_cmap('Blues'))
cmap.set_under('red')
cmap.set_over('magenta')

###############################################################################
### Define title etc. for plot
###############################################################################
print(f'Define title etc. for plot')
tmp_titlesubstr = [
        f'Averaged saturation of surface layer [-]',
        f'Period: {cftime.datetime.strftime(dates_D1[0], "%Y-%m-%d")} to {cftime.datetime.strftime(dates_D1[-2], "%Y-%m-%d")}',
        ]
kwargs_imshow2PDiff = {
        'title': None,
        #'title': '\n'.join(tmp_titlesubstr),
        'infostr': True,
        'v1_name': 'ERA5',
        'v2_name': 'HiCam',
        'var_vmax': 1,
        'var_vmin': 0,
        'diff_vmax': 0.4,
        'diff_vmin': -0.4,
        'var_cmap': cmap,
        'diff_cmap': mpl.cm.get_cmap('coolwarm'),
        'saveFile': f'./examples_Compare2DatasetsMap.pdf',
        #'saveFile': f'./Diff_{varName}_{cftime.datetime.strftime(dates_D1[0], "%Y%m%d%H%M")}.pdf',
        #'dpi': 100,
        'figsize': (10, 4),
        }

###############################################################################
### Start plotting
###############################################################################
print(f'Start plotting')
PlotLib.plot_imshow2PDiff(v1=var_D1_ave, v2=var_D2_ave,
        **kwargs_imshow2PDiff)
