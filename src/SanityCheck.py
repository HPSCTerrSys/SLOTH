import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import netCDF4 as nc
import argparse
import sys
import os
import copy

import nc_time_axis
import cftime

def get_PlotMinMaxMid_Percentil(data2plot, lower=5, upper=95):
    vmin = np.nanpercentile(data2plot, lower)
    vmax = np.nanpercentile(data2plot, upper)
    #vmid = np.nanpercentile(data2plot, 50)
    vmid = (vmax+vmin) / 2.
    '''
    tmp_scalFac = 1.5
    vmin = tmp_mean - (tmp_scalFac*tmp_std)
    vmax = tmp_mean + (tmp_scalFac*tmp_std)
    '''
    return vmin, vmax, vmid

parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
parser.add_argument('--infiles', '-i', type=str, nargs='+', required=True,
                    help='full file-path')
parser.add_argument('--varname', '-v', type=str,
                    help='variable name to check')
parser.add_argument('--realVarName', '-r', type=str, default=None,
                    help='real variable name to show')

args        = parser.parse_args()
infiles     = args.infiles
varname     = args.varname
realVarName = args.realVarName if args.realVarName is not None else args.varname

lowerP = 2
upperP = 100-lowerP

for infile in infiles:
    filename = os.path.splitext(f'{infile}')[0]
    #with nc.Dataset(f'{infile}', 'r') as nc_file:
    nc_file = nc.Dataset(f'{infile}', 'r')
    #print(f'nc_file: \n{nc_file}')
    nc_time = nc_file.variables['time']
    dates = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
    nc_var = nc_file.variables[f'{varname}']
    #print(f'nc_var: \n{nc_var}')
    #nc_var = nc_var[:15,...]
    nt, ny, nx = nc_var.shape
    #print(f'var.shape: {var.shape}')
    #sys.exit(0)

    print('###################################################################')
    print('###################################################################')
    var       = nc_var[...]
    var_mask  = var.mask
    if not var_mask.any():
        var_mask = np.zeros(var.shape, dtype=bool)
    var       = var.filled(fill_value=np.nan)
    var_min_T = np.nanmin(var, axis=0)
    var_max_T = np.nanmax(var, axis=0)
    var_sum_T = np.nansum(var, axis=0)
    var_sum_T[var_mask[0]] = np.nan

    #clevs = [0, 2, 5, 10, 15, 20, 30, 40, 60, 80,
    #         100, 140, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500]
    #cmap_data = [(1.0, 1.0, 1.0),
    #         (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
    #         (0.0, 1.0, 1.0),
    #         (0.0, 0.8784313797950745, 0.501960813999176),
    #         (0.0, 0.7529411911964417, 0.0),
    #         (0.501960813999176, 0.8784313797950745, 0.0),
    #         (1.0, 1.0, 0.0),
    #         (1.0, 0.6274510025978088, 0.0),
    #         (1.0, 0.0, 0.0),
    #         (1.0, 0.125490203499794, 0.501960813999176),
    #         (0.9411764740943909, 0.250980406999588, 1.0),
    #         (0.501960813999176, 0.125490203499794, 1.0),
    #         (0.250980406999588, 0.250980406999588, 1.0),
    #         (0.125490203499794, 0.125490203499794, 0.501960813999176),
    #         (0.125490203499794, 0.125490203499794, 0.125490203499794),
    #         (0.501960813999176, 0.501960813999176, 0.501960813999176),
    #         (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
    #         (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
    #         (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
    #         (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
    #         (0.4000000059604645, 0.20000000298023224, 0.0)]
    #cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    #norm = mcolors.BoundaryNorm(clevs, cmap.N)
    cmap = copy.copy(mpl.cm.get_cmap('Spectral'))
    #cmap = copy.copy(mpl.cm.get_cmap('viridis'))
    cmap.set_under('black')
    cmap.set_over('magenta')
  
    fig = plt.figure(figsize=(11,8), dpi=100)
    tmp_title_str = [
            f'Sanity-Check for {cftime.datetime.strftime(dates[0], "%Y-%m-%d")} to {cftime.datetime.strftime(dates[-1], "%Y-%m-%d")}',
            f'Var: {realVarName} in [mm]',
        ]
    fig.suptitle('\n'.join(tmp_title_str))
    gs  = fig.add_gridspec(nrows=2,ncols=2)
    #gs.update(wspace=0.3, hspace=0.4)
    min_ax = fig.add_subplot(gs[0,0])
    max_ax = fig.add_subplot(gs[0,1])
    sum_ax = fig.add_subplot(gs[1,0])
    his_ax = fig.add_subplot(gs[1,1])
    
    min_ax.set_title(f'{realVarName} min')
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(var_min_T, lower=lowerP, upper=upperP)
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
    tmp_infostr = [
            f'min: {np.nanmin(var_min_T):.2e}',
            f'max: {np.nanmax(var_min_T):.2e}',
            f'mean: {np.nanmean(var_min_T):.2e}',
            f'std: {np.nanstd(var_min_T):.2e}',
            f'q_{lowerP}: {np.nanpercentile(var_min_T, lowerP):.2e}',
            f'q_50: {np.nanpercentile(var_min_T, 50):.2e}',
            f'q_{upperP}: {np.nanpercentile(var_min_T, upperP):.2e}',
            ]
    min_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                 verticalalignment='top', transform=min_ax.transAxes,
                 fontsize=8)
    img_min = min_ax.imshow(var_min_T, origin='lower', vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_min, ax=min_ax, extend='both')
    

    max_ax.set_title(f'{realVarName} max')
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(var_max_T, lower=lowerP, upper=upperP)
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
    tmp_infostr = [
            f'min: {np.nanmin(var_max_T):.2e}',
            f'max: {np.nanmax(var_max_T):.2e}',
            f'mean: {np.nanmean(var_max_T):.2e}',
            f'std: {np.nanstd(var_max_T):.2e}',
            f'q_5{lowerP}: {np.nanpercentile(var_max_T, lowerP):.2e}',
            f'q_50: {np.nanpercentile(var_max_T, 50):.2e}',
            f'q_{upperP}: {np.nanpercentile(var_max_T, upperP):.2e}',
            ]
    max_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                 verticalalignment='top', transform=max_ax.transAxes,
                 fontsize=8)
    img_max = max_ax.imshow(var_max_T, origin='lower', vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_max, ax=max_ax, extend='both')


    sum_ax.set_title(f'{realVarName} sum')
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(var_sum_T, lower=lowerP, upper=upperP)
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
    tmp_infostr = [
            f'min: {np.nanmin(var_sum_T,):.2e}',
            f'max: {np.nanmax(var_sum_T,):.2e}',
            f'mean: {np.nanmean(var_sum_T,):.2e}',
            f'std: {np.nanstd(var_sum_T,):.2e}',
            f'q_{lowerP}: {np.nanpercentile(var_sum_T, lowerP):.2e}',
            f'q_50: {np.nanpercentile(var_sum_T, 50):.2e}',
            f'q_{upperP}: {np.nanpercentile(var_sum_T, upperP):.2e}',
            ]
    sum_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                 verticalalignment='top', transform=sum_ax.transAxes,
                 fontsize=8)
    img_sum = sum_ax.imshow(var_sum_T, origin='lower', vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_sum, ax=sum_ax, extend='both')

    his_ax.set_title(f'{realVarName} sum - distribution')
    #his_ax.set_xlabel(f'{realVarName} [mm / {nt} days]')
    his_ax.set_ylabel(f'# of occurrence')
    hist_data = var_sum_T[~np.isnan(var_sum_T)].flatten()
    range_min = np.nanpercentile(hist_data, lowerP)
    #range_min = np.nanmin(hist_data)
    range_max = np.nanpercentile(hist_data, upperP)
    #range_max = np.nanmax(hist_data)
    #range_max = threshold
    bins = int((range_max-range_min) // 10)
    print(f'BINS: {bins}')
    his_ax.hist(hist_data, log=True, range=(range_min, range_max), bins=bins)

    plt.savefig(f'{filename}_SanityCheck.png')

    # proper close netCDF file
    nc_file.close()
