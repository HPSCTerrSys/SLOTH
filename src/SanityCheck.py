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

def plot_SanityCheck(data, kind='sum', data_mask=None, figname='./SanityCheck', 
    lowerP=2, upperP=98, **kwargs):
    if not isinstance(data, np.ndarray):
        print(f'data is of type {type(data)} but <class "numpy.ndarray"> is required!')
        return None
    if data.ndim != 3:
        print(f'data is of dimension {data.ndim} but dimension 3 is required!')
        return None

    # get / handle kwargs
    fig_title       = kwargs.pop('fig_title', None)
    minax_title     = kwargs.pop('minax_title', None)
    maxax_title     = kwargs.pop('maxax_title', None)
    kinax_title     = kwargs.pop('kinax_title', None)
    hisax_title     = kwargs.pop('hisax_title', None)


    data_min_T = np.nanmin(data, axis=0)
    data_max_T = np.nanmax(data, axis=0)
    if kind=='sum':
        data_kin_T = np.nansum(data, axis=0)
        data_kin_T[data_mask[0]] = np.nan
    elif kind=='mean':
        data_kin_T = np.nanmean(data, axis=0)
        data_kin_T[data_mask[0]] = np.nan

    cmap = copy.copy(mpl.cm.get_cmap('Spectral'))
    cmap.set_under('black')
    cmap.set_over('magenta')

    fig = plt.figure(figsize=(11,8), dpi=100)
    fig.suptitle(fig_title)
    gs  = fig.add_gridspec(nrows=2,ncols=2)
    #gs.update(wspace=0.3, hspace=0.4)
    min_ax = fig.add_subplot(gs[0,0])
    max_ax = fig.add_subplot(gs[0,1])
    kin_ax = fig.add_subplot(gs[1,0])
    his_ax = fig.add_subplot(gs[1,1])

    min_ax.set_title(minax_title)
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(data_min_T, lower=lowerP, upper=upperP)
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
    tmp_infostr = [
            f'min: {np.nanmin(data_min_T):.2e}',
            f'max: {np.nanmax(data_min_T):.2e}',
            f'mean: {np.nanmean(data_min_T):.2e}',
            f'std: {np.nanstd(data_min_T):.2e}',
            f'q_{lowerP}: {np.nanpercentile(data_min_T, lowerP):.2e}',
            f'q_50: {np.nanpercentile(data_min_T, 50):.2e}',
            f'q_{upperP}: {np.nanpercentile(data_min_T, upperP):.2e}',
            ]
    min_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                 verticalalignment='top', transform=min_ax.transAxes,
                 fontsize=8)
    img_min = min_ax.imshow(data_min_T, origin='lower', vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_min, ax=min_ax, extend='both')


    max_ax.set_title(maxax_title)
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(data_max_T, lower=lowerP, upper=upperP)
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
    tmp_infostr = [
            f'min: {np.nanmin(data_max_T):.2e}',
            f'max: {np.nanmax(data_max_T):.2e}',
            f'mean: {np.nanmean(data_max_T):.2e}',
            f'std: {np.nanstd(data_max_T):.2e}',
            f'q_5{lowerP}: {np.nanpercentile(data_max_T, lowerP):.2e}',
            f'q_50: {np.nanpercentile(data_max_T, 50):.2e}',
            f'q_{upperP}: {np.nanpercentile(data_max_T, upperP):.2e}',
            ]
    max_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                 verticalalignment='top', transform=max_ax.transAxes,
                 fontsize=8)
    img_max = max_ax.imshow(data_max_T, origin='lower', vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_max, ax=max_ax, extend='both')


    kin_ax.set_title(kinax_title)
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(data_kin_T, lower=lowerP, upper=upperP)
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
    tmp_infostr = [
            f'min: {np.nanmin(data_kin_T,):.2e}',
            f'max: {np.nanmax(data_kin_T,):.2e}',
            f'mean: {np.nanmean(data_kin_T,):.2e}',
            f'std: {np.nanstd(data_kin_T,):.2e}',
            f'q_{lowerP}: {np.nanpercentile(data_kin_T, lowerP):.2e}',
            f'q_50: {np.nanpercentile(data_kin_T, 50):.2e}',
            f'q_{upperP}: {np.nanpercentile(data_kin_T, upperP):.2e}',
            ]
    kin_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                 verticalalignment='top', transform=kin_ax.transAxes,
                 fontsize=8)
    img_kin = kin_ax.imshow(data_kin_T, origin='lower', vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_kin, ax=kin_ax, extend='both')

    his_ax.set_title(hisax_title)
    his_ax.set_ylabel(f'# of occurrence')
    hist_data = data_kin_T[~np.isnan(data_kin_T)].flatten()
    range_min = np.nanpercentile(hist_data, lowerP)
    #range_min = np.nanmin(hist_data)
    range_max = np.nanpercentile(hist_data, upperP)
    #range_max = np.nanmax(hist_data)
    #range_max = threshold
    bins = int(abs((range_max-range_min) // 10))    
    print(f'BINS: {bins}')
    his_ax.hist(hist_data, log=True, range=(range_min, range_max), bins=bins)

    plt.savefig(f'{figname}')
    plt.close('all')

    pass

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--infile', '-i', type=str, nargs='+', required=True,
                    help='full file-path')
    parser.add_argument('--varname', '-v', type=str, default=None,
                    help='variable name to check')

    args        = parser.parse_args()
    infile      = args.infile
    in_varname  = args.varname

    lowerP = 2
    upperP = 100-lowerP

    filename = os.path.splitext(f'{infile}')[0]
    varname = in_varname if in_varname is not None else filename.split('/')[-1].replace('_ts', '')

    nc_file = nc.Dataset(f'{infile}', 'r')
    nc_time = nc_file.variables['time']
    dates = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
    nc_var = nc_file.variables[f'{varname}']
    nt, ny, nx = nc_var.shape
    nc_var_unit = nc_var.units
    var  = nc_var[...]
    var_mask  = var.mask
    if not var_mask.any():
        var_mask  = np.zeros(var.shape, dtype=bool)
        var       = var.filled(fill_value=np.nan)


    tmp_title_str = [
        f'Sanity-Check for {cftime.datetime.strftime(dates[0], "%Y-%m-%d")} to {cftime.datetime.strftime(dates[-1], "%Y-%m-%d")}',
        f'Var: {varname} in [{nc_var_unit}]',
    ]
    fig_title    = '\n'.join(tmp_title_str)
    figname      = f'{filename}_SanityCheck.png'
    minax_title  = f'{varname} min'
    maxax_title  = f'{varname} max'
    kinax_title  = f'{varname} sum'
    hisax_title  = f'{varname} sum - distribution'

    plot_SanityCheck(data=var, data_mask=var_mask, kind='mean', figname=figname,
        fig_title=fig_title, minax_title=minax_title, 
        maxax_title=maxax_title, kinax_title=kinax_title,
        hisax_title=hisax_title)
    
    # proper close netCDF file
    nc_file.close()

    # for infile in infiles:
    #     print(f'processing {infile}')
    #     filename = os.path.splitext(f'{infile}')[0]
    #     varname = in_varname if in_varname is not None else filename.split('/')[-1].replace('_ts', '')
    #     try:
    #         nc_file = nc.Dataset(f'{infile}', 'r')
    #         nc_time = nc_file.variables['time']
    #         dates = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
    #         nc_var = nc_file.variables[f'{varname}']
    #         nt, ny, nx = nc_var.shape
    #         nc_var_unit = nc_var.units
    #         var  = nc_var[...]
    #         var_mask  = var.mask
    #         if not var_mask.any():
    #             var_mask  = np.zeros(var.shape, dtype=bool)
    #             var       = var.filled(fill_value=np.nan)


    #         tmp_title_str = [
    #             f'Sanity-Check for {cftime.datetime.strftime(dates[0], "%Y-%m-%d")} to {cftime.datetime.strftime(dates[-1], "%Y-%m-%d")}',
    #             f'Var: {varname} in [{nc_var_unit}]',
    #         ]
    #         fig_title    = '\n'.join(tmp_title_str)
    #         figname      = f'{filename}_SanityCheck.png'
    #         minax_title  = f'{varname} min'
    #         maxax_title  = f'{varname} max'
    #         kinax_title  = f'{varname} sum'
    #         hisax_title  = f'{varname} sum - distribution'

    #         plot_SanityCheck(data=var, data_mask=var_mask, kind='mean', figname=figname,
    #             fig_title=fig_title, minax_title=minax_title, 
    #             maxax_title=maxax_title, kinax_title=kinax_title,
    #             hisax_title=hisax_title)
            
    #         # proper close netCDF file
    #         nc_file.close()
    #     except:
    #         print(f'Unexpected error: {sys.exc_info()}')
    #         #raise
    #         print(f'Own ERROR while reading {filename} --> CONTINUE!')
    #         continue
