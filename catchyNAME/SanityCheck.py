"""
author: Nikals WAGNER
e-mail: n.wagner@fz-juelich.de
version: 2021-03-05

Description:
SanityCheck.py is aimed to quickly proof if data show valid results or not.

Often it is required, to quick inspect data (e.g. modeloutput) to validate if
the modelresults are reasonable or not. No nice plots are required, simple
some stats showing as much information as possible about the inspected data
without getting confused.
This script is aimed to handle whatever data it gets and is adjusting
at some upper and lower percentile (default 2% adn 98%) of the passed data
to prevent outliers destroing the colorbar or histogram (check number with the
colorbar! Those do change from plot to plot!).
One usecase of this script could be to run automatically after each individual 
simulation of a long-term climate simulation to easily check if the simulation
results stay reasonable or if some problmes (time shift, instable simualtion, ) 
show up.
Another usecase would be to run manually for some specific data-set, using this
script standalone passing data via command line arguments.

Check examples/ to see how this works. 
"""
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

import cftime

def get_PlotMinMaxMid_Percentil(data2plot, lower=5, upper=95):
    vmin = np.nanpercentile(data2plot, lower)
    vmax = np.nanpercentile(data2plot, upper)
    vmid = (vmax+vmin) / 2.

    return vmin, vmax, vmid

def get_infostr(data, lowerP=2, upperP=98):
    tmp_infostr = [
        f'min: {np.nanmin(data):.2e}',
        f'max: {np.nanmax(data):.2e}',
        f'mean: {np.nanmean(data):.2e}',
        f'std: {np.nanstd(data):.2e}',
        f'q_{lowerP}: {np.nanpercentile(data, lowerP):.2e}',
        f'q_50: {np.nanpercentile(data, 50):.2e}',
        f'q_{upperP}: {np.nanpercentile(data, upperP):.2e}',
        ]
    return '\n'.join(tmp_infostr)

def plot_SanityCheck_3D(data, data_mask=None, kind='sum', 
    figname='./SanityCheck_3D.pdf', fig_title=None, minax_title=None, 
    maxax_title=None, kinax_title=None, hisax_title=None,
    lowerP=2, upperP=98, interactive=False, cmapName='Spectral'):

    ###########################################################################
    #### Small check if data fits requirements
    ###########################################################################
    if not isinstance(data, np.ndarray):
        print(f'data is of type {type(data)} but <class "numpy.ndarray"> is required!')
        return None
    if data.ndim != 3:
        print(f'data is of dimension {data.ndim} but dimension 3 is required!')
        return None


    ###########################################################################
    #### Calculate Min, Max, and Kind (sum or mean) for data
    ###########################################################################
    data_min_T = np.nanmin(data, axis=0)
    data_max_T = np.nanmax(data, axis=0)
    if kind=='sum':
        data_kin_T = np.nansum(data, axis=0)
    elif kind=='mean':
        data_kin_T = np.nanmean(data, axis=0)
    
    if data_mask is not None:
        data_kin_T[data_mask[0]] = np.nan

    ###########################################################################
    #### Defining cmap extend (values below min and above max)
    ###########################################################################
    cmap = copy.copy(mpl.cm.get_cmap(f'{cmapName}'))
    cmap.set_under('black')
    cmap.set_over('magenta')

    ###########################################################################
    #### Create figure
    ###########################################################################
    fig = plt.figure(figsize=(11,8))
    fig.suptitle(fig_title)
    gs  = fig.add_gridspec(nrows=2,ncols=2)
    #gs.update(wspace=0.3, hspace=0.4)
    min_ax = fig.add_subplot(gs[0,0])
    max_ax = fig.add_subplot(gs[0,1])
    kin_ax = fig.add_subplot(gs[1,0])
    his_ax = fig.add_subplot(gs[1,1])

    ###########################################################################
    #### Handling min-plot
    ###########################################################################
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
    tmp_infostr = get_infostr(data_min_T, lowerP=lowerP, upperP=upperP)
    t = min_ax.text(0.01, 0.99, tmp_infostr,
                 verticalalignment='top', transform=min_ax.transAxes,
                 fontsize=8)
    t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='grey'))
    img_min = min_ax.imshow(data_min_T, origin='lower', #vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_min, ax=min_ax, extend='both')

    ###########################################################################
    #### Handling max-plot
    ###########################################################################
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
    tmp_infostr = get_infostr(data_max_T, lowerP=lowerP, upperP=upperP)
    t = max_ax.text(0.01, 0.99, tmp_infostr,
                 verticalalignment='top', transform=max_ax.transAxes,
                 fontsize=8)
    t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='grey'))
    img_max = max_ax.imshow(data_max_T, origin='lower', #vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_max, ax=max_ax, extend='both')


    ###########################################################################
    #### Handling kind-plot (mean or sum)
    ###########################################################################
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
    tmp_infostr = get_infostr(data_kin_T, lowerP=lowerP, upperP=upperP)
    t = kin_ax.text(0.01, 0.99, tmp_infostr,
                 verticalalignment='top', transform=kin_ax.transAxes,
                 fontsize=8)
    t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='grey'))
    img_kin = kin_ax.imshow(data_kin_T, origin='lower', #vmin=tmp_vmin, vmax=tmp_vmax,
                            interpolation='none', cmap=cmap, norm=tmp_norm)
    fig.colorbar(img_kin, ax=kin_ax, extend='both')

    ###########################################################################
    #### Handling histogram of kind-plot
    ###########################################################################
    his_ax.set_title(hisax_title)
    his_ax.set_ylabel(f'# of occurrence')
    hist_data = data_kin_T[~np.isnan(data_kin_T)].flatten()
    range_min = np.nanpercentile(hist_data, lowerP)
    #range_min = np.nanmin(hist_data)
    range_max = np.nanpercentile(hist_data, upperP)
    #range_max = np.nanmax(hist_data)
    #range_max = threshold
    #bins = int(abs((range_max-range_min) // 10))    
    bins = 200   
    print(f'BINS: {bins}')
    his_ax.hist(hist_data, log=True, range=(range_min, range_max), bins=bins)

    ###########################################################################
    #### If interactive --> show plot now
    ###########################################################################
    if interactive is True:
        plt.show()
    else:
        fig.savefig(f'{figname}')

    plt.close('all')

def plot_SanityCheck_2D(data, data_mask=None, figname='./SanityCheck.pdf',
    lowerP=2, upperP=98, interactive=False, **kwargs):
    if not isinstance(data, np.ndarray):
        print(f'data is of type {type(data)} but <class "numpy.ndarray"> is required!')
        return None
    if data.ndim != 2:
        print(f'data is of dimension {data.ndim} but dimension 2 is required!')
        return None

    # get / handle kwargs
    fig_title       = kwargs.pop('fig_title', None)
    kinax_title     = kwargs.pop('kinax_title', None)
    kinax_vmin      = kwargs.pop('kinax_vmin', None)
    kinax_vmax      = kwargs.pop('kinax_vmax', None)
    hisax_title     = kwargs.pop('hisax_title', None)

    data_kin_T = data
    if data_mask is not None:
        data_kin_T[data_mask] = np.nan

    cmap = copy.copy(mpl.cm.get_cmap('Spectral'))
    cmap.set_under('black')
    cmap.set_over('magenta')

    fig = plt.figure(figsize=(11,8))#, dpi=100)
    fig.suptitle(fig_title)
    gs  = fig.add_gridspec(nrows=1,ncols=2)
    #gs.update(wspace=0.3, hspace=0.4)
    kin_ax = fig.add_subplot(gs[0,0])
    his_ax = fig.add_subplot(gs[0,1])

    kin_ax.set_title(kinax_title)
    tmp_vmin, tmp_vmax, tmp_vmid = get_PlotMinMaxMid_Percentil(data_kin_T, lower=lowerP, upper=upperP)
    tmp_vmin = tmp_vmin if kinax_vmin is None else kinax_vmin
    tmp_vmax = tmp_vmax if kinax_vmax is None else kinax_vmax
    if tmp_vmin == tmp_vmax:
        tmp_norm = None
        tmp_vmin = None
        tmp_vmax = None
    elif np.sign(tmp_vmin)*np.sign(tmp_vmax) < 0.0:
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0, vmax=tmp_vmax)
    else:
        print(f'handling tmp_vmin: {tmp_vmin}; tmp_vmid: {tmp_vmid}; tmp_vmax: {tmp_vmax}')
        #tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=tmp_vmid, vmax=tmp_vmax)
        tmp_norm = mcolors.TwoSlopeNorm(vmin=tmp_vmin, vcenter=0.005, vmax=tmp_vmax)
    tmp_infostr = get_infostr(data_kin_T)
    kin_ax.text(0.01, 0.99, tmp_infostr,
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
    #bins = int(abs((range_max-range_min) // 10))
    bins = 200
    print(f'BINS: {bins}')
    his_ax.hist(hist_data, log=True, range=(range_min, range_max), bins=bins)

    if interactive is True:
        plt.show()
    else:
        fig.savefig(f'{figname}')
    plt.close('all')

if __name__=='__main__':

    def parseSlices(args):
        out = tuple(( int(item) if item != 'None' else slice(None) for item in args))
        if out:
            return out
        else:
            return None

    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--infile', '-i', type=str, required=True,
                    help='full file-path')
    parser.add_argument('--varname', '-v', type=str, required=True,
                    help='variable name to check')
    parser.add_argument('--maskValueLower', '-l', type=float, default=None,
            help='in case netCDF file does not proper set fillValue / maskValue (default:None)')
    parser.add_argument('--maskValueUpper', '-u', type=float, default=None,
            help='in case netCDF file does not proper set fillValue / maskValue (default: None)')
    parser.add_argument('--ND', '-d', type=str, default='3D',
            help='2D or 3D data, hmm? jsut set --ND 2D or --ND 3D (default: 3D)')
    parser.add_argument('--Slices', type=str, nargs='+',  default=[0, -1],
    help='In case you pass data with dim not matching ND pass indexing here (default: 0 -1)')

    args           = parser.parse_args()
    infile         = args.infile
    varname        = args.varname
    maskValueLower = args.maskValueLower
    maskValueUpper = args.maskValueUpper
    ND             = args.ND
    Slices         = parseSlices(args.Slices)



    lowerP = 2
    upperP = 100-lowerP

    filename = os.path.splitext(f'{infile}')[0]

    nc_file = nc.Dataset(f'{infile}', 'r')
    #nc_time = nc_file.variables['time']
    #dates = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
    nc_var = nc_file.variables[f'{varname}']
    #nt, ny, nx = nc_var.shape
    #nc_var_unit = nc_var.units
    if ND == '2D' and nc_var.ndim == 2:
        var  = nc_var[...]
    elif ND == '3D' and nc_var.ndim == 3:
        var  = nc_var[...]
    else:
        print(f'inputfile dim ({nc_var.ndim}) does not match ND={ND}')
        print(f'-- do some fancy stuf below or EXIT')
        print(f'Slices: {Slices}')
        var = nc_var.__getitem__(Slices)
        print(f'nc_var.ndim: {nc_var.ndim}')
        print(f'var.shape: {var.shape}')
    
    var_mask  = var.mask
    if not var_mask.any():
        var_mask  = np.zeros(var.shape, dtype=bool)
        var       = var.filled(fill_value=np.nan)
    if maskValueLower is not None:
        var[var <= maskValueLower] = np.nan
    if maskValueUpper is not None:
        var[var >= maskValueUpper] = np.nan


    tmp_title_str = [
        f'Sanity-Check for {filename.split("/")[-1]}',
        f'Var: {varname} -- Slices: {Slices}',
    ]
    fig_title    = '\n'.join(tmp_title_str)
    figname      = f'{filename}_SanityCheck_t{Slices[0]}.pdf'
    minax_title  = f'{varname} min'
    maxax_title  = f'{varname} max'
    kinax_title  = f'{varname} sum'
    hisax_title  = f'{varname} sum - distribution'

    if ND == '3D':
        plot_SanityCheck_3D(data=var, data_mask=var_mask, kind='mean', figname=figname,
            fig_title=fig_title, minax_title=minax_title, 
            maxax_title=maxax_title, kinax_title=kinax_title,
            hisax_title=hisax_title)
    elif ND == '2D':
        plot_SanityCheck_2D(data=var, data_mask=var_mask, figname=figname,
            fig_title=fig_title, kinax_title=kinax_title,
            hisax_title=hisax_title, kinax_vmin=0, kinax_vmax=0.01,
            interactive=False)
    
    # proper close netCDF file
    nc_file.close()
