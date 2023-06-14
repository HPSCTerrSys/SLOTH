#!/usr/bin/env python
"""
SanityCheck.py is a script aimed at quickly validating the validity of data by 
providing simple statistical information.

Often, it is necessary to inspect data, such as model output, to determine if 
the results are reasonable or not. Instead of elaborate plots, the script 
focuses on presenting relevant statistics to avoid confusion while conveying 
as much information about the data as possible.

The script can handle various types of data and automatically adjusts the 
colorbar according to upper and lower percentiles to account for outliers. This 
prevents outliers from distorting the colorbar or histogram. Therefore it is 
important to check the values on on each colorbar, as they may vary from plot 
to plot.

One use case for this script is to run it automatically after each individual 
simulation in a long-term climate simulation. This allows for easy verification 
of the simulation results' reasonableness, detecting problems like time shifts 
or instability.

Another use case is to run the script manually on specific data sets by passing 
them as command line arguments.

For examples of how to use this script, refer to the 'examples/' directory.
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

def get_PlotMinMaxMid_Percentil(data, lower=5, upper=95):
    """
    Calculate the minimum, maximum, and midpoint values based on percentiles of the input data.

    Parameters
    ----------
    data : numpy.ma.MaskedArray
		Input data array.
    lower : int, optional
		Lower percentile value. Default is 5.
    upper : int, optional
		Upper percentile value. Default is 95.

    Returns
    -------
    tuple : 
		A tuple containing the minimum, maximum, and midpoint values.

    Examples
    --------
    >>> import numpy as np
    >>> import numpy.ma as ma
    >>> data = ma.masked_array([1.2, 3.4, 5.6, 7.8, 9.0], mask=[False, False, True, False, False])
    >>> min_val, max_val, mid_val = get_PlotMinMaxMid_Percentil(data, lower=10, upper=90)
    >>> print("Minimum value:", min_val)
    Minimum value: 1.86
    >>> print("Maximum value:", max_val)
    Maximum value: 8.64
    >>> print("Midpoint value:", mid_val)
    Midpoint value: 5.25

    """
    # Compress data to remove masked values, which are not taken into account.
    # Note that ma.compressed() is returning a 1D array! So pay attantion 
    # when to use!
    data_compressed = data.compressed()
    vmin = np.percentile(data_compressed, lower)
    vmax = np.percentile(data_compressed, upper)
    vmid = (vmax+vmin) / 2.

    return vmin, vmax, vmid

def get_infostr(data, lowerP=2, upperP=98):
    # Compress data to remove masked values, which are not taken into account
    data_compressed = data.compressed()
    tmp_infostr = [
        f'min: {np.min(data_compressed):.2e}',
        f'max: {np.max(data_compressed):.2e}',
        f'mean: {np.mean(data_compressed):.2e}',
        f'std: {np.std(data_compressed):.2e}',
        f'q_{lowerP}: {np.percentile(data_compressed, lowerP):.2e}',
        f'q_50: {np.percentile(data_compressed, 50):.2e}',
        f'q_{upperP}: {np.percentile(data_compressed, upperP):.2e}',
        ]
    return '\n'.join(tmp_infostr)

def plot_SanityCheck(data, kind='mean', 
    figname='./SanityCheck_3D.pdf', fig_title=None,
    lowerP=2, upperP=98, interactive=False, cmapName='Spectral'):
    """
    Plot a sanity check for given data.
    
    The Sanity Plot is a plot consisting of 4 sub-plots, trying to visualise 
    some important statistics in a compact way, aimed to determine if the data 
    set inspected is plausible or not (sanity).

    Parameters
    ----------
    data : numpy.ma.MaskedArray
        3D array of data.
    kind : {'sum', 'mean'}, optional
        Calculation type for the data statistics. Default is 'mean'.
    figname : str, optional
        File name to save the plot. Default is './SanityCheck_3D.pdf'.
    fig_title : str, optional
        Title of the plot.
    lowerP : int, optional
        Lower percentile value for plot limits. Default is 2.
    upperP : int, optional
        Upper percentile value for plot limits. Default is 98.
    interactive : bool, optional
        If True, display the plot interactively. If False, save the plot to 
        'figname'. Default is False.
    cmapName : str, optional
        Name of the colormap. Default is 'Spectral'.

    Returns
    -------
    None

    Notes
    -----
    - The 'data' input must be a 3D numpy masked array (t, y, x).
    - The 'kind' parameter specifies whether to calculate the sum or mean of the data.
    - The function generates a plot with subplots for the minimum, maximum, kind (sum or mean), and histogram of the data.
    - The colormap normalization is determined based on the percentiles of the data.
    - The plot can be displayed interactively or saved to a file.
    - If the 'interactive' parameter is set to True, the plot is displayed using plt.show(). If False, the plot is saved to the file specified by 'figname'.

    """

    ###########################################################################
    #### Small check if data fits requirements
    ###########################################################################
    if not isinstance(data, np.ma.MaskedArray):
        print(f'data is of type {type(data)} but <class "numpy.ma.core.MaskedArray"> is required!')
        return None
    if data.ndim != 3:
        print(f'data is of dimension {data.ndim} but dimension 3 is required!')
        return None

    ###########################################################################
    #### Calculate Min, Max, and Kind (sum or mean) for data
    ###########################################################################
    data_min_T = np.ma.min(data, axis=0)
    data_max_T = np.ma.max(data, axis=0)
    if kind=='sum':
        data_kin_T = np.ma.sum(data, axis=0)
    elif kind=='mean':
        data_kin_T = np.ma.mean(data, axis=0)
    
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
    min_ax.set_title('min')
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
    max_ax.set_title('max')
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
    kin_ax.set_title(kind)
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
    his_ax.set_title(f'histogram of {kind}')
    his_ax.set_ylabel(f'# of occurrence')
    hist_data = data_kin_T[~np.isnan(data_kin_T)].flatten()
    range_min = np.nanpercentile(hist_data, lowerP)
    range_max = np.nanpercentile(hist_data, upperP)
    bins = 200   
    his_ax.hist(hist_data, log=True, range=(range_min, range_max), bins=bins)

    ###########################################################################
    #### If interactive --> show plot now
    ###########################################################################
    if interactive is True:
        plt.show()
    else:
        fig.savefig(f'{figname}')

    plt.close('all')

if __name__=='__main__':

    def parseSlices(args):
        if args is not None:
            out = tuple(( int(item) if item != 'None' else slice(None) for item in args))
            return out
        else:
            return slice(None)

    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--infile', '-i', type=str, required=True,
            help='full file-path')
    parser.add_argument('--varname', '-v', type=str, required=True,
            help='variable name to check')
    parser.add_argument('--maskValueLower', '-l', type=float, default=None,
            help='in case netCDF file does not proper set fillValue / maskValue (default:None)')
    parser.add_argument('--maskValueUpper', '-u', type=float, default=None,
            help='in case netCDF file does not proper set fillValue / maskValue (default: None)')
    parser.add_argument('--Slices', type=str, nargs='+',  default=None,
            help='In case you pass data with dim not suitable for this SanityCheck (default: 0 -1)')

    args           = parser.parse_args()
    infile         = args.infile
    varname        = args.varname
    maskValueLower = args.maskValueLower
    maskValueUpper = args.maskValueUpper
    Slices         = parseSlices(args.Slices)

    lowerP = 2
    upperP = 100-lowerP

    filename = os.path.splitext(f'{infile}')[0]

    nc_file = nc.Dataset(f'{infile}', 'r')
    nc_var = nc_file.variables[f'{varname}']
    var = nc_var.__getitem__(Slices)
    if nc_var.ndim == 2:
        print(f'DEBUG: input data is 2D --> expand one dim to run this script')
        var = var[np.newaxis,...]
    
    if maskValueLower is not None:
        var = np.ma.masked_where(var < maskValueLower, var)
    if maskValueUpper is not None:
        var = np.ma.masked_where(var > maskValueUpper, var)
    
    tmp_title_str = [
        f'Sanity-Check for {filename.split("/")[-1]}',
        f'Var: {varname} -- Slices: {Slices}',
    ]
    fig_title    = '\n'.join(tmp_title_str)
    figname      = f'{filename}_SanityCheck.pdf'
    minax_title  = f'{varname} min'
    maxax_title  = f'{varname} max'
    kinax_title  = f'{varname} mean'
    hisax_title  = f'{varname} mean - distribution'

    plot_SanityCheck(data=var, kind='mean', 
            figname=figname, fig_title=fig_title, minax_title=minax_title, 
            maxax_title=maxax_title, kinax_title=kinax_title,
            hisax_title=hisax_title)
    
    # proper close netCDF file
    nc_file.close()
