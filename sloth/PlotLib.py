import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import datetime as dt
import cftime
import sys
import copy

def plot_XY_2VarMean_TwinX(x, y1, y2, ax, **kwargs):
    ''' 
    perform a XY plot with dates at x-axis for two variables with
    seperate y-axis
    accepts **kwargs
    '''
    title     = kwargs.pop('title', None)
    infostr   = kwargs.pop('infostr', '')
    y1_label  = kwargs.pop('y1_label', 'y1_default')
    y1_ylabel = kwargs.pop('y1_ylabel', 'y1')
    y1_hlines = kwargs.pop('y1_hlines', None)
    y2_label  = kwargs.pop('y2_label', 'y2_default')
    y2_ylabel = kwargs.pop('y2_ylabel', 'y2')
    y2_hlines = kwargs.pop('y2_hlines', None)

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.set_title(f'{title}')
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_ylabel(fr'{y1_ylabel}')
    ax.plot(x, y1, lw=0.8, color='blue', label=fr'{y1_label}')
    if y1_hlines:
        [ax.axhline(hline, ls='--', lw=0.5, color='blue') for hline in y1_hlines]
    if infostr:
        ax.text(0.01, 0.99, fr'{infostr}',
                verticalalignment='top', transform=ax.transAxes)
    ax.legend(loc='lower left')
    
    ax_twin=ax.twinx()
    ax_twin.xaxis.set_major_locator(locator)
    ax_twin.xaxis.set_major_formatter(formatter)
    ax_twin.set_ylabel(fr'{y2_ylabel}')
    ax_twin.plot(x, y2, lw=0.8, color='grey', label=fr'{y2_label}')
    if y2_hlines:
        [ax_twin.axhline(hline, ls='--', lw=0.5, color='grey') for hline in y2_hlines]
    ax_twin.legend(loc='lower right')

def get_PlotMinMaxMid_Percentil(data, lower=2, upper=98):
    if isinstance(data, np.ma.MaskedArray):
        # Compress data to remove masked values, which are not taken into
        # account by np.percentile
        data = data.compressed()
    vmin = np.percentile(data, lower)
    vmax = np.percentile(data, upper)
    vmid = (vmax+vmin) / 2.
    '''
    tmp_scalFac = 1.5
    vmin = tmp_mean - (tmp_scalFac*tmp_std)
    vmax = tmp_mean + (tmp_scalFac*tmp_std)
    '''
    return vmin, vmax, vmid


def get_infostr(data, lowerP=2, upperP=98):
    if isinstance(data, np.ma.MaskedArray):
        # Compress data to remove masked values, which are not taken into
        # account by np.percentile
        data = data.compressed()
    tmp_infostr = [
        f'min: {np.min(data):.2e}',
        f'max: {np.max(data):.2e}',
        f'mean: {np.mean(data):.2e}',
        f'std: {np.std(data):.2e}',
        f'q_{lowerP}: {np.percentile(data, lowerP):.2e}',
        f'q_50: {np.percentile(data, 50):.2e}',
        f'q_{upperP}: {np.percentile(data, upperP):.2e}',
        f'2D sum: {np.sum(data):.2e}',
        ]
    return '\n'.join(tmp_infostr)

def plot_imshow2PDiff(v1, v2, **kwargs):
    '''
    generates 3 imshows in a rwo
    plot 2 vars and difference of both
    '''
    title       = kwargs.pop('title', None)
    infostr     = kwargs.pop('infostr', False)
    v1_name     = kwargs.pop('v1_name', 'Variable 1')
    v2_name     = kwargs.pop('v2_name', 'Variable 2')
    var_vmin    = kwargs.pop('var_vmin', get_PlotMinMaxMid_Percentil(v1)[0])
    var_vmax    = kwargs.pop('var_vmax', get_PlotMinMaxMid_Percentil(v1)[1])
    diff_vmin   = kwargs.pop('diff_vmin', None)
    diff_vmax   = kwargs.pop('diff_vmax', None)
    # a bit more special default handling
    if 'var_cmap' in kwargs:
        var_cmap    = kwargs.pop('var_cmap')
    else:
        var_cmap = copy.copy(mpl.cm.get_cmap('Reds'))
        var_cmap.set_under('blue')
        var_cmap.set_over('magenta')
    diff_cmap   = kwargs.pop('diff_cmap', plt.get_cmap('coolwarm'))
    saveFile    = kwargs.pop('saveFile', None)
    dpi         = kwargs.pop('dpi', 100)
    figsize     = kwargs.pop('figsize', (8.27, 11.69))

    if var_vmin == var_vmax:
        var_norm = None
    elif np.sign(var_vmin)*np.sign(var_vmax) < 0.0:
        var_norm = mcolors.TwoSlopeNorm(vmin=var_vmin, vcenter=0, vmax=var_vmax)
    else:
        var_norm = mcolors.TwoSlopeNorm(vmin=var_vmin, vcenter=(var_vmax+var_vmin)/2., vmax=var_vmax)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    if title is not None:
        fig.suptitle(f'{title}')
    gs  = fig.add_gridspec(nrows=1,ncols=11)
    gs.update(wspace=0.3, hspace=0.4)
    var1_ax = fig.add_subplot(gs[0,0:3])
    var2_ax = fig.add_subplot(gs[0,3:6])
    var_cax = inset_axes(var2_ax,
               width="5%",  # width = 5% of parent_bbox width
               height="100%",  # height : 340% good for a (4x4) Grid
               loc='lower left',
               bbox_to_anchor=(1.05, 0., 1, 1),
               bbox_transform=var2_ax.transAxes,
               borderpad=0,
               )
    diff_ax  = fig.add_subplot(gs[0,7:10])
    diff_cax = inset_axes(diff_ax,
               width="5%",  # width = 5% of parent_bbox width
               height="100%",  # height : 340% good for a (4x4) Grid
               loc='lower left',
               bbox_to_anchor=(1.05, 0., 1, 1),
               bbox_transform=diff_ax.transAxes,
               borderpad=0,
               )

    var1_ax.set_title(fr'{v1_name}')
    if infostr:
        var1_ax.text(0.01, 0.99, get_infostr(v1),
                     verticalalignment='top', transform=var1_ax.transAxes,
                     fontsize=8)
    im = var1_ax.imshow(v1, origin='lower', norm=var_norm,# vmin=var_vmin, vmax=var_vmax,
            cmap=var_cmap, interpolation='none')

    var2_ax.set_title(fr'{v2_name}')
    var2_ax.get_yaxis().set_visible(False)
    if infostr:
        var2_ax.text(0.01, 0.99, get_infostr(v2),
                     verticalalignment='top', transform=var2_ax.transAxes,
                     fontsize=8)
    im = var2_ax.imshow(v2, origin='lower', norm=var_norm,#vmin=var_vmin, vmax=var_vmax,
            cmap=var_cmap, interpolation='none')
    cb = plt.colorbar(im, cax=var_cax, orientation='vertical',
            extend='both')

    diff_ax.set_title(fr'Diff: {v1_name} - {v2_name}')
    diff_ax.get_yaxis().set_visible(False)
    if infostr:
        diff_ax.text(0.01, 0.99, get_infostr(v1-v2),
                     verticalalignment='top', transform=diff_ax.transAxes,
                     fontsize=8)
    im = diff_ax.imshow(v1 - v2, origin='lower',
            vmin=diff_vmin, vmax=diff_vmax, 
            cmap=diff_cmap, interpolation='none')
    cb = plt.colorbar(im, cax=diff_cax, orientation='vertical',
            extend='both')

    if saveFile is not None:
        plt.savefig(f'{saveFile}', bbox_inches='tight', pad_inches=0)
    return fig, var1_ax, var2_ax, diff_ax

def plot_ClimateYearMonth(climateData, **kwargs):
    '''
    generates 3 imshows in a rwo
    plot 2 vars and difference of both
    '''
    title       = kwargs.pop('title', None)
    infostr     = kwargs.pop('infostr', False)
    var_vmin    = kwargs.pop('var_vmin', get_PlotMinMaxMid_Percentil(climateData)[0])
    var_vmax    = kwargs.pop('var_vmax', get_PlotMinMaxMid_Percentil(climateData)[1])
    # a bit more special default handling
    if 'var_cmap' in kwargs:
        var_cmap    = kwargs.pop('var_cmap')
    else:
        var_cmap = copy.copy(mpl.cm.get_cmap('Reds'))
        var_cmap.set_under('blue')
        var_cmap.set_over('magenta')
    saveFile    = kwargs.pop('saveFile', './DefaultSaveName_plot_ClimateYearMonth.png')
    dpi         = kwargs.pop('dpi', 300)
    figsize     = kwargs.pop('figsize', (8.5, 10.5))

    titles      = ['Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


    fig = plt.figure(figsize=figsize, dpi=dpi)
    if title is not None:
        fig.suptitle(f'{title}')
    gs  = fig.add_gridspec(nrows=3,ncols=4)
    gs.update(wspace=0.1, hspace=0.4)

    for month in range(climateData.shape[0]):
        month_data = climateData[month]
        ax = fig.add_subplot(gs[month//4,month%4])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title(titles[month], fontsize=7)

        if infostr:
            t = ax.text(0.03, 0.97, get_infostr(month_data),
                    verticalalignment='top', transform=ax.transAxes,
                    fontsize=3)
            t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='grey'))
        im = ax.imshow(month_data, origin='lower', vmin=var_vmin, vmax=var_vmax,
                cmap=var_cmap, interpolation='none')

    cbar_ax = fig.add_axes([0.95, 0.15, 0.03, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    fig.savefig(f'{saveFile}', bbox_inches='tight', pad_inches=0)

def plot_MappedSubAreas(mapper, fit_name='NotSet', search_rad=3, save_dir='../data'):
    for idx, ID in enumerate(mapper.ObsIDs):
        print(f'--- plotting ObsID {ID}')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        rawX = mapper.MapXIdx_raw[idx]
        rawY = mapper.MapYIdx_raw[idx]
        data2plot = mapper.SimMeanQ[rawY-search_rad:rawY+search_rad+1, rawX-search_rad:rawX+search_rad+1].copy()
        data2plot /= np.nanmax(mapper.SimMeanQ[rawY-search_rad:rawY+search_rad+1, rawX-search_rad:rawX+search_rad+1])
        im = ax.imshow(data2plot, origin='lower')
        fig.colorbar(im, ax=ax)
        # raw data is always the centre = search_rad
        ax.scatter(search_rad ,search_rad, c='red', marker='x', label='raw')
        x_sliced = mapper.MapXIdx_fit[idx] - ( rawX - search_rad)
        y_sliced = mapper.MapYIdx_fit[idx] - ( rawY - search_rad)
        ax.scatter(x_sliced, y_sliced, c='red', marker='o', label='best')
        ax.legend()
        title_strs = [
                   f'GRDC-ID: {ID}',
                   f'fit-routine used: {fit_name}',
                   f'ObsMeanQ: {mapper.ObsMeanQ[idx]:.2f} m^3/s',
                   f'SimMeanQ raw: {mapper.SimMeanQ[mapper.MapYIdx_raw[idx], mapper.MapXIdx_raw[idx]]:.2f} m^3/s',
                   f'SimMeanQ fit: {mapper.SimMeanQ[mapper.MapYIdx_fit[idx], mapper.MapXIdx_fit[idx]]:.2f} m^3/s',
                   # f'ObsMeanArea / SimMeanArea: {mapper.ObsMeanArea[idx] / mapper.SimMeanArea[idx]:.2f} m^2/m^2'
                   ]
        ax.set_title('\n'.join(title_strs))
        fig.savefig(f'{save_dir}/MappedSubAreas_{fit_name}_{ID}.png', bbox_inches='tight', pad_inches=0)
        plt.close('all')

