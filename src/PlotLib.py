import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import datetime as dt
import sys

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

def plot_imshow2PDiff(v1, v2, **kwargs):
    '''
    generates 3 imshows in a rwo
    plot 2 vars and difference of both
    '''
    title       = kwargs.pop('title', None)
    infostr     = kwargs.pop('infostr', False)
    v1_name     = kwargs.pop('v1_name', 'Variable 1')
    v2_name     = kwargs.pop('v2_name', 'Variable 2')
    var_vmin    = kwargs.pop('var_vmin', np.nanmin(v1))
    var_vmax    = kwargs.pop('var_vmax', np.nanmax(v1))
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
    saveFile    = kwargs.pop('saveFile', './DefaultSaveName_plot_imshow2pDiff.png')
    dpi         = kwargs.pop('dpi', 100)
    figsize     = kwargs.pop('figsize', (8.27, 11.69))

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
        tmp_infostr = [
                f'min: {np.nanmin(v1):.2e}',
                f'max: {np.nanmax(v1):.2e}',
                f'mean: {np.nanmean(v1):.2e}',
                f'median: {np.nanmedian(v1):.2e}',
                ]
        var1_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                     verticalalignment='top', transform=var1_ax.transAxes,
                     fontsize=8)
    im = var1_ax.imshow(v1, origin='lower', vmin=var_vmin, vmax=var_vmax,
            cmap=var_cmap, interpolation='none')

    var2_ax.set_title(fr'{v2_name}')
    var2_ax.get_yaxis().set_visible(False)
    if infostr:
        tmp_infostr = [
                f'min: {np.nanmin(v2):.2e}',
                f'max: {np.nanmax(v2):.2e}',
                f'mean: {np.nanmean(v2):.2e}',
                f'median: {np.nanmedian(v2):.2e}',
                ]
        var2_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                     verticalalignment='top', transform=var2_ax.transAxes,
                     fontsize=8)
    im = var2_ax.imshow(v2, origin='lower', vmin=var_vmin, vmax=var_vmax,
            cmap=var_cmap, interpolation='none')
    cb = plt.colorbar(im, cax=var_cax, orientation='vertical',
            extend='both')

    diff_ax.set_title(fr'Diff: {v1_name} - {v2_name}')
    diff_ax.get_yaxis().set_visible(False)
    if infostr:
        tmp_infostr = [
                f'min: {np.nanmin(v1-v2):.2e}',
                f'max: {np.nanmax(v1-v2):.2e}',
                f'mean: {np.nanmean(v1-v2):.2e}',
                f'median: {np.nanmedian(v1-v2):.2e}',
                ]
        diff_ax.text(0.01, 0.99, '\n'.join(tmp_infostr),
                     verticalalignment='top', transform=diff_ax.transAxes,
                     fontsize=8)
    im = diff_ax.imshow(v1 - v2, origin='lower',
            vmin=diff_vmin, vmax=diff_vmax, 
            cmap=diff_cmap, interpolation='none')
    cb = plt.colorbar(im, cax=diff_cax, orientation='vertical',
            extend='both')

    plt.savefig(f'{saveFile}', bbox_inches='tight', pad_inches=0)

