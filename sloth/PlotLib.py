import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import datetime as dt
import cftime
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

def get_PlotMinMaxMid_Percentil(data2plot, lower=2, upper=98):
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
        var1_ax.text(0.01, 0.99, get_infostr(v1),
                     verticalalignment='top', transform=var1_ax.transAxes,
                     fontsize=8)
    im = var1_ax.imshow(v1, origin='lower', vmin=var_vmin, vmax=var_vmax,
            cmap=var_cmap, interpolation='none')

    var2_ax.set_title(fr'{v2_name}')
    var2_ax.get_yaxis().set_visible(False)
    if infostr:
        var2_ax.text(0.01, 0.99, get_infostr(v2),
                     verticalalignment='top', transform=var2_ax.transAxes,
                     fontsize=8)
    im = var2_ax.imshow(v2, origin='lower', vmin=var_vmin, vmax=var_vmax,
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

    plt.savefig(f'{saveFile}', bbox_inches='tight', pad_inches=0)

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

    plt.savefig(f'{saveFile}', bbox_inches='tight', pad_inches=0)

def plot_HeatWaveInvest(abs_intensity_array, dailyTime_hw, varName, pixel_x, pixel_y, hwYear,
    rel_intensity_array, dailyMean_hw, clima_smooth, clima_90p,
    clima, vmin, vmax):
# def plot_HeatWaveInvest(abs_intensity_array, dailyTime_hw, varName, pixel_x, pixel_y, hwYear,
#     rel_intensity_array, dailyMean_hw, clima_smooth, date_start_clim, date_final_clim, clima_90p,
#     clima, vmin, vmax):
    ###############################################################################
    #### Plot stuff
    ###############################################################################
    fig, ax = plt.subplots(figsize=(8,6))
    # makes a grid on the background
    ax.grid()
    dailyMean_hw_celsius=dailyMean_hw[:,pixel_y,pixel_x]-273.15        
    ax.plot(dailyMean_hw_celsius, color='black', label='daily '+str(varName), linewidth=2)

    ################### filtering the noise in daily climatological data ##########
    #b, a = signal.butter(3, 0.1, btype='lowpass', analog=False) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
    #low_passed = signal.filtfilt(b, a, noisy_signal)
    ##b, a = ellip(4, 0.01, 120, 0.125)  # Filter to be applied.
    #clima_smooth=signal.filtfilt(b, a, clima[:,pixel_y,pixel_x], padlen=30)
    clima_smooth_celsius=clima_smooth-273.15
    #clima_smooth=savgol_filter(clima[:,pixel_y,pixel_x], 11, 3) #savitzky_golay(clima[:,pixel_y,pixel_x], 51, 3) # window size 51, polynomial order 3
    ax.plot(clima_smooth_celsius, color='green', label=f'daily clim. {varName} (smoothed)',  linewidth=2)
    ax.plot(clima[:,pixel_y,pixel_x]-273.15, color='green', linestyle='dashed', label=f'daily clim. {varName} (original)',  linewidth=2)

    clima_90p_plot_celsius=clima_90p[:,pixel_y,pixel_x]-273.15
    #clima_smooth=savgol_filter(clima[:,pixel_y,pixel_x], 11, 3) #savitzky_golay(clima[:,pixel_y,pixel_x], 51, 3) # window size 51, polynomial order 3
    ax.plot(clima_90p_plot_celsius, color='darkviolet', linewidth=2, label=f'daily 5-day window 90th percentile of reference climatology')


    labels_idx = np.arange(0,dailyTime_hw.shape[0], 15) #NoI) #15
    ax.set_xticks(labels_idx)
    labels = [cftime.datetime.strftime(item, '%Y-%m-%d') for item in dailyTime_hw[labels_idx] ]
    ax.set_xticklabels(labels)
    x = np.arange(dailyTime_hw.shape[0])
      
    #ax.legend()  
    #ax.legend(loc='lower left',fancybox=True, shadow=True) # ncol=4)
    ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

    plt.title(f'TSMP {varName} [X={pixel_x}, Y={pixel_y} heat waves investigation: {hwYear}',fontsize=12,fontweight='bold')                
    plt.ylabel(f'Daily mean temperatures ({varName}), K', fontsize=11)
    plt.ylim(vmin, vmax)
    #plt.xlabel('Days',fontsize=11)  

    ax.fill_between(x, clima_smooth_celsius, dailyMean_hw_celsius, where=dailyMean_hw_celsius>clima_smooth_celsius, facecolor='lightcoral', interpolate=True, alpha=1.0) #lightsalmon
    ax.fill_between(x, clima_smooth_celsius, dailyMean_hw_celsius, where=dailyMean_hw_celsius<clima_smooth_celsius, facecolor='deepskyblue', interpolate=True, alpha=1.0) 
    ax.fill_between(x, clima_90p_plot_celsius, dailyMean_hw_celsius, where=dailyMean_hw_celsius>clima_90p_plot_celsius, facecolor='darkred', interpolate=True, alpha=1.0)

    plt.savefig('examples_DetectHeatwaves.pdf', dpi=380)
    plt.close('all')

    # fig, ax = plt.subplots(figsize=(8,3))
    # # makes a grid on the background
    # ax.grid()
    # ax.plot(abs_intensity_array, color='blue', label='Absolute intensity', linewidth=2)

    # labels_idx = np.arange(0,dailyTime_hw.shape[0], 15) #NoI) #15
    # ax.set_xticks(labels_idx)
    # labels = [cftime.datetime.strftime(item, '%Y-%m-%d') for item in dailyTime_hw[labels_idx] ]
    # ax.set_xticklabels(labels)
    # x = np.arange(dailyTime_hw.shape[0])
      
    # #ax.legend()  
    # #ax.legend(loc='lower left',fancybox=True, shadow=True) # ncol=4)
    # ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

    # plt.title(f'TSMP {varName} [X={pixel_x} Y={pixel_y}] heat waves investigation: {hwYear}',fontsize=12,fontweight='bold')                
    # plt.ylabel('Abs. intensity of heat wave, K', fontsize=11)
    # plt.savefig(f'Absolute intensity_{varName}_diff_hwYear{hwYear}.png', dpi=380) 


    # ###############################################################################
    # #### Plot stuff +1
    # ###############################################################################
    # fig, ax = plt.subplots(figsize=(8,3))
    # # makes a grid on the background
    # ax.grid()
    # ax.plot(rel_intensity_array, color='blue', linestyle = 'dashed', label='Relative intensity', linewidth=2)

    # labels_idx = np.arange(0,dailyTime_hw.shape[0], 15) #NoI) #15
    # ax.set_xticks(labels_idx)
    # labels = [cftime.datetime.strftime(item, '%Y-%m-%d') for item in dailyTime_hw[labels_idx] ]
    # ax.set_xticklabels(labels)
    # x = np.arange(dailyTime_hw.shape[0])
      
    # #ax.legend()  
    # #ax.legend(loc='lower left',fancybox=True, shadow=True) # ncol=4)
    # ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

    # plt.title('TSMP '+str(varName) +' [X=' + str(pixel_x) +', Y=' + str(pixel_y) +']' +' heat waves investigation: '+str(hwYear),fontsize=12,fontweight='bold')                
    # plt.ylabel('Rel. intensity of heat wave', fontsize=11)
    # plt.savefig(f'Relative intensity_{varName}_diff_hwYear{hwYear}.png', dpi=380) 