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
import ParFlow_IO as pio
import VanG
import pars_ParFlowTCL as ppfl

"""
Example script to show how to use different methods of this repo.

Test-Case: compare tow datasets, as e.g. modifed run vs referenz run.

"""

def mappIndicator(ParFlowNamelist, IndicatorFile):
    ###############################################################################
    ### parse ParFlow-Namelist to get indicators values
    ###############################################################################
    PFNamelistDict = ppfl.pars_ParFlowNamelist(f'{ParFlowNamelist}')

    nz = int(PFNamelistDict['ComputationalGrid.NZ'])
    tcl_dz_keys = [f'Cell.{i}.dzScale.Value' for i in range(nz)]
    dz_mult = [float(PFNamelistDict[tcl_dz_key]) for tcl_dz_key in tcl_dz_keys]

    GeomInputs = PFNamelistDict['GeomInput.indinput.GeomNames']

    IndicatorInput = {GeomInput:int(PFNamelistDict[f'GeomInput.{GeomInput}.Value']) for GeomInput in GeomInputs}
    # print(f'IndicatorInput: {IndicatorInput}')
    vanGA_GeomNames = PFNamelistDict['Phase.Saturation.GeomNames']
    # print(f'vanGA_GeomNames: {vanGA_GeomNames}')
    tcl_vanGA_keys = [ ]
    vanG_a = {vanGA_GeomName:float(PFNamelistDict[f'Geom.{vanGA_GeomName}.Saturation.Alpha']) for vanGA_GeomName in vanGA_GeomNames}
    # print(f'vanG_a: {vanG_a}')
    vanG_n = {vanGA_GeomName:float(PFNamelistDict[f'Geom.{vanGA_GeomName}.Saturation.N']) for vanGA_GeomName in vanGA_GeomNames}
    vanG_sres = {vanGA_GeomName:float(PFNamelistDict[f'Geom.{vanGA_GeomName}.Saturation.SRes']) for vanGA_GeomName in vanGA_GeomNames}

    ###############################################################################
    ### prepare arrays to return
    ###############################################################################
    Indi3D = pio.read_pfb(f'{IndicatorFile}')
    shape3D= Indi3D.shape
    print(f'shape3D: {shape3D}')
    alpha  = np.full(shape3D,vanG_a['domain'])
    nvg    = np.full(shape3D,vanG_n['domain'])
    sres   = np.full(shape3D,vanG_sres['domain'])

    ###############################################################################
    ### mapp indicators to VanGenuchten parameter according to ParFlow-Namelist
    ###############################################################################
    for GeomName in vanGA_GeomNames:
        if GeomName == 'domain':
            continue
        alpha    = np.where(Indi3D == IndicatorInput[GeomName], vanG_a[GeomName], alpha)
        nvg      = np.where(Indi3D == IndicatorInput[GeomName], vanG_n[GeomName],   nvg)
        sres     = np.where(Indi3D == IndicatorInput[GeomName], vanG_sres[GeomName],  sres)

    return alpha, nvg, sres

###############################################################################
### Define some paths, filenames, etc and generate list of all files
###############################################################################
print(f'Define some paths, filenames, etc and generate list of all files')

rootDir = f'/p/scratch/cjibg35/tsmpforecast/development/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow'
rootdir_D1 = f'{rootDir}/postpro_skip'
rootdir_D2 = f'{rootDir}/postpro_Lin'
ParFlowNamelist = f'{rootDir}/ctrl/ParFlowCLM_DE05.tcl'
IndicatorFile   = f'{rootDir}/run/DE-0055_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_allv.pfb'

lvl = -1
varName = f'pressure'
filename = f'{varName}.nc'
season = 'DJF' # DJF, MAM, JJA, SO

files_D1 = {
        ## december, january, february
        'DJF':[f'{rootdir_D1}/1960_12/{filename}',
               f'{rootdir_D1}/1961_01/{filename}',
               f'{rootdir_D1}/1961_02/{filename}',
              ],
        ## march, apil, may
        'MAM':[f'{rootdir_D1}/1961_03/{filename}',
               f'{rootdir_D1}/1961_04/{filename}',
               f'{rootdir_D1}/1961_05/{filename}',
              ],
        ## june, july, august
        'JJA':[f'{rootdir_D1}/1961_06/{filename}',
               f'{rootdir_D1}/1961_07/{filename}',
               f'{rootdir_D1}/1961_08/{filename}',
              ],
        ## september, october
        'SO':[f'{rootdir_D1}/1961_09/{filename}',
              f'{rootdir_D1}/1961_10/{filename}',
             ],
        }
files_D1 = files_D1[season]
files_D2 = {
        ## december, january, february
        'DJF':[f'{rootdir_D2}/1960_12/{filename}',
               f'{rootdir_D2}/1961_01/{filename}',
               f'{rootdir_D2}/1961_02/{filename}',
              ],
        ## march, apil, may
        'MAM':[f'{rootdir_D2}/1961_03/{filename}',
               f'{rootdir_D2}/1961_04/{filename}',
               f'{rootdir_D2}/1961_05/{filename}',
              ],
        ## june, july, august
        'JJA':[f'{rootdir_D2}/1961_06/{filename}',
               f'{rootdir_D2}/1961_07/{filename}',
               f'{rootdir_D2}/1961_08/{filename}',
              ],
        ## september, october
        'SO':[f'{rootdir_D2}/1961_09/{filename}',
              f'{rootdir_D2}/1961_10/{filename}',
             ],
        }
files_D2 = files_D2[season]


###############################################################################
### read in data and save as ndarray
###############################################################################
print(f'read in data and save as ndarray')
for idx in range(len(files_D1)):
    print(f'idx: {idx}')
    with nc.Dataset(f'{files_D1[idx]}', 'r') as nc_file_D1, nc.Dataset(f'{files_D2[idx]}', 'r') as nc_file_D2:
        tmp_var1  = nc_file_D1.variables[varName][:,lvl,...]
        tmp_mask1 = np.ma.getmaskarray(tmp_var1)
        tmp_var1  = tmp_var1.filled(fill_value=np.nan)
        tmp_var1[tmp_var1 <= -1e+38] = np.nan
        
        tmp_var2  = nc_file_D2.variables[varName][:,lvl,...]
        tmp_mask2 = np.ma.getmaskarray(tmp_var2)
        tmp_var2  = tmp_var2.filled(fill_value=np.nan)
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
alpha, nvg, sres = mappIndicator(ParFlowNamelist=ParFlowNamelist, IndicatorFile=IndicatorFile)
satur_D1         = VanG.vanGenuchten(refP=var_D1, sSat=1.0, sRes=sres[lvl], nVanG=nvg[lvl], aVanG=alpha[lvl])
satur_D2         = VanG.vanGenuchten(refP=var_D2, sSat=1.0, sRes=sres[lvl], nVanG=nvg[lvl], aVanG=alpha[lvl])
print(f'satur_D1.shape: {satur_D1.shape}')
print(f'satur_D1 max: {np.nanmax(satur_D1)}')
print(f'satur_D2 min: {np.nanmin(satur_D2)}')
print(f'satur_D1 has inf: {np.isinf(satur_D1).any()}')
print(f'satur_D2 has inf: {np.isinf(satur_D2).any()}')
var_D1_ave = np.nanmean(satur_D1, axis=0)
var_D2_ave = np.nanmean(satur_D2, axis=0)

###############################################################################
### Define a special colormap
### seel also: https://matplotlib.org/stable/gallery/color/colormap_reference.html
###############################################################################
print(f'Define a special colormap')
cmap = copy.copy(mpl.cm.get_cmap('Blues'))
cmap.set_under('red')
cmap.set_over('magenta')
cmap_diff = mpl.cm.get_cmap('coolwarm_r')
cmap_diff.set_under('red')
cmap_diff.set_over('magenta')

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
        'title': '\n'.join(tmp_titlesubstr),
        'infostr': True,
        'v1_name': 'skip',
        'v2_name': 'Lin',
        'var_vmax': 1,
        'var_vmin': 0,
        'diff_vmax': 0.3,
        'diff_vmin': -0.3,
        'var_cmap': cmap,
        'diff_cmap': cmap_diff,
        'saveFile': f'./DE05_compare_{season}.pdf',
        #'saveFile': f'./Diff_{varName}_{cftime.datetime.strftime(dates_D1[0], "%Y%m%d%H%M")}.pdf',
        #'dpi': 400,
        'figsize': (10, 4),
        }

###############################################################################
### Start plotting
###############################################################################
print(f'Start plotting')
PlotLib.plot_imshow2PDiff(v1=var_D1_ave, v2=var_D2_ave,
        **kwargs_imshow2PDiff)
