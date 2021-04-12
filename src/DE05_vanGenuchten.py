import numpy as np
import heat as ht
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys

src_path='../src'
sys.path.append(src_path)
import ParFlow_IO as pio
import VanG
import pars_ParFlowTCL as ppfl

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

    return alpha, nvg, sres, Indi3D

rootDir = f'/p/scratch/cjibg35/tsmpforecast/development/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow'
ParFlowNamelist = f'{rootDir}/ctrl/ParFlowCLM_DE05.tcl'
IndicatorFile   = f'{rootDir}/run/DE-0055_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_allv.pfb'

lvl = -1
alpha, nvg, sres, Indi3D = mappIndicator(ParFlowNamelist=ParFlowNamelist, IndicatorFile=IndicatorFile)


fig = plt.figure(figsize=(6,6), dpi=100)
gs  = fig.add_gridspec(nrows=10,ncols=10)
gs.update(wspace=0.3, hspace=0.4)

cmap = mpl.cm.get_cmap('gist_rainbow')
cmap_indi3d = mpl.cm.get_cmap('tab20')

var1_ax = fig.add_subplot(gs[0:4,0:4])
var2_ax = fig.add_subplot(gs[0:4,5:9])
var3_ax = fig.add_subplot(gs[5:9,0:4])
var4_ax = fig.add_subplot(gs[5:9,5:9])
var1_cax = inset_axes(var1_ax,
           width="5%",  # width = 5% of parent_bbox width
           height="100%",  # height : 340% good for a (4x4) Grid
           loc='lower left',
           bbox_to_anchor=(1.05, 0., 1, 1),
           bbox_transform=var1_ax.transAxes,
           borderpad=0,
           )
var2_cax = inset_axes(var2_ax,
           width="5%",  # width = 5% of parent_bbox width
           height="100%",  # height : 340% good for a (4x4) Grid
           loc='lower left',
           bbox_to_anchor=(1.05, 0., 1, 1),
           bbox_transform=var2_ax.transAxes,
           borderpad=0,
           )
var3_cax = inset_axes(var3_ax,
           width="5%",  # width = 5% of parent_bbox width
           height="100%",  # height : 340% good for a (4x4) Grid
           loc='lower left',
           bbox_to_anchor=(1.05, 0., 1, 1),
           bbox_transform=var3_ax.transAxes,
           borderpad=0,
           )
var4_cax = inset_axes(var4_ax,
           width="5%",  # width = 5% of parent_bbox width
           height="100%",  # height : 340% good for a (4x4) Grid
           loc='lower left',
           bbox_to_anchor=(1.05, 0., 1, 1),
           bbox_transform=var4_ax.transAxes,
           borderpad=0,
           )

var1_ax.set_title(fr'alpha')
im = var1_ax.imshow(alpha[lvl], origin='lower',
        cmap=cmap, interpolation='none')
cb = plt.colorbar(im, cax=var1_cax, orientation='vertical',
        extend='both')

var2_ax.set_title(fr'nvg')
var2_ax.get_yaxis().set_visible(False)
im = var2_ax.imshow(nvg[lvl], origin='lower',
        cmap=cmap, interpolation='none')
cb = plt.colorbar(im, cax=var2_cax, orientation='vertical',
        extend='both')

var3_ax.set_title(fr'sres')
var3_ax.get_yaxis().set_visible(False)
im = var3_ax.imshow(sres[lvl], origin='lower',
        cmap=cmap, interpolation='none')
cb = plt.colorbar(im, cax=var3_cax, orientation='vertical',
        extend='both')

var4_ax.set_title(fr'Indi3D')
var4_ax.get_yaxis().set_visible(False)
im = var4_ax.imshow(Indi3D[lvl], origin='lower',
        cmap=cmap_indi3d, interpolation='none')
cb = plt.colorbar(im, cax=var4_cax, orientation='vertical',
        extend='both')

plt.savefig(f'vanGenuchten_lvl{lvl}.pdf', bbox_inches='tight', pad_inches=0)
