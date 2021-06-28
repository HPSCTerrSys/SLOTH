import numpy as np
import netCDF4 as nc
import glob
import os
import sys
import heat as ht

src_path='../src/'
sys.path.append(src_path)
# ANalysisTool does contain the function 'get_intervalSlice()'
import ANalysisTool as ANT
import SanityCheck 
src_path='../extern/ana_parflow-diagnostics_pythonheat'
sys.path.append(src_path)
# ANalysisTool does contain the function 'get_intervalSlice()'
import Diagnostics
import IO as htio


###############################################################################
### Define some paths, filenames, etc
###############################################################################
# NEEDED FOR HeAT
split = None

dataRootDir   = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples'
datasetName   = 'ERA5Climat_EUR11_ECMWF-ERA5_analysis_FZJ-IBG3'
procType      = 'postpro'
dataDate      = '1980_02'
pressVarName  = 'press'
etVarName     = 'et'
pressFiles    = [
          f'{dataRootDir}/{datasetName}/{procType}/1984_01/{pressVarName}.nc',
          f'{dataRootDir}/{datasetName}/{procType}/1984_02/{pressVarName}.nc'
          ]
etFiles    = [
          f'{dataRootDir}/{datasetName}/{procType}/1984_01/{etVarName}.nc',
          f'{dataRootDir}/{datasetName}/{procType}/1984_02/{etVarName}.nc'
          ]

maskFile      = f'{dataRootDir}/{datasetName}/{procType}/1980_01/mask.nc'
# specific_storage
sstoragefile  = f'{dataRootDir}/{datasetName}/tmp_static/specific_storage.pfb'
# porosity 
# should be possible to get this from indicator files and namelist...!
porofile      = f'{dataRootDir}/{datasetName}/tmp_static/porosity.pfb'
permFileZ     = f'{dataRootDir}/{datasetName}/tmp_static/perm_z.pfb'
permFileY     = f'{dataRootDir}/{datasetName}/tmp_static/perm_y.pfb'
permFileX     = f'{dataRootDir}/{datasetName}/tmp_static/perm_x.pfb'
pflname       = f'{dataRootDir}/{datasetName}/namelists/coup_oas.tcl'
indicatorfile = f'{dataRootDir}/{datasetName}/geo/parflow/ParFlow_SOIL_INDICATOR3_from_EUR03_x1592y1544z15_on_EUR11_x436y424z15.pfb'
slopeFileX = f'{dataRootDir}/{datasetName}/geo/parflow/slopex.pfb'
slopeFileY = f'{dataRootDir}/{datasetName}/geo/parflow/slopey.pfb'

press = []
for pressFile in pressFiles:
    print(f'ht.load_netcdf')
    tmp_press = ht.load_netcdf(pressFile, split=split, variable=pressVarName)
    press.append(tmp_press)
press = ht.concatenate(press, axis=0)
print(f'press.shape: {press.shape}')

et = []
for etFile in etFiles:
    print(f'ht.load_netcdf')
    tmp_et = ht.load_netcdf(etFile, split=split, variable=etVarName)
    et.append(tmp_et)
et = ht.concatenate(et, axis=0)
print(f'et.shape: {et.shape}')

mask = ht.load_netcdf(maskFile, split=split, variable='mask')[0]
# strange mask, but should be correct for sample data-set
# 99999 is land and other is lake / sea
# needed by Diagnostics:
# 0 = pixel to mask (sea)
# 1 = pixel to keep (land)
mask = ht.where(mask==99999,1,0)
print(f'mask.shape: {mask.shape}')

sstorage = htio.read_pfb(sstoragefile, split=split)
print(f'sstorage.shape: {sstorage.shape}')
# SanityCheck.plot_SanityCheck_3D(data=sstorage.numpy(),
#     kind='mean', data_mask=mask, 
#     fig_title='specific storage (z,y,x)', minax_title='min', maxax_title='max', 
#     kinax_title='mean', cmapName='Spectral')
# sys.exit()

poro = htio.read_pfb(porofile, split=split)
print(f'poro.shape: {poro.shape}')

permz = htio.read_pfb(permFileZ, split=split)
permy = htio.read_pfb(permFileY, split=split)
permx = htio.read_pfb(permFileX, split=split)
print(f'permz.shape: {permz.shape}')
print(f'permy.shape: {permy.shape}')
print(f'permx.shape: {permx.shape}')

slopex = htio.read_pfb(slopeFileX, split=split)
slopey = htio.read_pfb(slopeFileY, split=split)
print(f'slopex.shape: {slopex.shape}')
print(f'slopey.shape: {slopey.shape}')



print(f'About reading pfl-namelist')
indicatorMap = ANT.toolBox.mappIndicator(ParFlowNamelist=pflname, IndicatorFile=indicatorfile)
alpha   = ht.array(indicatorMap['alpha'], is_split=0, comm=ht.MPI_WORLD)
nvg     = ht.array(indicatorMap['nvg'], is_split=0, comm=ht.MPI_WORLD)
sres    = ht.array(indicatorMap['sres'], is_split=0, comm=ht.MPI_WORLD)
dzmult = ht.array(indicatorMap['dz_mult'], is_split=0, comm=ht.MPI_WORLD)
dz      = indicatorMap['dz']
dy      = indicatorMap['dy']
dx      = indicatorMap['dx']
nz      = indicatorMap['nz']
ny      = indicatorMap['ny']
nx      = indicatorMap['nx']

perm = ht.zeros((3,nz,ny,nx),split=split)
perm[0]=permz
perm[1]=permy
perm[2]=permx

# 3h dump
# or
# 0.25 dt intern ParFlow?
dt       = ht.float64(3)
mannings = ht.full(alpha.shape, 5.5e-5, split=split)

diag = Diagnostics.Diagnostics(Mask=mask, Perm=perm,
    Poro=poro,
    Sstorage=sstorage,
    Ssat=1., Sres=sres, Nvg=nvg, Alpha=alpha,
    Mannings=mannings, Slopex=slopex, Slopey=slopey,
    Dx=dx, Dy=dy, Dz=dz, 
    Dzmult=dzmult, 
    Nx=nx, Ny=ny, Nz=nz,
    Terrainfollowing=True, Split=split)

print(f'about to calculate satur and subsurfaceStorage')
satur = ht.full(press.shape,0,split=split)
subsurfaceStorage = ht.full(press.shape,0,split=split)
for t in range(press.shape[0]):
    print(f'handling t: {t}')
    satur[t], krel = diag.VanGenuchten(press[t])

    satur[t] = ht.where(mask!=0,satur[t],0.0)
    press[t] = ht.where(mask!=0,press[t],99999.0)

    # SanityCheck.plot_SanityCheck_3D(data=var, 
    #     # below is optional
    #     data_mask=var_mask, kind='mean', figname=figname,
    #     lowerP=2, upperP=98, interactive=False,
    #     fig_title=fig_title, minax_title=minax_title, maxax_title=maxax_title, 
    #     kinax_title=kinax_title, hisax_title=hisax_title, cmapName='Spectral_r')

    #subsurfaceStorage[t] = diag._SubsurfaceStorage(press[t], satur[t])

    # #Calculate relative saturation and relative hydraulic conductivity
    # dummy,krel = diag.VanGenuchten(press)

    # #Work with the ParFlow output for now, ...
    # satur = io.read_pfb(path + name + '.out.satur.'+ ('{:05d}'.format(t)) + '.pfb',split=split)
    # satur = ht.where(mask==1.0,satur,0.0)
    # #..., because something is wrong with the generated van Genuchten fields
    # #print(ht.sum((satur-satur_)*mask))

    #Obtain pressure at the land surface
    top_layer_press = diag.TopLayerPressure(press[t])

    #Returns an unmasked 3D field of subsurface storage, (L^3)
    subsurface_storage=diag.SubsurfaceStorage(press[t],satur[t])

    #Returns an unmasked 2D field of surface storage, (L^3)
    surface_storage = diag.SurfaceStorage(top_layer_press)

    #Calculate subsurface flow in all 6 directions for each grid cell (L^3/T)
    flowleft,flowright,flowfront,flowback,flowbottom,flowtop = diag.SubsurfaceFlow(press[t],krel)

    #Calculate overland flow (L^2/T)
    oflowx,oflowy = diag.OverlandFlow(top_layer_press)    
    oflowx = ht.where(mask[0]!=0,oflowx,0)
    oflowy = ht.where(mask[0]!=0,oflowy,0)

    #Calculate net overland flow for each top layer cell (L^3/T)
    net_overland_flow = diag._NetLateralOverlandFlow(oflowx,oflowy)

    #Column balance
    if t > 0:
      #Read source/sink values coming from CLM (m/h)
      sink = et[t]
      sink = ht.where(mask!=0,sink,0.0)
      #Convert into (L^3)
      for k in range (nz):
          sink[k,:,:] *= dz * dzmult[k] * dy * dx * dt


      #Source/sink integrated over column (L^3)
      sourcesink = ht.sum(sink, axis=0) 

      #Change in subsurface storage for each cell (L^3)
      dstorage_cell = old_subsurface_storage - subsurface_storage

      #Change in surface storage for each surface cell (L^3)
      dsurface_storage_cell = old_surface_storage - surface_storage

      #Surface balance (L^3)
      balance_surface = dsurface_storage_cell - dt * net_overland_flow

      #Divergence of the flux for each cell
      divergence_cell = dt * (flowleft-flowright+flowfront-flowback+flowbottom-flowtop)

      #Balance for each cell
      balance_cell = dstorage_cell + divergence_cell

      #Change in storage for each column
      dstorage_column = ht.sum(dstorage_cell*mask,axis=0)

      #Divergence of the flux for each column
      divergence_column = ht.sum(divergence_cell*mask,axis=0)

      #Balance for each column
      balance_column  = ht.sum(balance_cell*mask,axis=0)
      balance_column += balance_surface
      balance_column += sourcesink

      #Mass balance over full domain without flux at the top boundary
      print('Time step:',t, ', dstorage:',ht.sum(dstorage_column))
      print('Time step:',t, ', divergence:',ht.sum(divergence_column))
      print('Time step:',t, ', dsurface_storage:',ht.sum(dsurface_storage_cell))
      print('Time step:',t, ', netoverlandflow:',ht.sum(net_overland_flow))
      print('Time step:',t, ', surface_balance:',ht.sum(balance_surface))
      print('Time step:',t, ', source/sink:',ht.sum(sourcesink))
      print('Time step:',t, ', total balance:',ht.sum(balance_column))

    #New becomes old in the ensuing time step
    old_subsurface_storage = subsurface_storage
    old_surface_storage = surface_storage


# PLUS store as netCDF
