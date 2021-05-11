import heat as ht
import sys 
diagnostics_path='/p/home/jusers/gordeev1/juwels/ana_parflow-diagnostics_pythonheat'
sys.path.append(diagnostics_path)
import numpy as np
import os
import math
from Diagnostics import Diagnostics
import IO as io
import matplotlib.pyplot as plt


def read_pfb_batch(folder,prefix,*fnames,mask=None,default_val=None):
    """
    DESCRIPTION
               reads/masks .pfb files - input variables for Diagnostics class             
    INPUT
               folder - path of folder with pfb files
               prefix - model/month name e.g. cordex0.11_1980_06
               mask   - ht.array with 1 in active and 0 in inactive regions
               default_val - value assigned in inactive areas
               fnames - names of files to read, e.g. '.out.mask.pfb','.out.porosity.pfb'
    OUTPUT
               package - list containing ht.arrays, corresponding to *fnames 
    """
    package=list()
    if mask is not None:
        for fname in fnames:
            unmasked=io.read_pfb(folder+prefix+fname,split=None)
            masked=ht.where(mask==1.0,unmasked,default_val)
            package.append(masked)
        
        return package

    else:
        for fname in fnames:
            var=io.read_pfb(folder+prefix+fname,split=None)
            package.append(var)
        
        return package


def iterate(press_unpacked,satur_unpacked,mask,diag,nt):
    """
    DESCRIPTION
               iterate over monthly timesteps, apply Diagnostics class methods on every timestep
    INPUT
               press/satur_unpacked - 4D (nt,nz,ny,nx) arrays with ERA5 WU simulation output
               mask                 - ht.array with 1 in active and 0 in inactive regions
               diag                 - Diagnostics class object
               nt                   - int,number of timesteps within current month
    OUTPUT
               subsurf,toplayer,surf - lists with nt elements, each element is 2D/3D heat array,2D/3D heat array,diagnostics method output. 
    """
    subsurf=[]
    toplayer=[]
    surf=[]
    for t in range(nt):
        print("proscessing timestep {}".format(t))
        # select 3D DNDarray of shape (nz,ny,nx)
        # pressure
        press = press_unpacked[t,:,:,:]
        press = ht.where(mask>0.0,press,99999.0)

        # satur
        satur = satur_unpacked[t,:,:,:]
        satur = ht.where(mask==1.0,satur,0.0)

        #Returns an unmasked 3D field of subsurface storage, (L^3)
        subsurface_storage=diag.SubsurfaceStorage(press,satur)
        subsurf.append(subsurface_storage)
        
        #Obtain pressure at the land surface
        toplayer_press = diag.TopLayerPressure(press)
        toplayer.append(toplayer_press)
        #Returns an unmasked 2D field of surface storage, (L^3)
        surface_storage = diag.SurfaceStorage(toplayer_press)
        surf.append(surface_storage)

    return subsurf,toplayer,surf


def htarr_plot(heatarr,variable,timestep,zlayer=None,vmin=None):
    """
    DESCRIPTION
               imshow 2D/3D heat arrays
    INPUT
               heatarr - 2D/3D heat array of shape (nz,ny,nx) or (ny,nx)
               variable - string,variable name to show on the plot
               timestep - int,timestep  to show on the plot
    OUTPUT
               imshow plot
    """

    if heatarr.ndim == 3:
        # extract data from xarray dataset
        extr_surface = np.array(heatarr[zlayer],dtype=float)

    elif heatarr.ndim ==2:
        # extract data from xarray dataset
        extr_surface =np.array(heatarr,dtype=float)

    else:
        print("Error, false dimensionality of input heatarr")
        return None

    # plotting itself
    fig=plt.figure(figsize=(10,8))
    plt.imshow(extr_surface[::-1,:],interpolation="nearest",vmin=vmin)
    plt.title("{} timestep={}, zlayer={}".format(variable,timestep,zlayer))
    plt.colorbar()
    plt.show()


