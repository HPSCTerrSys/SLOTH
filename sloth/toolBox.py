""" toolBox - submodule of cathyNAME

author: Nikals WAGNER
e-mail: n.wagner@fz-juelich.de
version: 2021-05-11

Description:
tooBox.py is aimed to hold standalone functions used / developed / found
within the realm of cathyNAME development.
[...]

Check examples/ to see how this works. 
"""
import sys
import os
import numpy as np
import netCDF4 as nc
import datetime as dt
from . import pars_ParFlowTCL as ppfl
from . import ParFlow_IO as pio
import matplotlib.pyplot as plt

def calc_catchment(slopex, slopey, x, y):
    """ This function calculates the catchment related to a given pixel
    based on x- and y-slope files

    This algorithem is 

    INPUT: 
    slopex: 2D ndarray
        slopes in x-direction
    slopey: 2D ndarray
        slopes in y-direction
    x: int
        index in x-direction to calulate catchment from 
    y: int
        index in y-direction to calulate catchment from

    RETURN:
    catchment: 2D ndarray 
        ndarray of the same size as slopex/y. 0 = not part of catchment; 1 = part of catchment
    """
    dims = slopey.shape
    nx = dims[1]
    # convert slopes to: 'in which 1D index I do drain'
    drainx = np.zeros_like(slopex)
    drainx = np.where(slopex<0, 1, drainx)
    drainx = np.where(slopex>0, -1, drainx)
    drainx = np.where(slopex==0, -999, drainx)
    drainy = np.zeros_like(slopey)
    drainy = np.where(slopey<0, nx, drainy)
    drainy = np.where(slopey>0, -nx, drainy)
    drainy = np.where(slopey==0, -999, drainy)
    # plt.imshow(drainx)
    # plt.show()
    # FlatDrainX
    fdx = drainx.ravel(order='C')
    # FlatDrainY
    fdy = drainy.ravel(order='C')
    # FlatCatchment
    fc = np.zeros_like(fdx)
    # flat_idx order='C': idx = (y * nx) + x
    start = (y * nx ) + x

    openEnds = [start]
    while openEnds:
        # get one open end
        step = openEnds.pop()
        # print(f'poped: {step}')
        # print(f'len openEnds: {len(openEnds)}')
        # mark open end as catchment
        fc[step] = 1
        # GET surrounding Pixel
        # Nort = step - nx; East = step +1
        # West = step - 1 ;South = step + nx
        # NEWS = [step-nx, step+1, step-1, step+nx]
        NS = [step-nx, step+nx]
        EW = [step+1, step-1]
        # print(f'NS: {NS}')
        # print(f'EW: {EW}')
        # CHECK if surrounding drain to step (D2S) and are NOT catchment already
        try:
            NSD2S = [ idx for idx in  NS if (idx + fdy[idx] == step and not fc[idx] ) ] 
            EWD2S = [ idx for idx in  EW if (idx + fdx[idx] == step and not fc[idx] ) ] 
        except IndexError:
            print('FEHLER')
            continue
        # print(f'NSD2S: {NSD2S}')
        # print(f'EWD2S: {EWD2S}')
        D2S = NSD2S + EWD2S
        # print(f'D2S: {D2S}')

        # add all found pixes to openEnds
        openEnds += D2S

    # np.save('./catchment_mask', fc.reshape(dims))
    return fc.reshape(dims)

def calc_catchment_old_nonzeroslope(slopex, slopey, x, y):
    cmask = np.zeros_like(slopex)
    cdim = cmask.shape
    cmask[y,x] = 1
    loop = 0
    max_loop = 10000
    new_p = 1
    while (new_p != 0) and (loop < max_loop):
        print('start loop {} (found {} new points in last loop)'.format(loop, new_p))
        new_p = 0
        for y in range(1, cdim[0]-1):
            for x in range(1, cdim[1]-1):
                if(cmask[y,x] == 0):
                    slopdir = np.abs(slopex[y,x]) >= np.abs(slopey[y,x])
                    if slopdir:#RUNOFF IN Y-DIRECTION
                        if( cmask[y,x+1] == 1 and slopex[y,x] <= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                        elif( cmask[y,x-1] == 1 and slopex[y,x] >= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                    else:
                        if( cmask[y+1,x] == 1 and slopey[y,x] <= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                        elif( cmask[y-1,x] == 1 and slopey[y,x] >= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue

                else:
                    continue
        loop += 1
    np.save('./catchment_mask_old_nonzeroslope', cmask)

def calc_catchment_old_majorslope(slopex, slopey, x, y):
    cmask = np.zeros_like(slopex)
    cdim = cmask.shape
    cmask[y,x] = 1
    loop = 0
    max_loop = 10000
    new_p = 1
    while (new_p != 0) and (loop < max_loop):
        print('start loop {} (found {} new points in last loop)'.format(loop, new_p))
        new_p = 0
        for y in range(1, cdim[0]-1):
            for x in range(1, cdim[1]-1):
                if(cmask[y,x] == 0):
                    slopdir = np.abs(slopex[y,x]) >= np.abs(slopey[y,x])
                    if slopdir:#RUNOFF IN Y-DIRECTION
                        if( cmask[y,x+1] == 1 and slopex[y,x] <= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                        elif( cmask[y,x-1] == 1 and slopex[y,x] >= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                    else:
                        if( cmask[y+1,x] == 1 and slopey[y,x] <= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                        elif( cmask[y-1,x] == 1 and slopey[y,x] >= 0 ):
                            cmask[y,x] = 1
                            new_p += 1
                            continue
                else:
                    continue
        loop += 1
    np.save('./catchment_mask_old_majorslope', cmask)

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

def get_intervalSlice(dates, sliceInterval='month'):
    ''' This functions calculates interval slices of a given time-series

    This function calculates interval slices of a given time-series, 
    aiming to mimic pythons default slice notation `array[start:stop:step]`
    to also enable multi dimensional slicing. 
    The calculation takes the models dumpintervall into account, wherefore 
    this function is working for (nearly) every time-resolution.
    The calculation is based on datetime-objects and therefore really 
    calculates the slice of individual interval, even if the time-series 
    does not start or end at the first or last of a given interval.

    Use case:
    You got a time-series of hourly data points and want to calculate the
    monthly mean. This functions returns the slices needed to divide your
    data into monthly peaces.

    Return value:
    -------------
    Slices: list
        A list with length of found intervals with related slices

    This function assumes:
    ----------------------
    -) The passed time-series covers at least on intervall!
    -) At least hourly steps / dumpIntervals! 
    -) No seconds in dates - model output is at least on full minutes!

    Parameters
    ----------
    dates : NDarray 
        the time-series as datetime-object. 
    sliceInterval: str
        defining the interval. Supported are: 'day', 'month'
    '''

    # Check is passed sliceInterval is supported and exit if not.
    supportedSliceInterval = ['day', 'month']
    if sliceInterval not in supportedSliceInterval:
        print(f'ERROR: sliceInterval "{sliceInterval}" is not supported.')
        print(f'---    supported values are: {supportedSliceInterval}')
        print(f'---    EXIT program')
        return False

    print(f'########################################################################')
    print(f'#### calculate dumpInterval and check if equal for all data-points')
    print(f'########################################################################')
    # Calculating dumpInterval
    tmp_dumpInterval = np.diff(dates, n=1)
    # Check if dumpInterval equal for all data-points
    # and if so define set dumpInterval 
    if not np.all(tmp_dumpInterval == tmp_dumpInterval[0]):
        print('ERROR: calculated dumpInterval is not equal for all data-points')
        # In case of error: break function
        return False
    else:
        dumpInterval = tmp_dumpInterval[0]
        print(f'DONE --> DumpInterval: {dumpInterval}')        
    
    # Finding the first time-step of a interval, by looping over the time-series and 
    # check for each time-step if this is the first of the interval (break loop if found).
    print(f'########################################################################')
    print(f'#### Finding first of month')
    print(f'########################################################################')
    Offset = 0
    while True:
        # verify first date is first of interval.
        # idea explained with sliceInterval='month': 
        # first step - dumpInterval (or 0.5*dumpInterval) is first of month 00UTC
        # By checking 1*dumpInterval and 0.5*dumpInterval we do take into account, that
        # the time-axis of averaged model output might got shifted in between the 
        # time-bounds.
        print(f'Offset: {Offset}')
        tmp_first = dates[Offset]
        if sliceInterval == 'day':
            # using midnight as reference
            tmp_first_ref = tmp_first.replace(hour=0, minute=0, second=0)
        elif sliceInterval == 'month':
            # using the first of current month at midnight as reference
            tmp_first_ref = tmp_first.replace(day=1, hour=0, minute=0, second=0)

        # First time-step of dates is already first of a interval
        if (tmp_first == tmp_first_ref):
            print(f'check step {Offset} is first of a month at midnight')
            break
        #if (tmp_first - dumpInterval) == tmp_first_ref:
        #    print(f'check step {Offset} is first of a month at midnight')
        #    break
        # First time-step is not the exact first of a interval but set in
        # middle of time-bounds and therefore halfe a dumpInteral away from
        # first date of a interval
        elif (tmp_first - (0.5*dumpInterval)) == tmp_first_ref:
            print(f'check step {Offset} is first of a month at midnight')
            break
        else:
            print(f'ERROR: step {Offset} is not first step of month at midnight!')
            print(f'step: {tmp_first}')
            Offset += 1
            # hard break of loop if no 'beginning' is found after 50 time-steps
            if Offset > 50:
                break
            continue

    # Calculating the slices by looping over all time-steps and check if the 
    # next time-step belong to the next interval.
    # NWR 20210422:
    # there should be some clever and short solution for below loop ...!
    # now that we know the dumpInterval and first step is first of interval...
    print(f'########################################################################')
    print(f'#### getting month series / slices')
    print(f'########################################################################')
    t_lower = Offset
    Slices = []
    for t in range(Offset, dates.shape[0]):
        # I decided to go for the solution checking current month and next month
        # to catch the case if the dateset contains one month only!
        if sliceInterval == 'day':
            currInterval = dates[t].day
            nextInterval = (dates[t] + dumpInterval).day
        elif sliceInterval == 'month':
            currInterval = dates[t].month
            nextInterval = (dates[t] + dumpInterval).month
        if nextInterval != currInterval:
            Slices.append(slice(t_lower,t+1,None))
            t_lower = t+1

    return Slices 


def createNetCDF(fileName, domain=None, nz=None, author=None, 
    description=None, source=None, contact=None, institution=None, 
    history=None, timeCalendar=None, timeUnit=None, NBOUNDCUT=0):

    #######################################################################
    #### Default definitions for domain
    #######################################################################
    if domain == 'DE05':
        # https://icg4geo.icg.kfa-juelich.de/SoftwareTools/prepro_parflowclm_de05_static_fields/blob/master/doc/content/Grid.rst
        ny = nx = 2000     # [-]
        dy = dx = 0.0055   # [deg]
        # Rotated pole geographical coordinates
        rpol_Y  = 39.25    # [deg]
        rpol_X  = -162.0   # [deg]
        # South-West Corner of domain in rotated coordinates
        SWC_Y = -5.38725   # [deg]
        SWC_X = -10.82725  # [deg]
    elif domain == 'EU11_TSMP':
        ny = 432           # [-]
        nx = 444           # [-]
        dy = dx = 0.110    # [deg]
        # Rotated pole geographical coordinates
        rpol_Y  = 39.25    # [deg]
        rpol_X  = -162.0   # [deg]
        # South-West Corner of domain in rotated coordinates
        SWC_Y = -24.4720   # [deg]
        SWC_X = -29.4713   # [deg]
    elif domain == 'EU11': #################LPo, 25.06.2021!!!!!!!!! check if it is correct!!!!!!!!!!
        ny = 412           # [-]
        nx = 424           # [-]
        dy = dx = 0.110    # [deg]
        # Rotated pole geographical coordinates
        rpol_Y  = 39.25    # [deg]
        rpol_X  = -162.0   # [deg]
        # South-West Corner of domain in rotated coordinates
        SWC_Y = -23.375   # [deg] rlat
        SWC_X = -28.375   # [deg] rlon
    else:
        print(f'ERROR: passed domain is not supported. domain={domain} --> Exit')
        return False

    #######################################################################
    #### Checking is time- and / or z-axis is used
    #######################################################################
    withTime = False
    if timeUnit is not None and timeCalendar is not None:
        withTime = True
    else:
        print('NOT creating time-axis')
        print(f'--  timeUnit = {timeUnit}; timeCalendar = {timeCalendar}') 

    withZlvl = False
    if nz is not None:
        withZlvl = True
    else:
        print('NOT creating z-axis')

    #######################################################################
    #### Create netCDF file (overwrite if exist)
    #######################################################################
    # If no file-extension is passed, add '.nc' as default
    fileRoot, fileExt = os.path.splitext(fileName)
    if not fileExt:
       fileExt = '.nc'
    fileName = f'{fileRoot}{fileExt}'

    # Create netCDF file
    nc_file = nc.Dataset(f'{fileName}', 'w', format='NETCDF4')
    # Add basic information
    nc_file.author      = f'{author}'
    nc_file.contact     = f'{contact}'
    nc_file.institution = f'{institution}' 
    nc_file.description = f'{description}'
    nc_file.history     = f'{history}'
    nc_file.source      = f'{source}'

    # Create dimensions
    # Take into account to 'cut' pixel at domain border (NBOUNDCUT)
    drlon = nc_file.createDimension('rlon',nx-2*NBOUNDCUT)
    drlat = nc_file.createDimension('rlat',ny-2*NBOUNDCUT)
    if withZlvl:
        dlvl = nc_file.createDimension('lvl',nz)
    dtime = nc_file.createDimension('time',None)

    rlon = nc_file.createVariable('rlon', 'f4', ('rlon',),
                                zlib=True)
    rlon.standard_name = "grid_longitude"
    rlon.long_name = "rotated longitude"
    rlon.units = "degrees"
    rlon.axis = "X"
    # Take into account to 'cut' pixel at domain border (NBOUNDCUT)
    rlon_values = np.array([SWC_X + (i*dx) for i in range(NBOUNDCUT, nx-NBOUNDCUT)])
    rlon[...] = rlon_values[...]

    rlat = nc_file.createVariable('rlat', 'f4', ('rlat',),
                                    zlib=True)
    rlat.standard_name = "grid_latitude"
    rlat.long_name = "rotated latitude"
    rlat.units = "degrees"
    rlat.axis = "Y"
    # Take into account to 'cut' pixel at domain border (NBOUNDCUT)
    rlat_values = np.array([SWC_Y + (i*dy) for i in range(NBOUNDCUT, ny-NBOUNDCUT)])
    rlat[...] = rlat_values[...]

    if withZlvl:
        lvl = nc_file.createVariable('lvl', 'f4', ('lvl',),
                      zlib=True)
        lvl.standard_name = "level"
        lvl.long_name = "ParFlow layers"
        lvl.units = "-"
        lvl.axis = "Z"
        lvl_values = np.arange(nz)
        lvl[...] = lvl_values[...]

    if withTime:
        ncTime = nc_file.createVariable('time', 'i2', ('time',), zlib=True)
        ncTime.units = f'{timeUnit}'
        ncTime.calendar = f'{timeCalendar}'

    # Create grid-mapping for rotated-pole grid
    rotated_pole = nc_file.createVariable('rotated_pole', 'i2', zlib=True)
    rotated_pole.long_name = "coordinates of the rotated North Pole"
    rotated_pole.grid_mapping_name = "rotated_latitude_longitude"
    rotated_pole.grid_north_pole_latitude = rpol_Y
    rotated_pole.grid_north_pole_longitude = rpol_X

    # Close netCDF file for save and return
    nc_file.close()
    return fileName

def mappIndicator(ParFlowNamelist, IndicatorFile):
    ###############################################################################
    ### parse ParFlow-Namelist to get indicators values
    ###############################################################################
    PFNamelistDict = ppfl.pars_ParFlowNamelist(f'{ParFlowNamelist}')

    nz = int(PFNamelistDict['ComputationalGrid.NZ'])
    ny = int(PFNamelistDict['ComputationalGrid.NY'])
    nx = int(PFNamelistDict['ComputationalGrid.NX'])
    tcl_dz_keys = [f'Cell.{i}.dzScale.Value' for i in range(nz)]
    dz_mult = [float(PFNamelistDict[tcl_dz_key]) for tcl_dz_key in tcl_dz_keys]
    dz_mult = np.asarray(dz_mult)

    dx = float(PFNamelistDict['ComputationalGrid.DX'])
    dy = float(PFNamelistDict['ComputationalGrid.DY'])
    dz = float(PFNamelistDict['ComputationalGrid.DZ'])

    # indi_input indinput, etc is individual for each namelist...!
    # maybe read all GeomInputs, stre as list, loop over all and merge.
    # if one is crashing, skip
    GeomInputNames = PFNamelistDict['GeomInput.Names']
    for GeomInputName in GeomInputNames:
        try:
            GeomInputs = PFNamelistDict[f'GeomInput.{GeomInputName}.GeomNames']

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
                nvg      = np.where(Indi3D == IndicatorInput[GeomName], vanG_n[GeomName], nvg)
                sres     = np.where(Indi3D == IndicatorInput[GeomName], vanG_sres[GeomName], sres)
        except KeyError:
            continue

    outDict = {}
    outDict['alpha']   = alpha
    outDict['nvg']     = nvg
    outDict['sres']    = sres
    outDict['dz_mult'] = dz_mult
    outDict['dz']      = dz
    outDict['dy']      = dy
    outDict['dx']      = dx
    outDict['nz']      = nz
    outDict['ny']      = ny
    outDict['nx']      = nx

    return outDict
        
if __name__ == '__main__':
    print('Im there!')
