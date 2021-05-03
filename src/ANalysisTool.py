import sys
import os
import numpy as np
import netCDF4 as nc
import datetime as dt
import pars_ParFlowTCL as ppfl
import ParFlow_IO as pio

class toolBox:
    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    # As this is a simple collection of functions no constructor should be 
    # needed.
    def __init__(self):
        pass

    def get_intervalSlice(dates, sliceInterval='month'):
    # def get_intervalSlice(dates, dumpInterval, sliceInterval='month'):
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
            print(f'DONE')        
        
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

            if (tmp_first - dumpInterval) == tmp_first_ref:
                print(f'check step {Offset} is first of a month at midnight')
                break
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
        # now that we know the dumpInterval and first step is first of month...
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
        elif domain == 'EU11':
            ny = 432           # [-]
            nx = 444           # [-]
            dy = dx = 0.110    # [deg]
            # Rotated pole geographical coordinates
            rpol_Y  = 39.25    # [deg]
            rpol_X  = -162.0   # [deg]
            # South-West Corner of domain in rotated coordinates
            SWC_Y = -24.4720   # [deg]
            SWC_X = -29.4713   # [deg]
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
        GeomInputs = PFNamelistDict['GeomInput.indi_input.GeomNames']

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