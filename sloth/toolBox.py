""" toolBox - submodule of SLOTH

author: Nikals WAGNER
e-mail: n.wagner@fz-juelich.de
version: 2022-05-04

Description:
[...]
"""
import sys
import os
import glob
import numpy as np
import datetime
import matplotlib.pyplot as plt
from calendar import monthrange
from . import IO as io
from scipy import ndimage as nd

import sloth.slothHelper as slothHelper


def calc_catchment(slopex, slopey, x, y):
    """ This function calculates the catchment related to a given outlet
    pixel based on x- and y-slope files

    This algorithem is quiet simple, but I found no faster one, at least
    not without using other external libs.
    The idea is to start with a list of pixels belonging to the catchment, for 
    example with a list with a single entry -- the given outlet pixel. 
    Than loop over this list, pop (!) the current pixel and 
    i)  mark the current pixel as catchment and
    ii) check all sorounding pixels, if those does drain into the current pixel
    and does not belong to the catchment already.
    All found sorrounding pixels, which does drain into the current one and
    does not belong to the catchment already are appended to the list the loop
    is itterating over. Than the next itteration is started.
    This way the algorithem does find every pixel belonging to the catchment of
    the passes initial pixel(s) -- e.g. the outlet pixel.

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
        # Get one open end
        step = openEnds.pop()
        # Mark open end as catchment
        fc[step] = 1
        # GET surrounding Pixel
        # Nort = step - nx; East = step +1
        # West = step - 1 ;South = step + nx
        # NEWS = [step-nx, step+1, step-1, step+nx]
        NS = [step-nx, step+nx]
        EW = [step+1, step-1]
        # CHECK if surrounding drain to step (D2S) and are NOT catchment already.
        # Step is the current handled open end.
        try:
            NSD2S = [ idx for idx in  NS if (idx + fdy[idx] == step and not fc[idx] ) ] 
            EWD2S = [ idx for idx in  EW if (idx + fdx[idx] == step and not fc[idx] ) ] 
        except IndexError:
            print('FEHLER')
            continue
        D2S = NSD2S + EWD2S
        # add all found pixes to openEnds
        openEnds += D2S

    return fc.reshape(dims)

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

def spher_dist_v1(lon1, lat1, lon2, lat2, Rearth=6373):
    """ calculate the spherical / haversine distance

    Source: https://www.kompf.de/gps/distcalc.html
    This function is supposed to proper handle different shaped coords
    latX and lonX is supposed to be passed in rad

    return 2D ndarray
    """
    term1 = np.sin(lat1) * np.sin(lat2)
    term2 = np.cos(lat1) * np.cos(lat2)
    term3 = np.cos(lon2 - lon1)
    return Rearth * np.arccos(term1+term2*term3)

def find_nearest_Index_2D(point, coord2D):
    dist = np.abs(coord2D - point)
    idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
    return idx[0], idx[1]

def plusOneMonth(currDate):
    num_days = monthrange(currDate.year, currDate.month)[1]
    return currDate + datetime.timedelta(hours=24*num_days)

def trunc(values, decs=0):
    """ truncates a passed float value by given floating point digit

    This funciton does truncate a given float value or floting ndarray
    by the precision defined by `decs`.

    Example: trunc(values=2.12345, decs=3) --> 2.123
    
    Input value:
    ------------
    values: ndarray
        A ndarray of any dimension (also scalar).
    decs: int
        A integer value defining the truncation precision.

    Return value:
    -------------
    values: ndarray
        A ndarray of same type as input

    """
    return np.trunc(values*10**decs)/(10**decs)

def fill(data, invalid=None, transkargs={}):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell.
    Source: https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array

    Input:
        data:       numpy array of any dimension
        invalid:    a binary array of same shape as 'data'.
                    data value are replaced where invalid is True
                    If None (default), use: invalid  = np.isnan(data)
        transkargs: optional arguments one want to pass to
                    distance_transform_edt()

    Output:
        Return a filled array.
    """
    # check if data and invalid shape is equal, as otherwise filling is not
    # possible
    if data.shape != invalid.shape:
        print(f'ERROR: data.shape {data.shape} != invalid.shape (invalid.shape)')
    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid,
                                    return_distances=False,
                                    return_indices=True,
                                    **transkargs)
    return data[tuple(ind)]

def get_prudenceMask(lat2D, lon2D, prudName):
    """ return a prudance mask

    Return a boolean mask-array (True = masked, False = not masked) based on
    a passed set of longitude and latitude values and the name of the prudence
    region.
    The shape of the mask-array is set equal to the shape of input lat2D.
    Source: http://prudence.dmi.dk/public/publications/PSICC/Christensen&Christensen.pdf p.38

    INPUT:
      lat2D:    ndarray
                2D latitude information for each pixel
      lon2D:    ndarray
                2D longitude information for each pixel
      prudName: str
                Short name of prudence region
    RETURN
      prudMask: ndarray
                Ndarray of dtype boolean of the same shape as lat2D.
                True = masked; False = not masked
    """
    if (prudName=='BI'):
        prudMask = np.where((lat2D < 50.0) | (lat2D > 59.0)  | (lon2D < -10.0) | (lon2D >  2.0), True, False)
    elif (prudName=='IP'):
        prudMask = np.where((lat2D < 36.0) | (lat2D > 44.0)  | (lon2D < -10.0) | (lon2D >  3.0), True, False)
    elif (prudName=='FR'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 50.0)  | (lon2D < -5.0) | (lon2D >  5.0), True, False)
    elif (prudName=='ME'):
        prudMask = np.where((lat2D < 48.0) | (lat2D > 55.0)  | (lon2D < 2.0) | (lon2D >  16.0), True, False)
    elif (prudName=='SC'):
        prudMask = np.where((lat2D < 55.0) | (lat2D > 70.0)  | (lon2D < 5.0) | (lon2D >  30.0), True, False)
    elif (prudName=='AL'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 48.0)  | (lon2D < 5.0) | (lon2D >  15.0), True, False)
    elif (prudName=='MD'):
        prudMask = np.where((lat2D < 36.0) | (lat2D > 44.0)  | (lon2D < 3.0) | (lon2D >  25.0), True, False)
    elif (prudName=='EA'):
        prudMask = np.where((lat2D < 44.0) | (lat2D > 55.0)  | (lon2D < 16.0) | (lon2D >  30.0), True, False)
    else:
        print(f'prudance region {prudName} not found --> EXIT')
        eys.exit(1)

    return prudMask

def stampLLSM(data, invalid, LLSM, LLSMThreshold=2):
    """ stamps a LLSM to passed data

    Some times its needed to force different land datasets to use the same Land
    Lake Sea Mask (LLSM). One example is using different data sets to provide a
    model with soilindicators and landcover data. If the origin of both
    datsets is different, most propably the used LLSM is different too. This,
    for exmple could lead to the situation along the coastline that one dataset
    is indicating water (ocean) while the other is indicating land.
    This inconsistency has to be fixe, especially within the ealm of coupled
    models.

    To achive this, this function is
    i) masking invalid pixels within the passed data set
    ii) interpolating remaining data over the maske regions
    iii) setting pixel to water value according to a passed LLSM

    This way guaranthees each land pixel is holding land informtion and each
    water pixel is marked as water.

    INPUT:
      data: ndarray
            A 2D nd array holding the datset to stamp the LLSM
      invalid: scalar
               A scalar witht the value indicating invalid pixel
      LLSM: ndarray
            A ndarray of the same shape as data holding the LLSM
            (Land=2 Lake=1 Sea=0)

    RETURN:
      out: ndarray MASKED!
           A ndarray of the same shape as data, with stamped LLSM
    """

    # some checks to verify we can progress
    if not isinstance(data, np.ndarray):
        print(f'data is of type {type(data)} but <class "numpy.ndarray"> is required!')
        sys.exit(1)
    if not data.ndim == 2:
        print(f'data is of dimension {data.ndim} but dimension 2 is required!')
        sys.exit(1)

    # i)
    out = np.ma.masked_where(data==invalid, data)
    # ii)
    out = fill(data=out, invalid=out.mask)
    # iii)
    out = np.ma.masked_where((LLSM < LLSMThreshold), out)

    return out


def mapDataRange_lin(X, y_min=0, y_max=1):
    """ Map src data range linear to other data range.

    Mapping a data range to another data range could be quiet usefull. One 
    example is to normalize a data range (map [x,y] --> [0,1]) to better compare
    to other data ranges.

    This function is linear mapping arbitrarry source data ranges (X) to 
    arbitrary target data ranges (Y). 
    An intermediat step is used by first transform both (X and Y) into data 
    ranges starting from zero (X' and Y'), as those data range can be easily 
    mapped with
    (I)   y' = y'_max / x'_max * x'
    whereby 
    (II)  y' = y - y_min   AND   x' = x - x_min
    Putting together (I) and (II) does lead to
    (III) y = (y_max-y_min) / (x_max-x_min) * (x-x_min) + y_min

    Input values:
    -------------
    X: ndarray
        Source data range to remap
    y_min: scalar
        Min. value of target data range
    y_max: scalar
        Max. value of target data range
    x_min: scalar
        Min. value of source data range, if this should not be calculated based
        on X.
    x_max: scalar
        Max. value of source data range, if this should not be calculated based
        on X.

    Return value:
    -------------
    return: ndarray
        The target data range
    """
    x_min = np.min(X)
    x_max = np.max(X)

    Y = (y_max-y_min) / (x_max-x_min) * (X-x_min) + y_min

    return Y

if __name__ == '__main__':
    print('Im there!')
