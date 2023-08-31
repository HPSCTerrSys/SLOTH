""" toolBox - submodule of SLOTH
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
import geopandas as gpd
import xarray as xr
import pandas as pd
from shapely.geometry import Polygon

import sloth.slothHelper as slothHelper


def calc_catchment(slopex, slopey, x, y):
    """ Calculate the catchment area associated with a given outlet pixel 
    using x- and y-slope data.

    This function implements a simple yet efficient algorithm to determine the 
    catchment area by starting with a list of pixels that initially contains 
    only the outlet pixel(s). It iteratively processes the pixels in the list, 
    marking them as part of the catchment and identifying the surrounding pixels 
    that drain into the current pixel and are not already included in the 
    catchment. Foun pixels are added to the list and the algorithm continues 
    until all pixels belonging to the catchment area are discovered.

    Parameters
    ----------
    slopex : ndarray
        2D slopes in x-direction
    slopey : ndarray
        2D slopes in y-direction
    x : int
        index in x-direction to calulate catchment from 
    y : int
        index in y-direction to calulate catchment from

    Returns
    -------
    catchment : ndarray
        2D ndarray of the same size as slopex/y. 0 = not part of catchment; 1 = part of catchment
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

    Parameters
    ----------
    dates : NDarray 
        the time-series as datetime-object. 
    sliceInterval : str
        defining the interval. Supported are: 'day', 'month'

    Returns
    -------
    Slices : list
        A list with length of found intervals with related slices

    Notes 
    -----
    This function assumes:   
    i) The passed time-series covers at least on intervall!
    ii) At least hourly steps / dumpIntervals! 
    iii) No seconds in dates - model output is at least on full minutes!

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

    This function calculates the real distance (in [km]) between two points on
    earth given in lat/lon coordinates. 

    Parameters
    ----------
    lon1 : ndarray
        Longitude value in [rad] of first point in [rad]. Could be any dim
    lat1 : ndarray
        Latitude value in [rad] of first point in [rad]. Could be any dim
    lon2 : ndarray
        Longitude value in [rad] of second point in [rad]. Could be any dim
    lat2 : ndarray
        Latitude value in [rad] of second point in [rad]. Could be any dim
    Rearth : int or float 
        The earth radius in [km].

    Returns
    -------
    ndarray
        The distance is returned in [km]

    Notes
    -----
    Source: https://www.kompf.de/gps/distcalc.html

    """
    term1 = np.sin(lat1) * np.sin(lat2)
    term2 = np.cos(lat1) * np.cos(lat2)
    term3 = np.cos(lon2 - lon1)
    return Rearth * np.arccos(term1+term2*term3)

def find_nearest_Index_2D(point, coord2D):
    """
    Find the nearest index in a 2D array.

    Given a 2D array `coord2D` and a specific `point` value, this function returns
    the index of the 2D array that holds the closest values to the `point`.

    Parameters
    ----------
    point : scalar
        The value for which to find the nearest index in the 2D array.
    coord2D : ndarray
        The 2D array in which to search for the nearest index.

    Returns
    -------
    (int, int)
        The index in the first and second dimensions of `coord2D` that holds the
        closest values to the `point`.

    Notes
    -----
    This function calculates the absolute differences between `coord2D` and `point`
    using `np.abs()`. It then determines the index of the minimum value in the
    differences array using `np.argmin()`. The resulting flattened index is converted
    into a tuple of indices representing the original shape of `coord2D` using
    `np.unravel_index()`.

    Examples
    --------
    >>> import numpy as np
    >>> arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> find_nearest_Index_2D(5, arr)
    (1, 1)

    >>> arr = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
    >>> find_nearest_Index_2D(0.25, arr)
    (0, 1)

    """
    dist = np.abs(coord2D - point)
    idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
    return idx[0], idx[1]

def plusOneMonth(currDate):
    """
    Return the passed date incremented by one month.

    This function provides a simple way to calculate a date that is one month ahead
    of the given `currDate`. Unlike the bash command-line tool `date`, which can
    handle this calculation easily, Python's `datetime` objects are based on hours
    or seconds, requiring special treatment for adding one month due to varying
    month lengths. This function serves as a wrapper to handle this task in a clean
    manner, facilitating usage in other code snippets.

    Parameters
    ----------
    currDate : datetime
        An arbitrary date.

    Returns
    -------
    datetime
        The passed date incremented by one month.

    Notes
    -----
    This function determines the number of days in the month of `currDate` using
    `monthrange()` from the `calendar` module. It then returns `currDate` incremented
    by a timedelta calculated as `24 * num_days` hours. This ensures that the returned
    date correctly reflects one month ahead while accounting for varying month lengths.

    Examples
    --------
    >>> import datetime
    >>> currDate = datetime.datetime(2023, 5, 15)
    >>> plusOneMonth(currDate)
    datetime.datetime(2023, 6, 15, 0, 0)

    >>> currDate = datetime.datetime(2023, 12, 31)
    >>> plusOneMonth(currDate)
    datetime.datetime(2024, 1, 31, 0, 0)

    """
    num_days = monthrange(currDate.year, currDate.month)[1]
    return currDate + datetime.timedelta(hours=24*num_days)

def trunc(values, decs=0):
    """
    Truncate a passed float value or floating ndarray.

    This function truncates a given float value or floating ndarray by a 
    specified truncation precision (decimal digits).

    Parameters
    ----------
    values : ndarray
        A ndarray of any dimension (including scalar) containing float values to be truncated.    
    decs : int, optional
        An integer value defining the truncation precision. It specifies the number of decimal places to keep (default is 0).

    Returns
    -------
    ndarray
        A ndarray of the same type as the input, with the specified number of decimal places truncated from the original values.

    Examples
    --------
    >>> import numpy as np
    >>> values = np.array([2.12345, 3.45678, 4.56789])
    >>> trunc(values, decs=3)
    array([2.123, 3.456, 4.567])

    >>> value = 2.98765
    >>> trunc(value, decs=2)
    2.98

    """
    return np.trunc(values*10**decs)/(10**decs)

def fill(data, invalid=None, transkargs={}):
    """
    Fill invalid data points with nearest neighbor interpolation.

    This function replaces invalid data points in the input array (`data`) with the
    value of the nearest valid data point. Invalid data points are indicated by the
    `invalid` array or, if not provided, by NaN (Not a Number) values in `data`.

    Parameters
    ----------
    data : ndarray
        An array containing the data to be filled.
    invalid : ndarray, optional
        A binary array of the same shape as `data`. Data values are replaced where
        `invalid` is True. If not provided, invalid data points are identified using
        `np.isnan(data)` (default).
    transkargs : dict, optional
        Additional arguments to pass to the `distance_transform_edt()` function.
        Default is an empty dictionary.

    Returns
    -------
    ndarray
        A filled ndarray, where invalid data points have been replaced by the value
        of the nearest valid data point.

    Notes
    -----
    This function uses nearest neighbor interpolation to fill invalid data points in
    the input array. It calculates the Euclidean distance transform of the `invalid`
    array to find the indices of the nearest valid data points. The original `data`
    array is then updated with the values from the nearest valid data points.

    If the shapes of `data` and `invalid` are not equal, an error is raised, as filling
    is not possible. 
    
    See also: https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array

    Example
    -------
    >>> import numpy as np
    >>> data = np.array([1.0, 2.0, np.nan, 4.0, np.nan])
    >>> fill(data)
    array([1., 2., 4., 4., 4.])

    >>> invalid = np.isnan(data)
    >>> fill(data, invalid)
    array([1., 2., 4., 4., 4.])

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
    """
    Return a Prudence mask based on latitude and longitude values.

    This function generates a boolean mask array where True represents masked areas
    and False represents non-masked areas. The mask is determined based on a set of
    latitude and longitude values and the name of the Prudence region.

    Parameters
    ----------
    lat2D : ndarray
        2D array containing latitude information for each pixel.
    lon2D : ndarray
        2D array containing longitude information for each pixel.
    prudName : str
        Short name of the Prudence region.

    Returns
    -------
    prudMask : ndarray
        Boolean ndarray of the same shape as `lat2D`, indicating masked and
        non-masked areas.
        True = masked; False = not masked.

    Notes
    -----
    The Prudence mask is created based on specific latitude and longitude ranges
    for each Prudence region. The function checks the `prudName` parameter and
    generates the corresponding mask using NumPy's `np.where()` function.

    The available Prudence region names and their corresponding latitude and
    longitude ranges are as follows (True is masked, Fals is not masked):

    - 'BI': Latitude: <50.0 or >59.0, Longitude: <-10.0 or >2.0
    - 'IP': Latitude: <36.0 or >44.0, Longitude: <-10.0 or >3.0
    - 'FR': Latitude: <44.0 or >50.0, Longitude: <-5.0 or >5.0
    - 'ME': Latitude: <48.0 or >55.0, Longitude: <2.0 or >16.0
    - 'SC': Latitude: <55.0 or >70.0, Longitude: <5.0 or >30.0
    - 'AL': Latitude: <44.0 or >48.0, Longitude: <5.0 or >15.0
    - 'MD': Latitude: <36.0 or >44.0, Longitude: <3.0 or >25.0
    - 'EA': Latitude: <44.0 or >55.0, Longitude: <16.0 or >30.0

    If an unsupported Prudence region name is provided, an error message is printed,
    and the program exits.

    For reference see also: http://prudence.dmi.dk/public/publications/PSICC/Christensen&Christensen.pdf p.38

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

    Some times its needed to force different datasets to use the same 
    Land-Lake-Sea-Mask (LLSM). One example is combining different data sets to 
    force a model with. If the datasets are taken fomr different sources, most 
    propably the used LLSM is different too. This, for exmple, could lead to the 
    situation along the coastline that one dataset is indicating water (ocean) 
    while the other is indicating land. This inconsistency has to be fixe, 
    especially within the realm of coupled models.


    To achive this, this function is\n
    i)   masking invalid pixels within the passed data set  
    ii)  interpolate remaining data over the maske regions  
    iii) set pixel to water value according to a passed LLSM  

    This guaranthees each land pixel is holding land informtion and each
    water pixel is marked as water.

    Parameters
    ----------
    data : ndarray
        A 2D nd array holding the datset to stamp the LLSM
    invalid : scalar
        A scalar witht the value indicating invalid pixel
    LLSM : ndarray
        A ndarray of the same shape as data holding the LLSM
        (Land=2 Lake=1 Sea=0)

    Returns
    -------
    ndarray
        A (MASKED!) ndarray of the same shape as data, with stamped LLSM

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


def mapDataRange_lin(X, y_min=0, y_max=1,
        x_min=None, x_max=None,
        cutMinMax=False):
    """Map a source data range linearly to a target data range.

    Perform a linear mapping of a source data range (X) to a target data range (Y). The
    function calculates the mapping using the formula: y = (y_max - y_min) / (x_max - x_min) * (x - x_min) + y_min.

    Parameters
    ----------
    X : ndarray
        Source data range to remap.
    y_min : scalar, optional
        Minimum value of the target data range (default is 0).
    y_max : scalar, optional
        Maximum value of the target data range (default is 1).
    x_min : scalar, optional
        Minimum value of the source data range. If not provided, it is calculated based on X.
    x_max : scalar, optional
        Maximum value of the source data range. If not provided, it is calculated based on X.
    cutMinMax : bool, optional
        If True, cut the mapped values outside the target data range to the minimum and maximum values.
    
    Returns
    -------
    ndarray
        The target data range after linear mapping.

    Notes
    -----
    This function is useful for mapping a source data range to a different target data range. It can be
    used for various purposes, such as normalizing data ranges to facilitate comparison with other ranges.
    An intermediat step is used by first transform both (X and Y) into data 
    ranges starting from zero (X' and Y'), as those data range can be easily 
    mapped with

    i)   y' = y'_max / x'_max * x'
    ii)  y' = y - y_min   AND   x' = x - x_min
    iii) y = (y_max-y_min) / (x_max-x_min) * (x-x_min) + y_min

    Examples
    --------
    >>> X = np.array([0, 1, 2, 3, 4, 5])
    >>> mapped = mapDataRange_lin(X, y_min=10, y_max=20)
    >>> print(mapped)
    [10. 12. 14. 16. 18. 20.]

    >>> X = np.array([100, 200, 300, 400, 500])
    >>> mapped = mapDataRange_lin(X, y_min=-1, y_max=1, x_min=100, x_max=500)
    >>> print(mapped)
    [-1.  -0.5  0.   0.5  1. ]

    """
    # Calculate x_min and x_max if not passed
    if x_min is None:
        x_min = np.min(X)
    if x_max is None:
        x_max = np.max(X)

    # Map X to new data range
    Y = (y_max-y_min) / (x_max-x_min) * (X-x_min) + y_min

    # Cut Y for min and max values if wanted
    if cutMinMax:
        Y = np.where(Y<=y_min, y_min, Y)
        Y = np.where(Y>=y_max, y_max, Y)

    return Y

def intersection_calculations(df_data, corners, area_of_interest, Name_area, crs_utm, nr_yr, nr_entries, save_dir):
    """ Calculate spatial mean values for a shapefile in interest.

    This function calculates spatial mean values (example for a specific 
    region or province) taking into account how much (ratio) of the 
    gridpoint is actually being intersected inside our region of interest.
    
    Parameters
    -------------
    df_data: dataframe
         dataframe containing infromation about each gridpoint
    
    corners: dataframe (can also be a nectdf or any other data type)
         containing infromation on the four longitudes and latitudes 
         that surrounds each gridpoint in the source dataframe
    
    area_of_interest: a shapefile or geodataframe
    	shapefile of the area of interest

    Name_area: the field name in the shapefile, that the dissolve will be based on
     
    crs_utm: projected coordinate reference system (utm)
    
    nr_yr: variable
    	number of years of interest
    
    nr_entries:variable
    	number of hours/days/or months etc..
    
    save_dir: path
    	path for saving the output 
     """
    

    # create a geodataframe from the dataframe that we want to work with
    df_data = pd.read_csv(df_data)
    print("creating a geodataframe of the given dataframe")
    gdf_data = gpd.GeoDataFrame(df_data, geometry=gpd.points_from_xy(df_data.lon,
                                                                     df_data.lat))
    ds_corners = xr.open_dataset(corners)

    # convert the netcdf to dataframe in order to be able 
    # to convert it to geopandas and create polygons
    df_corners = ds_corners.to_dataframe()
    df_corners = df_corners.reset_index()
    df_corners = df_corners.iloc[:, 4:]

    # create a list of the four pairs (lon,lat) that creates the polygons
    print("creating a list of four pairs (lon,lat) to create polygons")
    list_poly = []
    polygon_geom = list(zip(df_corners.grid_corner_lon, df_corners.grid_corner_lat))
    for i in range(0, len(polygon_geom), 4):
        list_poly.append(Polygon(polygon_geom[i:i + 4]))

    # create a geodataframe from the polygon list that was created
    print("creating a geodataframe with the polygons")
    gdf_polygons = gpd.GeoDataFrame(geometry=list_poly)
    # for some reason the polygons were doubled so a drop_duplicates() was used
    gdf_polygons = gdf_polygons.drop_duplicates() 

    # if the two goedataframes are do not have a coordinate systems
    # we have to set it first before going further setting the crs

    gdf_data_wgs = gdf_data.set_crs('EPSG:4326') # setting the crs to WGS84
    # in order to calculate the area later, it is needed to convert 
    # the geodataframe to UTM
    gdf_data_utm32N = gdf_data_wgs.to_crs(crs_utm) 

    gdf_polygons_wgs = gdf_polygons.set_crs('EPSG:4326')  # setting the crs to WGS84
    gdf_polygons_utm32N = gdf_polygons_wgs.to_crs(crs_utm)
    print("Geodataframes are created")

    # calculate the area of the polygons.
    # it is important for calculating the ratio (weight) of the area of 
    # each polygon interscted
    gdf_polygons_utm32N['area'] = gdf_polygons_utm32N.area

    # after calculating the area, in come cases (when calculating for whole country)
    # the area of polygons will differ for example in when comparing cities along 
    # the meridian lines, there will be slight difference in the area

    # joining thw two geodataframe using sjoin, in order to join the dataframe 
    # that has the values with its corresponding polygons
    print("Joining the two geodataframes")
    join_within_right_gdf_utm32N = gdf_data_utm32N.sjoin(gdf_polygons_utm32N, how="right", predicate="within")
    
    shapefile_gdf = gpd.read_file(area_of_interest)
    print(shapefile_gdf.crs)   # to check wether it has the same crs as the geodataframe that we created
    # if not convert it to utm (to be able to calculate the area)

    print("intersecting")
    overlay_gdf = gpd.overlay(join_within_right_gdf_utm32N, shapefile_gdf, how='intersection')
    overlay_gdf['area_intersected'] = overlay_gdf.area
    # the weight calculated below represents how much area from each polygon
    # is inside the shapefile after intesection
    overlay_gdf['weight'] = overlay_gdf['area_intersected']/overlay_gdf['area'] 


    # after this step we should have a geodataframe for the intersected polygons 
    # that are inside the shapefile of interest the geodataframe will have
    # the lon,lat, with the corresponding values for each date (each date in a column), 
    # area of the polygon before intersection, area after intersection,
    # and the ratio (area_intersected/area) and it can also include the names of the 
    # regions (in order to be able to dissolve it into mean values for each region)

    # in the next step, a dataframe was created to add the values of each 
    # polygon (gridpoint) after calculating its weight inside 
    # the intersection for each date
    nr_yr = nr_yr
    nr_entries = nr_entries # ex: number of dates :hourly, daily, monthly..etc
    nr_column = nr_yr*nr_entries
    df_data_inter = overlay_gdf.iloc[:,3:nr_column+2].multiply(overlay_gdf['weight'], axis="index")

    # copy the geometry column (important for creating a geodataframe) and other 
    # relevant data to the new dataframe
    # like the name of the regions (important for calculating the mean values for each region)
    Name_area = Name_area
    df_data_inter[Name_area] = overlay_gdf[Name_area] # this sould be adapted depending on the shapefile available
    df_data_inter['geometry'] = overlay_gdf['geometry']
    geometry = df_data_inter['geometry']
    gdf_data_inter = gpd.GeoDataFrame(df_data_inter, crs="epsg:25832", geometry=geometry)

    print("calculating the mean value for each region of interest")
    dissolve_gdf = gdf_data_inter.dissolve(by=Name_area, aggfunc='mean', as_index=True, level=None, sort=True, observed=False, dropna=True)
    dissolve_gdf_wgs = dissolve_gdf.to_crs('EPSG:4326') # convert back to wgs 1984
    
    print("saving as a shapefile")
    dissolve_gdf_wgs.to_file(save_dir)

    

if __name__ == '__main__':
    print('Im there!')
