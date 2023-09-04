import numpy as np
import netCDF4 as nc
import os
from struct import pack, unpack

import sloth.slothHelper as slothHelper
import sloth.coordTrafo

################################################################################
################################ PFB ###########################################
################################################################################
def write_packed(f, fmt, val):
    f.write(pack(fmt, val))

def create_pfb(filename, var, delta=(1, 1, 1), subgrids=(1, 1, 1)):
    """
    Create a ParFlow PFB file from a 3D array of variable data.

    Parameters
    ----------
    filename : str
        Name of the PFB file to be created.
    var : ndarray
        3D array of variable data to be stored in the PFB file.
    delta : tuple of float, optional
        Grid spacing values in the x, y, and z directions. Default is (1, 1, 1).
    subgrids : tuple of int, optional
        Number of subgrids in the x, y, and z directions. Default is (1, 1, 1).

    Examples
    --------
    >>> import numpy as np
    >>> data = np.random.rand(10, 20, 30)
    >>> create_pfb('output.pfb', data, delta=(0.5, 0.5, 0.5))
    ...
    This creates a PFB file named 'output.pfb' from a 3D array 'data' with 
    custom delta settings.

    """
    nz, ny, nx = var.shape
    dz, dy, dx = delta
    sz, sy, sx = subgrids

    filepfb = open(filename, 'wb')

    # Write start indices of global domain in x, y, z direction
    write_packed(filepfb, '>d', 0)
    write_packed(filepfb, '>d', 0)
    write_packed(filepfb, '>d', 0)

    # Write number of global gridpoints in x, y, z direction
    write_packed(filepfb, '>i', nx)
    write_packed(filepfb, '>i', ny)
    write_packed(filepfb, '>i', nz)

    # Write delta x, delta y and delta z
    write_packed(filepfb, '>d', dx)
    write_packed(filepfb, '>d', dy)
    write_packed(filepfb, '>d', dz)

    nSubGrid = np.prod(subgrids)

    nnx = int(nx / sx)
    nny = int(ny / sy)
    nnz = int(nz / sz)

    # Write the subgrid grid ID
    write_packed(filepfb, '>i', nSubGrid)

    for iz in np.arange(sz)*nnz:
        for iy in np.arange(sy)*nny:
            for ix in np.arange(sx)*nnx:
                #print(ix,iy,iz, nnx,nny,nnz)
                # Write start indices in x, y, z direction
                write_packed(filepfb, '>i', int(ix))
                write_packed(filepfb, '>i', int(iy))
                write_packed(filepfb, '>i', int(iz))

                # Write number of grid points in x, y and z direction for this 
                # subgrid
                write_packed(filepfb, '>i', nnx)
                write_packed(filepfb, '>i', nny)
                write_packed(filepfb, '>i', nnz)

                # Write the relative(to global) grid refinement in this subgrid
                # 0=same resolution as global
                write_packed(filepfb, '>i', 0)
                write_packed(filepfb, '>i', 0)
                write_packed(filepfb, '>i', 0)

                # Assuming the data is stored in 3D array called varArray of 
                # global size nz*ny*nx
                fmt = ">%dd" % (nnz*nny*nnx)
                filepfb.write(pack(fmt, *var[iz:iz+nnz,
                                               iy:iy+nny,
                                               ix:ix+nnx].flatten()))

    filepfb.close()

def read_pfb(filename):
    """
    Read a ParFlow PFB file and return the data as a numpy ndarray.

    Parameters
    ----------
    filename : str
        Name of the PFB file to be read.

    Returns
    -------
    data : ndarray
        3D array containing the data read from the PFB file.

    Examples
    --------
    >>> data = read_pfb('input.pfb')
    >>> print(data.shape)
    (10, 20, 30)
    ...
    This reads a PFB file named 'input.pfb' and returns the data as a 3D array.

    """
    with open(filename, "rb") as f:
        # read meta informations of datafile
        meta_inf = np.fromfile(f, dtype='>f8', count = 3)
        x1 = meta_inf[0]
        y1 = meta_inf[1]
        z1 = meta_inf[2]


        meta_inf = np.fromfile(f, dtype='>i4', count = 3)
        nx = meta_inf[0]
        ny = meta_inf[1]
        nz = meta_inf[2]
        nn = nx * ny * nz


        meta_inf = np.fromfile(f, dtype='>f8', count = 3)
        dx = meta_inf[0]
        dy = meta_inf[1]
        dz = meta_inf[2]


        meta_inf = np.fromfile(f, dtype='>i4', count = 1)
        nsubgrid = meta_inf[0]

        data =  np.ndarray(shape=(nz,ny,nx), dtype='>f8')

        for s in range(nsubgrid):
            meta_inf = np.fromfile(f, dtype='>i4', count = 9)
            ix = meta_inf[0]
            iy = meta_inf[1]
            iz = meta_inf[2]
            # print("---{0} Start Index (X,Y,Z):".format(s+1), ix, iy, iz)

            nx = meta_inf[3]
            ny = meta_inf[4]
            nz = meta_inf[5]
            nn = nx*ny*nz
            # print("---{0} Dimensions (X,Y,Z):".format(s+1), nx, ny, nz)

            rx = meta_inf[6]
            ry = meta_inf[7]
            rz = meta_inf[8]
            # print("---{0} Offsets (X,Y,Z):".format(s+1), rx, ry, rz)

            tmp_data = np.fromfile(f, dtype='>f8', count=nn).reshape((nz,ny,nx))

            data[iz:iz+nz, iy:iy+ny, ix:ix+nx] = tmp_data

    return data

def read_pfbMetaData(filename):
    with open(filename, "rb") as f:
        # read meta informations of datafile
        meta_inf = np.fromfile(f, dtype='>f8', count = 3)
        x1 = meta_inf[0]
        y1 = meta_inf[1]
        z1 = meta_inf[2]
        print(f'x1: {x1}; y1: {y1}; z1: {z1}')


        meta_inf = np.fromfile(f, dtype='>i4', count = 3)
        nx = meta_inf[0]
        ny = meta_inf[1]
        nz = meta_inf[2]
        nn = nx * ny * nz
        print(f'nx: {nx}; ny: {ny}; nz: {nz}; nn: {nn}')


        meta_inf = np.fromfile(f, dtype='>f8', count = 3)
        dx = meta_inf[0]
        dy = meta_inf[1]
        dz = meta_inf[2]
        print(f'dx: {dx}; dy: {dy}; dz: {dz}')


        meta_inf = np.fromfile(f, dtype='>i4', count = 1)
        nsubgrid = meta_inf[0]
        print(f'nsubgrid: {nsubgrid}')

    return None

################################################################################
############################# netCDF ###########################################
################################################################################

def createNetCDF(fileName, domain=None, nz=None, calcLatLon=False,
    author=None,
    description=None, source=None, contact=None, institution=None,
    history=None, timeCalendar=None, timeUnit=None, NBOUNDCUT=0):
    """
    Create a NetCDF file with typical metadata and dimensions.

    Parameters
    ----------
    fileName : str
        Name of the output NetCDF file.
    domain : str or None, optional
        Path to the domain definition file or a valid CORDEX/Griddes domain name.
        If None, default domain definitions will be used.
    nz : int or None, optional
        Number of vertical levels. If None, no z-axis will be created.
    calcLatLon : bool, optional
        Flag indicating whether to calculate latitude and longitude values.
        If True, lat and lon variables will be created.
    author : str or None, optional
        Name of the author.
    description : str or None, optional
        Description of the NetCDF file.
    source : str or None, optional
        Source of the data.
    contact : str or None, optional
        Contact information.
    institution : str or None, optional
        Institution associated with the data.
    history : str or None, optional
        History information.
    timeCalendar : str or None, optional
        Calendar type for the time-axis.
    timeUnit : str or None, optional
        Unit of measurement for the time-axis.
    NBOUNDCUT : int, optional
        Number of pixels to cut at the domain border.

    Returns
    -------
    fileName : str
        Name of the created NetCDF file.

    Example
    -------
    >>> # Create a NetCDF file with default domain definitions and no z-axis
    >>> createNetCDF("output.nc")
    
    >>> # Create a NetCDF file with a specific domain and 10 vertical levels
    >>> createNetCDF("output.nc", domain="my_domain.txt", nz=10, calcLatLon=True)

    """

    #######################################################################
    #### Get domain definitions
    #######################################################################
    availableCORDEXDomains  = slothHelper.get_listOfCordexGrids()
    availableGriddesDomains = slothHelper.get_listOfGriddes()
    # Check if 'domain' is pointing to a domain path
    if os.path.exists(domain):
        domainDef = slothHelper.get_griddesDomDef(domain)
    # Check if 'domain' is a official CORDEX name pattern
    elif domain in availableCORDEXDomains:
        domainDef = slothHelper.get_cordexDomDef(domain)
    # Check if 'domain' is provided by SLOTH
    elif domain in availableGriddesDomains:
        # Read griddes file from configs dir
        # Configs is located under `sloth/` (os.path.dirname(__file__))
        griddesFileName = f'{domain}_griddes.txt'
        griddesFile     = f'{os.path.dirname(__file__)}/configs/{griddesFileName}'
        if not os.path.isfile(griddesFile):
            print(f'ERROR: There is no griddes file with name {griddesFile} --> EXIT')
            sys.exit(1)
        domainDef = slothHelper.get_griddesDomDef(griddesFile)
    else:
        print(f'ERROR: passed domain is not supported. domain={domain} --> Exit')
        return False

    nx     = domainDef['Nlon']
    ny     = domainDef['Nlat']
    dx     = domainDef['dlon']
    dy     = domainDef['dlat']
    rpol_X = domainDef['NPlon']
    rpol_Y = domainDef['NPlat']
    SWC_X  = domainDef['SWlon']
    SWC_Y  = domainDef['SWlat']

    #######################################################################
    #### Checking if time- and / or z-axis is used
    #######################################################################
    withTime = False
    if timeUnit is not None and timeCalendar is not None:
        withTime = True
    else:
        #print('NOT creating time-axis')
        #print(f'--  timeUnit = {timeUnit}; timeCalendar = {timeCalendar}')
        pass

    withZlvl = False
    if nz is not None:
        withZlvl = True
    else:
        #print('NOT creating z-axis')
        pass

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

    if calcLatLon:
        rlon2D, rlat2D = np.meshgrid(rlon_values, rlat_values)
        lat2D, lon2D = sloth.coordTrafo.undo_grid_rotation(
                rlat = rlat2D, rlon = rlon2D, 
                np_lat = rpol_Y, np_lon = rpol_X)

        lat = nc_file.createVariable('lat', 'f4', ('rlat','rlon'),
                                    zlib=True)
        lat.standard_name = "latitude"
        lat.long_name = "latitude"
        lat.units = "degrees_north"
        lat.coordinates = "lon lat"
        lat.grid_mapping = "rotated_pole"
        lat[...] = lat2D[...]

        lon = nc_file.createVariable('lon', 'f4', ('rlat','rlon'),
                                    zlib=True)
        lon.standard_name = "longitude"
        lon.long_name = "longitude"
        lon.units = "degrees_east"
        lon.coordinates = "lon lat"
        lon.grid_mapping = "rotated_pole"
        lon[...] = lon2D[...]
    else:
        #print(f'-- no lat lon values used')
        pass

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
        ncTime = nc_file.createVariable('time', 'f8', ('time',), zlib=True)
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

def readSa(file):
    """
    Reads data from a file in ASCI format and returns a NumPy array.

    Parameters
    ----------
    file : str
        The file path to read the data from.

    Returns
    -------
    numpy.ndarray
        A NumPy array containing the data read from the file.

    Example
    -------
    >>> data = readSa('data.txt')
    >>> print(data)
    [[1.2 3.4]
     [5.6 7.8]]

    """
    with open(file, 'r') as f:
        header = f.readline()
        nx, ny, nz = (int(item) for item in header.split(' '))

        data = np.genfromtxt(f, dtype=float)
        data = data.reshape((nz, ny, nx))

        return data

def writeSa(file, data):
    """
    Writes data to a file in ASCI format.

    Parameters
    ----------
    file : str
        The file path to write the data to.
    data : numpy.ndarray
        The NumPy array containing the data to be written.

    Returns
    -------
    None

    Example
    -------
    >>> data = np.array([[1.2, 3.4], [5.6, 7.8]])
    >>> writeSa('output.txt', data)

    """
    nz, ny, nx = data.shape
    with open(file, 'w') as f:
        f.write(f'{nx} {ny} {nz}\n')
        # Below should be more easy with a flatten array...
        # But how to flatt? C or F order? 
        # If knowen and tested change below.
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f'{data[k,j,i]}\n')
