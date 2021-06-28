""" Example script to show how to write data in netCDF format.

netCDF is a standalone file format, meaning every meta-information needed to 
work with the data inside the netCDF file (units, descriptions, source, author, 
coordinates, etc), has to be provided with the netCDF files itself.
To achieve this one has to provide those meta-information already while 
writing / creating the netCDF file.

"""
import netCDF4 as nc
import numpy as np
import datetime as dt

###############################################################################
#### Collecting meta-information
###############################################################################
# To get some 'real' coordinates the following line is using the grid-definition
# of the DE05 Domain, which can be found here:
# https://icg4geo.icg.kfa-juelich.de/SoftwareTools/prepro_parflowclm_de05_static_fields/blob/master/doc/content/Grid.rst
nx = ny = 2000    # number of pixel in x and y direction
dx = dy = 0.0055  # spatial resolution in [deg]
# As DE05 is defined in rotated coordinates, the following coordinates are
# given in rotated coordinates as well.
llc_X = -10.82725 # lower left corner of grid in [deg] (r_longitude)
llc_Y = -5.38725  # lower left corner of grid in [deg] (r_latitude)
urc_X = 0.16725   # upper right corner of grid in [deg] (r_longitude)
urc_Y = 5.60725   # upper right corner of grid in [deg] (r_latitude)
# To define rotated coordinates in relation to geographical coordinates we need
# to know the geographical coordinates of the rotated north pole
rotatedpole_X = -162.0 # rotated pole in geographical coordinates in [deg] (longitude)
rotatedpole_Y = 39.25  # rotated pole in geographical coordinates in [deg] (latitude)

saveFile = '../data/examples_Write2NetCDF_long.nc'

###############################################################################
#### Data to store in netCDF file
###############################################################################
# create example data to store with the netCDF file:
# using 5 time steps, and nx, ny of DE05 grid
np.random.seed(42)
data = np.random.rand(5,ny,nx)

###############################################################################
#### Create netCDF file and fill with attributes
###############################################################################
# Open the file defined in 'saveFile' with mode='w' to write/create the file 
# and overwrites in case file already exist. 
nc_file = nc.Dataset(saveFile, 'w', format='NETCDF4')
# Add some global attributes to provide the user with some basic information 
# what is store with the netCDF file while inspection with e.g. ncdump.
nc_file.author      = 'Niklas WAGNER'
nc_file.contact     = 'n.wagner@fz-juelich.de'
nc_file.institution = 'FZJ - IBG-3'
nc_file.description = 'Write a short description of your data to ship with the netCDF files!'
nc_file.history     = f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}'
nc_file.source      = 'add source here'

# Create dimensions
# Dimensions are mandatory for netCDF files, as every variable has to be 
# assigned to the correct related dimension.
# We want to store 3D data in time, Y, X
dtime = nc_file.createDimension('time',None)
drlat = nc_file.createDimension('rlat',ny)
drlon = nc_file.createDimension('rlon',nx)

# Create coordinate variables. 
# For each dimension one should create a variable with the same name as the 
# dimension. With this variable one should store the correct coordinate values,
# as a dimension in netCDF does not contain any real value. 
rlon = nc_file.createVariable('rlon', 'f4', ('rlon',))
# If the variable is created, add some description to is:
# 'standard_name' is the name of the variable according e.g. cf-convention.
rlon.standard_name = "grid_longitude"
# A longer name, if 'standard_name' is some short-form
rlon.long_name     = "rotated longitude"
# Units of the variable
rlon.units         = "degrees"
# 
rlon.axis          = "X"
# The values of the variable, in this case the rotated coordinates .
# 'urc_X+0.0001' is used here, as np.arange(start, stop, step), does exclude 
# 'stop' value, which I need to include. 
rlon_values        = np.arange(llc_X, urc_X+0.0001, dx)
rlon[...]          = rlon_values[...]

rlat = nc_file.createVariable('rlat', 'f4', ('rlat',))
rlat.standard_name = "grid_latitude"
rlat.long_name     = "rotated latitude"
rlat.units         = "degrees"
rlat.axis          = "Y"
rlat_values        = np.arange(llc_Y, urc_Y+0.0001, dy)
rlat[...]          = rlat_values[...]

# This example is using a grid defined on rotated coordinates, wherefore
# it is needed to define a variable called 'rotated_pole' and fill with 
# the correct values.
# This variable is needed to provide third-party tools with all information
# about the coordinates used, e.g. to correctly remap the data with CDO. 
rotated_pole = nc_file.createVariable('rotated_pole', 'i2')
rotated_pole.long_name = "coordinates of the rotated North Pole"
# 'grid_mapping_name' has to follow cf-convention
rotated_pole.grid_mapping_name = "rotated_latitude_longitude"
rotated_pole.grid_north_pole_longitude = rotatedpole_X
rotated_pole.grid_north_pole_latitude = rotatedpole_Y

# Create the actual variable we want to store the data at.
# Name of the variable: 'TestData'
# dtype: 'f4' = float32
# dimension: ('time', 'rlat', 'rlon',)
#   The dimension have to exist.
#   The dimension have to fit the shape of the data
# fill__value: -9999
#   Missing values are masked with the value defined here
# zlib: True
#   This enables data compression to save storage
ncVar = nc_file.createVariable('TestData', 'f4', ('time', 'rlat', 'rlon',),
                                fill_value=-9999,
                                zlib=True)
ncVar.standard_name = 'test_name'
ncVar.long_name = 'variable to test writing netCDF'
ncVar.units ='-'
# See above while creating 'rotated_pole' variable.
# Variables need to refer to the 'rotated_pole' variable.
ncVar.grid_mapping = 'rotated_pole'

# Create time variable.
# Time values are stored quiet easy with netCDF - by simple integers 'i2'.
# Those integers values do represent steps since a give reference date, which
# is defined with the time.units.
# E.g:
# units = 'days since 1979-01-01 00:00:00'
# time  = 0, 1, 2, 3, 4
# do represent:
# '1979-01-01 00:00:00', '1979-01-02 00:00:00', '1979-01-03 00:00:00',
# '1979-01-04 00:00:00', '1979-01-05 00:00:00'
# This works e.g. with: 'seconds since ', 'hours since ', 'moths since', etc
# With the correct 'calendar' attribute this also takes different calendars 
# into account (e.g. leap and noleap).
# For full information see:
# https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch04s04.html
ncTime = nc_file.createVariable('time', 'i2', ('time',))
ncTime.standard_name = 'time'
ncTime.units = 'days since 1979-01-01 00:00:00'
ncTime.calendar = '365_day'

###############################################################################
#### Fill netCDF file with data
###############################################################################
# While creating the variable we want to store the data to we assigned this to 
# the name 'ncVar', which therefore can be treated as a kind of usual array, 
# and we can assign the data generated previously ([...] means every dimension:
ncVar[...] = data[...]
# The same for the time, where we do use a fictive time-axis:
ncTime[...] = np.arange(data.shape[0])

# Finlay close the newly created netCDF file to save made changes.
# Do inspect the file with 'ncdump -h' / 'ncview' to see if everything is 
# correct.
nc_file.close()
