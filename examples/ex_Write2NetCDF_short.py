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
import sys

sloth_path='../'
sys.path.append(sloth_path)
import sloth


###############################################################################
### Define some paths, filenames, etc
###############################################################################
saveFile = '../data/examples_Write2NetCDF_short.nc'

###############################################################################
#### Data to store in netCDF file
###############################################################################
# create example data to store with the netCDF file:
# using 5 time steps, and nx, ny of DE05 grid
np.random.seed(42)
data = np.random.rand(5,2000,2000)

###############################################################################
#### Create netCDF file and fill with basic attributes
###############################################################################
# For mor detailed information about how createNetCDF() does work, see
# sloth/toolBox.py --> createNetCDF()
netCDFFileName = sloth.toolBox.createNetCDF(saveFile, domain='DE05', 
	author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
	institution='FZJ - IBG-3', history=f'Created: {dt.datetime.now().strftime("%Y-%m-%d %H:%M")}',
	description='Write a short description of your data to ship with the netCDF files!',
	source='add source here')

###############################################################################
#### Create the actual variable we want to store the data at.
###############################################################################
with nc.Dataset(netCDFFileName, 'a') as nc_file:
	# Name of the variable: 'TestData'
	ncVar = nc_file.createVariable('TestData', 'f4', ('time', 'rlat', 'rlon',),
	                                fill_value=-9999,
	                                zlib=True)
	ncVar.standard_name = 'test_name'
	ncVar.long_name = 'variable to test writing netCDF'
	ncVar.units ='-'
	ncVar.grid_mapping = 'rotated_pole'

	ncTime = nc_file.createVariable('time', 'i2', ('time',))
	ncTime.standard_name = 'time'
	ncTime.units = 'days since 1979-01-01 00:00:00'
	ncTime.calendar = '365_day'

	ncVar[...] = data[...]
	ncTime[...] = np.arange(data.shape[0])
