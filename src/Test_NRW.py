import numpy as np
import sys
import os
import glob

import VAlidationTool as vat
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc

import sys
sys.path.append("../extern/transform_rotated_pole")
import rotate
sys.path.append("../extern/parflow")
from ParFlow_IO import read_pfb


# nx, ny = (200, 200)
# x = np.linspace(-5, 5, nx)
# y = np.linspace(15, 45, ny)
# SimLons, SimLats = np.meshgrid(x, y)
# plt.imshow(SimLons, origin='lower')
# plt.show()
# plt.imshow(SimLats, origin='lower')
# plt.show()
# SimMeanQ = mpl.image.imread('../data/SimEG.png')
# plt.imshow(SimMeanQ)
# plt.show()
# print(SimLons)
# print(SimLats)
# sys.exit(1)

# ObsX_idx 	= [28,64,149,48,105,80]
# ObsY_idx 	= [185,150,81,152,117,134]
# ObsMeanQ 	= np.asarray([0,0,0,0,0,0])
# ObsLons 	= SimLons[ObsY_idx,ObsX_idx] 
# ObsLats 	= SimLats[ObsY_idx,ObsX_idx]
# print(ObsLons)
# print(ObsLats)
# # num_points=10
# #ObsLons=np.random.randint(-5, 5,num_points)
# #ObsLats=np.random.randint(15, 45,num_points)
# # print(ObsLons)
# # print(ObsLats)
# # sys.exit(1)

# ObsIDs = np.arange(ObsLons.shape[0])

# new_obj = vat.mapper()
# new_obj.SimLons 	= SimLons
# new_obj.SimLats 	= SimLats
# new_obj.ObsLons 	= ObsLons
# new_obj.ObsLats 	= ObsLats
# new_obj.ObsIDs 		= ObsIDs
# new_obj.SimMeanQ 	= SimMeanQ
# new_obj.ObsMeanQ	= ObsMeanQ
 

# new_obj.MapBestQ()

# # new_obj.MapRaw()
# # print(new_obj.MapYIdx_raw)
# # print(new_obj.MapXIdx_raw)

# # new_obj.writeMap2File('./FirstTestFile.csv')

# # new_obj.ObsLats = ObsLats

# # print(new_obj.MapXIdx_raw)
# # print(new_obj.MapYIdx_raw)

####################################################################
# REAL TEST CASE
####################################################################

# nc_file_path = '/home/n.wagner/development/VAlidationTool_data/CarinaClim'
# files = sorted(glob.glob(f'{nc_file_path}/pgw*.nc'))
# for file in files:
# 	nc_file = nc.Dataset(file)
# 	pgw = nc_file['pgw'][...]
# 	rlat = nc_file['rlat'][...]
# 	rlon = nc_file['rlon'][...]
# 	rotated_pole = nc_file['rotated_pole']
# 	# np_lat = getattr(rotated_pole, 'grid_north_pole_latitude')
# 	# np_lon = getattr(rotated_pole, 'grid_north_pole_longitude')
# 	np_lat = 39.25 
# 	np_lon = -162.0
# 	SimLats, SimLons = rotate.undo_grid_rotation(rlat,rlon,np_lat,np_lon)
# 	nc_file.close()

# read in pressure head
Toplayerpress_file = '../data/Toplayerpress.npy'

pfb_file_path = '/home/n.wagner/development/VAlidationTool_data/ERA5Clim'
if os.path.isfile(Toplayerpress_file):
	print('read already existing file')
	Toplayerpress = np.load(Toplayerpress_file)
else:
	files = sorted(glob.glob(f'{pfb_file_path}/*.out.press.*.pfb'))
	Toplayerpress = None
	for file in files:
		tmp_pressure_head = read_pfb(file)
		tmp_pressure_head = tmp_pressure_head[-1]
		print(tmp_pressure_head.shape)
		if not Toplayerpress is None:
			Toplayerpress = np.append(Toplayerpress, tmp_pressure_head[np.newaxis,...], axis=0)
		else:
			Toplayerpress = tmp_pressure_head
			# add 'time-axis'
			Toplayerpress = Toplayerpress[np.newaxis, ...]

	np.save(Toplayerpress_file, Toplayerpress)

mask = read_pfb(f'{pfb_file_path}/cordex0.11_1980_08.out.mask.pfb')
Toplayerpress = np.where(mask[-1]==99999, Toplayerpress, np.nan)


slopex = read_pfb(f'{pfb_file_path}/slopex.pfb')
slopey = read_pfb(f'{pfb_file_path}/slopey.pfb')
mannings = 5.5E-5

shape2D = (slopey.shape[0], slopex.shape[1])
# flowx= np.zeros(shape2D)
# flowy= np.zeros(shape2D)
# dirx = np.where(slopex > 0.0, -1.0, 1.0)  
# diry = np.where(slopey > 0.0, -1.0, 1.0)

Ponding = np.where(Toplayerpress>0,Toplayerpress,0.0)
#We need only the positive pressure values and set the rest to zero, which results in zero overland flow
flowx = np.absolute(slopex[0,:,:])**(1./2.)/mannings * Ponding**(5./3.)
flowy = np.absolute(slopey[0,:,:])**(1./2.)/mannings * Ponding**(5./3.)
# flowx = dirx * (np.absolute(slopex[0,:,:]))**(1./2.)/mannings * Ponding**(5./3.)
# flowy = diry * (np.absolute(slopey[0,:,:]))**(1./2.)/mannings * Ponding**(5./3.)

# Flow x Edgelengh / 3600 = m^3 / sec
total_discharge = (flowx * 12500) + (flowy * 12500) / 3600
SimMeanQ = np.nanmean(total_discharge, axis=0)

nc_file = '/home/n.wagner/development/VAlidationTool_data/ERA5Clim/lffd1980010100c.nc'
nc_dataset = nc.Dataset(nc_file)
SimLats = nc_dataset['lat'][4:-4, 4:-4]
SimLons = nc_dataset['lon'][4:-4, 4:-4]


# NEW DATASET
# get related files
file_path = '/home/n.wagner/development/VAlidationTool_data/hydroobs_discharge_GRDC/o.data/GRDC'
files = sorted(glob.glob(f'{file_path}/*.day'))
#indexHeader, indexList = vat.toolBox.create_GRDCIndex(files)

# 1) index the dataset 
#    there is no need for a return value...
indexHeader_v1, indexList_v1 = vat.toolBox.read_GRDCIndex(indexFile='../data/index_GRDC_USER.csv')
# print(f'check if headers are identical: {indexHeader == indexHeader_v1}')
# print(f'check if lists are identical: {indexList == indexList_v1}')

# 2) read a Index file if already existing (want to save those locally)
GRDCHeader, GRDCStations = vat.toolBox.filter_GRDCIndex(indexFile='../data/index_GRDC_USER.csv', key='Country', value='DE')
# GRDCHeader, GRDCStations = vat.toolBox.filter_GRDCIndex(indexFile='../data/index_GRDC_USER.csv', key='GRDC-No', value='6335060')
print(f'check if headers are identical: {GRDCHeader == indexHeader_v1}')
print(f'check if lists are identical: {GRDCStations == indexList_v1}')

print(indexHeader_v1)

ObsGRDC = vat.toolBox.read_GRDCFiles(GRDCIndexObj=(GRDCHeader, GRDCStations))

# print(ObsGRDC)
print(ObsGRDC['6335060']['data'])

ObsLats 	= []
ObsLons 	= []
ObsIDs  	= []
ObsMeanQ 	= []
for key in ObsGRDC.keys():
	ObsLats.append(ObsGRDC[key]['lat'])
	ObsLons.append(ObsGRDC[key]['lon'])
	ObsIDs.append(key)
	ObsMeanQ.append(np.nanmean(ObsGRDC[key]['data']))
ObsLats 	= np.asarray(ObsLats, dtype=float)
ObsLons 	= np.asarray(ObsLons, dtype=float)
ObsIDs  	= np.asarray(ObsIDs, dtype=int)
ObsMeanQ 	= np.asarray(ObsMeanQ)

BestFit = vat.mapper()
BestFit.SimLons 	= SimLons
BestFit.SimLats 	= SimLats
BestFit.ObsLons 	= ObsLons
BestFit.ObsLats 	= ObsLats
BestFit.ObsIDs 		= ObsIDs
BestFit.SimMeanQ 	= SimMeanQ
BestFit.ObsMeanQ	= ObsMeanQ

HighestQ = vat.mapper()
HighestQ.SimLons 	= SimLons
HighestQ.SimLats 	= SimLats
HighestQ.ObsLons 	= ObsLons
HighestQ.ObsLats 	= ObsLats
HighestQ.ObsIDs 	= ObsIDs
HighestQ.SimMeanQ 	= SimMeanQ
HighestQ.ObsMeanQ	= ObsMeanQ


HighestQ.MapHighQ()
BestFit.MapBestQ()

# print(f'SimLons: {BestFit.SimLons}')
# print(f'SimLats: {BestFit.SimLats}')
# print(f'ObsLons: {BestFit.ObsLons}')
# print(f'ObsLats: {BestFit.ObsLats}')

# print(f'ObsIDs: {BestFit.ObsIDs}')

# print(f'SimMeanQ: {BestFit.SimMeanQ}')
# print(f'ObsMeanQ: {BestFit.ObsMeanQ}')

# print(f'BestFit.MapXIdx_raw: {BestFit.MapXIdx_raw}')
# print(f'BestFit.MapYIdx_raw: {BestFit.MapYIdx_raw}')
# print(f'BestFit.MapXIdx_fit: {BestFit.MapXIdx_fit}')
# print(f'BestFit.MapYIdx_fit: {BestFit.MapYIdx_fit}')

# search_rad = 1
# for idx, ID in enumerate(BestFit.ObsIDs):
# 	fig = plt.figure()
# 	ax = fig.add_subplot(1, 1, 1)

# 	rawX = BestFit.MapXIdx_raw[idx]
# 	rawY = BestFit.MapYIdx_raw[idx] 
# 	data2plot = BestFit.SimMeanQ[rawY-search_rad:rawY+search_rad+1, rawX-search_rad:rawX+search_rad+1]
# 	data2plot /= np.nanmax(BestFit.SimMeanQ[rawY-search_rad:rawY+search_rad+1, rawX-search_rad:rawX+search_rad+1])
# 	im = ax.imshow(data2plot)
# 	fig.colorbar(im, ax=ax)
# 	# raw data is always the centre = search_rad
# 	ax.scatter(search_rad ,search_rad, c='red', marker='x', label='raw')
# 	x_sliced = BestFit.MapXIdx_fit[idx] - ( rawX - search_rad)
# 	y_sliced = BestFit.MapYIdx_fit[idx] - ( rawY - search_rad)
# 	print(x_sliced, y_sliced)
# 	ax.scatter(x_sliced, y_sliced, c='red', marker='o', label='best')
# 	x_sliced = HighestQ.MapXIdx_fit[idx] - ( rawX - search_rad)
# 	y_sliced = HighestQ.MapYIdx_fit[idx] - ( rawY - search_rad)
# 	print(x_sliced, y_sliced)
# 	ax.scatter(x_sliced, y_sliced, c='green', marker='o', label='highest')
# 	ax.legend()
# 	ax.set_title(f'GRDC-ID: {ID}')
# 	fig.savefig(f'../data/{ID}.png')
# 	plt.close('all')