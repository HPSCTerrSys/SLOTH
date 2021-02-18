import numpy as np
import sys
import os
import glob
import time

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


# pfb_file_path = '/home/n.wagner/development/VAlidationTool_data/LPOClim'
pfb_file_path = '/home/n.wagner/development/VAlidationTool_data/ERA5Clim'
# mask = read_pfb(f'{pfb_file_path}/cordex0.11_1960_04.out.mask.pfb')
mask = read_pfb(f'{pfb_file_path}/cordex0.11_1980_08.out.mask.pfb')
# read in pressure head
Toplayerpress_file = '../data/Toplayerpress.npy'
if os.path.isfile(Toplayerpress_file):
	print('read already existing file')
	Toplayerpress = np.load(Toplayerpress_file)
else:
	files = sorted(glob.glob(f'{pfb_file_path}/*.out.press.*.pfb'))
	Toplayerpress = None
	for file in files:
		tmp_pressure_head = read_pfb(file)
		tmp_pressure_head = tmp_pressure_head[-1]
		tmp_pressure_head = np.where(mask[-1]==99999, tmp_pressure_head, 0.0)
		print(f'max press: {np.nanmax(tmp_pressure_head)}')
		print(f'mean press: {np.nanmean(tmp_pressure_head)}')
		print(tmp_pressure_head.shape)
		if not Toplayerpress is None:
			Toplayerpress = np.append(Toplayerpress, tmp_pressure_head[np.newaxis,...], axis=0)
		else:
			Toplayerpress = tmp_pressure_head
			# add 'time-axis'
			Toplayerpress = Toplayerpress[np.newaxis, ...]


	np.save(Toplayerpress_file, Toplayerpress)

# Toplayerpress = np.where(mask[-1]==99999, Toplayerpress, 0.0)


slopex = read_pfb(f'{pfb_file_path}/slopex.pfb')
slopey = read_pfb(f'{pfb_file_path}/slopey.pfb')
mannings = 5.5E-5

# 329, 166 are taken from BIBIS GRDC dataset for ID 6842900
# vat.toolBox.calc_catchment(slopex[0], slopey[0], 329, 166)
# sys.exit()

# shape2D = (slopey.shape[0], slopex.shape[1])
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
total_discharge =  ( (flowx * 12500) + (flowy * 12500) ) / 3600

# for idx in range(total_discharge.shape[0]):
# 	plt.imshow(total_discharge[idx])
# 	plt.show()
# sys.exit()
SimMeanQ = np.nanmean(total_discharge, axis=0)


nc_file = '/home/n.wagner/development/VAlidationTool_data/ERA5Clim/lffd1980010100c.nc'
nc_dataset = nc.Dataset(nc_file)
SimLats = nc_dataset['lat'][4:-4, 4:-4]
SimLons = nc_dataset['lon'][4:-4, 4:-4]


# NEW DATASET
# get related files
file_path = '/home/n.wagner/development/VAlidationTool_data/hydroobs_discharge_GRDC/o.data/GRDC'
files = sorted(glob.glob(f'{file_path}/*.day'))
# indexHeader, indexList = vat.toolBox.create_GRDCIndex(files)

# 1) index the dataset 
#    there is no need for a return value...
# indexHeader_v1, indexList_v1 = vat.toolBox.read_GRDCIndex(indexFile='../data/index_GRDC_USER.csv')

# 2) read a Index file if already existing (want to save those locally)
GRDC_DE = vat.GRDCdataset(GRDCfiles=files, GRDCindexFile='../data/index_GRDC_USER.csv')
GRDC_DE.create_indexFile()
GRDC_DE.filter_index(key='Country', value='DE')
GRDC_DE.filter_index_date(start='1980-08', end='1980-08')
# GRDC_DE.dump_index()
GRDC_DE.read_files(start='1980-08-01', end='1980-08-31')
print(f'GRDC_DE.id.shape: {GRDC_DE.id.shape}')
print(f'GRDC_DE.data.shape: {GRDC_DE.data.shape}')
print(f'GRDC_DE.lats.shape: {GRDC_DE.lats.shape}')
print(f'GRDC_DE.lons.shape: {GRDC_DE.lons.shape}')
print(f'GRDC_DE.time.shape: {GRDC_DE.time.shape}')
print(f'GRDC_DE.meanArea.shape: {GRDC_DE.meanArea.shape}')

"""
BestFit = vat.mapper()
BestFit.SimLons 	= SimLons
BestFit.SimLats 	= SimLats
BestFit.ObsLons 	= GRDC_DE.lons
BestFit.ObsLats 	= GRDC_DE.lats
BestFit.ObsIDs 		= GRDC_DE.id
BestFit.SimMeanQ 	= SimMeanQ
BestFit.ObsMeanQ	= np.nanmean(GRDC_DE.data, axis=1)

HighestQ = vat.mapper()
HighestQ.SimLons 	= SimLons
HighestQ.SimLats 	= SimLats
HighestQ.ObsLons 	= GRDC_DE.lons
HighestQ.ObsLats 	= GRDC_DE.lats
HighestQ.ObsIDs 	= GRDC_DE.id
HighestQ.SimMeanQ 	= SimMeanQ
HighestQ.ObsMeanQ	= np.nanmean(GRDC_DE.data, axis=1)

BeasArea = vat.mapper()
BeasArea.SimLons 	= SimLons
BeasArea.SimLats 	= SimLats
BeasArea.ObsLons 	= GRDC_DE.lons
BeasArea.ObsLats 	= GRDC_DE.lats
BeasArea.ObsIDs 	= GRDC_DE.id
BeasArea.SimMeanQ 	= SimMeanQ
BeasArea.ObsMeanQ	= np.nanmean(GRDC_DE.data, axis=1)
BeasArea.ObsMeanArea= GRDC_DE.meanArea

HighestQ.MapHighQ()
BestFit.MapBestQ()
BeasArea.MapBestCatchment(slopex=slopex[0], slopey=slopey[0])#, meanArea=BeasArea.ObsMeanArea)

print(f'SimLons: {BestFit.SimLons}')
print(f'SimLats: {BestFit.SimLats}')
print(f'ObsLons: {BestFit.ObsLons}')
print(f'ObsLats: {BestFit.ObsLats}')

print(f'ObsIDs: {BestFit.ObsIDs}')

print(f'SimMeanQ: {BestFit.SimMeanQ}')
print(f'ObsMeanQ: {BestFit.ObsMeanQ}')

print(f'BestFit.MapXIdx_raw: {BestFit.MapXIdx_raw}')
print(f'BestFit.MapYIdx_raw: {BestFit.MapYIdx_raw}')
print(f'BestFit.MapXIdx_fit: {BestFit.MapXIdx_fit}')
print(f'BestFit.MapYIdx_fit: {BestFit.MapYIdx_fit}')

vat.toolBox.plot_MappedSubAreas(HighestQ)
"""

Mapper = vat.mapper(SimLons=SimLons, SimLats=SimLats,
	                ObsLons=GRDC_DE.lons, ObsLats=GRDC_DE.lats,
	                ObsIDs=GRDC_DE.id, 
	                SimMeanQ=SimMeanQ, ObsMeanQ=np.nanmean(GRDC_DE.data, axis=1))

Mapper.MapHighQ()
MapHighQ_X = Mapper.MapXIdx_fit
MapHighQ_Y = Mapper.MapYIdx_fit
Mapper.MapBestQ()
MapBestQ_X = Mapper.MapXIdx_fit
MapBestQ_Y = Mapper.MapYIdx_fit
Mapper.ObsMeanArea = GRDC_DE.meanArea
Mapper.MapBestCatchment(slopex=slopex[0], slopey=slopey[0])
MapBestCatchment_X = Mapper.MapXIdx_fit
MapBestCatchment_Y = Mapper.MapYIdx_fit

save_dir='../data'
for idx, ID in enumerate(Mapper.ObsIDs):
	print(f'--- plotting ObsID {ID}')
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	rawX 	= Mapper.MapXIdx_raw[idx]
	rawY 	= Mapper.MapYIdx_raw[idx] 
	HighQX 	= MapHighQ_X[idx]
	HighQY	= MapHighQ_Y[idx]
	BestQX 	= MapBestQ_X[idx]
	BestQY 	= MapBestQ_Y[idx]
	BestCaX = MapBestCatchment_X[idx]
	BestCaY = MapBestCatchment_Y[idx]

	plot_data = GRDC_DE.data[idx,:]
	GRDC_mean = np.nanmean(plot_data)
	ax.plot(plot_data, label='GRDC')

	plot_data = total_discharge[::8,rawY, rawX]
	raw_mean = np.nanmean(plot_data)
	if (raw_mean // GRDC_mean) < 100:
		ax.plot(plot_data, label='SimRaw')

	plot_data = total_discharge[::8,BestQY, BestQX]
	BestQ_mean = np.nanmean(plot_data)
	if (BestQ_mean // GRDC_mean) < 100:
		ax.plot(plot_data, label='BestQ')

	plot_data = total_discharge[::8,BestCaY, BestCaX]
	BestCA_mean = np.nanmean(plot_data)
	if (BestCA_mean // GRDC_mean) < 100:
		ax.plot(plot_data, label='BestCA')

	plot_data = total_discharge[::8,HighQY, HighQX]
	HighQ_mean = np.nanmean(plot_data)
	if (HighQ_mean // GRDC_mean) < 100:
		ax.plot(plot_data, label='HightQ')

	title_strs = [
	           f'GRDC-ID: {ID}',
	           # f'fit-routine used: {fit_name}',
	           # f'ObsMeanQ: {Mapper.ObsMeanQ[idx]:.2f} m^3/s',
	           # f'SimMeanQ raw: {Mapper.SimMeanQ[Mapper.MapYIdx_raw[idx], Mapper.MapXIdx_raw[idx]]:.2f} m^3/s',
	           # f'SimMeanQ fit: {Mapper.SimMeanQ[Mapper.MapYIdx_fit[idx], Mapper.MapXIdx_fit[idx]]:.2f} m^3/s',
	           # f'ObsMeanArea / SimMeanArea: {Mapper.ObsMeanArea[idx] / Mapper.SimMeanArea[idx]:.2f} m^2/m^2'
	           ]
	ax.set_title('\n'.join(title_strs))
	ax.legend()

	fig.savefig(f'{save_dir}/LinePlot_{ID}.png', bbox_inches='tight', pad_inches=0)
	plt.close('all')