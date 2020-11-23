import numpy as np
import csv
import sys
import datetime as dt

# REMOVE
import matplotlib.pyplot as plt

class mapper:

	###########################################################################
	############################# Definition ##################################
	###########################################################################
	def __init__(self, SimLons=None, SimLats=None, 
		               ObsLons=None, ObsLats=None, 
		               ObsIDs=None, 
		               SimMeanQ=None, ObsMeanQ=None):
	    # input
		self.SimLons 		= SimLons
		self.SimLats 		= SimLats
		self.ObsLons 		= ObsLons
		self.ObsLats 		= ObsLats

		self.ObsIDs 		= ObsIDs

		self.SimMeanQ 		= SimMeanQ
		self.ObsMeanQ 		= ObsMeanQ

		# output
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def SimLons(self):
		return self.__SimLons
	@SimLons.setter
	def SimLons(self, SimLons):
		""" SimLons is expected to be an 2D ndarray """
		if SimLons is None:
			print(f'Initialize SimLons as NoneType')
			self.__SimLons = None
			return None
		if not isinstance(SimLons, np.ndarray):
			print(f'SimLons is of type {type(SimLons)} but <class "numpy.ndarray"> is required!')
			self.__SimLons = None
			return None
		if not SimLons.ndim == 2:
			print(f'SimLons is of dimension {SimLons.ndim} but dimension 2 is required!')
			self.__SimLons = None
			return None
		self.__SimLons 		= SimLons
		# set mapped idx to None again, since the SimLons was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def SimLats(self):
		return self.__SimLats
	@SimLats.setter
	def SimLats(self, SimLats):
		""" SimLats is expected to be an 2D ndarray """
		if SimLats is None:
			print(f'Initialize SimLats as NoneType')
			self.__SimLats = None
			return None
		if not isinstance(SimLats, np.ndarray):
			print(f'SimLats is of type {type(SimLats)} but <class "numpy.ndarray"> is required!')
			self.__SimLats = None
			return None
		if not SimLats.ndim == 2:
			print(f'SimLats is of dimension {SimLats.ndim} but dimension 2 is required!')
			self.__SimLats = None
			return None
		self.__SimLats 		= SimLats
		# set mapped idx to None again, since the SimLats was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def ObsLons(self):
		return self.__ObsLons
	@ObsLons.setter
	def ObsLons(self, ObsLons):
		""" ObsLons is expected to be an 2D ndarray """
		if ObsLons is None:
			print(f'Initialize ObsLons as NoneType')
			self.__ObsLons = None
			return None
		if not isinstance(ObsLons, np.ndarray):
			print(f'ObsLons is of type {type(ObsLons)} but <class "numpy.ndarray"> is required!')
			self.__ObsLons = None
			return None
		if not ObsLons.ndim == 1:
			print(f'ObsLons is of dimension {ObsLons.ndim} but dimension 1 is required!')
			self.__ObsLons = None
			return None
		self.__ObsLons 		= ObsLons
		# set mapped idx to None again, since the ObsLons was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def ObsLats(self):
		return self.__ObsLats
	@ObsLats.setter
	def ObsLats(self, ObsLats):
		""" ObsLats is expected to be an 2D ndarray """
		if ObsLats is None:
			print(f'Initialize ObsLats as NoneType')
			self.__ObsLats = None
			return None
		if not isinstance(ObsLats, np.ndarray):
			print(f'ObsLats is of type {type(ObsLats)} but <class "numpy.ndarray"> is required!')
			self.__ObsLats = None
			return None
		if not ObsLats.ndim == 1:
			print(f'ObsLats is of dimension {ObsLats.ndim} but dimension 1 is required!')
			self.__ObsLats = None
			return None
		self.__ObsLats 		= ObsLats
		# set mapped idx to None again, since the ObsLats was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def ObsIDs(self):
		return self.__ObsIDs
	@ObsIDs.setter
	def ObsIDs(self, ObsIDs):
		""" ObsIDs is expected to be an 2D ndarray """
		if ObsIDs is None:
			print(f'Initialize ObsIDs as NoneType')
			self.__ObsIDs = None
			return None
		if not isinstance(ObsIDs, np.ndarray):
			print(f'ObsIDs is of type {type(ObsIDs)} but <class "numpy.ndarray"> is required!')
			self.__ObsIDs = None
			return None
		if not ObsIDs.ndim == 1:
			print(f'ObsIDs is of dimension {ObsIDs.ndim} but dimension 1 is required!')
			self.__ObsIDs = None
			return None
		self.__ObsIDs 		= ObsIDs
		# set mapped idx to None again, since the ObsIDs was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def SimMeanQ(self):
		return self.__SimMeanQ
	@SimMeanQ.setter
	def SimMeanQ(self, SimMeanQ):
		""" SimMeanQ is expected to be an 2D ndarray """
		if SimMeanQ is None:
			print(f'Initialize SimMeanQ as NoneType')
			self.__SimMeanQ = None
			return None
		if not isinstance(SimMeanQ, np.ndarray):
			print(f'SimMeanQ is of type {type(SimMeanQ)} but <class "numpy.ndarray"> is required!')
			self.__SimMeanQ = None
			return None
		if not SimMeanQ.ndim == 2:
			print(f'SimMeanQ is of dimension {SimMeanQ.ndim} but dimension 2 is required!')
			self.__SimMeanQ = None
			return None
		self.__SimMeanQ 	= SimMeanQ
		# set mapped idx to None again, since the SimMeanQ was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	@property
	def ObsMeanQ(self):
		return self.__ObsMeanQ
	@ObsMeanQ.setter
	def ObsMeanQ(self, ObsMeanQ):
		""" ObsMeanQ is expected to be an 2D ndarray """
		if ObsMeanQ is None:
			print(f'Initialize ObsMeanQ as NoneType')
			self.__ObsMeanQ = None
			return None
		if not isinstance(ObsMeanQ, np.ndarray):
			print(f'ObsMeanQ is of type {type(ObsMeanQ)} but <class "numpy.ndarray"> is required!')
			self.__ObsMeanQ = None
			return None
		if not ObsMeanQ.ndim == 1:
			print(f'ObsMeanQ is of dimension {ObsMeanQ.ndim} but dimension 2 is required!')
			self.__ObsMeanQ = None
			return None
		self.__ObsMeanQ 	= ObsMeanQ
		# set mapped idx to None again, since the ObsMeanQ was changed!
		self.MapXIdx_raw	= None
		self.MapYIdx_raw	= None
		self.MapXIdx_fit	= None
		self.MapYIdx_fit	= None

	"""
	One can think about additional functions like:
	setGridFromDef
	readNetCDF
	[...]
	or similar
	"""

	###########################################################################
	########################## Auxiliary tools ################################
	###########################################################################

	def spher_dist_v1(self, lon1, lat1, lon2, lat2, Rearth=6373):
		""" calculate the spherical / haversine distance

		Source: https://www.kompf.de/gps/distcalc.html
		This function is supposed to proper handle different shaped coords
		latX and lonX is supposed to be passed in rad

		return 2D ndarray
		"""
		term1 = np.sin(lat1) * np.sin(lat2)
		term2 = np.cos(lat1) * np.cos(lat2)
		term3 = np.cos(lon2 - lon1)
		# tmp_bool = ( (-1 <= (term1+term2*term3)) & ((term1+term2*term3) >= 1) & ((term1+term2*term3) == 0) & (np.isnan(term1+term2*term3)) ).all()
		# print(f'arccos: {tmp_bool}')
		return Rearth * np.arccos(term1+term2*term3)

	def spher_dist_v2(self, lon1, lat1, lon2, lat2, Rearth=6373):
		""" calculate the spherical / haversine distance

		Source: WSH
		This function is supposed to proper handle different shaped coords
		latX and lonX is supposed to be passed in rad

		return 2D ndarray
		"""
		dlon = lon2 - lon1
		dlat = lat2 - lat1
		a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
		c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))  
		return Rearth * c


   	###########################################################################
	########################### Core functions ################################
	########################################################################### 

	def checkt4MapRaw(self):
		"""
		This is a separate function to keep MapRaw readable
		"""
		# are all variables defined?
		if self.SimLons is None:
			print('self.SimLons is not defined yet, but required by function self.MapRaw()!')
			return False
		if self.SimLats is None:
			print('self.SimLats is not defined yet, but required by function self.MapRaw()!')
			return False
		if self.ObsLons is None:
			print('self.ObsLons is not defined yet, but required by function self.MapRaw()!')
			return False
		if self.ObsLats is None:
			print('self.ObsLats is not defined yet, but required by function self.MapRaw()!')
			return False
		if self.ObsIDs is None:
			print('self.ObsIDs is not defined yet, but required by function self.MapRaw()!')
			return False
		# are the variables / shapes reasonable?
		if self.SimLons.shape != self.SimLats.shape:
			print(f'The shape of self.SimLons {self.SimLons.shape} is not equal the shape of self.SimLats {self.SimLats.shape}!')
			return False
		if self.ObsLons.shape != self.ObsLats.shape:
			print(f'The shape of self.ObsLons {self.ObsLons.shape} is not equal the shape of self.ObsLats {self.ObsLats.shape}!')
			return False
		if self.ObsLons.shape != self.ObsIDs.shape:
			print(f'The shape of self.ObsLons / self.ObsLats {self.ObsLons.shape} is not equal the shape of self.ObsIDs {self.ObsIDs.shape}!')
			return False

		return True

	def MapRaw(self):
		"""Map the passed OBS data on the SimGrid
		"""
		#check if all needed data are already defined
		if not self.checkt4MapRaw():
			print('checkt4MapRaw() failed --> self.MapRaw() canceled!')
			return None

		tmp_MapXIdx_raw = []
		tmp_MapYIdx_raw = []
		# loop over ObsCoords via ObsIDs
		for idx, ObsID in enumerate(self.ObsIDs):
			ObsLon = self.ObsLons[idx]
			ObsLat = self.ObsLats[idx]
			# calculate distance between Obs and SimGrid on Earth surface
			dist = self.spher_dist_v1(np.deg2rad(self.SimLons),
				                      np.deg2rad(self.SimLats),
				                      np.deg2rad(ObsLon),
				                      np.deg2rad(ObsLat))
			# find index of smallest distance
			mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
			# temporally store found idx in list
			tmp_MapYIdx_raw.append(mapped_idx[0])
			tmp_MapXIdx_raw.append(mapped_idx[1])

		# store found idx in class variable
		self.MapYIdx_raw = np.array(tmp_MapYIdx_raw)
		self.MapXIdx_raw = np.array(tmp_MapXIdx_raw)

	def check4MapBestQ(self):
		"""
		This is a separate function to keep MapBestQ readable

		Because self.MapRaw() is called as part of self.MapBestQ(), 
		this function only checks additional requirements
		"""
		# are all variables defined?
		if self.SimMeanQ is None:
			print('self.SimMeanQ is not defined yet, but required by function self.MapBestQ()!')
			return False
		if self.ObsMeanQ is None:
			print('self.ObsMeanQ is not defined yet, but required by function self.MapBestQ()!')
			return False
		if self.SimMeanQ.shape != self.SimLons.shape:
			print(f'The shape of self.SimMeanQ {self.SimMeanQ.shape} is not equal the shape of self.SimLons / self.SimLats {self.SimLons.shape}!')
			return False
		if self.ObsMeanQ.shape != self.ObsLons.shape:
			print(f'The shape of self.ObsMeanQ {self.ObsMeanQ.shape} is not equal the shape of self.ObsLons / self.ObsLats {self.ObsLons.shape}!')
			return False

		return True

	def MapBestQ(self, search_rad=1):
		"""Map the passed OBS on the SimGrid by choosing best fitting Q

		First MapRaw, than adjust according to Q
		"""
		self.MapRaw()
		#check if all needed data are already defined
		if not self.check4MapBestQ():
			print('check4MapBestQ() failed --> self.MapBestQ() canceled!')
			return None

		tmp_MapXIdx_fit = []
		tmp_MapYIdx_fit = []
		#for i,j in zip(self.MapXIdx_raw, self.MapYIdx_raw):
		for idx, ObsID in enumerate(self.ObsIDs):
			y = self.MapYIdx_raw[idx]
			x = self.MapXIdx_raw[idx]
			# extract sub area of self.SimMeanQ according to current ObsID
			# and search radius search_rad (+1 cause last is exclusive)
			sub_SimMeanQ = self.SimMeanQ[y-search_rad:y+search_rad+1, x-search_rad:x+search_rad+1]
			sub_ObsMeanQ = self.ObsMeanQ[idx]
			# distance between SimMeanQ and ObsMeanQ within search_rad
			dist = np.abs(sub_SimMeanQ - sub_ObsMeanQ)
			mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)

			tmp_realYidx = (y - search_rad) + mapped_idx[0]
			tmp_realXidx = (x - search_rad) + mapped_idx[1]

			# print(f'sub_SimMeanQ: {sub_SimMeanQ}')
			# print(f'sub_ObsMeanQ: {sub_ObsMeanQ}')
			# print(f'mapped_idx[0]: {mapped_idx[0]}')
			# print(f'mapped_idx[1]: {mapped_idx[1]}')

			# plt.imshow(self.SimMeanQ, vmax=2.08e4)
			# plt.scatter(x, y, c='red', marker='x')
			# plt.scatter(tmp_realXidx, tmp_realYidx, c='red', marker='x')
			# plt.colorbar()
			# plt.show()
			# plt.imshow(sub_SimMeanQ)
			# plt.scatter(mapped_idx[1], mapped_idx[0], c='red', marker='x')
			# plt.colorbar()
			# plt.show()

			tmp_MapYIdx_fit.append(tmp_realYidx)
			tmp_MapXIdx_fit.append(tmp_realXidx)

		self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
		self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)

	def MapHighQ(self, search_rad=1):
		"""Map the passed OBS on the SimGrid by choosing the highes Q value

		First MapRaw, than adjust according to Q
		"""
		self.MapRaw()
		#check if all needed data are already defined
		if not self.check4MapBestQ():
			print('check4MapBestQ() failed --> self.MapBestQ() canceled!')
			return None

		tmp_MapXIdx_fit = []
		tmp_MapYIdx_fit = []
		#for i,j in zip(self.MapXIdx_raw, self.MapYIdx_raw):
		for idx, ObsID in enumerate(self.ObsIDs):
			y = self.MapYIdx_raw[idx]
			x = self.MapXIdx_raw[idx]
			# extract sub area of self.SimMeanQ according to current ObsID
			# and search radius search_rad (+1 cause last is exclusive)
			sub_SimMeanQ = self.SimMeanQ[y-search_rad:y+search_rad+1, x-search_rad:x+search_rad+1]
			sub_ObsMeanQ = self.ObsMeanQ[idx]
			# get index of max Q in search_rad
			mapped_idx = np.unravel_index(np.argmax(sub_SimMeanQ, axis=None),sub_SimMeanQ.shape)

			tmp_realYidx = (y - search_rad) + mapped_idx[0]
			tmp_realXidx = (x - search_rad) + mapped_idx[1]

			# print(f'sub_SimMeanQ: {sub_SimMeanQ}')
			# print(f'sub_ObsMeanQ: {sub_ObsMeanQ}')
			# print(f'mapped_idx[0]: {mapped_idx[0]}')
			# print(f'mapped_idx[1]: {mapped_idx[1]}')

			# plt.imshow(self.SimMeanQ, vmax=2.08e4)
			# plt.scatter(x, y, c='red', marker='x')
			# plt.scatter(tmp_realXidx, tmp_realYidx, c='red', marker='x')
			# plt.colorbar()
			# plt.show()
			# plt.imshow(sub_SimMeanQ)
			# plt.scatter(mapped_idx[1], mapped_idx[0], c='red', marker='x')
			# plt.colorbar()
			# plt.show()

			tmp_MapYIdx_fit.append(tmp_realYidx)
			tmp_MapXIdx_fit.append(tmp_realXidx)

		self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
		self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)

	def writeMap2File(self, file):
		""" write the mapped coordinates to a given file

		currently only csv-format
		"""
		with open(file,'w', newline='') as outFile:
			writer = csv.writer(outFile, delimiter=',')
			header=['ObsID', 
			        'ObsLon', 'ObsLat', 
			        'MapXIdx_raw', 'MapYIdx_raw',
			        'related SimLon', 'related SimLat']
			writer.writerow(header)
			for idx, ObsID in enumerate(self.ObsIDs):
				row = [ObsID, 
				       self.ObsLons[idx], self.ObsLats[idx], 
				       self.MapXIdx_raw[idx], self.MapYIdx_raw[idx],
				       self.SimLons[self.MapXIdx_raw[idx], self.MapYIdx_raw[idx]], self.SimLats[self.MapXIdx_raw[idx], self.MapYIdx_raw[idx]]]
				writer.writerow(row)


class toolBox:
	###########################################################################
	############################# Definition ##################################
	###########################################################################
	def __init__(self):
		# self.files 				= files
		# self.metaDataKeywords 	= metaDataKeywords
		# self.header_lines 		= header_lines
		# self.skipp_lines 		= skipp_lines
		pass

	def create_GRDCIndex(files, 
		                 header_lines=1, meta_lines=40,
		                 outputFile='../data/index_GRDC_USER.csv'):
		# 'Time series' is handled special what is what is why keywords are handled special
		keywords = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 'Time series']
		header_out  = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 'Date start', 'Date end', 'File']
		list_out 	= []
		# delete index file if already exist to not append same lsit at the end of file
		with open(outputFile, "w") as fout:
			wr = csv.writer(fout)
			# write also location where files was found for later usage
			wr.writerow(header_out)

		for single_file in files:
			write2csv = {}

			with open(single_file, 'r', encoding="utf8", errors='ignore') as f:
				metadata = [next(f) for x in range(meta_lines)]
				for line in metadata:
					tmp_line = line
					tmp_line = tmp_line.rstrip("\n")
					tmp_line = ' '.join(tmp_line.split())
					tmp_line = tmp_line.split(':')
					for key in keywords:
						if any(key in x for x in tmp_line):
							tmp_data2write = (tmp_line[-1].strip())
							if not key == 'Time series':
								write2csv[key] = tmp_data2write
							else:
								tmp_data2write = tmp_data2write.split(' - ')
								write2csv['Date start'] = tmp_data2write[0]
								write2csv['Date end'] = tmp_data2write[1]

				# write also location where files was found for later usage
				write2csv['File'] = f'{single_file}'


			with open(outputFile, "a") as fout:
				wr = csv.writer(fout)
				row_out = [write2csv[key] for key in write2csv.keys() ]
				list_out.append(row_out)
				wr.writerow(row_out)

		return header_out, list_out

	def read_GRDCIndex(indexFile):
		with open(indexFile, "r", encoding="utf8", errors='ignore') as fin:
			reader = csv.reader(fin)
			header = next(reader, None)
			data = [list(row) for row in reader]

		return header, data


	def filter_GRDCIndex(indexFile, key, value):
		with open(indexFile, "r", encoding="utf8", errors='ignore') as fin:
			reader = csv.reader(fin)
			header = next(reader, None)

			# return only those where key == value
			key_idx = header.index(key)
			data = [list(row) for row in reader if list(row)[key_idx] == value]

			return header, data

	def read_GRDCFiles(GRDCIndexObj, metaLines=40, delimiter=';'):
		"""reads the date and 'original' data only!
		"""
		indexHeader = GRDCIndexObj[0]
		indexList 	= GRDCIndexObj[1]
		file_idx 	= indexHeader.index('File')
		id_idx 		= indexHeader.index('GRDC-No')
		lat_idx 	= indexHeader.index('Latitude')
		lon_idx 	= indexHeader.index('Longitude')

		# station_header 	= ['YYYY-MM-DD', 'Original']
		# time_idx = station_header.index('YYYY-MM-DD')
		# data_idx = station_header.index('Original')

		out = {}
		for station in indexList:
			with open(station[file_idx], "r", encoding="utf8", errors='ignore') as f:
				# skip meta data
				_ = [next(f) for x in range(metaLines)]
				# set file pointer
				reader = csv.reader(f, delimiter=delimiter)
				# read header
				tmp_header = next(reader, None)
				tmp_header = [ entry.strip() for entry in tmp_header ]
				# get needed /  wanted index
				time_idx = tmp_header.index('YYYY-MM-DD')
				print('time_idx:', time_idx)
				data_idx = tmp_header.index('Original')
				print('data_idx:', data_idx)
				# loop over all data
				tmp_time = []
				tmp_data = []
				for row in reader:
					# print(row)
					tmp_time.append(dt.datetime.strptime(row[time_idx], '%Y-%m-%d'))
					tmp_data.append(row[data_idx].strip())
					# print(tmp_time)
				tmp_time = np.asarray(tmp_time)
				tmp_data = np.asarray(tmp_data, dtype=float)
				tmp_data[tmp_data==-999] = np.nan
				tmp_data[tmp_data==-99] = np.nan

				out[station[id_idx]] = {'time':tmp_time, 'data':tmp_data, 
			                            'lat':station[lat_idx] ,'lon':station[lon_idx]}

		return out


	def extract_Stations(indexFile, key, value):
		with open(indexFile, "r") as fin:
			reader = csv.DictReader(fin, delimiter=",")

			tmp_IDs 	= []
			tmp_Lats 	= []
			tmp_Lons 	= []
			for row in reader:
				if row[key] == value:
					tmp_IDs.append(row['GRDC-No'])
					tmp_Lats.append(row['Latitude'])
					tmp_Lons.append(row['Longitude'])

		return np.asarray(tmp_IDs), np.asarray(tmp_Lats), np.asarray(tmp_Lons)

if __name__ == '__main__':
	print('Im there!')