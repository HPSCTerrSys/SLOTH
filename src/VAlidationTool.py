import numpy as np
import csv
import sys
import os
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
                       SimMeanQ=None, ObsMeanQ=None,
                       SimMeanArea=None, ObsMeanArea=None):
        # input
        self.SimLons        = SimLons
        self.SimLats        = SimLats
        self.ObsLons        = ObsLons
        self.ObsLats        = ObsLats

        self.ObsIDs         = ObsIDs

        self.SimMeanQ       = SimMeanQ
        self.ObsMeanQ       = ObsMeanQ

        # output
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

        self.ObsMeanArea    = None
        self.SimMeanArea    = None

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
        self.__SimLons      = SimLons
        # set mapped idx to None again, since the SimLons was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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
        self.__SimLats      = SimLats
        # set mapped idx to None again, since the SimLats was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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
        self.__ObsLons      = ObsLons
        # set mapped idx to None again, since the ObsLons was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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
        self.__ObsLats      = ObsLats
        # set mapped idx to None again, since the ObsLats was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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
        self.__ObsIDs       = ObsIDs
        # set mapped idx to None again, since the ObsIDs was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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
        self.__SimMeanQ     = SimMeanQ
        # set mapped idx to None again, since the SimMeanQ was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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
        self.__ObsMeanQ     = ObsMeanQ
        # set mapped idx to None again, since the ObsMeanQ was changed!
        self.MapXIdx_raw    = None
        self.MapYIdx_raw    = None
        self.MapXIdx_fit    = None
        self.MapYIdx_fit    = None

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

            tmp_MapYIdx_fit.append(tmp_realYidx)
            tmp_MapXIdx_fit.append(tmp_realXidx)

        self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
        self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)

    def MapBestCatchment(self, search_rad=1, slopey=None, slopex=None):#, meanArea=None):
        """Map the passed OBS on the SimGrid by choosing that pixel within the search_rad
        which related catchment area fits best.

        First MapRaw, than adjust according to Q
        """
        self.MapRaw()
        #check if all needed data are already defined
        if not self.check4MapBestQ():
            print('check4MapBestQ() failed --> self.MapBestQ() canceled!')
            return None

        tmp_MapXIdx_fit = []
        tmp_MapYIdx_fit = []
        tmp_SimManArea = []
        #for i,j in zip(self.MapXIdx_raw, self.MapYIdx_raw):
        for idx, ObsID in enumerate(self.ObsIDs):
            y           = self.MapYIdx_raw[idx]
            x           = self.MapXIdx_raw[idx]
            meanArea    = self.ObsMeanArea[idx]

            tmp_catchmentAreaList = []
            tmp_catchmentXList = []
            tmp_catchmentYList = []
            for x_inc in range(-search_rad, search_rad+1):
                for y_inc in range(-search_rad, search_rad+1):
                    tmp_x = x + x_inc
                    tmp_y = y + y_inc
                    tmp_catchmentMask = toolBox.calc_catchment(slopex, slopey, tmp_x, tmp_y)
                    tmp_catchmentMask *= 12.5*12.5
                    tmp_catchmentArea= np.nansum(tmp_catchmentMask)
                    tmp_catchmentAreaList.append(tmp_catchmentArea)
                    tmp_catchmentXList.append(tmp_x)
                    tmp_catchmentYList.append(tmp_y)

            tmp_catchmentArea = np.asarray(tmp_catchmentAreaList)
            tmp_catchmentX = np.asarray(tmp_catchmentXList)
            tmp_catchmentY = np.asarray(tmp_catchmentYList)
            # now tmp_catchmentArea is a numpy.ndarray holding the catchmentsize of the pixel
            # related to X and Y values stored in tmp_catchmentX and tmp_catchmentY
            # X and Y are absolute pixels and NOT of a subset.

            dist = np.abs( tmp_catchmentArea - meanArea )
            mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
            # print(f'mapped_idx best catchment {mapped_idx}')

            tmp_best_fitting_area_x = tmp_catchmentX[mapped_idx[0]]
            tmp_best_fitting_area_y = tmp_catchmentY[mapped_idx[0]]
            tmp_best_fitting_area   = tmp_catchmentArea[mapped_idx[0]]

            tmp_MapYIdx_fit.append(tmp_best_fitting_area_y)
            tmp_MapXIdx_fit.append(tmp_best_fitting_area_x)
            tmp_SimManArea.append(tmp_best_fitting_area)

        self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
        self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)
        self.SimMeanArea = np.array(tmp_SimManArea)

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


class coreDataset():
    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    def __init__(self, data=None):
        self.data               = data
        self.dataDim            = None
        pass

    @property
    def data(self):
        return self.__data
    @data.setter
    def data(self, data):
        """ data is expected to be an XD ndarray """
        if data is None:
            print(f'Initialize data as NoneType')
            self.__data = None
            return None
        if not isinstance(data, np.ndarray):
            print(f'data is of type {type(data)} but <class "numpy.ndarray"> is required!')
            self.__data = None
            return None

        self.__data         = data
        self.dataDim        = data.ndim


class GRDCdataset(coreDataset):
    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    def __init__(self, data=None, GRDCfiles=None, GRDCindexFile=None,
                       GRDCindexObj=None):
        # print('before')
        super().__init__(data)
        # print('after')
        self.GRDCfiles          = GRDCfiles if GRDCfiles is not None else None
        self.GRDCindexFile      = GRDCindexFile if GRDCindexFile is not None else './index_GRDC_DEFAULT.csv'
        self.GRDCindexObj       = None 

        self.id                 = None
        self.data               = None
        self.lats               = None
        self.lons               = None
        self.time               = None
        self.meanArea           = None

        # 'Time series' is handled special what is what is why keywords are handled special
        self.GRDCkeywords       = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 'Catchment area', 'Time series']
        self.GRDCheader_out     = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 'Catchment area', 'Date start', 'Date end', 'File']

        # force create index-file if GRDCindexFile not passed but GRDCfiles.
        if GRDCfiles is not None and GRDCindexFile is None:
            self.create_indexFile(force=True)

    @property
    def GRDCindexObj(self):
        return self.__GRDCindexObj
    @GRDCindexObj.setter
    def GRDCindexObj(self, GRDCindexObj):
        """ GRDCindexObj is expected to be an XD ndarray """
        if GRDCindexObj is None:
            print(f'Initialize GRDCindexObj as NoneType')
            self.__GRDCindexObj = None
            return None
        if not isinstance(GRDCindexObj, tuple):
            print(f'GRDCindexObj is of type {type(GRDCindexObj)} but <tuple> is required!')
            self.__data = None
            return None
        self.__GRDCindexObj         = GRDCindexObj

        # unset read variables if related index is changed
        # This way I want to prevent from changing the index
        # band getting confused with strange data, belong to another index
        self.id         = None
        self.data       = None
        self.lats       = None
        self.lons       = None
        self.time       = None
        self.meanArea   = None

    ###########################################################################
    ########################## Auxiliary tools ################################
    ###########################################################################

    def create_indexFile(self, files=None, 
                         keywords=None, header_out=None,
                         header_lines=1, meta_lines=40,
                         force=False):

        if os.path.isfile(self.GRDCindexFile) and force == False: 
            print(f'self.GRDCindexFile ({self.GRDCindexFile}) already exist -- read in the existing file instead. Use "force=True" to overwrite')
            self.GRDCindexObj = self.read_indexFile(indexFile=self.GRDCindexFile)
            return None

        files       = files if not files is None else self.GRDCfiles
        keywords    = keywords if not keywords is None else self.GRDCkeywords
        header_out  = header_out if not header_out is None else self.GRDCheader_out
        # 'Time series' is handled special what is what is why keywords are handled special
        # keywords = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 'Catchment area', 'Time series']
        # header_out  = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 'Catchment area', 'Date start', 'Date end', 'File']
        list_out    = []
        # delete index file if already exist to not append same list at the end of file
        with open(self.GRDCindexFile, "w") as fout:
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


            with open(self.GRDCindexFile, "a") as fout:
                wr = csv.writer(fout)
                row_out = [write2csv[key] for key in write2csv.keys() ]
                list_out.append(row_out)
                wr.writerow(row_out)
        self.GRDCindexObj = (header_out, list_out)

    def read_indexFile(self, indexFile=None):
        indexFile = indexFile if not indexFile is None else self.GRDCindexFiled
        with open(indexFile, "r", encoding="utf8", errors='ignore') as fin:
            reader = csv.reader(fin)
            header = next(reader, None)
            data = [list(row) for row in reader]

        self.GRDCindexObj = (header, data)
        return header, data

    def filter_index(self, key, value, indexObj=None, store=True):
        """ this function is filtering an index from indexObj
        """
        # use self.GRDCindexObj is nothing is passed
        indexObj = indexObj if not indexObj is None else self.GRDCindexObj

        header  = indexObj[0]
        key_idx = header.index(key)
        data    = indexObj[1]
        data_filtered = [row for row in data if row[key_idx] == value]

        # only store result in current object if wanted.
        # true is default but to keep the option to say no
        if store:
            self.GRDCindexObj = (header, data_filtered)

        return (header, data_filtered)

    def filter_index_date(self, start, end, indexObj=None, store=True, form='%Y-%m'):
        """ this function is filtering an index from indexObj
        """
        # use self.GRDCindexObj is nothing is passed
        indexObj = indexObj if not indexObj is None else self.GRDCindexObj
        if not isinstance(start, dt.datetime):
            #print(f'start: {start}')
            try:
                start = dt.datetime.strptime(start, form)
            except ValueError:
                print('####################################################################')
                print('ERROR in filter_index_date()')
                print(f'-- passed start-date {start} does not match passed date-form {form}')
                print('-- Forget to pass a form at all?')
                print('####################################################################')
                raise

        if not isinstance(end, dt.datetime):
            #print(f'end: {end}')
            end = dt.datetime.strptime(end, form)

        header      = indexObj[0]
        start_idx   = header.index('Date start')
        end_idx     = header.index('Date end')
        data        = indexObj[1]

        try: 
            tmp_form = '%Y-%m'
            data_filtered = [row for row in data if dt.datetime.strptime(row[start_idx], tmp_form) <= start]
            data_filtered = [row for row in data_filtered if dt.datetime.strptime(row[end_idx], tmp_form) >= end]
        except ValueError:
            tmp_form = '%Y'
            data_filtered = [row for row in data if dt.datetime.strptime(row[start_idx], tmp_form) <= start]
            data_filtered = [row for row in data_filtered if dt.datetime.strptime(row[end_idx], tmp_form) >= end]

        # only store result in current object if wanted.
        # true is default but to keep the option to say no
        if store:
            self.GRDCindexObj = (header, data_filtered)

        return (header, data_filtered)

    def dump_index(self, keys2dump=['GRDC-No', 'River', 'Country',
                                    'Date start', 'Date end']):
        header  = self.GRDCindexObj[0]
        data    = self.GRDCindexObj[1]

        key_idx = [header.index(idx) for idx in keys2dump]

        out_header  = [header[idx] for idx in key_idx]
        print('\t'.join(out_header))
        
        for n, row in enumerate(data):
            tmp_out = [row[idx] for idx in key_idx]
            # truncate long str
            tmp_out = [(item[:10] + '..') if len(item) > 10 else item for item in tmp_out]

            print(f'{n:03}', '\t'.join(tmp_out))

    def read_files(self, start, end, indexObj=None,
                         metaLines=40, delimiter=';',
                         form='%Y-%m-%d', dischargeKey='Calculated'):
        """reads the date and 'original' data only!
        """
        if not isinstance(start, dt.datetime):
            start = dt.datetime.strptime(start, form)
        if not isinstance(end, dt.datetime):
            end = dt.datetime.strptime(end, form)

        indexObj    = indexObj if not indexObj is None else self.GRDCindexObj
        indexHeader = indexObj[0]
        indexList   = indexObj[1]
        file_idx    = indexHeader.index('File')
        id_idx      = indexHeader.index('GRDC-No')
        lat_idx     = indexHeader.index('Latitude')
        lon_idx     = indexHeader.index('Longitude')
        meanArea    = indexHeader.index('Catchment area')
        start_idx   = indexHeader.index('Date start')
        end_idx     = indexHeader.index('Date end')

        # station_header    = ['YYYY-MM-DD', 'Original']
        # time_idx = station_header.index('YYYY-MM-DD')
        # data_idx = station_header.index('Original')

        out = {}
        tmp_out_id   = []
        tmp_out_data = []
        tmp_out_lats = []
        tmp_out_lons = []
        tmp_out_time = []
        tmp_out_meanArea = []
        for station in indexList:
            #print(f'start reading GRDC_no: {station[id_idx]}')

            # catch different GRDC time-stamps
            try:
                start_station = dt.datetime.strptime(station[start_idx], '%Y')
                # calculating month between start and start_staion
                sliceStart = (start.year - start_station.year) * 12 + start.month - start_station.month
                #sliceEnd   = end   - dt.datetime.strptime(station[end_idx], '%Y')
            except ValueError:
                sliceStart = (start - dt.datetime.strptime(station[start_idx], '%Y-%m')).days
                #sliceEnd   = end   - dt.datetime.strptime(station[end_idx], '%Y-%m')
            #print(sliceStart)

            #print(f'sliceStart: {sliceStart}')
            with open(station[file_idx], "r", encoding="utf8", errors='ignore') as f:
                # skip meta data
                _ = [next(f) for x in range(metaLines)]
                # set file pointer
                reader = csv.reader(f, delimiter=delimiter)
                # read header
                tmp_header = next(reader, None)

                # slice not needed data
                _ = [next(f) for x in range(sliceStart)]

                tmp_header = [ entry.strip() for entry in tmp_header ]
                # get needed /  wanted index
                time_idx = tmp_header.index('YYYY-MM-DD')
                # print('time_idx:', time_idx)
                data_idx = tmp_header.index(f'{dischargeKey}')
                # print('data_idx:', data_idx)
                # loop over all data
                tmp_time = []
                tmp_data = []
                for row in reader:
                    #print(f'row: {row}')
                    # print(row)
                    try:
                        tmp = dt.datetime.strptime(row[time_idx], '%Y-%m-%d')
                    except ValueError:
                        # monthly files does contain date format as: YYYY-mm-00
                        # but datetime cannot handle day=00, so cut it 
                        tmp_date = row[time_idx].split('-')
                        tmp_date = '-'.join(tmp_date[:-1])
                        tmp = dt.datetime.strptime(tmp_date, '%Y-%m')
                        end = end.replace(day=1)
                    tmp_time.append(tmp)
                    tmp_data.append(row[data_idx].strip())
                    #print(f'tmp: {tmp} vs end: {end}')
                    if tmp >= end:
                        break
                    # print(tmp_time)
                #print(f'-- found {len(tmp_time)} data-points')
                # tmp_time = np.asarray(tmp_time)
                # tmp_data = np.asarray(tmp_data, dtype=float)
                ## mask missing values
                #tmp_data[tmp_data==-999] = np.nan
                #tmp_data[tmp_data==-99] = np.nan

                tmp_out_id.append(station[id_idx])
                tmp_out_data.append(np.asarray(tmp_data, dtype=float))
                tmp_out_lats.append(station[lat_idx])
                tmp_out_lons.append(station[lon_idx])
                tmp_out_time.append(np.asarray(tmp_time, dtype=object))
                tmp_out_meanArea.append(station[meanArea])

        self.id         = np.asarray(tmp_out_id, dtype=int)
        self.data       = np.asarray(tmp_out_data)
        self.lats       = np.asarray(tmp_out_lats, dtype=float)
        self.lons       = np.asarray(tmp_out_lons, dtype=float)
        self.time       = np.asarray(tmp_out_time)
        self.meanArea   = np.asarray(tmp_out_meanArea, dtype=float)
                # out[station[id_idx]] = {'time':tmp_time, 'data':tmp_data, 
             #                            'lat':station[lat_idx] ,'lon':station[lon_idx],
             #                            'meanArea':station[meanArea]}

        return out



class toolBox:
    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    def __init__(self):
        pass

    def calc_catchment(slopex, slopey, x, y):
        """ This function calculates the catchment related to a given pixel
        based on x- and y-slope files

        This algorithem is 

        INPUT: 
            slopex  - slopes in x-direction [2D ndarray]
            slopey  - slopes in y-direction [2D ndarray]
            x       - index in x-direction to calulate catchment from [int]
            y       - index in y-direction to calulate catchment from [int]

        RETURN:
            catchment   - 2D ndarray of the same size as slopex/y. 0=not part of catchment; 1=part of catchment
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
            im = ax.imshow(data2plot)
            fig.colorbar(im, ax=ax)
            # raw data is always the centre = search_rad
            ax.scatter(search_rad ,search_rad, c='red', marker='x', label='raw')
            x_sliced = mapper.MapXIdx_fit[idx] - ( rawX - search_rad)
            y_sliced = mapper.MapYIdx_fit[idx] - ( rawY - search_rad)
            ax.scatter(x_sliced, y_sliced, c='red', marker='o', label='best')
            title_strs = [
                       f'GRDC-ID: {ID}',
                       f'fit-routine used: {fit_name}',
                       f'ObsMeanQ: {mapper.ObsMeanQ[idx]:.2f} m^3/s',
                       f'SimMeanQ raw: {mapper.SimMeanQ[mapper.MapYIdx_raw[idx], mapper.MapXIdx_raw[idx]]:.2f} m^3/s',
                       f'SimMeanQ fit: {mapper.SimMeanQ[mapper.MapYIdx_fit[idx], mapper.MapXIdx_fit[idx]]:.2f} m^3/s',
                       # f'ObsMeanArea / SimMeanArea: {mapper.ObsMeanArea[idx] / mapper.SimMeanArea[idx]:.2f} m^2/m^2'
                       ]
            ax.set_title('\n'.join(title_strs))
            fig.savefig(f'{save_dir}/Mapplot_{ID}.png', bbox_inches='tight', pad_inches=0)
            plt.close('all')


if __name__ == '__main__':
    print('Im there!')