import numpy as np
import csv
import sys
import os
from . import toolBox

class mapper:
    ''' [...]
        Due to several reasons (e.g. coarse resolutions) the river corridors
        within TSMP (as all other models) does not have to match the real corridor
        for all pixels along the river. 
        [...]
    '''

    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    def __init__(self, SimLons=None, SimLats=None, 
                       ObsLons=None, ObsLats=None, 
                       ObsIDs=None, 
                       SimMeanQ=None, ObsMeanQ=None,
                       SimMeanArea=None, ObsMeanArea=None):
        ''' Default constructor of mapper-class

        Object-variables:
        -----------------
        SimLons: 2D ndarray
            Two dimensional ndarray containing the lon value for each point of 
            SimGrid
        SimLats: 2D ndarray
            Two dimensional ndarray containing the lat value for each point of 
            SimGrid
        ObsLons: 1D ndarray
            One dimensional ndarray containing the lon values for each individual 
            GRDC station stored with the object
        ObsLats: 1D ndarray
            One dimensional ndarray containing the lat values for each individual 
            GRDC station stored with the object
        ObsID: 1D ndarray
            One dimensional ndarray containing the station ID for each individual 
            GRDC station stored with the object
        SimMeanQ: 2D ndarray
            Two dimensional ndarray containing the mean discharge for each point
            of SimGrid (mean Q for entire period or ref period) 
        ObsMeanQ: 1D ndarray
            One dimensional ndarray containing the mean discharge for each individual 
            GRDC station stored with the object
        SimMeanArea: 1D ndarray
            One dimensional ndarray containing the best fitting catchment area
            for each individual GRDC station stored with the object (calculation 
            based on Sim data)
        ObsMeanArea: 1D ndarray
            One dimensional ndarray containing the mean catchment area of each individual 
            GRDC station stored with the object
        '''
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
            print(f'ObsMeanQ is of dimension {ObsMeanQ.ndim} but dimension 1 is required!')
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
        ''' This is a separate function to keep MapRaw() functions readable.

        This function checks if the passed data fulfills some basic 
        requirements as e.g. if all needed data are defined and if 
        SimMeanQ and SimLons are of same shape.
        This is basically to avoid trivial errors while using this class.

        Return value:
        -------------
        __: boolean 
            True if check is passed, False if some errors are detected.

        '''
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
        ''' This functions straight forward maps OBS on SimGrid.

        This function maps OBS data according to its Lat/Lon values 
        to SimGrid, by calculating the 'real' distance between OBS location
        and all points in SimGrid. The point with the smallest distance is 
        choose as the correct location of OBS on SimGrid. 
        The 'real' distance is hereby the distance in [m] between OBS and SimGrind 
        calculated on a sphere.

        Return value:
        -------------
        No return value!
        This function sets / updates the object variables self.MapYIdx_fit, 
        and self.MapXIdx_fit directly

        '''

        #check if all needed data are already defined
        if not self.checkt4MapRaw():
            print('checkt4MapRaw() failed --> self.MapRaw() canceled!')
            return None

        # Create empty lists for results
        tmp_MapXIdx_raw = []
        tmp_MapYIdx_raw = []
        # Loop over all ObsIDs stored with the object
        for idx, ObsID in enumerate(self.ObsIDs):
            ObsLon = self.ObsLons[idx]
            ObsLat = self.ObsLats[idx]
            # Calculate distance between Obs and SimGrid on Earth surface
            dist = self.spher_dist_v1(np.deg2rad(self.SimLons),
                                      np.deg2rad(self.SimLats),
                                      np.deg2rad(ObsLon),
                                      np.deg2rad(ObsLat))
            # Find index of smallest distance
            mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
            # Temporally store found idx in list
            tmp_MapYIdx_raw.append(mapped_idx[0])
            tmp_MapXIdx_raw.append(mapped_idx[1])

        # update object-variables with found information
        self.MapYIdx_raw = np.array(tmp_MapYIdx_raw)
        self.MapXIdx_raw = np.array(tmp_MapXIdx_raw)

    def check4MapXXX(self):
        ''' This is a separate function to keep MapXXX functions readable.

        This function checks if the passed data fulfills some basic 
        requirements as e.g. SimMeanQ and SimLons are of same shape.
        This is basically to avoid trivial errors while using this class.
        Because self.MapRaw() is called as part of all self.MapBestXXX() 
        functions, this function only checks additional requirements

        Return value:
        -------------
        __: boolean 
            True if check is passed, False if some errors are detected.

        '''
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
        ''' This functions maps OBS on SimGrid by choosing that pixel with best 
        fitting discharge (Q) within given radius.

        The best fitting Q is found by calculating the discharge of
        each pixel within a given radius around the origin pixel.
        The origin pixel is hereby defined by MapRaw().
        That pixel with Q closes to GRDC data is than set.

        Return value:
        -------------
        No return value!
        This function sets / updates the object variables self.MapYIdx_fit, 
        and self.MapXIdx_fit directly

        Parameters
        ----------
        search_rad : int 
            defining the radius around the origin pixel to search for best fitting Q. 
        '''

        # First MapRaw(), than adjust according to best fitting Q
        self.MapRaw()
        #check if all needed data are already defined
        if not self.check4MapXXX():
            print('check4MapXXX() failed --> self.MapBestQ() canceled!')
            return None

        # Create empty lists for results
        tmp_MapXIdx_fit = []
        tmp_MapYIdx_fit = []
        # Loop over all ObsIDs stored with the object
        for idx, ObsID in enumerate(self.ObsIDs):
            # Set the original pixel around which to search
            y = self.MapYIdx_raw[idx]
            x = self.MapXIdx_raw[idx]
            # extract subset of self.SimMeanQ around origin pixel
            # and search radius search_rad (+1 cause last is exclusive)
            sub_SimMeanQ = self.SimMeanQ[y-search_rad:y+search_rad+1, x-search_rad:x+search_rad+1]
            sub_ObsMeanQ = self.ObsMeanQ[idx]
            # Calculate difference between SimMeanQ and ObsMeanQ (GRDC) within subset
            dist = np.abs(sub_SimMeanQ - sub_ObsMeanQ)
            # Get index of best fitting Q (min(dist)) within subset
            mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)

            # Convert subset-index to global index and
            # store (x|y) of best fitting Q in results list
            tmp_realYidx = (y - search_rad) + mapped_idx[0]
            tmp_realXidx = (x - search_rad) + mapped_idx[1]
            tmp_MapYIdx_fit.append(tmp_realYidx)
            tmp_MapXIdx_fit.append(tmp_realXidx)

        # update object-variables with found information
        self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
        self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)

    def MapHighQ(self, search_rad=1):
        ''' This functions maps OBS on SimGrid by choosing that pixel with highest 
        discharge (Q) within given radius.

        The highest discharge (Q) is simply found by applying np.argmax() to 
        the subset (search_rad around origin pixel) of SimMeanQ.
        The origin pixel is hereby defined by MapRaw()
        SimMeanQ is part of the objects-variables

        Return value:
        -------------
        No return value!
        This function sets / updates the object variables self.MapYIdx_fit, 
        and self.MapXIdx_fit directly

        Parameters
        ----------
        search_rad : int 
            defining the radius around the origin pixel to search for highest Q. 
        '''
        # First MapRaw(), than adjust according to highest Q
        self.MapRaw()
        #check if all needed data are already defined
        if not self.check4MapXXX():
            print('check4MapXXX() failed --> self.MapBestQ() canceled!')
            return None

        # Create empty lists for results
        tmp_MapXIdx_fit = []
        tmp_MapYIdx_fit = []
        # Loop over all ObsIDs stored with the object
        for idx, ObsID in enumerate(self.ObsIDs):
            # Set the original pixel around which to search
            y = self.MapYIdx_raw[idx]
            x = self.MapXIdx_raw[idx]
            # extract sub area of self.SimMeanQ according to current ObsID
            # and search radius search_rad (+1 because last is exclusive)
            sub_SimMeanQ = self.SimMeanQ[y-search_rad:y+search_rad+1, x-search_rad:x+search_rad+1]
            # sub_ObsMeanQ = self.ObsMeanQ[idx] # REMOVE
            # Get index of Q within subset
            mapped_idx = np.unravel_index(np.argmax(sub_SimMeanQ, axis=None),sub_SimMeanQ.shape)

            # Convert subset-index to global index and
            # store (x|y) of highest Q in results list
            tmp_realYidx = (y - search_rad) + mapped_idx[0]
            tmp_realXidx = (x - search_rad) + mapped_idx[1]
            tmp_MapYIdx_fit.append(tmp_realYidx)
            tmp_MapXIdx_fit.append(tmp_realXidx)

        # update object-variables with found information
        self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
        self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)

    def MapBestCatchment(self, search_rad=1, dx=12500., dy=12500., slopey=None, slopex=None):
        ''' This functions maps OBS on SimGrid by choosing that pixel which related 
        catchment area fits best to GRDC.

        The best fitting catchment is found by calculating the catchment of
        each pixel within a given radius around the origin pixel.
        The origin pixel is hereby defined by MapRaw().
        That pixel with catchment area closes to GRDC data is than set.

        Return value:
        -------------
        No return value!
        This function sets / updates the object variables self.MapYIdx_fit, 
        self.MapXIdx_fit, and self.SimMeanArea directly

        Parameters
        ----------
        search_rad : int 
            defining the radius around the origin pixel to search for best fitting catchment. 
        dy: float
            defining the y-resolution of slope-grid
        dx: float
            defining the x-resolution of slope-grid
        slopey: 2D ndarray
            defining the ParFlow slopes in y-direction used to calculate the catchment
        slopex: 2D ndarray
            defining the ParFlow slopes in x-direction used to calculate the catchment
        '''

        # First MapRaw(), than adjust according to best fitting catchment-size
        self.MapRaw()
        #check if all needed data are already defined
        if not self.check4MapXXX():
            print('check4MapXXX() failed --> self.MapBestQ() canceled!')
            return None

        # Create empty lists for results
        tmp_MapXIdx_fit = []
        tmp_MapYIdx_fit = []
        tmp_SimManArea = []
        # Loop over all ObsIDs stored with the object
        for idx, ObsID in enumerate(self.ObsIDs):
            # Set the original pixel around which to search
            y        = self.MapYIdx_raw[idx]
            x        = self.MapXIdx_raw[idx]
            # Set catchment area defined by GRDC to compare with
            meanArea = self.ObsMeanArea[idx]

            # Create lists for temporary results
            tmp_catchmentAreaList = []
            tmp_catchmentXList = []
            tmp_catchmentYList = []
            # Loop over all pixel within given search-radius around the 
            # origin pixel, including the original pixel.
            for x_inc in range(-search_rad, search_rad+1):
                for y_inc in range(-search_rad, search_rad+1):
                    tmp_x = x + x_inc
                    tmp_y = y + y_inc
                    # find catchment for (tmp_x|tmp_y) 
                    # for more information on calculation of catchment see:
                    # catcyNAME/toolBox.py --> calc_catchment()
                    tmp_catchmentMask = toolBox.calc_catchment(slopex, slopey, tmp_x, tmp_y)
                    # Calculate area of catchment by multiplying with dx and dy
                    # Than change units from [m^2] to [km^2] to stay compatible to GRDC
                    tmp_catchmentMask *= dx*dy
                    tmp_catchmentMask += 1./(1000.*1000.)
                    tmp_catchmentArea= np.nansum(tmp_catchmentMask)
                    # Add information for current inspected pixel in 
                    # temporary results list
                    tmp_catchmentAreaList.append(tmp_catchmentArea)
                    tmp_catchmentXList.append(tmp_x)
                    tmp_catchmentYList.append(tmp_y)

            # Convert temporary results list to ndarray
            tmp_catchmentArea = np.asarray(tmp_catchmentAreaList)
            tmp_catchmentX = np.asarray(tmp_catchmentXList)
            tmp_catchmentY = np.asarray(tmp_catchmentYList)
            # Now tmp_catchmentArea is a numpy.ndarray holding the catchment-size
            # of the pixel related to X and Y values stored in tmp_catchmentX 
            # and tmp_catchmentY;
            # X and Y are absolute pixels of passed slopex and slopey and NOT of 
            # a subset.

            # Calculate difference between calculated catchment areas and GRDC data
            dist = np.abs( tmp_catchmentArea - meanArea )
            # Get index of best fitting catchment (min(dist)) within temporary 
            # result list
            mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
            # Store information of best fitting pixel with results list 
            tmp_best_fitting_area_x = tmp_catchmentX[mapped_idx[0]]
            tmp_best_fitting_area_y = tmp_catchmentY[mapped_idx[0]]
            tmp_best_fitting_area   = tmp_catchmentArea[mapped_idx[0]]
            tmp_MapYIdx_fit.append(tmp_best_fitting_area_y)
            tmp_MapXIdx_fit.append(tmp_best_fitting_area_x)
            tmp_SimManArea.append(tmp_best_fitting_area)

        # update object-variables with found information
        self.MapYIdx_fit = np.array(tmp_MapYIdx_fit)
        self.MapXIdx_fit = np.array(tmp_MapXIdx_fit)
        self.SimMeanArea = np.array(tmp_SimManArea)

    def writeMap2File(self, file):
        ''' Write mapped coordinates to a given file.

        To save mapped coordinates for later use or better comparison, 
        one can write / dump those to a file.
        Currently this function does only write data t CSV-format.
        Currently already existing files are overwritten.

        Return value:
        -------------
        No return value!

        Parameters
        ----------
        file : str 
            defining the file-path to which data should be written.

        '''
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
