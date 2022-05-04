import numpy as np
import csv
import sys
import os
import datetime as dt


class GRDCdataset():
    ''' Class aimed to easy handle GRDC datasets.

    This class (GRDCdataset()) is aimed to easily handle GRDC datasets, 
    which are usually stored in individual files per station and therefore
    are not really easy to inspect / to get an overview.

    Therefore this class offers functions to:
    -) read in entire GRC datasets
    -) create and store index-files in CSV-format of the dataset
    -) filter the dataset according to different keys
    -) read in (filtered) datasets in ndarrays 

    The overall principle of this class is to store a index-file like
    variable with each instance of this class -- GRDCindexObj.
    [...]
    '''
    ###########################################################################
    ############################# Definition ##################################
    ###########################################################################
    def __init__(self, GRDCfiles=None, GRDCindexFile=None, GRDCstationsFile=None):
        ''' Default constructor of GRDCdataset-class

        Object-variables:
        -----------------
        GRDCfiles: list
            List of all GRDC-files belonging to the dataset of interesst.
        GRDCindexFile: str
            The full path to the index file. for both: to store index file at,
            or to read index file from.
        
        '''
        self.GRDCfiles          = GRDCfiles
        self.GRDCindexFile      = GRDCindexFile if GRDCindexFile is not None else './index_GRDC_DEFAULT.csv'
        self.GRDCstationsFile    = GRDCstationsFile

        self.GRDCstationsDict   = None
        self.GRDCindexObj       = None
        self.id                 = None
        self.data               = None
        self.lats               = None
        self.lons               = None
        self.time               = None
        self.meanArea           = None

        # 'Time series' is handled special what is what is why keywords are handled special
        # GRDCkeywords: list
        #    list of keywords which the preamble of individual GRDC files should 
        #    be searched
        # GRDCheader_out: list
        #    list of keywords which should be used as header information in 
        #    GRDCindexFile.
        self.GRDCkeywords       = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 
                'Catchment area', 'Field delimiter', 'Meta Lines', 'Time series']
        self.GRDCheader_out     = ['GRDC-No', 'River', 'Station', 'Country', 'Latitude', 'Longitude', 
                'Catchment area', 'Field delimiter', 'Meta Lines', 'Date start', 'Date end', 'File']

        # force create index-file if GRDCindexFile not passed but GRDCfiles.
        if GRDCfiles is not None and GRDCindexFile is None:
            self.create_indexFile(force=True)
        if GRDCstationsFile is not None and self.GRDCstationsDict is None:
            self.read_GRDCstationsFile(delimiter=',')

    @property
    def GRDCindexFile(self):
        return self.__GRDCindexFile
    @GRDCindexFile.setter
    def GRDCindexFile(self, GRDCindexFile):
        """ handle setting of GRDCindexFile"""
        if GRDCindexFile is None:
            print(f'Initialize GRDCindexFile as NoneType')
            self.__GRDCindexFile = None
            return None
        if not isinstance(GRDCindexFile, str):
            print(f'GRDCindexFile is of type {type(GRDCindexFile)} but <class "str"> is required!')
            self.__GRDCindexFile = None
            return None
        self.__GRDCindexFile = GRDCindexFile
        # unset read variables if GRDCindexFile is changed
        # This way I want to prevent getting inconsistent datasets.
        self.id         = None
        self.data       = None
        self.lats       = None
        self.lons       = None
        self.time       = None
        self.meanArea   = None

    @property
    def GRDCfiles(self):
        return self.__GRDCfiles
    @GRDCfiles.setter
    def GRDCfiles(self, GRDCfiles):
        """ handle setting GRDCfiles"""
        if GRDCfiles is None:
            print(f'Initialize GRDCfiles as NoneType')
            self.__GRDCfiles = None
            return None
        if not isinstance(GRDCfiles, list):
            print(f'GRDCfiles is of type {type(GRDCfiles)} but <list> is required!')
            self.__GRDCfiles = None
            return None
        self.__GRDCfiles = GRDCfiles
        # unset read variables if GRDCfiles is changed
        # This way I want to prevent getting inconsistent datasets.
        self.id         = None
        self.data       = None
        self.lats       = None
        self.lons       = None
        self.time       = None
        self.meanArea   = None

    @property
    def GRDCstationsFile(self):
        return self.__GRDCstationsFile
    @GRDCstationsFile.setter
    def GRDCstationsFile(self, GRDCstationsFile):
        """ handle setting GRDCstationsFile"""
        if GRDCstationsFile is None:
            print(f'Initialize GRDCstationsFile as NoneType')
            self.__GRDCstationsFile = None
            return None
        if not isinstance(GRDCstationsFile, str):
            print(f'GRDCstationsFile is of type {type(GRDCfiles)} but <str> is required!')
            self.__GRDCstationsFile = None
            return None
        self.__GRDCstationsFile = GRDCstationsFile
        # unset read variables if GRDCfiles is changed
        # This way I want to prevent getting inconsistent datasets.
        self.id         = None
        self.data       = None
        self.lats       = None
        self.lons       = None
        self.time       = None
        self.meanArea   = None

    ###########################################################################
    ########################## Auxiliary tools ################################
    ###########################################################################
    def create_indexFile(self, force=True):
        ''' This function indicates a GRDC data set

        GRDC data-sets are usually stored in individual files per station, what
        makes those hard to browse for humans, especially if one wants to get 
        an overview of the entire dataset.
        This functions indicates all GRDC files passed to related 
        object-constructor, and dumps those to a indexFile in CSV-format. 
        In addition the same index is kept as an object-variable to easy filter
        the data-set in later use. 

        Return value:
        -------------
        No return value!
            This function sets / updates the object variables GRDCindexObj
            directly and does store information to CSV-file

        Parameters
        ----------
        force: boolean
            True:  force to iterate through GRDCfiles again, even if 
                   GRDCindexFile does already exist
            False: force to read GRDCindexFile if already exist

        '''

        # Check if GRDCindexFile does already exist
        # Check if the user does not force to indicate GRDC data-set again
        # If both matches: read in existing GRDCindexFile instead of 
        # indicating again
        if os.path.isfile(self.GRDCindexFile) and force == False: 
            print(f'self.GRDCindexFile ({self.GRDCindexFile}) already exist -- read in the existing file instead indicatin all files again.') 
            print(f'-- Use keyword "force=True" to overwrite')
            self.read_indexFile()
            # Return None to exit function here
            return None

        # Create list to store results (individual station ID, name, area etc.)
        list_out    = []
        # Open GRDCindexFile and force overwriting to not append same 
        # list at the end of a already existing file
        with open(self.GRDCindexFile, "w") as fout:
            wr = csv.writer(fout)
            # Write header to file
            wr.writerow(self.GRDCheader_out)

        # Loop over all individual GRDC-files
        for single_file in self.GRDCfiles:
            # Create dict for tmp-results
            write2csv = {}
            # Open individual GRDC file
            with open(single_file, 'r', encoding="utf8", errors='ignore') as f:
                # Save only metadata lines 
                # -> better for performance and storage compared
                #    to save entire file
                # Metadata lies does start with "#". 
                # -> end loop is first char is not "#"
                # -> count number of metalines for later read function
                metaLines = 0
                # Loop over all metadata-linesa
                for line in f:
                    # Apply some filter
                    tmp_line = line
                    tmp_line = tmp_line.rstrip("\n")
                    tmp_line = ' '.join(tmp_line.split())
                    tmp_line = tmp_line.split(':')
                    # Check is line is metadata-line
                    if tmp_line[0][0] != '#':
                        # line does not start with "#" -> no metadata-line
                        # breake the for loop
                        break
                    metaLines += 1
                    # Go through all GRDCkeywords and check if those are
                    # part of metadata line. If so store related information
                    '''
                    E.g metadata contains:
                    # Station:               KVEMO NATANEBI
                    # Country:               GE
                    # Latitude (dec. °):       41.929444
                    # Longitude (de. °):       41.803889
                    than write2csv could look like:
                    write2csv['Station']   = KVEMO NATANEBI
                    write2csv['Country']   = GE
                    write2csv['Latitude']  = 41.929444
                    write2csv['Longitude'] = 41.803889
                    etc.
                    '''
                    for key in self.GRDCkeywords:
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
                write2csv['Meta Lines'] = metaLines
                # special treatment for 'Catchment area'
                # This seems to be empty in newer datasets but is shiped with 
                # a GEOjson file. So if 'Catchment area' is empty in individual
                # files, mark as NaN and read in later from extra file
                if write2csv["Catchment area"] == '':
                    write2csv["Catchment area"] = 'NaN'

            # Write found metadata to GRDCindexFile (append at end)
            with open(self.GRDCindexFile, "a") as fout:
                wr = csv.writer(fout)
                row_out = [write2csv[key] for key in self.GRDCheader_out ]
                #row_out = [write2csv[key] for key in write2csv.keys() ]
                # remember found metadata also for GRDCindexObj
                list_out.append(row_out)
                wr.writerow(row_out)

        # set / update GRDCindexObj with found information 
        self.GRDCindexObj = (self.GRDCheader_out, list_out)

    def read_indexFile(self):
        ''' This function reads a GRDCindexFile as created by create_indexFile()

        As we can write GRDCindexFiles with create_indexFile(), we also need a
        function to read in those files.

        Return value:
        -------------
        No return value!
            This function sets / updates the object variables GRDCindexObj
            directly

        Parameters
        ----------
        No parameters needed!

        '''
        # Open GRDCindexFile and read in with csv-lib straight forward.
        with open(self.GRDCindexFile, "r", encoding="utf8", errors='ignore') as fin:
            reader  = csv.reader(fin)
            header  = next(reader, None)
            entries = [list(row) for row in reader]

        self.GRDCindexObj = (header, entries)

    def check_filter_index(self, key, value, operant, indexObj):
        ''' This is a separate function to keep filter_index() functions readable.

        This function checks if the passed data fulfills some basic 
        requirements as e.g. if the passed operant is supported. 
        This is basically to avoid trivial errors while using this class.

        Return value:
        -------------
        __: boolean 
            True if check is passed, False if some errors are detected.

        '''
        # Check if passed operant is supported, while supported operants
        # are hard-coded below
        supportedOperants = ['=', '==', '<', '<=', '>', '>=', 'in']
        if operant not in supportedOperants:
            print(f'ERROR GRDCdataset.filter_index(): operant ({operant}) is not supported ({supportedOperants})')
            return False
        # Check if value is a list in case operant == 'in'
        if not ( (operant == 'in') == (isinstance(value, list)) ):
            print(f'ERROR GRDCdataset.filter_index(): value type ({type(value)}) is not supported by operant = "{operant}"')
            return False
        # Check if key is in header
        try:
            header  = indexObj[0]
            key_idx = header.index(key)
        except ValueError:
            print(f'ERROR GRDCdataset.filter_index(): key ({key}) is not found in indexObj-header ({indexObj[0]})')
            return False

        return True

    def filter_index(self, key, value, operant=None, indexObj=None, store=True):
        ''' This function filters a GRDCindexObj

        To access individual entries of a GRDC data-set one need to filter 
        the GRDCindexObj accordingly. This is archived with the following 
        function where one can pass a key and a value to filter out every
        entry not matching. E.g. key=country and value=DE filters out 
        every station outside of Germany, that the GRDCindexObj does hold
        stations inside Germany only. The key has to be part of the 
        GRDCindexFile-header and therefore of self.GRDCheader_out.

        As the usual behavior of this function is to update / overwrite 
        self.GRDCindexObj with the outcome of the filtering, one might 
        not be able to apply multiple key-value pairs, if those exclude 
        each other. E.g. if one want to filter for two station names 
        separately. Therefore this function does return an GRDCindexObj
        like tuple, which could be saved for intermediate results. An 
        example and more explanation need to be provided...

        Return value:
        -------------
        filtered GRDCindexObj: tuple
            Tuple of two lists, where first entry is the index header, and second entry is the index body
        __ : __
            This function sets / updates the object variables GRDCindexObj
            directly

        Parameters
        ----------
        key: list
            List of string of GRDCindexFile header to filter
        value: list
            List of values (str) which to keep after filtering
        operant: str
            String defining how to filter data. E.g 
            -) Finding all stations inside DE:
            key='Country', value='DE', operant='=='
            -) Finding all stations which catchment is bigger than X:
            key='Catchment area', value=X, operant='>'
        indexObj: tuple
            Tuple of two lists, where first entry is the index header, and second entry is the index body
            This parameter is reserved for more complex filtering AND NOT YET USED!
        store: boolean
            True: set / update self.GRDCindexObj with filtered index
            This parameter is reserved for more complex filtering AND NOT YET USED!

        '''
        # use self.GRDCindexObj if nothing is passed
        indexObj = indexObj if indexObj is not None else self.GRDCindexObj
        #check if passed arguments are correct / expected
        if not self.check_filter_index(key=key, value=value, operant=operant, indexObj=indexObj):
            print('ERROR: check_filter_index() failed --> self.filter_index() canceled!')
            sys.exit(1)

        header  = indexObj[0]
        key_idx = header.index(key)
        entries = indexObj[1]
        if operant in ['=', '==']:
            try:
                entries_filtered = [row for row in entries if float(row[key_idx]) == value]
            except ValueError:
                entries_filtered = [row for row in entries if row[key_idx] == value]
        elif operant == 'in':
            try:
                entries_filtered = [row for row in entries if float(row[key_idx]) in value]
            except ValueError:
                entries_filtered = [row for row in entries if row[key_idx] in value]
        elif operant == '<':
            entries_filtered = [row for row in entries if float(row[key_idx]) < value]
        elif operant == '<=':
            entries_filtered = [row for row in entries if float(row[key_idx]) <= value]
        elif operant == '>':
            entries_filtered = [row for row in entries if float(row[key_idx]) > value]
        elif operant == '>=':
            entries_filtered = [row for row in entries if float(row[key_idx]) >= value]
        else:
            print(f'GRDCdataset.filter_index(): Could not manage to filter the data ERROR')
        # only store result in current object if wanted.
        # true is default but to keep the option to say no
        if store:
            self.GRDCindexObj = (header, entries_filtered)

        return (header, entries_filtered)

    def filter_index_date(self, start, end, indexObj=None, store=True, form='%Y-%m'):
        ''' This function filters a GRDCindexObj for time-periods

        To access individual time periods of a GRDC data-set one need to filter 
        the GRDCindexObj accordingly. This is archived with the following
        function, where one can pass a start and end date to filter out
        every entry not matching the defined time-period. E.g. 
        start='1980-01-01' and end='1980-01-31' would filter out every entry 
        not matching this time-period, that only measurements between 
        1980-01-01 and 1980-01-31 are kept.
        This function acts different compared to filter_index only in the 
        way, that now dates have to be compared, which is different to 
        other keys.

        Return value:
        -------------
        filtered GRDCindexObj: tuple
            Tuple of two lists, where first entry is the index header, and second entry is the index body
        __ : __
            This function sets / updates the object variables GRDCindexObj
            directly

        Parameters
        ----------
        start: 
        end: 
        indexObj: tuple
            Tuple of two lists, where first entry is the index header, and second entry is the index body
            This parameter is reserved for more complex filtering
        store: boolean
            True: set / update self.GRDCindexObj with filtered index
            This parameter is reserved for more complex filtering
        form:

        '''
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
        ''' This function dumps / prints the GRDCindexObj to std-out

        To inspect the current GRDCindexObj it might be useful to easily
        inspect its content. Therefore this functions was written.
        This function prints out all keys of the GRDCindexObj the user 
        wants to, defined by 'keys2dump' to the std-out (terminal).

        Return value:
        -------------
        No return value!

        Parameters
        ----------
        keys2dump: list
            List of GRDCindexObj header keys (e.g. defined in GRDCheader_out)

        '''
        header  = self.GRDCindexObj[0]
        entries = self.GRDCindexObj[1]

        # Find Header-indexes of keys which should be printed
        key_idx = [header.index(idx) for idx in keys2dump]

        # Find header-values of keys which should be printed ...
        out_header  = [header[idx] for idx in key_idx]
        # ... and print those header-values
        print('\t'.join(out_header))
        
        # Print the body of GRDCindexObj belonging to header-values
        for n, row in enumerate(entries):
            tmp_out = [row[idx] for idx in key_idx]
            # Truncate output if string is too long
            tmp_out = ['{:10.10}'.format(item) for item in tmp_out] 
            # tmp_out = [(item[:10] + '..') if len(item) > 10 else item for item in tmp_out]

            print(f'{n:03}', '\t'.join(tmp_out))

    def read_files(self, start, end,
                form='%Y-%m-%d', dischargeKey='Calculated'):
        ''' This function reads GRDC data into a ndarray

        This is the key function of the GRDCdataset-class, which reads in 
        GRDC-data to ndarrays. 
        This function reads in all files related to a entry in GRDCindexObj.
        So to not read in the entire GRDC data-set one has to apply above
        filter-functions before, reducing the entries in GRDCindexObj.


        Return value:
        -------------
        No return value!
            Below object-variables are updated directly:
            self.id         = 1D ndarray of read in station IDs
            self.data       = 2D ndarray of read in discharge data [ID, #datapoints]
            self.lats       = 1D ndarray of read in Lat values 
            self.lons       = 1D ndarray of read in Lon values 
            self.time       = 2D ndarray of read in time-steps [ID, #datapoints]
            self.meanArea   = 1D ndarray of read in meanAreas (catchment)

        Parameters
        ----------
        start: str or datetime-object 
            First time step of time-perios to read in. This needs to be 
            passed as datetime-object or as string. If passed as string, 
            one need also to pass the correct 'form' to convert string to 
            datetie-object
        end: str or datetime-object
            Last time step of time-perios to read in. This needs to be 
            passed as datetime-object or as string. If passed as string, 
            one need also to pass the correct 'form' to convert string to 
            datetie-object
        form: str

            Format string to convert 'start' and 'end' to datetime-object if those are passed as str
        dischargeKey: str
            String: 'Calculated' or 'Original'

        '''
        if dischargeKey not in ['Calculated', 'Original']:
            print(f'ERROR:')
            print(f'dischargeKey revived ({dischargeKey}) is not supported!')
            return None
        if not isinstance(start, dt.datetime):
            start = dt.datetime.strptime(start, form)
        if not isinstance(end, dt.datetime):
            end = dt.datetime.strptime(end, form)

        indexHeader   = self.GRDCindexObj[0]
        indexList     = self.GRDCindexObj[1]
        file_idx      = indexHeader.index('File')
        id_idx        = indexHeader.index('GRDC-No')
        lat_idx       = indexHeader.index('Latitude')
        lon_idx       = indexHeader.index('Longitude')
        meanArea_idx  = indexHeader.index('Catchment area')
        delimiter_idx = indexHeader.index('Field delimiter')
        metaLines_idx = indexHeader.index('Meta Lines')
        start_idx     = indexHeader.index('Date start')
        end_idx       = indexHeader.index('Date end')

        tmp_out_id   = []
        tmp_out_data = []
        tmp_out_lats = []
        tmp_out_lons = []
        tmp_out_time = []
        tmp_out_meanArea = []
        #print(f'indexList: {indexList}')
        for station in indexList:
            #print(f'start reading GRDC_no: {station[id_idx]}')

            #print(f'sliceStart: {sliceStart}')
            with open(station[file_idx], "r", encoding="utf8", errors='ignore') as f:

                # handle file pointer via csv module
                reader = csv.reader(f, delimiter=station[delimiter_idx])

                # skip meta data
                _ = [next(f) for x in range(station[metaLines_idx])]

                # read header
                tmp_header = next(reader, None)
                tmp_header = [ entry.strip() for entry in tmp_header ]
                # get needed /  wanted index
                time_idx = tmp_header.index('YYYY-MM-DD')
                data_idx = tmp_header.index(f'{dischargeKey}')

                # loop over all data
                tmp_time = []
                tmp_data = []
                for row in reader:
                    # skip lines before start-date
                    tmp = self.read_GRDCdate(row[time_idx])
                    if tmp < start:
                        continue
                    # break for lines after end-date
                    elif tmp > end:
                        break

                    #print(f'row: {row}')
                    tmp_time.append(tmp)
                    tmp_data.append(row[data_idx].strip())
                    #print(f'tmp: {tmp} vs end: {end}')

                tmp_out_id.append(station[id_idx])
                tmp_out_data.append(np.asarray(tmp_data, dtype=float))
                tmp_out_lats.append(station[lat_idx])
                tmp_out_lons.append(station[lon_idx])
                #print(f'tmp_time: {tmp_time}')
                tmp_out_time.append(np.asarray(tmp_time, dtype=object))
                # special treatment for catchment area which might be store in
                # extra GRDCstationsFile...
                tmpMeanArea = station[meanArea_idx] if station[meanArea_idx] != 'NaN' else None
                if tmpMeanArea is None and self.GRDCstationsDict is not None:
                    tmpMeanArea = self.GRDCstationsDict[station[id_idx]]['area']
                tmp_out_meanArea.append(tmpMeanArea)

        self.id         = np.asarray(tmp_out_id, dtype=int)
        self.data       = np.asarray(tmp_out_data)
        self.lats       = np.asarray(tmp_out_lats, dtype=float)
        self.lons       = np.asarray(tmp_out_lons, dtype=float)
        self.time       = np.asarray(tmp_out_time)
        self.meanArea   = np.asarray(tmp_out_meanArea, dtype=float)

    def read_GRDCdate(self, date_str):
        """ read GRDC date string

        It is knowen that some GRDC datasets does contain dates in a form where
        the day is indicated with zero. E.g. 1925-01-00.
        day=00 cannot be handled by datetime, wherefore a special treatment is 
        needed and done below.
        """
        try:
            date_out = dt.datetime.strptime(date_str, '%Y-%m-%d')
        except ValueError:
            # monthly files could contain date format as: YYYY-mm-00
            # but datetime cannot handle day=00, so cut it
            tmp_date = date_str.split('-')
            tmp_date = '-'.join(tmp_date[:-1])
            date_out = dt.datetime.strptime(tmp_date, '%Y-%m')
            date_out = date_out.replace(day=1)
        return date_out

    def read_GRDCstationsFile(self, delimiter=','):
        """ read and store GRDC stations file as obj-dict

        Some metadata information are not necessarily stored with the 
        metadatalines of each individual station file, as for example the 
        catchment area. For those situations this function does read the
        GRDC stations file and store as a dict {grdc_no:row} to make those 
        data easily available.
        
        NOTE:
            The GRDC stations files has to be provided as csv file! 
            The originally download format might be .xlsx, which has to be
            converted to .csv. 
            I decided to do so to not get this function overloaded by many
            extern libs, but stay with basic libs (as far as possible).
        """
        with open(self.GRDCstationsFile, 'r') as f:
            dictReader = csv.DictReader(f)
            GRDCstationsDict = {row['grdc_no']:row for row in dictReader}
            GRDCstationsDict['fieldnames'] = dictReader.fieldnames
            self.GRDCstationsDict = GRDCstationsDict

