import numpy as np
import csv
import sys
import os
import datetime as dt

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
