import numpy as np
import sys
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import VAlidationTool as vat



# 1) define where your GRDC dataset is located
file_path = '../data/example_GRDC'
files     = sorted(glob.glob(f'{file_path}/*.mon'))
# daily files are also possible
#files     = sorted(glob.glob(f'{file_path}/*.day'))

# 2) initialize GRDC dataset-object
#    This instance of the VAT-class holds the data and provides some useful 
#    functions to process the data-set
GRDC_example = vat.GRDCdataset(GRDCfiles=files, GRDCindexFile='../data/example_GRDC/index_GRDC_USER.csv')

# 3) create a index files if not already exist
#    One of above mentioned useful functions is 'create_indexFile()'. 
#    This function create a easy to read overview of all via 'files' defined
#    GRDC files. 
#    The index-file is names according to 'GRDCindexFile' and written in csv
#    format, no matter what you name the file. 
#    Try out and open this file via any text editor or Excel.
#    NOTE: 
#    even if this is a useful function to get an overview, this is also
#    mandatory for all following functions of VAT, so run this function!
#    'force=True' simply forces the function to overwrite already existing
#    index-files. This is easier for the beginning.
GRDC_example.create_indexFile(force=True)

# 4) filter your data
#    As you can see in the index-file, full GRDC datasets contains many 
#    stations all for different time-periods.
#    So some filter functions are provided. In principle you can filter 
#    for any keyword you can find in the header of the index-file, but some
#    keys are more useful than other...
#    Via default settings the index stored with vat.GRDCdataset() is updated
#    with this filter-functions, meaning you can apply multiple filters to
#    your dataset which are all added.
#    To be some clear, you can first filter for a country:
GRDC_example.filter_index(key='Country', value='FR')
#    and after this filter for a time-period - BUT DIFFERENT AS YOU THINK
#    Note here the difference between 'filter-index' and 'filter_index_date':
#    While 'filter_index' filters according to key-values pair straight forward
#    'filter_index_date' does filter those stations not fully covering the 
#    provided time-period. 
GRDC_example.filter_index_date(start='1980-01', end='1980-12', form='%Y-%m')
#    At the end 'GRDC_example' does contain stations locate in France and
#    holding data for the period 1980-08 to 1980-12.

# 5) read the data in
#    Last but not least, you can read the GRDC data in.
#    Before this point you were only collection and filtering meta-data, to
#    keep storage small and performance fast.
#    Note:
#    As 'filter_index_date' do remove those stations without data in given 
#    time-period, here you can really define which time period you want to
#    store with you final data-set
GRDC_example.read_files(start='1980-08-01', end='1980-12-31', form='%Y-%m-%d')

# 6) inspect your data
#    The read in data 
print(f'GRDC_example.id.shape: {GRDC_example.id.shape} ({type(GRDC_example.id)})')
print(f'GRDC_example.id: {GRDC_example.id}')
print(f'GRDC_example.data.shape: {GRDC_example.data.shape} ({type(GRDC_example.data)})')
print(f'GRDC_example.data: {GRDC_example.data}')
print(f'GRDC_example.lats.shape: {GRDC_example.lats.shape} ({type(GRDC_example.lats)})')
print(f'GRDC_example.lats: {GRDC_example.lats}')
print(f'GRDC_example.lons.shape: {GRDC_example.lons.shape} ({type(GRDC_example.lons)})')
print(f'GRDC_example.lons: {GRDC_example.lons}')
print(f'GRDC_example.time.shape: {GRDC_example.time.shape} ({type(GRDC_example.time)})')
print(f'GRDC_example.time: {GRDC_example.time}')
print(f'GRDC_example.meanArea.shape: {GRDC_example.meanArea.shape} ({type(GRDC_example.meanArea)})')
print(f'GRDC_example.meanArea: {GRDC_example.meanArea}')
