import numpy as np
import sys
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc

src_path='../src/'
sys.path.append(src_path)
import VAlidationTool as vat



# 1) define where your GRDC dataset is located
file_path = '../data/example_GRDC'
files     = sorted(glob.glob(f'{file_path}/*.mon'))
# daily files are also possible
#files     = sorted(glob.glob(f'{file_path}/*.day'))

# 2) initialize GRDC dataset-object
#    Below step creates a GRDCdataset instance defined with VAT, holding 
#    the data and provides some useful functions to process the data-set.
GRDC_example = vat.GRDCdataset(GRDCfiles=files)

# 3) create a index file if not already exist
#    One of above mentioned useful functions is 'create_indexFile()'. 
#    This function create a easy to read overview of all via 'files' defined
#    GRDC files. The index-file is placed next to this examples_VAT4GRDC.py 
#    script and stored as csv-file. 
#    To get a feeling for the index file, open the index file once it is 
#    created with any text editor or Excel.
GRDC_example.create_indexFile(force=True)

# 4) filter your data
#    As you can see in the index-file, full GRDC datasets contains many 
#    stations all for different time-periods.
#    So some filter functions are provided. In principle you can filter 
#    for any keyword you can find in the header of the index-file.
#    Via default settings the internal index of GRDC_example is updated
#    with this filter-functions, meaning you can apply multiple filters to
#    your dataset which are all added.
#    To be some clear, you can e.g.first filter for a country:
GRDC_example.filter_index(key='Country', value='FR')
#    and after this filter for a time-period - BUT DIFFERENT AS YOU THINK
#    Note here the difference between 'filter-index' and 'filter_index_date':
#    While 'filter_index' filters according to key-values pair straight forward
#    'filter_index_date' does filter those stations not fully covering the 
#    provided time-period. 
GRDC_example.filter_index_date(start='1979-01', end='1980-12', form='%Y-%m')
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
GRDC_example.read_files(start='1979-12-01', end='1980-12-31', form='%Y-%m-%d')

# 6) inspect your data
#    The read in data 
print(f'############################################################################')
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
