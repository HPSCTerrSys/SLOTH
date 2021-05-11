import numpy as np
import sys
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc

sloth_path='../'
sys.path.append(sloth_path)
import sloth


# 1) define where your GRDC dataset is located
file_path = '/p/scratch/cslts/shared_data/tmp_TestDataSet/samples/GRDC'
# provide list of individual GRDC files
files     = sorted(glob.glob(f'{file_path}/*.mon'))
# daily files are also possible
#files     = sorted(glob.glob(f'{file_path}/*.day'))

# 2) initialize GRDC dataset-object
#    Below step creates a GRDCdataset instance defined with VAT, holding 
#    the data and provides some useful functions to process the data-set.
# For mor detailed information about how GRDCdataset() does work, see
# sloth/GRDCdataset.py --> GRDCdataset()
GRDC_example = sloth.GRDCdataset.GRDCdataset(GRDCfiles=files)

# 3) filter your data
#    GRDC datasets contains many stations all for different time-periods.
#    So some filter functions are provided. In principle you can filter 
# NWR 20210429
# ADD FUNCTION TO PRINT INDEX / FILTER KEYS?
#    for any keyword you can find in the header of the index-file.
#    Via default settings the internal index of GRDC_example is updated
#    with this filter-functions, meaning you can apply multiple filters to
#    your dataset which are all added up.
#    To be some clear, you can e.g.first filter for a country:
# For mor detailed information about how filter_index() does work, see
# sloth/GRDCdataset.py --> filter_index()
GRDC_example.filter_index(key='Country', value='DE')
#    and after this filter for a time-period - BUT DIFFERENT AS YOU THINK
#    Note here the difference between 'filter-index' and 'filter_index_date':
#    While 'filter_index' filters according to key-values pair straight forward
#    'filter_index_date' does filter out those stations not fully covering the 
#    provided time-period. 
# For mor detailed information about how filter_index_date() does work, see
# sloth/GRDCdataset.py --> filter_index_date()
GRDC_example.filter_index_date(start='1979-01', end='1980-12', form='%Y-%m')
#    At the end 'GRDC_example' does contain stations locate in Germany and
#    holding data for the period 1980-01 to 1980-12.

# 4) read the data in
#    Last but not least, you can read the GRDC data in.
#    Before this point you were only collection and filtering meta-data, to
#    keep storage small and performance fast.
#    Note:
#    As 'filter_index_date' do remove those stations without data in given 
#    time-period, here you can really define which time period you want to
#    store with you final data-set
# For mor detailed information about how read_files() does work, see
# sloth/GRDCdataset.py --> read_files()
GRDC_example.read_files(start='1980-01-01', end='1980-12-31', form='%Y-%m-%d')

# 5) inspect your data
#    The read in data 
print(f'############################################################################')
print(f'Found {GRDC_example.data.shape[0]} stations with {GRDC_example.data.shape[1]} data-points for applied filter.')
print(f'Possible accessible data:')
print(f'GRDC_example.id.shape: {GRDC_example.id.shape} ({type(GRDC_example.id)})')
# print(f'GRDC_example.id: {GRDC_example.id}')
print(f'GRDC_example.data.shape: {GRDC_example.data.shape} ({type(GRDC_example.data)})')
# print(f'GRDC_example.data: {GRDC_example.data}')
print(f'GRDC_example.lats.shape: {GRDC_example.lats.shape} ({type(GRDC_example.lats)})')
# print(f'GRDC_example.lats: {GRDC_example.lats}')
print(f'GRDC_example.lons.shape: {GRDC_example.lons.shape} ({type(GRDC_example.lons)})')
# print(f'GRDC_example.lons: {GRDC_example.lons}')
print(f'GRDC_example.time.shape: {GRDC_example.time.shape} ({type(GRDC_example.time)})')
# print(f'GRDC_example.time: {GRDC_example.time}')
print(f'GRDC_example.meanArea.shape: {GRDC_example.meanArea.shape} ({type(GRDC_example.meanArea)})')
# print(f'GRDC_example.meanArea: {GRDC_example.meanArea}')
print(f'############################################################################')