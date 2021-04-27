import heat as ht
import numpy as np
import netCDF4 as nc
from struct import pack, unpack
import sys

# Import different IO-routines for ParFlow output
src_path='../src/'
sys.path.append(src_path)
import ParFlow_IO as pio
ext_path='../extern/ana_parflow-diagnostics_pythonheat/'
sys.path.append(ext_path)
import IO as htIO

split = None
outRoot = f'../data'
rootTestDataSet = '/p/scratch/cslts/shared_data/tmp_TestDataSet'
inFile = f'{rootTestDataSet}/pfbReadAndWrite/ParFlowTest.nc'

print(f'reading in array')
htArray = ht.load_netcdf(inFile, split=split, variable='press')
print(f'htArray.shape: {htArray.shape:}')
npArray = htArray.numpy()
print(f'npArray.shape: {npArray.shape:}')

htIO.create_pfb(filename=f'{outRoot}/anaHeAT.pfb', dndarray=htArray, 
        delta=(1, 1, 1), subgrids=(1, 1, 1))
pio.create_pfb(filename=f'{outRoot}/orig.pfb', var=npArray, 
        delta=(1, 1, 1), subgrids=(1, 1, 1))

readHtArray = pio.read_pfb(filename=f'{outRoot}/anaHeAT.pfb')
readNpArray = pio.read_pfb(filename=f'{outRoot}/orig.pfb')

print(f'readHtArray and readNpArray equal?')
test1 = (readHtArray==readNpArray).all()
print(f'--> {test1}')

print(f'readHtArray and htArray equal?')
test2 = (readHtArray==htArray).all()
print(f'--> {test2}')
