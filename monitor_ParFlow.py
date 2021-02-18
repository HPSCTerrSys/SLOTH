import heat as ht
import sys 
src_path='./src/'
sys.path.append(src_path)
import numpy as np
import os
import math
from Diagnostics import Diagnostics
import IO as io
import matplotlib.pyplot as plt
import configparser as ConfigParser
# import user-defined functions
import monitor_functions as mf


#########################
# read config data 
#########################
config = ConfigParser.ConfigParser()
config_file_name='config/config_file.ini' 
config.read(config_file_name)

# folder containing gzip pfb simulation output for a particular year and month
pfb_folder = config['PATH']['pfb_folder'].strip('"')
# folder containing netcdf parflow simulation output
netcdfs_folder =  config['PATH']['netcdfs_folder'].strip('"')
# choose particular month
variable="satur"

#spatial parameters of the parflow simlulation
dx=int(config['GEOMETRY']['dx'])
dy=int(config['GEOMETRY']['dy'])
dz=int(config['GEOMETRY']['dz'])

nx=int(config['GEOMETRY']['nx'])
ny=int(config['GEOMETRY']['ny'])
nz=int(config['GEOMETRY']['nz'])

dzmult=config['GEOMETRY']['dzmult'].split(',')
dzmult=np.array(dzmult,dtype=float)
dzmult=ht.array(dzmult,dtype=ht.float64)

split=None
name = 'cordex0.11_1980_06'


#########################
# read static information
#########################
sstorage,mask,poro=mf.read_pfb_batch(pfb_folder,name,'.out.specific_storage.pfb','.out.mask.pfb','.out.porosity.pfb')
#mask file containes 3d file (15,424,436) with 0's for sea and 10^5 for land
#set the mask to one in active and zero in inactive regions
mask  = ht.where(mask>0.0, 1.0, mask)
permx,permy,permz=mf.read_pfb_batch(pfb_folder,name,'.out.perm_x.pfb','.out.perm_y.pfb','.out.perm_z.pfb',mask=mask,default_val=0.0)

mannings = None
slopex   = None
slopey   = None
terrainfollowing = False

# fill out permeability 4D tensor
perm = ht.zeros((3,nz,ny,nx),split=split)
perm[0]=permz
perm[1]=permy
perm[2]=permx

shape2D=(ny, nx)
shape3D=(nz, ny, nx)
#shape4D=(nt, nz, ny, nx)

#some other constant values
ssat  = ht.full(shape3D,1.0,dtype=ht.float64,split=None)
nvg   = ht.full(shape3D,2.0,dtype=ht.float64,split=None)
sres  = ht.full(shape3D,0.2,dtype=ht.float64,split=None)
alpha = ht.full(shape3D,1.0,dtype=ht.float64,split=None)


#########################
# apply Diagnostics class 
#########################
diag = Diagnostics(mask, perm, poro, sstorage, ssat, sres, nvg, alpha, mannings, slopex, slopey, dx, dy, dz, dzmult, nx, ny, nz, terrainfollowing, split)

# read parflow pressure output 
press_fpath = "1980_06/press.nc"
press_packed = io.read_nc4( {press_fpath:"press"} ,dir=netcdfs_folder)

# unpacked_data (nt,nz,ny,nx) array DNDarray
press_unpacked= press_packed[press_fpath]["press"]
nt =  press_unpacked.shape[0]

# read and unpack parflow saturation output 
satur_fpath = "1980_06/satur.nc"
satur_packed = io.read_nc4( {satur_fpath:"satur"} ,dir=netcdfs_folder)
satur_unpacked= satur_packed[satur_fpath]["satur"]

# iterate over monthly timesteps, apply Diagnostics class methods on every timestep
subsurf,toplayer,surf=mf.iterate(press_unpacked,satur_unpacked,mask,diag,nt)



