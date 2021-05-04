# LPo, 2021
# This script calculates discharge from press.nc (from postprocessed ParFlow output after TSMP simulations
# with MPI-ESM-LR r1i1p1 forcing data within Hi-CAM project

import netCDF4 as nc
import datetime as dt
import numpy as np
import gc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os
import matplotlib.colors as colors
from matplotlib.colors import LogNorm

#src_path='/p/project/cesmtst/poshyvailo1/ANALYSIS_postpro_genericVAlidationTool/src/'
#sys.path.append(src_path)
#import SanityCheck

split=None
year_start = 2005
year_final = 2005
month_start_def =1
month_final_def =1

tsmp_var = 'press' 

WORK_DIR = f'/p/scratch/cjibg35/poshyvailo1/HiCam-CORDEX_EUR-11_MPI-ESM-LR_histo_r1i1p1_FZJ-IBG3-TSMP120EC_v00aJuwelsCpuProdTt-1949_2005/'

#####################################################
dx = 12500.0 #[m]
dy = 12500.0 #[m]
dz = 2.0 #[m]
dzmult = [7.50, 7.50, 5.0, 5.0, 2.0, 0.50, 0.35, 0.25, 0.15, 0.10, 0.065, 0.035, 0.025, 0.015, 0.01] #from coup.tcl file of PRL
nx = 436
ny = 424
nz = 15
coeff_x = dx / 3600.0 #3600sec -- 1 hour; [m/h] -> [m/s]
coeff_y = dy / 3600.0 
terrainfollowing = True
dzmult_rearanged=[] #[m]
           
#####################################################
print(f'**********************************************************')
for year in range(year_start, year_final+1):
     if (year == 1949):
           month_start =12
           month_final = 12  
     else:
           month_start = month_start_def
           month_final = month_final_def
           
     for month in range(month_start, month_final+1):        
           print(f'Year of calculation: {year}')
           print(f'Month of calculation: {month}')
           
           path = WORK_DIR + 'postpro/'+ str(year) + '_' +str(month).zfill(2)+'/'
           outDir = WORK_DIR + 'ParFlow_var/'
           
           if( month==1 or month==3 or month==5 or month==7 or month==8 or month==10 or month==12):
                nt = 31*8
           elif(month==2):
                nt = 28*8
           else:
                nt =30*8                                  

           ###################################################################################################################################################
           ############## reading press from "postpro"            
           path_var = path + tsmp_var + '.nc'
           print(f'Variable from postpro-folder is: {tsmp_var}')
           print(f'Reading in the file: {path_var}')  

           #shape4D=(nt+1,nz, ny, nx)
           shape3D_t=(nt+1,ny,nx)
           wtd=np.empty(shape3D_t)
           wtd_temp=np.empty(shape3D_t)

           with nc.Dataset(str(path_var), 'r') as nc_file:
                nc_time = nc_file.variables['time']  
                dates = nc.num2date(nc_time[:],units=nc_time.units,calendar=nc_time.calendar)
                standard_name_time=nc_file.variables['time'].standard_name                              
                long_name_time=nc_file.variables['time'].long_name                                             
                 
                nc_rlat=nc_file.variables[str('rlat')] #424, Y
                units_rlat=nc_file.variables['rlat'].units
                axis_rlat=nc_file.variables['rlat'].axis
                long_name_rlat=nc_file.variables['rlat'].long_name  
                standard_name_rlat=nc_file.variables['rlat'].standard_name              
                                           
                nc_rlon=nc_file.variables[str('rlon')] #436, X
                units_rlon=nc_file.variables['rlon'].units
                axis_rlon=nc_file.variables['rlon'].axis
                long_name_rlon=nc_file.variables['rlon'].long_name  
                standard_name_rlon=nc_file.variables['rlon'].standard_name              
            
                tsmp_var_read = nc_file.variables[tsmp_var]                               
                tsmp_var_read_np=tsmp_var_read[:,:,:,:] 
                print(f'tsmp_var_read.shape {tsmp_var_read.shape}')
                
                var_mask  = np.where(abs(tsmp_var_read[:,:,:,:])<1.e3, tsmp_var_read, np.nan) #mask press values with abs. values >1.e3
                print(var_mask)
                
                ##########################
                for a in range(0,nz+1):
                     dzmult_rearanged.append(dz*sum(dzmult[a:nz]))
                     print(f'{dzmult[a:nz]}')
                     if (a==nz+1):
                          dzmult_rearanged.append(0)
                print(f'Depths of the boundaries of each layer: {dzmult_rearanged}')

                ##########################
                #Elevations of the middle points of each layer --  the distance from the bedrock to the middle points of each layer;
                #d is different for each pixel. For example, if we find P > 0 at the 1st-14th layer, then we use the elevation at the 14th layer as d.
                dzmult_rearanged_np=np.array(dzmult_rearanged, dtype=np.float32)
                d = np.nanmax(dzmult_rearanged) - (dzmult_rearanged_np[:-1] + dzmult_rearanged_np[1:]) / 2
                print(f'Elevations of the middle points of each layer: {d}')
                print(f'{np.nanmax(dzmult_rearanged)}')

                #########################
                positive_location = np.zeros(var_mask.shape)
                positive_location[var_mask>0] = 1 #replace value to "1" where press>0

                overall_index = np.argmin(positive_location,axis=1)-1 #np.argmin returns first index where postitive_location(var_mask)=0 (meaning, the first index where var_mask<0;
                #then "-1" gives the highest level where positive_location(var_mask)=1

                t, rlat, rlon = np.indices(wtd_temp.shape) #return an array representing the indices of a grid
                wtd_temp = dzmult_rearanged_np[0]-(var_mask[t,overall_index[t,rlat,rlon],rlat,rlon]+d[overall_index[t,rlat,rlon]])
   
                # mean over 1 month (all time steps)       
                wtd_mean=np.nanmean(wtd_temp, axis=0)

           #--------
           # plotting
           #--------
           fig, ax = plt.subplots()
           cax = ax.imshow(wtd_mean, cmap=cm.jet,origin='lower', norm=LogNorm(vmin=0.01, vmax=100), interpolation='none') # vmin=0, vmax=10)#,
           cbar = fig.colorbar(cax,extend='max')
           cbar.set_label('unit [m]', rotation=270, labelpad=15, y=0.5) #shift in y-direction: set y
           ax.set_title('Water table depth (WTD): '+ str(year) +'-' +str(month).zfill(2), fontweight='bold')
           plt.grid(linestyle='dotted')

           #-----------------------
           plt.savefig('WTD_'+str(tsmp_var)+'_'+str(year)+'_'+str(month)+'.png', bbox_inches='tight', pad_inches=0, dpi=380)   
           plt.show()    
                
           #file_nc.close()                                   
           gc.collect()                    
##################################################### THE END :-)
