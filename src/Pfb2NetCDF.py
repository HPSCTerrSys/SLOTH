import numpy as np
import netCDF4 as nc
import glob
import argparse
import sys
import os
import datetime as dt
import ParFlow_IO as pio

def Pfb2NetCDF(infiles, varname, outfile):
    """ Converts given ParFlow output (.pfb) into netCDF format(.nc)

    This function converts given .pfb file into .nc file, while no cf-convention
    is followed. Only 2D and 3D .pfb files are supported, while it is assumed
    that not time-axis is used (usual ParFlow behavior).
    Outputfile name is fix and equal to varname.

    Parameters
    ----------
    infiles : list of str
        Path to input .pfb file.
    varname : str
        Name of variable stored with .pfb file.
    outfile : str
        Path to output netCDF file.

    Returns
    -------
    bool
        allways True.

    """
    if not isinstance(infiles, list):
        print(f'infiles is of type {type(infiles)} but <class "list"> is required!')
        return None

    # handle default outfile behavior
    if outfile is None:
        outFileName = os.path.splitext(f'{infile[0]}')[0]
        outFile = f'{outFileName}.nc'
    else:
        outFile = f'{outfile}'
    
    tmp_data = []
    for infile in infiles:
        pfbFileName = f'{infile}'
        data = pio.read_pfb(pfbFileName)
        tmp_data.append(data)

    pfbFile = np.asarray(tmp_data)
    #pfbFile = pio.read_pfb(pfbFileName)
    
    outFile = f'{outfile}'
    ncfile = nc.Dataset(outFile,'w')

    if pfbFile.ndim == 2:
        ny, nx = pfbFile.shape
        yDim = ncfile.createDimension('y', ny)
        xDim = ncfile.createDimension('x', nx)
        ncVar = ncfile.createVariable(f'{varname}','f4',('y', 'x',))
    elif pfbFile.ndim == 3:
        nz, ny, nx = pfbFile.shape
        zDim = ncfile.createDimension('z', nz)
        yDim = ncfile.createDimension('y', ny)
        xDim = ncfile.createDimension('x', nx)
        ncVar = ncfile.createVariable(f'{varname}','f4',('z', 'y', 'x',))
    elif pfbFile.ndim == 4:
        nt, nz, ny, nx = pfbFile.shape
        tDim = ncfile.createDimension('time', nt)
        zDim = ncfile.createDimension('z', nz)
        yDim = ncfile.createDimension('y', ny)
        xDim = ncfile.createDimension('x', nx)
        ncVar = ncfile.createVariable(f'{varname}','f4',('time', 'z', 'y', 'x',))

    ncVar[...] = pfbFile[...]
    ncfile.close()
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--infiles', '-i', type=str, nargs='+', required=True,
                        help='full pfb-file-path')
    parser.add_argument('--varname', '-f', type=str, required=True,
                        help='variable name contained in pfb-file')
    parser.add_argument('--outfile', '-o', type=str, default=None,
                        help='path to output file, if not set, input filename is taken (.py-->.nc)')
    args = parser.parse_args()

    infiles        = args.infiles
    varname        = args.varname
    outfile        = args.outfile
    print(f'infiles: {infiles}')

    # Notice outfile vs outFile
    Pfb2NetCDF(infiles=infiles, varname=varname, outfile=outfile)

