import numpy as np
import netCDF4 as nc
import glob
import argparse
import sys
import datetime as dt
import ParFlow_IO as pio

def Pfb2NetCDF(infile, varname, outdir):
    """ Converts given ParFlow output (.pfb) into netCDF format(.nc)

    This function converts given .pfb file into .nc file, while no cf-convention
    is followed. Only 2D and 3D .pfb files are supported, while it is assumed
    that not time-axis is used (usual ParFlow behavior).
    Outputfile name is fix and equal to varname.

    Parameters
    ----------
    infile : str
        Path to input .pfb file.
    varname : str
        Name of variable stored with .pfb file.
    outdir : str
        Path to directory where netCDF file should be stored.

    Returns
    -------
    bool
        allways True.

    """

    pfbFileName = f'{infile}'
    pfbFile = pio.read_pfb(pfbFileName)
    outFile = f'{outdir}/{varname}.nc'
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

    ncVar[...] = pfbFile[...]
    ncfile.close()
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--infile', '-i', type=str, required=True,
                        help='full pfb-file-path')
    parser.add_argument('--varname', '-f', type=str, required=True,
                        help='variable name contained in pfb-file')
    parser.add_argument('--outdir', '-o', type=str, required=True,
                        help='path to dir where to store output')
    args = parser.parse_args()

    infile        = args.infile
    varname       = args.varname
    outdir        = args.outdir
    Pfb2NetCDF(infile=infile, varname=varname, outdir=outdir)

