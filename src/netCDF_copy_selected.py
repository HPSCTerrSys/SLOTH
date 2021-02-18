import numpy as np
import netCDF4 as nc
import argparse

""" Copy secected variables of given netCDF file.

netCDF does not allow to delete individual variables or dimension. Therefore one
need to copy every variable of one files to another except that one to delete.
This script does the job providing a white- and blacklist option.

Usage:
Use the help function: python SCRIPTNAME -h

Inspiration is taken from the stackoverflow post:
https://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
"""
def netCDF_copy_selected(in_file, out_file, blacklist, whitelist):
    with nc.Dataset(in_file) as src, nc.Dataset(out_file, "w") as dst:
        # in case whitelist is given ignore blaklist and fill blacklist with each
        # variable name in nc-file except of whitelist
        if whitelist:
            blacklist = []
            for key in src.variables:
                if key not in whitelist:
                    blacklist.append(key)
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name not in blacklist:
                x = dst.createVariable(name, variable.datatype, variable.dimensions, zlib=True)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--infile', '-i', type=str, required=True,
                        help='absolut path to in_netCDF file')
    parser.add_argument('--outfile', '-o', type=str, default='./out.nc',
                        help='absolut path to out_netCDF file')
    parser.add_argument('--blacklist', '-b', nargs='+', type=str, default=None,
                        help='which var not to copy')
    parser.add_argument('--whitelist', '-w', nargs='+', type=str, default=None,
                        help='which var to copy')
    args = parser.parse_args()

    in_file = args.infile
    out_file = args.outfile
    blacklist = args.blacklist
    whitelist = args.whitelist

    netCDF_copy_selected(in_file=in_file, out_file=out_file, 
            blacklist=blacklist, whitelist=whitelist)
