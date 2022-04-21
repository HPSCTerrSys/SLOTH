import glob
import os
import sys
import configparser


def get_listOfGriddes():
    """ return a lsit of all availabe griddes definition files
    """
    griddesFiles = glob.glob(f'{os.path.dirname(__file__)}/configs/*_griddes.txt')
    griddesDomainNames = [ item.split('/')[-1] for item in griddesFiles]
    griddesDomainNames = [ item.split('_')[0] for item in griddesDomainNames]
    return griddesDomainNames

def get_griddesDomDef(GridName, returnRaw=False):
    """ returns a individual domain definition from a cdo griddes file

    This function does return a domain defenition from a grides files created
    with CDO

    Input value:
    ------------
    GridName: str
        Name of the target grid

    Return value:
    -------------
    domainDefinition: dict
        A dict containing all parameters needed to describe the requested
        domain stored in individual keys.

    """
    # Read griddes file from configs dir
    # Configs is located under `sloth/` (os.path.dirname(__file__))
    griddesFileName = f'{GridName}_griddes.txt'
    griddesFile     = f'{os.path.dirname(__file__)}/configs/{griddesFileName}'
    if not os.path.isfile(griddesFile):
        print(f'ERROR: There is no griddes file with name {griddesFile} --> EXIT')
        sys.exit(1)

    fileDict = {}
    with open(griddesFile) as f:
        for line in f:
            # skip line if '#' is in as this is a comment
            if '#' in line: continue
            lineItems = line.split('=')
            lineKey   = lineItems[0].strip()
            lineValue = lineItems[1].strip()
            fileDict[lineKey] = lineValue

    # To simply return griddes content
    if returnRaw:
        return fileDict
    # To return formatted dict
    else:
        outDict = {}
        outDict['SWlon'] = float(fileDict['xfirst'])
        outDict['SWlat'] = float(fileDict['yfirst'])
        outDict['dlon']  = float(fileDict['xinc'])
        outDict['dlat']  = float(fileDict['yinc'])
        # Sometime values are like '-162.f', where the f needs to be removed 
        tmp_NPlon        = fileDict['grid_north_pole_longitude'].replace('f','')
        outDict['NPlon'] = float(tmp_NPlon)
        # Sometime values are like '-162.f', where the f needs to be removed 
        tmp_NPlat        = fileDict['grid_north_pole_latitude'].replace('f','')
        outDict['NPlat'] = float(tmp_NPlat)
        outDict['Nlon']  = int(fileDict['xsize'])
        outDict['Nlat']  = int(fileDict['ysize'])

        return outDict

def get_listOfCordexGrids():
    config = configparser.ConfigParser()
    config.read(f'{os.path.dirname(__file__)}/configs/CordexGrid.conf')
    return config.sections()

def get_cordexDomDef(GridName):
    """ returns a cordex domain definition

    This function does return a codex domain defenition in rotated coordinates
    as defined in this document:
    http://is-enes-data.github.io/cordex_archive_specifications.pdf
    The related table is stored under `sloth/configs/`

    Input value:
    ------------
    GridName: str
        A valid name of a cordex domain

    Return value:
    -------------
    domainDefinition: dict
        A dict containing all parameters needed to describe the requested
        cordex domain stored in individual keys.

    """
    # read CORDEX definition from config-file
    config = configparser.ConfigParser()
    config.read('./configs/CordexGrid.conf')
    domainDefinition = {}
    domainDefinition['SWlon'] = float(config[GridName]['West'])
    domainDefinition['SWlat'] = float(config[GridName]['South'])
    domainDefinition['dlon']  = float(config[GridName]['deg'])
    domainDefinition['dlat']  = float(config[GridName]['deg'])
    domainDefinition['NPlon'] = float(config[GridName]['NP_lon'])
    domainDefinition['NPlat'] = float(config[GridName]['NP_lat'])
    domainDefinition['Nlon']  = int(config[GridName]['N_lon'])
    domainDefinition['Nlat']  = int(config[GridName]['N_lat'])

    return domainDefinition

