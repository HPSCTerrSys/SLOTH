""" coordTrafo - submodule of SLOTH

This sub-module contains two inverse functions to switch between rotated and 
non-rotated coordinates. The original version of these functions was written 
by Tobias TESCH at FZJ-IBG3, but was migrated to SLOTH for easier handling and 
less dependencies.
"""
import numpy as np

def undo_grid_rotation(rlat,rlon,np_lat,np_lon):
    '''Undo the rotated pole grid transformation 
    
    Parameters
    ----------
    rlat, rlon : ndarray
        1d arrays of rotated latitude/longitude values (in degree)

    np_lat, np_lon: ndarray
        geographical latitude/longitude of rotated north pole (in degree) 

    Returns
    -------
    lat, lon : ndarray
        2d arrays of geographical latitude/longitude values (in degree) 

    Notes
    -----
    For reference see
    M. Baldauf et al. "Kurze Beschreibung des Lokal-Modells Kürzestfrist 
    COSMO-DE (LMK) und seiner Datenbanken auf dem Datenserver des DWD", pp. 
    21-22 https://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/cosmo_de/cosmo_de_dbbeschr_version_2_4_161124.pdf?__blob=publicationFile&v=4 

    '''
    np_lon = np.radians(np_lon)
    np_lat = np.radians(np_lat)
    rlon,rlat = np.meshgrid(np.radians(rlon),np.radians(rlat))

    lon = np.degrees(np_lon - np.arctan2(np.cos(rlat)*np.sin(rlon),np.sin(rlat)*np.cos(np_lat)-np.sin(np_lat)*np.cos(rlat)*np.cos(rlon)))
    lon[lon<-180] += 360		#To ensure values in the range [-180,180]
    lon[lon>180] -= 360
    lat = np.degrees(np.arcsin(np.sin(rlat)*np.sin(np_lat) + np.cos(rlat)*np.cos(rlon)*np.cos(np_lat))) 

    return lat,lon

def rotate_grid(lat, lon, np_lat, np_lon):
    ''' Transform to rotated pole grid 
    
    Parameters
    ----------
    lat, lon : ndarray
        1d arrays of geographical latitude/longitude values (in degree)
    np_lat, np_lon : double
        geographical latitude/longitude of rotated north pole (in degree) 

    Returns
    -------
    rlat, rlon : ndarray
        2d arrays of rotated latitude/longitude values (in degree) 

    Notes
    -----
    For reference see 
    M. Baldauf et al. "Kurze Beschreibung des Lokal-Modells Kürzestfrist 
    COSMO-DE (LMK) und seiner Datenbanken auf dem Datenserver des DWD", pp. 
    21-22 https://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/cosmo_de/cosmo_de_dbbeschr_version_2_4_161124.pdf?__blob=publicationFile&v=4

    '''
    np_lon = np.radians(np_lon)
    np_lat = np.radians(np_lat)
    lon, lat = np.meshgrid(np.radians(lon),np.radians(lat))

    rlon = np.degrees(np.arctan2(-np.cos(lat)*np.sin(lon - np_lon),
                                 -np.cos(lat)*np.sin(np_lat)*np.cos(lon-np_lon)+np.sin(lat)*np.cos(np_lat)
                      )) 
    # Ensure values in the range [-180,180]
    rlon[rlon<-180] += 360 
    rlon[rlon>180]  -= 360

    rlat = np.degrees(np.arcsin(np.sin(lat)*np.sin(np_lat)+np.cos(lat)*np.cos(np_lat)*np.cos(lon - np_lon)))

    return rlat, rlon


