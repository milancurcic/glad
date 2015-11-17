"""
glad.util

Provides utility functions that help manipulate GLAD data.
"""

def argmin_datetime(time0,time_array):
    """
    Returns the index of datetime array time_array element for which
    the datetime time0 is nearest.

    Parameters
    ----------
    time0 : datetime instance
        Time for which the index we seek.
    time_array: ndarray or list
        An array or list of datetime instances.

    Returns
    -------
    n : integer
        Index corresponding to the element of time_array for which
        the datetime instance time0 is nearest.
    """
    from matplotlib.dates import date2num
    from numpy import abs,argmin

    n = argmin(abs(date2num(time0)-date2num(time_array)))
    return n


def get_drifter_ids_from_ascii(filename):
    """
    Returns a list of all drifter ids that are available in the GLAD
    dataset

    Parameters
    ----------
    filename : string
        Path to GLAD ascii data file.

    Returns
    -------
    ids : list
        A list of integers that represent all drifter ids available in 
        the GLAD data file.
    """
    data = [line.rstrip() for line in open(filename,'r').readlines()]

    ids = []
    for line in data[5::100]:
        id = int(line[7:10])
        if not id in ids:
            ids.append(id)

    return ids



def haversine(lon1,lat1,lon2,lat2):
    """
    Returns the great circle distance between two points 
    on a spherical Earth.

    Parameters
    ----------
    lon1 : float
        Longitude (degrees) of first point.
    lat1 : float
        Latitude (degrees) of first point.
    lon2 : float
        Longitude (degrees) of second point.
    lat2 : float
        Latitude (degrees) of second point.

    Returns
    -------
    dist : float
        Distance in kilometers.
    """
    from math import radians,cos,sin,asin,sqrt

    R_Earth = 6.371E6

    # convert decimal degrees to radians 
    lon1,lat1,lon2,lat2 = map(radians,[lon1,lat1,lon2,lat2])

    # Haversine formula 
    dlon = lon2-lon1
    dlat = lat2-lat1
    a = sin(0.5*dlat)**2+cos(lat1)*cos(lat2)*sin(0.5*dlon)**2
    c = 2*asin(sqrt(a))

    # distance in kilometers
    dist = c*R_Earth*1E-3

    return dist


def write_all_to_netcdf(filename,nc_output_path):
    """
    Converts all the drifters from ascii to NetCDF. This is a 
    recommended first step.

    Parameters
    ----------
    filename : string
        Path to GLAD ascii data file.
    nc_output_path : string
        Path for NetCDF output files.   
    """
    from glad import GladDrifter

    glad_ids = get_drifter_ids_from_ascii(filename)

    for id in glad_ids:
        d = GladDrifter(id)
        print 'glad: Converting drifter '+str(id)
        d.read_from_ascii(filename)
        d.write_to_netcdf(nc_output_path
                          + '/glad_drifter_'+'%3.3i' % id+'.nc')
