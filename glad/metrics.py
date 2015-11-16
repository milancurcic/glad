"""
glad.metrics

Provides functions for calculation of some metrics useful for GLAD.
"""

def center_of_mass(drifters,time):
    """
    Computes the center of mass latitude and longitude given a list
    or array of GladDrifter instances.

    Parameters
    ----------
    drifters : list, ndarray
        A list or array containing any number of GladDrifter instances.
    time : datetime instance
        A time for which the center_of_mass is desired

    Returns
    -------
    clat, clon : tuple
        A tuple of center latitude and longitude
    """
    import numpy as np
    from glad.util import argmin_datetime

    lats = np.array([])
    lons = np.array([])

    for d in drifters:
        n = argmin_datetime(time,d.time)
        lons = np.append(lons,d.lon[n])
        lats = np.append(lats,d.lat[n])

    clat = np.mean(lats)
    clon = np.mean(lons)

    return clat,clon


def distance(d1,d2,time):
    """
    Calculates distance between drifters d1 and d2 at a given time.

    Parameters
    ----------
    d1 : GladDrifter instance
        Drifter 1.
    d2 : GladDrifter instance
        Drifter 2.
    time : datetime instance
        Time at which the distance between drifters is desired.

    Returns
    -------
    dist : float
        Distance in kilometers.
    """
    from glad.util import argmin_datetime,haversine

    for d in [d1,d2]:
        if time < d.time[0] or time > d.time[-1]:
            raise ValueError('Time not in valid range for '\
                            +'drifter '+str(d.id))

    n1 = argmin_datetime(d1.time,time)
    n2 = argmin_datetime(d2.time,time)

    dist = haversine(d1.lon[n1],d1.lat[n1],d2.lon[n2],d2.lat[n2])

    return dist


def absolute_dispersion(drifters,time,starttime):
    """
    Calculates absolute dispersion A^2, given desired current and 
    initial time.

    Parameters
    ----------
    drifters : list or ndarray
        A list or numpy array of GladDrifter instances.
    time : datetime instance
        Time at which to compute absolute dispersion.
    starttime : datetime instance
        Start time.
    
    Returns
    -------
    A2 : float
        Absolute dispersion in km^2.
    """
    dists = np.array([])

    for d in drifters:

        n1 = argmin_datetime(t,d.time)
        n0 = argmin_datetime(t0,d.time)

        dist_squared = np.append(dist_squared,\
            haversine(d.lon[n1],d.lat[n1],d.lon[n0],d.lat[n0])**2)

    A2 = np.mean(dist_squared)

    return A2


def cloud_dispersion(drifters,time):
    """
    Calculates cloud dispersion C^2 (dispersion relative to center of 
    mass) for a given time.

    Parameters
    ----------
    drifters : list or ndarray
        A list or numpy array of GladDrifter instances.
    time : datetime instance
        Time at which to compute cloud dispersion.
    
    Returns
    -------
    C2 : float
        Cloud dispersion in km^2.
    """
    clat,clon = center_of_mass(drifters,time)

    dist_squared = np.array([])
    for d in drifters:
        n0 = argmin_datetime(t,d.time)
        dists = np.append(dist_squared,\
            haversine(d.lon[n0],d.lat[n0],clon,clat)**2)

    C2 = np.mean(dist_squared)

    return C2




def relative_dispersion(drifters):
    """
    Relative dispersion D^2, averaged over all 2-pair combinations
    in the given drifter set. 

    Important note: Lengths (and thus time) must be equal for all drifters.

    Parameters
    ----------
    drifters : list or ndarray
        A list or numpy array of GladDrifter instances.
    
    Returns
    -------
    D2 : float
        Relative dispersion in km^2.

    """
    from itertools import combinations

    time = drifters[0].time
    dist_squared = []

    for p in combinations(drifters,2):

        pairdist = []

        for n in range(len(time)):
            pairdist.append(haversine(p[0].lon[n],p[0].lat[n],p[1].lon[n],p[1].lat[n])**2)

        dist_squared.append(pairdist)

    dist_squared = np.array(dist_squared)

    D2 = np.mean(dists,axis=0)

    return D2
