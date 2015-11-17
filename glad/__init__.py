#!/usr/bin/env python

"""
glad is a Python module that provides a base GladDrifter
class with methods to read the drifter data from ascii
file and read and write data from/to NetCDF files.
"""

class GladDrifter():
    """
    A GLAD drifter class with a method to write out 
    data into a NetCDF file.
    """
    def __init__(self,id):
        """
        GladDrifter constructor.

        Parameters
        ----------
        id : int
            GLAD drifter ID number, larger than 0 and smaller
            than 451. Some values are not available.
        """
        self.id = int(id)


    def has_time(self,time):
        """
        Checks if time is within GLAD drifter time window.

        Parameters
        ----------
        time : datetime instance
            Time to be checked.

        Returns
        -------
        valid : boolean
            True if self.time[0] <= time <= self.time[-1], otherwise
            False.
        """ 
        valid = self.time[0] <= time <= self.time[-1]
        
        return valid


    def read_from_ascii(self,filename):
        """
        Reads GLAD drifter data from ascii file.

        Parameters
        ----------
        filename : character string
            Path to GLAD data file name (typically GLAD_15min_filtered.dat)
        """ 
        import numpy as np
        from datetime import datetime

        data = [line.rstrip() for line in open(filename,'r').readlines()\
                if 'CARTHE_'+'%3.3i' % self.id in line]

        if len(data) == 0:
            raise ValueError('No GLAD drifter ''%3.3i' % self.id\
                                        +' found in '+filename)

        self.time = []
        self.lon = []
        self.lat = []
        self.u = []
        self.v = []
        self.pos_error = []
        self.vel_error = []

        for line in data:
            line = line.split()
            self.time.append(datetime.strptime(line[1]+line[2][:8],\
                                 '%Y-%m-%d%H:%M:%S'))
            self.lat.append(np.float(line[3]))
            self.lon.append(np.float(line[4]))
            self.pos_error.append(np.float(line[5]))
            self.u.append(np.float(line[6]))
            self.v.append(np.float(line[7]))
            self.vel_error.append(np.float(line[8]))

        self.time = np.array(self.time)
        self.lat = np.array(self.lat)
        self.lon = np.array(self.lon)
        self.pos_error = np.array(self.pos_error)
        self.u = np.array(self.u)
        self.v = np.array(self.v)
        self.vel_error = np.array(self.vel_error)

    
    def read_from_netcdf(self,filename):
        """
        Reads GLAD data from NetCDF file.
 
        Parameters
        ----------
        filename : string
            Path to NetCDF file from which to read the drifter.
        """
        from matplotlib.dates import num2date
        from netCDF4 import Dataset
        import numpy as np

        nc = Dataset(filename,'r',format='NETCDF3_CLASSIC')
        self.time = np.array([x.replace(tzinfo=None) for x \
                             in num2date(nc.variables['Time'][:])])
        self.lon = nc.variables['lon'][:]
        self.lat = nc.variables['lat'][:]
        self.u   = nc.variables['u'][:]
        self.v   = nc.variables['v'][:]
        self.pos_error = nc.variables['pos_error'][:]
        self.vel_error = nc.variables['vel_error'][:]
        nc.close()


    def write_to_netcdf(self,filename=None):
        """
        Writes data into NetCDF file of the name glad_drifter_id.nc.

        Parameters
        ----------
        filename : string, optional
            A relative path to the output filename. If not provided
            a default filename glad_drifter_id.nc is written.
        """
        from netCDF4 import Dataset
        from matplotlib.dates import date2num

        if not filename:
            filename = 'glad_drifter_'+'%3.3i' % self.id+'.nc'

        nc = Dataset(filename,'w',format='NETCDF3_CLASSIC')
        nc.createDimension('Time',size=0)
        nc.createVariable('Time','f8',dimensions=('Time'))[:] = date2num(self.time)
        nc.createVariable('lon','f8',dimensions=('Time'))[:] = self.lon
        nc.createVariable('lat','f8',dimensions=('Time'))[:] = self.lat
        nc.createVariable('pos_error','f8',dimensions=('Time'))[:] = self.pos_error
        nc.createVariable('u','f8',dimensions=('Time'))[:] = self.u
        nc.createVariable('v','f8',dimensions=('Time'))[:] = self.v
        nc.createVariable('vel_error','f8',dimensions=('Time'))[:] = self.vel_error
        nc.close()


    def travel_distance(self,starttime=None,endtime=None):
        """ 
        Returns a distance traveled by a drifter between 
        starttime and endtime, in kilometers.

        Parameters
        ----------
        starttime : datetime instance, optional
            Start time.
        endtime : datetime instance, optional
            End time.
        """
        from glad.util import argmin_datetime,haversine

        if not starttime:
            starttime = self.time[0]
        
        if not endtime:
            endtime = self.time[-1]
        
        # start and end indices on drifter's time coordinate
        n0 = argmin_datetime(starttime,self.time)
        n1 = argmin_datetime(endtime,self.time)

        #if d.time[n1] < endtime-timedelta(hours=1):
        #    return 0
        dist = 0
        for n in range(n0+1,n1):
            dist += haversine(self.lon[n],self.lat[n],\
                              self.lon[n-1],self.lat[n-1])

        return dist
