#!/usr/bin/env python

class gladDrifter():
    """
    A GLAD drifter class with a method to write out 
    data into a NetCDF file.
    """
    def __init__(self,id):
        """
        gladDrifter constructor.

        Parameters
        ----------
        id : int
            GLAD drifter ID number, larger than 0 and smaller
            than 451. Some values are not available.
        """
        self.id = int(id)


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
