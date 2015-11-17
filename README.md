# glad-0.1.0

Python library for processing and analysis of 
Grand LAgrangian Depoloyment (GLAD) drifter data.
Data was collected during the field experiment by the
[Consortium for Advanced Research on Transport 
of Hydrocarbon in the Environment (CARTHE)]
(http://carthe.org/).
CARTHE is funded by [Gulf of Mexico Research Initiative (GoMRI)]
(http://gulfresearchinitiative.org/).

## Getting started

### Downloading the GLAD data

This code is only useful with original GLAD data.
Download the data from the [GRIIDC repository](https://data.gulfresearchinitiative.org/data/R1.x134.073:0004).

### Reading and writing data

```python
>>> from glad import GladDrifter
>>> from glad.util import get_drifter_ids_from_ascii
>>> d = GladDrifter(1) # initialize the drifter with id 1
>>> d.read_from_ascii('data/GLAD_15min_filtered.dat') # populate the instance with data
>>> # Let's get all the GLAD drifter IDs available
>>> glad_ids = get_drifter_ids_from_ascii('data/GLAD_15min_filtered.dat')
>>> len(glad_ids) # How many drifters are there?
297
>>> # All the drifters can now be loaded through a loop
>>> drifters = []
>>> for id in glad_ids:
...     d = GladDrifter(id)
...     d.read_from_ascii('data/GLAD_15min_filtered.dat')
```
To work with many or all GLAD drifters, it can get slow to read them
from ascii file every time. Writing the data into NetCDF makes for much
faster reading on demand. It is recommended to process all drifters from 
ascii to NetCDF first:
```python
>>> from glad.util import write_all_to_netcdf
>>> write_all_to_netcdf('data/GLAD_15min_filtered.dat','data')
```
From hereon, we can load the data quickly from NetCDF files:
```python
>>> from glad.util import read_all_from_netcdf
>>> drifters = read_all_from_netcdf('data')
```

### Calculating some metrics

Distance (km) between 2 drifters at time t:
```python
>>> from glad.metrics import *
>>> import datetime
>>> t = datetime.datetime(2012,8,5)
>>> d0 = drifters[0]
>>> d1 = drifters[1]
>>> distance(d0,d1,t)
54.71736600279378
```
Absolute dispersion (km^2) for the same pair from August 1 to August 5:
```python
>>> t0 = datetime.datetime(2012,8,1)
>>> absolute_dispersion([d0,d1],t0,t)
414.97181526613241
```
Absolute dispersion for the whole dataset:
```python
>>> absolute_dispersion(drifters,t0,t)
4091.5763167954137
```
Cloud dispersion (km^2) at time `t`:
```python
>>> cloud_dispersion(drifters,t)
18561.313882889885
```
Relative dispersion (km^2) for the drifter pair at time `t`:
```python
>>> relative_dispersion([d0,d1],t)
2993.9901422836929
```
