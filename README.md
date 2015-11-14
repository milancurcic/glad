# glad

Python library for processing and analysis of 
Grand LAgrangian Depoloyment (GLAD) drifter data.
Data was collected during the field experiment by the
[Consortium for Advanced Research on Transport 
of Hydrocarbon in the Environment (CARTHE)]
(http://carthe.org/).
CARTHE is funded by [Gulf of Mexico Research Initiative (GoMRI)]
(http://gulfresearchinitiative.org/)

## Getting started

This code is only useful with original GLAD data.
Download the data from the [GRIIDC repository](https://data.gulfresearchinitiative.org/data/R1.x134.073:0004).

## Example usage

```python
>>> from glad import GladDrifter
>>> d = GladDrifter(1) # initialize the drifter with id 1
>>> d.read_from_ascii('data/GLAD_15min_filtered.dat') # populate the instance with data
>>> d.write_to_netcdf() # output to netcdf
```
