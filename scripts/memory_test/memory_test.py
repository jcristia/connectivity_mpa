# determine how much memory gets used when opendrift does the final open of the
# nc file at the end of the simulation
# see basemodel.py and io_netcdf.py

# it looks like the memory load is when it creates a numpy zero array that matches
# the dimensions of the nc file

import numpy as np
import netCDF4 as nc
import os
import tracemalloc

tracemalloc.start()

# test on this 17GB file
od_file = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\fvcom_results\cal03brcl_21_0003.nc'

# these things don't use much memory
infile = nc.Dataset(od_file, 'r')
steps_output = len(infile.dimensions['time'])
start_time = infile.variables['time'][0]
start_time = infile.variables['time'][:]

# I think this is the costly bit from opendrift io script
history = np.ma.array(
    np.zeros([333000, 360]),   # this is about the size of the nc file (particles * timestep outputs)
    dtype=np.dtype('float32'), mask=[True])
uvel = np.ma.array(np.ones([333000, 360]), dtype=np.dtype('float32'), mask=[True])
#history['uvel'] = uvel[list(range(0,333000)), 0:360]

# snapshot = tracemalloc.take_snapshot()
# top_stats = snapshot.statistics('lineno')
# print("[ Top 10 ]")
# for stat in top_stats[:10]:
#     print(stat)

current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage is {current / 1048576}MB; Peak was {peak / 1048576}MB")
tracemalloc.stop()

# tracemalloc output seems a bit variable, but for the above size, it looks like
# peak usage is 33GB
# keep in mind, Opendrift opens this array and then adds variables to it