# test how particles move at reader boundary
# https://github.com/OpenDrift/opendrift/blob/master/examples/example_reader_boundary.py

# environment: opendrift_mpaconn

######## Model

from opendrift.models.oceandrift import OceanDrift
o = OceanDrift(loglevel=20)


######## Readers

from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_shape

file_nep = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\processed\NEP36_1h_20110101_20110316.nc'
file_ssc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\salishsea\salishseacast\forcing\SalishSea_1h_20110101_20110316_opendrift.nc'
file_lnd = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\landmask_FINAL_wgs84.shp'

reader_nep = reader_netCDF_CF_unstructured.Reader(file_nep, latstep=0.01, lonstep=0.01, buffer=0.08) #0.01, see my measurments in excel file
reader_ssc = reader_netCDF_CF_unstructured.Reader(file_ssc, latstep=0.004, lonstep=0.004, buffer=0.08) #0.004, reducing this by 0.001 from seagrass chapter. It better matches the grid in the north.
# I reduced the buffer size from what is the default amount in the reader. Buffer size needs to be large enough to cover potential movement of a particle within one time step. I originally had the ssc one at 0.04, but I got errors that it was not large enough.
reader_lnd = reader_shape.Reader.from_shpfiles(file_lnd)

print(reader_nep)
print(reader_ssc)
print(reader_lnd)

reader_lnd.plot()
reader_ssc.plot()
reader_nep.plot()

o.add_reader([reader_lnd, reader_ssc, reader_nep])


######### Seed particles

import numpy as np
lats = np.linspace(48.4, 48.6, 50)
lons = np.linspace(-124.68, -125, 50)
lons, lats = np.meshgrid(lons, lats)

o.seed_elements(lons, lats, radius=0, number=2500, time=reader_ssc.start_time)

from datetime import datetime
from datetime import timedelta

# shp = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\simulations\mpa_20210106_testing\mpas_testing\mpas_testing_20210110.shp'
# time_step = timedelta(hours=4)
# num_steps = 1
# for i in range(num_steps):
#     o.seed_from_shapefile(shp, number=500, featurenum=5, time=reader_ssc.start_time + i*time_step)
# o.elements_scheduled

######### Configure and Run

#o.list_configspec()
o.set_config('general:use_auto_landmask', False)  # so so important if you want to use your own landmask
#o.set_config('drift:current_uncertainty', 0.22)
o.set_config('general:coastline_action', 'stranding')
o.set_config('drift:scheme', 'runge-kutta')

#o.run(end_time=reader_ssc.start_time + timedelta(days=1), time_step=60, time_step_output=1800, outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\simulations\mpa_20210106_testing\output_2.nc', export_variables=["age_seconds", "land_binary_mask"])
o.run(end_time=reader_ssc.start_time + timedelta(days=14), time_step=60, time_step_output=1800, outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\simulations\mpa_20210110_testing_boundary\output_2.nc', export_variables=["age_seconds", "land_binary_mask"])

######### Outputs

print(o)
o.plot()
o.animation(filename='output_2.mp4')

#####################

