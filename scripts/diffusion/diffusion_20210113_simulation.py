# compare horizontal_diffusion and current_uncertainty
# https://opendrift.github.io/gallery/example_horizontal_diffusion.html
# \MPA_connectivity\opendrift-1.4.2\examples\example_horizontal_diffusion.py

# environment: opendrift_mpaconn

######## Model and Readers

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_shape

file_nep = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\processed\NEP36_1h_20110101_20110316.nc'
file_ssc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\salishsea\salishseacast\forcing\SalishSea_1h_20110101_20110316_opendrift.nc'
file_lnd = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\landmask_FINAL_wgs84.shp'

reader_nep = reader_netCDF_CF_unstructured.Reader(file_nep, latstep=0.01, lonstep=0.01, buffer=0.08)
reader_ssc = reader_netCDF_CF_unstructured.Reader(file_ssc, latstep=0.004, lonstep=0.004, buffer=0.08)
reader_lnd = reader_shape.Reader.from_shpfiles(file_lnd)

print(reader_nep)
print(reader_ssc)
print(reader_lnd)


######## First run
o1 = OceanDrift(loglevel=20)
o1.add_reader([reader_lnd, reader_ssc, reader_nep])
from datetime import datetime
from datetime import timedelta
lon = -123.44
lat = 49.12
o1.seed_elements(lon, lat, number=1000, time=reader_ssc.start_time)
o1.set_config('general:use_auto_landmask', False)  # so so important if you want to use your own landmask
o1.set_config('drift:current_uncertainty', 0.22)
o1.set_config('general:coastline_action', 'stranding')
o1.set_config('drift:scheme', 'runge-kutta')
o1.run(
    end_time=reader_ssc.start_time + timedelta(days=1), 
    time_step=60, 
    time_step_output=1800, 
    outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\diffusion\output_1.nc', export_variables=["age_seconds", "land_binary_mask"])
o1.plot(filename='output_1.jpg')

######## Second run
o2 = OceanDrift(loglevel=20)
o2.add_reader([reader_lnd, reader_ssc, reader_nep])
from datetime import datetime
from datetime import timedelta
lon = -123.44
lat = 49.12
o2.seed_elements(lon, lat, number=1000, time=reader_ssc.start_time)
o2.set_config('general:use_auto_landmask', False)  # so so important if you want to use your own landmask
o2.set_config('drift:horizontal_diffusivity', 1.5)
o2.set_config('general:coastline_action', 'stranding')
o2.set_config('drift:scheme', 'runge-kutta')
o2.run(
    end_time=reader_ssc.start_time + timedelta(days=1), 
    time_step=60, 
    time_step_output=1800, 
    outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\diffusion\output_2.nc', export_variables=["age_seconds", "land_binary_mask"])
o2.plot(filename='output_2.jpg')

######### Compare

o2.animation(compare=[o1], filename='compare_1.mp4')

#####################

