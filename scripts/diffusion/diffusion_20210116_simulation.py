# compare normal current_uncertainty with my custom code for current uncertainty

# environment: opendrift_mpaconn

######## Model and Readers

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_shape

file_nep = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\processed\NEP36_1h_20110101_20110316.nc'
file_ssc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\salishsea\salishseacast\forcing\SalishSea_1h_20110101_20110316_opendrift.nc'
file_lnd = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\landmask_FINAL_wgs84.shp'

reader_nep = reader_netCDF_CF_unstructured.Reader(file_nep, latstep=0.01, lonstep=0.01, buffer=0.08, name='NEP')
reader_ssc = reader_netCDF_CF_unstructured.Reader(file_ssc, latstep=0.004, lonstep=0.004, buffer=0.08, name='SSC')
reader_lnd = reader_shape.Reader.from_shpfiles(file_lnd)


######## First run (set current_uncertainty)
o1 = OceanDrift(loglevel=20)
o1.add_reader([reader_lnd, reader_ssc, reader_nep])
from datetime import datetime
from datetime import timedelta
lon = -126.4
lat = 49.00
o1.seed_elements(lon, lat, number=1000, time=reader_ssc.start_time)
o1.set_config('general:use_auto_landmask', False)  # so so important if you want to use your own landmask
o1.set_config('drift:current_uncertainty', 0.81)
o1.set_config('general:coastline_action', 'stranding')
o1.set_config('drift:scheme', 'runge-kutta')
o1.run(
    end_time=reader_ssc.start_time + timedelta(days=1), 
    time_step=60, 
    time_step_output=1800, 
    outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\diffusion\output_3.nc', export_variables=["age_seconds", "land_binary_mask"])
o1.plot(filename='output_3.jpg')

######## Second run (rely on custom hardcoded values)
o2 = OceanDrift(loglevel=20)
o2.add_reader([reader_lnd, reader_ssc, reader_nep])
from datetime import datetime
from datetime import timedelta
lon = -126.4
lat = 49.00
o2.seed_elements(lon, lat, number=1000, time=reader_ssc.start_time)
o2.set_config('general:use_auto_landmask', False)  # so so important if you want to use your own landmask
o2.set_config('general:coastline_action', 'stranding')
o2.set_config('drift:scheme', 'runge-kutta')
o2.run(
    end_time=reader_ssc.start_time + timedelta(days=1), 
    time_step=60, 
    time_step_output=1800, 
    outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\diffusion\output_4.nc', export_variables=["age_seconds", "land_binary_mask"])
o2.plot(filename='output_4.jpg')

######### Compare

o2.animation(compare=[o1], filename='compare_2.mp4')

#####################

