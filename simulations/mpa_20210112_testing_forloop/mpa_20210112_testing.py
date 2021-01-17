# Opendrift particle tracking
# environment: opendrift_mpaconn

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_shape
from datetime import datetime
from datetime import timedelta
from osgeo import ogr
import numpy as np
import os

######## Set up directory structure

root = os.path.dirname(__file__)
outputs = os.path.join(root, 'outputs')
nc_out = os.path.join(outputs, 'nc')
np_out = os.path.join(outputs, 'npy')
logs = os.path.join(outputs, 'logs')
for d in [outputs, nc_out, np_out, logs]:
    if not os.path.exists(d):
        os.makedirs(d)


######## Readers

nc_nep = 'NEP36_1h_20110101_20110316.nc'
nc_ssc = 'SalishSea_1h_20110101_20110316_opendrift.nc'

base = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models'
base_nep = r'nep_nemo\processed'
base_ssc = r'salishsea\salishseacast\forcing'
base_lnd = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline'
shp_lnd = 'landmask_FINAL_wgs84.shp'

file_nep = os.path.join(base, base_nep, nc_nep)
file_ssc = os.path.join(base, base_ssc, nc_ssc)
file_lnd = os.path.join(base_lnd, shp_lnd)
reader_nep = reader_netCDF_CF_unstructured.Reader(
    file_nep, latstep=0.01, lonstep=0.01, buffer=0.08, name='NEP')
reader_ssc = reader_netCDF_CF_unstructured.Reader(
    file_ssc, latstep=0.004, lonstep=0.004, buffer=0.08, name='SSC')
reader_lnd = reader_shape.Reader.from_shpfiles(file_lnd)


######### Simulation for each shapefile

shp_directory = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpas_shp'
shp_files = os.listdir(shp_directory)
shapefiles = []
for file in shp_files:
    if file.endswith('.shp'):
        shapefiles.append(os.path.join(shp_directory, file))

# can break down list so that I can run subsets concurrently
# can still use the same root and the same folders
shapefiles = shapefiles[:3]

for shp in shapefiles:

    ######## Shapefile attributes
    shp = ogr.Open(shp)
    lyr = shp.GetLayer(0)
    for feature in lyr:
        uID = feature.GetField('uID_202011')
        particles = feature.GetField('part_num')
        break    
    part_fact = 0.1 # for testing
    particles = int(particles * part_fact)

    ######## Model
    o = OceanDrift(
        loglevel=0, 
        logfile= os.path.join(logs, 'log_{}.log'.format(uID)),
        seed=None
        )
    o.add_reader([reader_lnd, reader_ssc, reader_nep])

    ######### Seed particles
    time_step = timedelta(hours=4)
    num_steps = 6
    for i in range(num_steps):
        o.seed_from_shapefile(
            shp, 
            origin_marker=uID,
            number=particles, 
            time=reader_ssc.start_time + i*time_step
            )
    # starting coordinates, for use in biology script
    npy_lon = os.path.join(np_out, 'lon_{}.npy'.format(uID))
    npy_lat = os.path.join(np_out, 'lat_{}.npy'.format(uID))    
    np.save(npy_lon, o.elements_scheduled.lon)
    np.save(npy_lat, o.elements_scheduled.lat)

    ######### Configure and Run
    #o.list_configspec()
    o.set_config('general:use_auto_landmask', False)  # use custom landmask
    o.set_config('drift:current_uncertainty', 0) # using basemodel hardcoded values
    o.set_config('general:coastline_action', 'stranding')
    o.set_config('drift:scheme', 'runge-kutta')

    output_nc = os.path.join(nc_out, 'output_{}.nc'.format(uID))
    print("Simulation for uID {}".format(uID))
    o.run(
        #steps=2,
        end_time=reader_ssc.start_time + timedelta(days=5), 
        time_step=60, 
        time_step_output=1800, 
        outfile= output_nc, 
        export_variables=[
            'age_seconds', 
            'land_binary_mask',
            'origin_marker']
        )

    ######### Outputs
    print(o)
    #o.plot(filename='output_3.jpg')
    #o.animation()


#####################

