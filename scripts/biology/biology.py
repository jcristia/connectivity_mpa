# Script to add bio part of a biophysical model
# This script takes the output netcdf file from an Opendrift simulation and 
# modifies particle end trajectories by consdering precompetency period, 
# settlement, and mortality.

# John Cristiani
# University of British Columbia
# origin: 2019-04-19
# chapter 2 updates: 2021-01-15

# env: biology


import netCDF4 as nc
import numpy as np
from shapely.geometry import shape, Point, LineString, Polygon
import pandas as pd
import geopandas
import logging
import os
logging.basicConfig(level=logging.INFO)



###################
# paths and variables
# folder structure and filenames are important
# If there are any significant changes then check how filenames are
# concatenated at the bottom
###################

input_nc_dir = r'outputs/nc' # where nc files were output
input_npy_dir = r'outputs/npy' # where npy files were output
#### for cluster:
#input_shp_dir = r'shp_inputs/release_polys' # shpfiles used for opendrift release
#shp_og = r'shp_inputs/shp_og/mpas_.shp' # full release_polys dataset
#shp_og_buff = r'shp_inputs/shp_og_buff/mpas_.shp'  # buffered for checking settlement. 
# This is to account for possible slivers between coastline
#### for local:
input_shp_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpas_shp' # shpfiles used for opendrift release
shp_og = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpa_.shp' # full release_polys dataset
shp_og_buff = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpa_.shp'  # buffered for checking settlement
####
crs_input_shp = {'init' :'epsg:3005'} # BC Albers
output_shp_dir = r'outputs/shp' # connection lines

precomp = 0
# if I am using 'stranding' in opendrift, then I likely need at least a small 
# precompetency period. Otherwise everything just ends up settling at home. 
# Be careful setting this. It depends on the time_step and time_step_output 
# you used in the run. It is in units of the timestep output. If time step 
# output is 30 minutes, then a precomp of 2 is 1 hour.

# get these values from the simulation script
time_step_output = 0.5 # in hours. It will be in seconds in the opendrift script
interval_of_release = 4 # in hours, interval can't be less than time step output
# if no delayed release then just put same value as time_step_output
num_of_releases = 3 # if no delayed release then just put 1

# allow particles to settle?
settlement_apply = True

# mortality
mortality_rate = 0.15 # instantaneous daily rate
mort_period = 8 # after how many time_step_outputs to apply mortality rate 
# MAKE THIS A FACTOR OF 24
# The mortality rate will be scaled appropriately so that it still matches the 
# daily rate. This option is given because it seems uncessary to apply it at 
# every time step, but for some species with short PLDs, it will make sense to 
# apply it more often than once per day. If mortality rate is 0, the also set 
# this to 0.

# if it is an Opendrift backwards run
# !!! I noticed some inconsistent behavior with the shape of resulting arrays
# for each particle in the time dimension. If they are not all the same length
# this would create an issue for how I analyze backwards runs. I will need to
# look into this more if I need to do backwards runs.
# So, for now this feature is considered in development. 
backwards_run = False

# PLD limit
# I run opendrift simulations for what I expect the max PLD to be that I am 
# considering. Then in this script I can set a smaller PLD and see which 
# connections are made if I had only run the simulation up to a certain 
# timestep. I need to do PLDs all at once on each run of this script because 
# mortality is random and I want all PLDs done on one random selection of 
# particles instead of on different selections.
# Provide PLDs in a list in units of days
plds = [1, 3, 7, 21, 60]



###################
# BIOLOGY FUNCTIONS
###################

###################
# assign a polygon unique ID to each particle
###################

def get_particle_originPoly(
    origin_marker, traj, crs_input_shp, lat_np, lon_np, backwards_run):

    logging.info("Getting origin coordinates of each particle")

    # get starting lat/lon for each particle
    lons = np.load(lon_np)
    lats = np.load(lat_np)
    # reverse order for backwards runs (opendrift records them in reverse in 
    # the netcdf file)
    if backwards_run:
        lons = lons[::-1]
        lats = lats[::-1]
    
    # get origin_marker (uID)
    om = origin_marker[0][0]
    if np.ma.is_masked(om):
        om = origin_marker[0][origin_marker[0].mask == False][0]

    # add origin coords , trajid and uID to df
    df = pd.DataFrame()
    df['o_coords'] = list(zip(lons, lats))
    df['o_coords'] = df['o_coords'].apply(Point)
    df['traj_id'] = list(traj)
    df['origin_marker'] = om
    points = geopandas.GeoDataFrame(df, geometry='o_coords')
    points.crs = {'init' :'epsg:4326'}
    points = points.to_crs(crs_input_shp)
    origin = pd.DataFrame(data=points)
    
    return origin


###################
# determine precompetency and release intervals
###################

def calc_precomp(
    precomp, time_step_output, particles_per_release, interval_of_release, 
    num_of_releases, traj):

    # the timesteps when we release particles
    timesteps_with_release = []
    for release in range(num_of_releases):
        ts = (float(interval_of_release) / float(time_step_output)) * release
        timesteps_with_release.append(int(ts))

    # when the precompetency period ends for each group of particle releases
    precomp_end_timestep = []
    for release in timesteps_with_release:
        ts_e = release + precomp
        precomp_end_timestep.append(ts_e)

    # the range of time periods of the precomp period for each group of particles
    precomp_range = []
    for p in precomp_end_timestep:
        precomp_range.append([p-precomp, p])

    # the corresponding particle IDs for each release
    particle_range = []
    if num_of_releases == 1:
        particle_range = [[1, len(traj) + 1]]
    else:
        for release in range(1,num_of_releases+1):
            # Opendrift keeps particles in order that they are released. Hopefully this never changes.
            p_range = [1 + ((release-1) * particles_per_release),(release * particles_per_release) +1]
            particle_range.append(p_range)

    return timesteps_with_release, precomp_end_timestep, precomp_range, particle_range


###################
# settle a particle when it drifts over a release_polys patch
# account for precompetency period
###################

def settlement(settlement_apply, origin, shp_og_buff, timestep, status, lon, lat, traj, crs_input_shp, precomp, precomp_range, particle_range, mortality_rate):

    poly  = geopandas.GeoDataFrame.from_file(shp_og_buff)
    poly.crs = crs_input_shp
    dest_df = pd.DataFrame(columns=['d_coords','traj_id','dest_id','time_int'])
    #pd_i = 0

    if settlement_apply: # if this is false, then it will just join the blank dest_df to origin, and the get_destination_coords function will fill in the rest
        for i in range(1,len(timestep)):

            output_str = "settlement time step " + str(i) + " of " + str(len(timestep)-1)
            logging.info(output_str)
            # get traj ids for particles that are active or where they were active on the previous step (just stranded)
            # NOTE: special case I may need to fix in the future: when running backwards I had a particle that was 1 on the very first time step. However, since it always seems to mask them after they are 1, I could just select where == 1 and not worry about if the previous step was 0.
            t_strand = traj[np.where((status[:,[i]] == 1) & (status[:,[i-1]] == 0))[0]]
            t_active = traj[np.where(status[:,[i]] == 0)[0]]
        
            # if we already settled it on a previous interation of the for loop then remove it from the list so we don't check it again
            #for p in dest_df["traj_id"]:
            #    if p in t_strand:
            #        index = np.argwhere(t_strand==p)
            #        t_strand = np.delete(t_strand, index)
            #    if p in t_active:
            #        index = np.argwhere(t_active==p)
            #        t_active = np.delete(t_active, index)
            t_strand = np.setdiff1d(t_strand, dest_df.traj_id.values)
            t_active = np.setdiff1d(t_active, dest_df.traj_id.values)

            if precomp > 0: # remove from p_active ones that are in their precomp period
                for period in precomp_range:
                    if i in range(period[0],period[1]):
                        period_index = precomp_range.index(period)
                        # get particles that are still in their precomp period
                        p_in_precomp = range(particle_range[period_index][0],particle_range[period_index][1])
                        #for y in p_in_precomp:
                        #    if y in t_active:
                        #        index = np.argwhere(t_active==y)
                        #        t_active = np.delete(t_active, index)
                        t_active = np.setdiff1d(t_active, p_in_precomp)

            t = np.concatenate((t_strand,t_active))

            if len(t) == 0:
                continue

            #lons = []
            #lats = []
            #for par in t:
            #    index = par - 1  # t is the actual trajectory ID, the index for that value is 1 less
            #    lons.append(lon[index][i])
            #    lats.append(lat[index][i])
            lons = lon[t-1,i]
            lats = lat[t-1,i]

            df = pd.DataFrame()
            df['d_coords'] = list(zip(lons, lats))
            df['d_coords'] = df['d_coords'].apply(Point)
            df['traj_id'] = list(t)
            points = geopandas.GeoDataFrame(df, geometry='d_coords')
            points.crs = {'init' :'epsg:4326'}
            points = points.to_crs(crs_input_shp)
            pointInPolys = geopandas.tools.sjoin(points, poly, how='inner')
            #for row in pointInPolys.itertuples(index=False):
            #    dest_df.loc[pd_i] = [row[0],row[1],row[6],i]
            #    pd_i += 1
            pointInPolys = pointInPolys.rename(columns={'uID':'dest_id'})
            pointInPolys['time_int'] = i
            dest_df = dest_df.append(pointInPolys[['d_coords','traj_id','dest_id','time_int']], ignore_index=True)

            # HERE, I would need to
            # If it is not the last step
            # join to origin, create a df of ones that settled at origin patch
            # remove from df ones that are in t_strand - these are ok to settle at home
            # Then, in dest_df, remove the ones that are remaining in this df

    
    # join the two tables
    # The resulting data frame is the particles that settled in another patch
    # to get all particles including the ones that did not settle change to:  how='outer'
    logging.info("merging destination and origin dataframes")
    # need to coerce merge. traj_id must be numeric. The dest_df data types were all "object"
    # this was not a problem on windows, but when running on the cluster it woud give an error
    dest_df = dest_df.infer_objects()
    origin = origin.infer_objects()
    dest_df.traj_id = dest_df.traj_id.astype('float')
    origin.traj_id = origin.traj_id.astype('float')
    origin_dest = dest_df.merge(origin, on='traj_id', how='outer')

    return origin_dest

###################
# add the final destination coordinates to particles that did not settle on a patch
###################

def get_destination_coords(origin_dest, traj, lon, lat, timestep, crs_input_shp, status):

    logging.info("getting destination coordinates")
    lons_dest = []
    lats_dest = []
    time_steps = []
    for i in range(len(traj)):
        if np.ma.is_masked(lon[i][-1]): # if the last value is masked (but really just to check if any values are masked, similar to getParticleOriginPoly). If it is masked then it must have stranded, and therefore we can search by where it is 1.
            # changed this statement from '== 1' to '> 0'. I found that if a particle goes outside of the grid it gets coded as '2 - missing data'. Opendrift says that anything above 0 is considered deactivated, so the actual number doesn't matter (at least for my purposes, just that it is bigger than 0.
            j = np.where(status[i] > 0)[0][0]
            lo = lon[i][j]
            lons_dest.append(lo)
            la = lat[i][j]
            lats_dest.append(la)

            # old way
            #lons_dest.append(lon[i][lon[i].mask==False][-1]) # get the last coordinate where it is not masked
            #lats_dest.append(lat[i][lat[i].mask==False][-1])

            # timestep
            index = j
            time_steps.append(index)
            # old way
            # this is an optimized way to get the last time step
            # it flips the masking, then gets locations where masking is True, then I just take the last one
            #d = status[i].mask == False
            #index = np.where(d)
            #time_steps.append(index[0][-1])
        else: # otherwise just get the last coordinate
            lons_dest.append(lon[i][-1])
            lats_dest.append(lat[i][-1])
            time_steps.append(len(lon[i]))

    df = pd.DataFrame()
    df['Coordinates'] = list(zip(lons_dest, lats_dest))
    df['Coordinates'] = df['Coordinates'].apply(Point)
    df['traj_id'] = list(traj)
    df['time_step'] = time_steps
    points_dest = geopandas.GeoDataFrame(df, geometry='Coordinates')
    points_dest.crs = {'init' :'epsg:4326'}
    points_dest = points_dest.to_crs(crs_input_shp)

    # this takes a long time. Should consider how to optimize it.
    #df1 = origin_dest[origin_dest.dest_id.isnull()]
    #for row in df1.itertuples(index=True):
    #    traj_id = row[2]
    #    dcoords = points_dest.loc[points_dest['traj_id'] == traj_id, 'Coordinates'].iloc[0]
    #    time = points_dest.loc[points_dest['traj_id'] == traj_id, 'time_step'].iloc[0]
    #    origin_dest['d_coords'][row[0]] = dcoords
    #    origin_dest['time_int'][row[0]] = time

    # I was originally doing itertuples, which I found out is extremely slow and a big no-no.
    # stick with vectorization when possible
    # join, fill in values where null, remove columns
    logging.info("joining destination coordinates to dataframe")
    points_dest = points_dest.infer_objects()
    points_dest.traj_id = points_dest.traj_id.astype('float')
    origin_dest = origin_dest.merge(points_dest, on='traj_id')
    origin_dest['time_int'].loc[origin_dest['time_int'].isnull()] = origin_dest['time_step']
    origin_dest['d_coords'].loc[origin_dest['d_coords'].isnull()] = origin_dest['Coordinates']
    origin_dest = origin_dest.drop(['time_step', 'Coordinates'], axis=1)

    origin_dest = origin_dest.sort_values(by=['traj_id'])
    origin_dest = origin_dest.reset_index(drop=True)

    return origin_dest


###################
# calculate mortality
###################

def calc_mortality(mortality_rate, traj, timestep, origin_dest, time_step_output, mort_period, interval_of_release, num_of_releases):

    logging.info("calculating mortality")

    mortality_p = pd.DataFrame(columns=['traj_id','mortstep'])

    if mortality_rate > 0:
        timestep_days = time_step_output / 24   # proportion of a day for one timestep
        mort_timesteps = np.arange(0, len(timestep)-1, mort_period)   # timesteps to apply mortality
        # run all this by dad and Patrick
        # instantaneous mortality rate for the mortality application interval
        inst_rate = 1 - (math.exp(math.log(1-mortality_rate) * (timestep_days * mort_period)))
        # so as long as the interval that I calculate mortality stays the same throughout a simulation, I don't need to worry about a changing rate or how many new particles enter the system

        # need to not consider particles that are not released yet, so find the periods and particles ranges
        # timesteps when particles are released
        timesteps_with_release = []
        for release in range(num_of_releases):
            ts = (float(interval_of_release) / float(time_step_output)) * release
            timesteps_with_release.append(int(ts))
        timesteps_with_release = np.array(timesteps_with_release)
        # particles in that period
        particle_range = []
        if num_of_releases == 1:
            particle_range = [[1, len(traj) + 1]]
        else:
            for release in range(1,num_of_releases+1):
                # Opendrift keeps particles in order that they are released. Hopefully this never changes.
                p_range = [1 + ((release-1) * particles_per_release),(release * particles_per_release) +1]
                particle_range.append(p_range)
        particle_range = np.array(particle_range)

        mortality_p = pd.DataFrame(columns=['traj_id','mortstep'])

        for i in mort_timesteps[1:]:
            mortality_selection = traj[:]
            # remove ones that have not seeded yet
            if i <= np.max(timesteps_with_release):
                # get indices of all timesteps >= current time step
                periods_exempt = particle_range[np.where(timesteps_with_release[:] >= i)]
                # get particles that should not be considered for mortality
                p_remove = np.arange(np.min(periods_exempt), np.max(periods_exempt))
                # remove from mortality_selection
                mortality_selection = np.setdiff1d(mortality_selection, p_remove)

            # remove ones that have already been added to mortality_p
            mortality_selection = np.setdiff1d(mortality_selection, mortality_p['traj_id'].values)

            # remove ones that have settled before this time step
            # select from origin_dest where time_int is this timestep and where dest_id is not null (so actually settled somewhere, I can still killed stranded ones)
            # .values turns the selection into a numpy array
            p_settled = origin_dest['traj_id'].loc[(origin_dest["time_int"] <= i) & (origin_dest['dest_id'].notnull())].values
            mortality_selection = np.setdiff1d(mortality_selection, p_settled)

            # select random particles based on inst_rate
            num_to_kill = int(len(mortality_selection) * inst_rate)
            mortality_selection = np.random.choice(mortality_selection, num_to_kill, replace=False)

            # append this selection to mortality_p with the timestep that they were killed
            df = pd.DataFrame({'traj_id':mortality_selection, 'mortstep':i})
            mortality_p = mortality_p.append(df, ignore_index=True, sort=True)


    # need to coerce merge. traj_id must be numeric. The dest_df data types were all "object"
    # this was not a problem on windows, but when running on the cluster it woud give an error
    mortality_p = mortality_p.infer_objects()
    origin_dest = origin_dest.infer_objects()
    mortality_p.traj_id = mortality_p.traj_id.astype('float')
    origin_dest.traj_id = origin_dest.traj_id.astype('float')
    # join to origin_dest
    origin_dest_mort = origin_dest.merge(mortality_p, on='traj_id', how='outer')
    # we still want the join to happen even if mortality_p is empty. It will just make the column NaN which we later turn to -1.
    # Since I will be doing a lot of runs and writing additional scripts to analyze the data, I think it is important that all data has the same columns, even if I didn't apply a mortality rate to it.

    return origin_dest_mort, mortality_p

###################
# add in starting time interval
# added 2020-01-13 so that we know the full time period of particles that may not have been released at the first time step
# This allows us calculate different connections for different PLDs
###################

def start_time_int(origin_dest_mort, timesteps_with_release, particle_range, traj):
    
    logging.info("adding in particle start time")

    df = pd.DataFrame()
    df['traj_id'] = list(traj)

    # starting time step of each particle
    # this could probably be optimized
    time_int_start = []
    for t in range(len(timesteps_with_release)):
        for particle in range(particle_range[t][0],particle_range[t][1]):
            time_int_start.append(timesteps_with_release[t])
    
    df['time_int_s'] = time_int_start # name shortened for shapefile

    # need to coerce merge. traj_id must be numeric. The dest_df data types were all "object"
    # this was not a problem on windows, but when running on the cluster it woud give an error
    df = df.infer_objects()
    origin_dest_mort = origin_dest_mort.infer_objects()
    df.traj_id = df.traj_id.astype('float')
    origin_dest_mort.traj_id = origin_dest_mort.traj_id.astype('float')   
    origin_dest_mort = origin_dest_mort.merge(df, on='traj_id')

    return origin_dest_mort








###################
# OUTPUTS
###################

#### output destinaiton points to shapefile ####
def out_shp_dest_points(origin_dest_mort, crs_input_shp, shp_out, date_start):
    logging.info("writing points to shapefile")
    # can only have one geometry column
    # remove origin spatial column since for origin I am just concernced about origin poly ID
    od = origin_dest_mort.drop(['o_coords'], axis=1)
    od = geopandas.GeoDataFrame(od, geometry='d_coords')
    od.crs = crs_input_shp
    od = od.fillna(-1) # fill NaN with -1, otherwise NaN gets turned to 0 on export. This could be confusing when analyzing the data
    od['date_start'] = date_start
    od.to_file(filename=shp_out, driver='ESRI Shapefile')

#### create connection lines ####
def connection_lines(shp_out, shp_og, crs_input_shp, conn_lines_out, date_start, pld_int, pld):

    logging.info("writing connection lines to shapefile")
    od = geopandas.read_file(shp_out)
    sg = geopandas.read_file(shp_og)

    ### on od, select particles where time_int_s minus time_int is less than or equal to PLD
    od_pld = od[(od.time_int - od.time_int_s <= pld_int)]

    # get each unique combination of originID and destID and get count of particles that survived
    od_unique = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID','dest_id']).size().reset_index(name='Freq')
    # how to read this:
    # first we select from od the ones that settled and survived
    # then we groupby unique combinations of uID and dest_id
    # then we get the count of those unique combinations
    # this normally makes uID the index and doesn't have a column name for count (the series we created), so we reset index and give the count a column name

    # df of time interval where first settlement occurred
    df_time_int = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID','dest_id'])['time_int'].min().reset_index(name='time_int')
    
    # set up for creating self connection lines. Size of circle lines based on amount settled and average area of all patches.
    def CircleCoords(xLeft, yCenter, r, n): # credit: MGET. Also, see circle_coords.xlsx for explanation of equation.
        return [(xLeft + r - math.cos(2*math.pi/n*x)*r, math.sin(2*math.pi/n*x)*r + yCenter) for x in range(n+1)]
    # min and max quantities used for normalization
    quantity_min = od_unique[od_unique.dest_id == od_unique.uID].Freq.min()
    quantity_max = od_unique[od_unique.dest_id == od_unique.uID].Freq.max()
    # get average area
    area_mean = sg.area.mean()
    # get radius of a circle with this area
    radius = math.sqrt(area_mean/math.pi)
    
    # for each unique combinaton create line from centroids
    connection_lines = pd.DataFrame(columns=['from_id','to_id','quantity','totalori','prob','time_int', 'line', 'pld'])
    conn_i = 0
    for row in od_unique.itertuples(index=False):
        # get total amount of particles released from patch
        total = od.uID[od.uID ==  row[0]].value_counts().values[0]

        # time interval where first settlement occurred
        time_int = df_time_int[(df_time_int.uID == row[0]) & (df_time_int.dest_id == row[1])]['time_int'].values[0]
    
        # get centroid of from and to patches
        centroid_origin = sg[sg.uID == row[0]].centroid
        centroid_dest = sg[sg.uID == row[1]].centroid

        if row[0] != row[1]:
            geom_line = LineString([centroid_origin.tolist()[0], centroid_dest.tolist()[0]])
        else:
            # normalize the quantites to 0.5 - 1 range (or I can do 0-1 but then the smallest one won't show up)
            #quantity_norm = 0.5 * (row[2] - quantity_min) / float(quantity_max - quantity_min) + 0.5
            quantity_norm = (row[2] - quantity_min) / float(quantity_max - quantity_min)
            radius_adj = radius * quantity_norm
            geom_line = LineString(CircleCoords(centroid_origin.x.tolist()[0], centroid_origin.y.tolist()[0], radius_adj, 90))
    
        connection_lines.loc[conn_i] = [row[0],row[1],float(row[2]),float(total),row[2]/float(total), time_int,geom_line, pld]
        conn_i += 1
    
    connection_lines['date_start'] = date_start   
    connection_lines = geopandas.GeoDataFrame(connection_lines, geometry='line')
    connection_lines.crs = crs_input_shp
    connection_lines.to_file(filename=conn_lines_out, driver='ESRI Shapefile')

#### output patch centroids to shapefile (for use in network analysis) ####
def out_shp_patch_centroids(shp_og, patch_centroids_out, crs_input_shp, date_start):
    sg = geopandas.read_file(shp_og)
    # copy poly to new GeoDataFrame
    points = sg.copy()
    # change the geometry
    points.geometry = points['geometry'].centroid
    # same crs
    points.crs = crs_input_shp
    points['date_start'] = date_start
    points.to_file(filename=patch_centroids_out, driver='ESRI Shapefile')





###################
# RUN biology
# run all functions, even if you aren't applying settlement and/or mortality
###################

# create output shp folder
if not os.path.exists(output_shp_dir):
    os.makedirs(output_shp_dir)

# get release shapefiles
sg_files = os.listdir(input_shp_dir)
shapefiles = []
for file in sg_files:
    if file.endswith('.shp'):
        shapefiles.append(os.path.join(input_shp_dir, file))

# reduce that range if necessary
shapefiles = shapefiles[:3]

# run biology functions for each release shapefile
for shp in shapefiles:

    # get base name
    base = os.path.splitext(os.path.basename(shp))[0]

    # output from Opendrift
    nc_output = os.path.join(input_nc_dir, 'output_' + base + '.nc')

    # the lat/lon numpy files of starting coords saved from the opendrift run
    lat_np = os.path.join(input_npy_dir, 'lat_' + base + '.npy')
    lon_np = os.path.join(input_npy_dir, 'lon_' + base + '.npy')

    release_polys = shp

    dataset = nc.Dataset(nc_output, "r+")
    lon = dataset.variables["lon"]
    lat = dataset.variables["lat"]
    traj = dataset.variables["trajectory"]
    status = dataset.variables["status"]
    print(status)
    timestep = dataset.variables["time"]
    origin_marker = dataset.variables['origin_marker']
    age_seconds = dataset.variables['age_seconds']
    # this will give the actual 30 SECOND interval that it was deactivated
    date_start = dataset.time_coverage_start

    particles_per_release = int(len(traj) / num_of_releases)

    origin = get_particle_originPoly(
        origin_marker, traj, crs_input_shp, lat_np, lon_np, backwards_run
        )

    if precomp == 0:
        timesteps_with_release = None
        precomp_end_timestep = None
        precomp_range = None
        particle_range = None
    else:
        timesteps_with_release, precomp_end_timestep, \
        precomp_range, particle_range = calc_precomp(
            precomp, time_step_output, particles_per_release, 
            interval_of_release, num_of_releases, traj
            )

    origin_dest = settlement(
        settlement_apply, origin, shp_og_buff, timestep, status, lon, lat, traj, 
        crs_input_shp, precomp, precomp_range, particle_range, mortality_rate
        )

    origin_dest = get_destination_coords(
        origin_dest, traj, lon, lat, timestep, crs_input_shp, status
        )

    origin_dest_mort, mortality_p = calc_mortality(
        mortality_rate, traj, timestep, origin_dest, time_step_output, 
        mort_period, interval_of_release, num_of_releases
        )

    origin_dest_mort = start_time_int(
        origin_dest_mort, timesteps_with_release, particle_range, traj
        )

    ### outputs

    shp_out = os.path.join(output_shp_dir, 'dest_biology_pts_' + base + '.shp')
    out_shp_dest_points(origin_dest_mort, crs_input_shp, shp_out, date_start)

    if settlement_apply:
        for pld in plds:
            # check that pld is not longer than length of timestep
            pld_int = (pld * 24) / time_step_output
            if pld_int > len(timestep):
                logging.error("PLD provided is greater than length of timestep")
                break
            conn_lines_out = os.path.join(output_shp_dir, 'connectivity_' + base + '_pld' + str(pld) + '.shp')
            connection_lines(
                shp_out, shp_og, crs_input_shp, conn_lines_out, date_start, pld_int, pld
                )

patch_centroids_out = os.path.join(output_shp_dir, 'patch_centroids.shp')
out_shp_patch_centroids(
    shp_og, patch_centroids_out, crs_input_shp, date_start
    )