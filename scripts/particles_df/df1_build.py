# COMBINE all the pieces together


import arcpy
import pandas as pd
import numpy as np
import pyarrow.feather as feather

output_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particles_df\output'
# Distance data prep gdb:
# This gdb was created when caluclating the distance of particles. It also 
# contains a lot of useful general info that I don't want to repeat.
dist_dataprep_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb'
dist_summary_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_summary.gdb'
fetch_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch3_post.gdb'
mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'
ocean_id = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\ocean_groupings_glmm_mpas'

#######################################
# Common pieces, outside for-loop

# get list of mpas to exclude
to_include = []
with arcpy.da.SearchCursor(mpas_exclude, ['uID_20201124', 'exclude']) as cursor:
    for row in cursor:
        if not row[1] == 1:
            to_include.append(row[0])

# get mpas as df
field_names = [i.name for i in arcpy.ListFields(mpas)]
cursor = arcpy.da.SearchCursor(mpas, field_names)
df_mpas = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_mpas = df_mpas[df_mpas.uID_20201124.isin(to_include)] # remove mpas to exlude
keep_fields = ['uID_20201124', 'provstate', 'Shape_Area']
drop_fields = []
for col in df_mpas.columns:
    if col not in keep_fields:
        drop_fields.append(col)
df_mpas = df_mpas.drop(drop_fields, axis=1)
df_mpas.loc[df_mpas.provstate=='British Columbia', 'provstate'] = 'BC'  # rename to BC
df_mpas = df_mpas.rename(columns={'uID_20201124':'mpa_part_id_orig'})

# Get origin XY, read as in-memory
origin_pts = os.path.join(dist_dataprep_gdb, 'origin_coordinates')
arcpy.env.outputCoordinateSystem = mpas
o_pts = arcpy.CopyFeatures_management(origin_pts, 'in_memory/origpts')
arcpy.AddXY_management(o_pts)
field_names = [i.name for i in arcpy.ListFields(o_pts)]
cursor = arcpy.da.SearchCursor(o_pts, field_names)
df_opts = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_opts = df_opts.drop(['OBJECTID', 'Shape'], axis=1)
df_opts = df_opts.rename(columns={'POINT_X':'origin_x', 'POINT_Y':'origin_y'})

# Get exposure(fetch) value
fetch_pts = os.path.join(fetch_gdb, 'origincoords_fetch')
field_names = [i.name for i in arcpy.ListFields(fetch_pts)]
cursor = arcpy.da.SearchCursor(fetch_pts, field_names)
df_fetch = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_fetch = df_fetch.drop(['OBJECTID', 'uID_ref_origin'], axis=1)

# Get ocean grouping (somewhat arbitrary grouping by region to include as a random effect)
field_names = [i.name for i in arcpy.ListFields(ocean_id)]
cursor = arcpy.da.SearchCursor(ocean_id, field_names)
df_ocid = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_ocid = df_ocid.drop(['OBJECTID', 'uID_original', 'partID', 'group_name', 'mpa_area'], axis=1)
df_ocid = df_ocid.rename(columns={'uID_20201124':'mpa_part_id_orig'})


#######################################
# For loop, by period-pld

plds = ['1','3','7','10','21','30','40','60']
periods = ['1101','1105','1108','1401','1405','1408','1701','1705','1708']
prev_length = 0

for period in periods:
    for pld in plds:

        print(f'processing {period}_{pld}')
        fc = os.path.join(dist_dataprep_gdb, f'destpts_{period}_pld{pld}')

        # get XY coords while dataset is still in fc format
        dest_pts = arcpy.CopyFeatures_management(fc, 'in_memory/destpts')
        arcpy.AddXY_management(dest_pts)
        field_names = [i.name for i in arcpy.ListFields(dest_pts)]
        cursor = arcpy.da.SearchCursor(dest_pts, field_names)
        df_dpts = pd.DataFrame(data=[row for row in cursor], columns=field_names)
        df_dpts = df_dpts.drop(['OBJECTID', 'Shape', 'time_int_s', 'mortstep', 'time_int', 'dest_id'], axis=1)
        df_dpts = df_dpts.rename(columns={'POINT_X':'dest_x', 'POINT_Y':'dest_y', 'uID_part':'mpa_part_id_orig'})
        arcpy.Delete_management(dest_pts)

        # Add fields
        df_dpts['pld'] = int(pld)
        df_dpts['year'] = pd.to_numeric(df_dpts.date_start.str[:4])
        df_dpts['month'] = pd.to_numeric(df_dpts.date_start.str[5:7])
        df_dpts['period_id'] = period
        df_dpts = df_dpts.drop(['date_start'], axis=1)

        # join destination points to origin points
        # The order is key here. By joining to df_opts and doing an inner or
        # left join, I am only keeping particles that are from mpas that should
        # not be excluded.
        df_od = df_opts.merge(df_dpts, how='inner', on='traj_id_u')
        df_od = df_od.drop(['traj_id_y'], axis=1)
        df_od = df_od.rename(columns={'traj_id_x':'traj_id'})

        # join mpas to get provstate
        df_od = df_od.merge(df_mpas, how='inner', on='mpa_part_id_orig')
        df_od['mpa_area'] = df_od.Shape_Area.astype(int)
        df_od = df_od.drop(['Shape_Area'], axis=1)

        # join fetch
        df_od = df_od.merge(df_fetch, how='inner', on='traj_id_u')

        # join ocean grouping id
        df_od = df_od.merge(df_ocid, how='inner', on='mpa_part_id_orig')

        # read in particle distances and join
        dist_fc = os.path.join(dist_summary_gdb, f'distance_{period}_pld{pld}')
        field_names = [i.name for i in arcpy.ListFields(dist_fc)]
        cursor = arcpy.da.SearchCursor(dist_fc, field_names)
        df_dist = pd.DataFrame(data=[row for row in cursor], columns=field_names)        
        df_dist = df_dist.drop(['OBJECTID', 'uID_ref_dest', 'uID_ref_origin'], axis=1)
        # there will be fewer particles in this dataset since it is only ones that settled/stranded
        df_od = df_od.merge(df_dist, how='left', on='traj_id_u')
        df_od = df_od.rename(columns={'PathCost':'distance'})
        # add 1/0 field for successful settlement. Only ones with a distance
        # value successfully stranded on the coastline. It's only in this script
        # where I did the intersect with the entire coastline and 30m bathy. The
        # timing columns (i.e. dest_id, time_int) are not helpeful in this case
        # since they only relate to mpa settlement.
        df_od['settle'] = np.where(~df_od.distance.isnull(), 1, None)

        # convert all floats to ints to cut down on size
        df_od['mpa_part_id_orig'] = df_od.mpa_part_id_orig.astype(int)
        df_od['origin_x'] = df_od.origin_x.astype(int)
        df_od['origin_y'] = df_od.origin_y.astype(int)
        df_od['dest_x'] = df_od.dest_x.astype(int)
        df_od['dest_y'] = df_od.dest_y.astype(int)

        # add master unique ID
        df_od = df_od.reset_index(drop=True)
        df_od = df_od.rename(columns={'traj_id_u':'traj_id_series'})
        df_od = df_od.drop(['traj_id'], axis=1)
        df_od['traj_id_overall'] = df_od.index + 1 + prev_length
        prev_length = prev_length + len(df_od)

        # output to arrow file
        out_file = os.path.join(output_dir, f'particles_{period}_{pld}.feather')
        feather.write_feather(df_od, out_file)


################################
