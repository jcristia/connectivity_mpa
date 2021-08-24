# associate distances between reference points to all particles


import arcpy
import pandas as pd
import numpy as np


prep_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb'
euclines_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_euclines.gdb'
summary_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_summary.gdb'



# read all euc lines into 1 dataframe
arcpy.env.workspace = euclines_gdb
df_lines = pd.DataFrame(columns=['origin_id', 'DestID', 'PathCost'])
fcs_euc = arcpy.ListFeatureClasses('euc*')
for fc in fcs_euc:
    field_names = ['origin_id', 'DestID', 'PathCost']
    cursor = arcpy.da.SearchCursor(fc, field_names)
    df_distances = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    df_lines = df_lines.append(df_distances, ignore_index=True)


# read near_ORIGIN in as dataframe. This links every particle's starting 
# position to a reference point.
arcpy.env.workspace = prep_gdb
field_names = ['traj_id_u', 'uID_ref_origin']
cursor = arcpy.da.SearchCursor('near_ORIGIN', field_names)
df_orig_ref = pd.DataFrame(data=[row for row in cursor], columns=field_names)


# for each particles point feature class, get distances
arcpy.env.workspace = prep_gdb
fcs_particles = arcpy.ListFeatureClasses('near_1*')
for fc in fcs_particles:
    
    # read in fc as df
    field_names = ['traj_id_u', 'uID_ref_dest']
    cursor = arcpy.da.SearchCursor(fc, field_names)
    df_particles = pd.DataFrame(data=[row for row in cursor], columns=field_names)

    # get origin reference point
    df_particles = df_particles.merge(df_orig_ref, 'inner', on='traj_id_u')

    # get least cost path
    df_particles = df_particles.merge(
        df_lines, 
        'inner', 
        left_on=['uID_ref_origin', 'uID_ref_dest'],
        right_on=['origin_id', 'DestID']
        )

    # save to table
    df_particles = df_particles.drop(['origin_id', 'DestID'], axis=1)
    date_pld = fc.split('_', 1)[1]
    x = np.array(np.rec.fromrecords(df_particles.values))
    names = df_particles.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    arcpy.da.NumPyArrayToTable(x, os.path.join(summary_gdb, f'distance_{date_pld}'))

