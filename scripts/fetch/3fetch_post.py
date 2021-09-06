# relate the fetch results of the reference points back to the particles


import arcpy
import pandas as pd
import numpy as np

out_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch3_post.gdb'
fetch_pts = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch2_measure.gdb\FetchPoints_0905_1305'
origin_pts = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch1_dataprep.gdb\origincoords_near'

arcpy.env.workspace = out_gdb


# read in fetch pts
field_names = [i.name for i in arcpy.ListFields(fetch_pts)]
cursor = arcpy.da.SearchCursor(fetch_pts, field_names)
df_fpts = pd.DataFrame(data=[row for row in cursor], columns=field_names)
# calculate total fetch length per point
fetch_fields = []
for col in df_fpts.columns:
    if col.startswith('bearing'):
        fetch_fields.append(col)
df_fpts['total_exposure'] = df_fpts[fetch_fields].sum(axis=1)


# read in origin coordinates
field_names = [i.name for i in arcpy.ListFields(origin_pts)]
cursor = arcpy.da.SearchCursor(origin_pts, field_names)
df_opts = pd.DataFrame(data=[row for row in cursor], columns=field_names)


# join fetch pts to origin coordinates and clean up
df_merge = df_opts.merge(df_fpts, how='left', on='uID_ref_origin')
drop_fields = []
for col in df_merge.columns:
    if col not in ['uID_ref_origin', 'traj_id_u', 'total_exposure']:
        drop_fields.append(col)
df_merge = df_merge.drop(drop_fields, axis=1)


# write to arc table
x = np.array(np.rec.fromrecords(df_merge.values))
names = df_merge.dtypes.index.tolist()
x.dtype.names = tuple(names)
arcpy.da.NumPyArrayToTable(x, os.path.join(arcpy.env.workspace, f'origincoords_fetch'))