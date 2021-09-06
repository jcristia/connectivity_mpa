# create reference points within MPAs to calculate fetch
# relate each particles origin point to one of these reference points

import arcpy
import pandas as pd

mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'
origin_pts = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb\origin_coordinates'
land = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_FINAL'
states_general = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_boundaries.gdb\states_general'
landmask = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_mb_nep'
out_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch1_dataprep.gdb'

arcpy.env.workspace = out_gdb


# read in mpa uIDs to exclude
to_include = []
with arcpy.da.SearchCursor(mpas_exclude, ['uID_20201124', 'exclude']) as cursor:
    for row in cursor:
        if not row[1] == 1:
            to_include.append(row[0])
# select from MPAs and copy to gdb
tmp_sel = arcpy.SelectLayerByAttribute_management(mpas, 'NEW_SELECTION', '"uID_20201124" IN {}'.format(str(tuple(to_include))))
arcpy.CopyFeatures_management(tmp_sel, 'mpas_1exclude')
# erase land from MPAs
arcpy.Erase_analysis('mpas_1exclude', land, 'mpas_2erase')


# create fishnet of points over entire area
arcpy.CreateFishnet_management(
    'refcoords_1fishnet',
    '288809.7027 -34801.6535',
    '288809.7027 1630803.6947',
    500, 500,
    None, None,
    None,
    'LABELS',
    'mpas_2erase',
    'POLYGON')
# clip points to MPAs
arcpy.Clip_analysis('refcoords_1fishnet_label', 'mpas_2erase', 'refcoords_2clip')
# clean up
arcpy.Delete_management('refcoords_1fishnet')
arcpy.Delete_management('refcoords_1fishnet_label')
arcpy.Delete_management('mpas_1exclude')


# there are some very small MPAs where we don't get any points inside.
# For these, make points along the boundary and merge the two datasets together.
arcpy.GeneratePointsAlongLines_management(
    'mpas_2erase',
    'refcoords_3boundarypts',
    'DISTANCE',
    Distance=500,
    Include_End_Points='END_POINTS'
)
# Merge
arcpy.Merge_management(['refcoords_2clip', 'refcoords_3boundarypts'], 'refcoords_4merge')

# Remove all extra fields
fields = arcpy.ListFields('refcoords_4merge')
for field in fields:
    if field.name not in ['OBJECTID', 'Shape']:
        arcpy.DeleteField_management('refcoords_4merge', field.name)
# add field for a reference pt uID
arcpy.AddField_management('refcoords_4merge', 'uID_ref_origin', 'LONG')
arcpy.CalculateField_management('refcoords_4merge', 'uID_ref_origin', '!OBJECTID!')


#### do near analysis of origin points to fishnet reference points

# first, remove those in excluded MPAs
trajidu_exclude = []
dest_pts = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb\destpts_1101_pld1'
    # dest pts datasets are the only ones that link the traj_id_u to the mpa part id
with arcpy.da.SearchCursor(dest_pts, ['uID_part', 'traj_id_u']) as cursor:
    for row in cursor:
        if row[0] not in to_include:
            trajidu_exclude.append(row[1])
print('removing particles that are in excluded MPAs')
# THIS TAKES FOREVER, but it works.
arcpy.MakeFeatureLayer_management(origin_pts, 'temp_lyr')
with arcpy.da.UpdateCursor('temp_lyr', ['traj_id_u']) as cursor:
    for row in cursor:
        if row[0] in trajidu_exclude:
            cursor.deleteRow()
print('copying output')
arcpy.env.outputCoordinateSystem = 'refcoords_4merge'
arcpy.CopyFeatures_management('temp_lyr', 'origincoords_near')
arcpy.Delete_management('temp_lyr')
# Near analysis
arcpy.Near_analysis('origincoords_near', 'refcoords_4merge', 5000, 'NO_LOCATION', 'NO_ANGLE', 'PLANAR')
arcpy.AddField_management('origincoords_near', 'uID_ref_origin', 'LONG')
arcpy.CalculateField_management('origincoords_near', 'uID_ref_origin', '!NEAR_FID!')


# remove any ref points that aren't referred to in origincoords_near
arcpy.CopyFeatures_management('refcoords_4merge', 'refcoords_5remove')
# get unique ref IDs from origincoords_near in list
field_names = [i.name for i in arcpy.ListFields('origincoords_near')]
cursor = arcpy.da.SearchCursor('origincoords_near', field_names)
df_near = pd.DataFrame(data=[row for row in cursor], columns=field_names)
uid_ref_origin_unique = list(df_near.uID_ref_origin.unique())
with arcpy.da.UpdateCursor('refcoords_5remove', ['uID_ref_origin']) as cursor:
    for row in cursor:
        if row[0] not in uid_ref_origin_unique:
            cursor.deleteRow()


# so I now have two datasets for the following scripts:
# refcoords_5remove: these will be the reference points that I will calculate
# fetch for.
# origincoords_near: these are the original starting points for every particle,
# there is an ID that links each particle to a ref point in refcoords_5remove.


# create land dataset for use in fetch script
# In this case, I need a land dataset that includes Alaska
states = ['Alaska', 'California', 'Oregon']
sel_states = arcpy.SelectLayerByAttribute_management(states_general, 'NEW_SELECTION', '"STATE_NAME" IN {}'.format(tuple(states)))
arcpy.Erase_analysis(sel_states, landmask, 'land_1erase')
arcpy.Merge_management([land, 'land_1erase'], 'land_2merge')
arcpy.Delete_management('land_1erase')