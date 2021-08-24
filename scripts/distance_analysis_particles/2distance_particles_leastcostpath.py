# Calculate overwater distance with barriers between origin and destination for 
# all particles (using reference points created in previous script).

import arcpy
import pandas as pd


prep_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb'
euclines_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_euclines.gdb'
summary_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_summary.gdb'
land_raster = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb\land_r_500'


############################################
# build table of each unique combination of origin ref and dest ref so that I
# don't have to create millions of unecessary lines

df_combos = pd.DataFrame(columns=['uID_ref_origin', 'uID_ref_dest'])

# make df out of near_ORIGIN
arcpy.env.workspace = prep_gdb
field_names = [i.name for i in arcpy.ListFields('near_ORIGIN')]
cursor = arcpy.da.SearchCursor('near_ORIGIN', field_names)
df_orig = pd.DataFrame(data=[row for row in cursor], columns=field_names)

# for each near dest fc, make df
fcs_dest = arcpy.ListFeatureClasses('near_1*')
for fc in fcs_dest:
    field_names = [i.name for i in arcpy.ListFields(fc)]
    cursor = arcpy.da.SearchCursor(fc, field_names)
    df_dest = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    # join df_orig to df_dest, keeping only records from dest
    df_join = df_dest.merge(df_orig, 'inner', 'traj_id_u')
    # append back to df_combos and group by unqiue combination of origin and dest
    # so that the dataframe doesn't grow too large
    df_combos = df_combos.append(df_join)
    df_combos = df_combos.groupby(['uID_ref_origin', 'uID_ref_dest']).size().reset_index().rename(columns={0:'count'})
    df_combos = df_combos.drop(['count'], axis=1)



##############################################
# for each origin point in df_combos, do the least cost path analysis

origin_list = df_combos.uID_ref_origin.unique()
cent_orig = os.path.join(prep_gdb, 'CENTROIDS_origin')
cent_dest = os.path.join(prep_gdb, 'CENTROIDS_dest')
land_r_500 = os.path.join(prep_gdb, 'land_r_500')
arcpy.env.snapRaster = land_r_500
arcpy.env.extent = land_r_500
arcpy.env.workspace = euclines_gdb
for origin_pt in origin_list:

    print(f'processing origin point {origin_pt}')

    # get list of destination reference points for that origin point
    dest_list = list(df_combos[df_combos.uID_ref_origin==origin_pt].uID_ref_dest.unique())

    # from CENTROIDS_origin, select that point and convert to raster
    centroid_sel = arcpy.SelectLayerByAttribute_management(cent_orig, 'NEW_SELECTION', f'uID_ref_origin={origin_pt}')
    origin = 'centroid_origin_{}'.format(str(origin_pt))
    arcpy.FeatureToRaster_conversion(centroid_sel, 'uID_ref_origin', origin, 500)       

    # Euclidean distance with barriers
    in_rast = origin
    max_dist = ''
    cellsize = land_r_500
    outDirectionRaster = '' # you would use this one if you don't have barriers
    distance_method = 'PLANAR'
    in_barrier_data = land_r_500
    outBackDirectionRaster = 'eucbackdirect' # this direction raster is for when you have barriers
    outEucDistance = arcpy.sa.EucDistance(in_rast, max_dist, cellsize, outDirectionRaster, distance_method, in_barrier_data, outBackDirectionRaster)
    outEucDistance.save('eucdistance')

    # Cost Path as Polyline
    # create destination layer, select from CENTROIDS_dest with dest_list
    cent_dest_sel = arcpy.SelectLayerByAttribute_management(cent_dest, 'NEW_SELECTION', f'"uID_ref_dest" IN {str(tuple(dest_list))}')
    dest = f'centroid_dest_{origin_pt}'
    arcpy.FeatureToRaster_conversion(cent_dest_sel, 'uID_ref_dest', dest, 500)
    inputDestinationLayer = dest
    inputCostLayer = 'eucdistance'
    inputDirectionLayer = 'eucbackdirect'
    outLines = f'euc_lines_{origin_pt}'
    pathType = 'EACH_CELL'
    destfield = ''
    arcpy.sa.CostPathAsPolyline(inputDestinationLayer, inputCostLayer, inputDirectionLayer, outLines, pathType, destfield)

    # add origin_uid attribute
    arcpy.AddField_management(outLines, 'origin_id', 'SHORT')
    with arcpy.da.UpdateCursor(outLines, ['origin_id']) as cursor:
        for row in cursor:
            row[0] = int(origin_pt)
            cursor.updateRow(row)

    arcpy.Delete_management(origin)
    arcpy.Delete_management('eucbackdirect')
    arcpy.Delete_management('eucdistance')
    arcpy.Delete_management(dest)

