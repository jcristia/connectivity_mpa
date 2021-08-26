# prep data for overwater distance analysis between origin and destination of
# all particles that settle

# I should have put all this in functions. Oh well.

# The challenge is...
# Some MPAs are very large and it doesn't make sense to just represent it as one
# point in a distance analysis, but I also can't do a distance analysis between 
# 200 million particles. However, I don't need this amount of precision anyways.
# Therefore, I am creating points to represent general areas, and then getting
# distances between them.
# The major challenge with this though is that a lot of points strand a bit
# inland, so when creating points, some of them get made on land. This script
# handles that. All points are within 1km of a reference point, so when
# presenting results, I need to remember that distance precision is 1km.

# Outputs:
# Two point datasets that represent reference origin and destination points for
# particles. I will use these datasets in the distance analysis.
# There is also a 'near' dataset for the origin and destination of ALL particles.

# Also to note:
# I am only calculating distance for particles that strand/settle.
# My analysis is looking at distance by PLD, but only ones that successfully
# recruit. The whole point is to see the relationships of the distance traveled 
# by SUCCESSFUL recruits to PLD. Those lost to sea we would assume would travel
# long distances, but the whole point is to see how physical processes/structure
# influence the distance traveled by successful recruits. 


import arcpy
from arcpy.sa import *
import numpy as np
import pandas as pd


dest_fcs = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts\COMBINED.gdb'
npys_path = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts\sim1101\outputs\npy' # for starting coords. They are the same for every time period.
land = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_FINAL'
out_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis_particles\spatial\dist_dataprep.gdb'
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'
mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'
contours = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\B07_bathy_merge'

arcpy.env.workspace = out_gdb

# make land raster and polygon to use for snapping and for distance analysis
if not arcpy.Exists('land_r_500'):
    arcpy.CopyFeatures_management(land, 'land_fc')
    arcpy.AddField_management('land_fc', 'land', 'SHORT')
    arcpy.CalculateField_management('land_fc', 'land', 1)
    arcpy.FeatureToRaster_conversion('land_fc', 'land', 'land_r_500', '500')
    arcpy.RasterToPolygon_conversion('land_r_500', 'land_r_500_poly', 'NO_SIMPLIFY')


# get origin coordinates from npys
if not arcpy.Exists('origin_coordinates'):
    sections = 9
    df = pd.DataFrame(columns=['traj_id_u', 'traj_id', 'o_lon', 'o_lat'])
    for section in range(1,sections+1):
        lons = np.load(os.path.join(npys_path, f'lon_{str(section)}.npy'))
        lats = np.load(os.path.join(npys_path, f'lat_{str(section)}.npy'))
        df_temp = pd.DataFrame()
        df_temp['o_lon'] = list(lons)
        df_temp['o_lat'] = list(lats)
        df_temp['traj_id'] = df_temp.index + 1
        df = df.append(df_temp, ignore_index=True)        
    df['traj_id_u'] = df.index + 1
    # output to fc
    x = np.array(np.rec.fromrecords(df.values))
    names = df.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    sr = arcpy.SpatialReference(4326)
    arcpy.da.NumPyArrayToFeatureClass(x, os.path.join(out_gdb, 'origin_coordinates'), ('o_lon', 'o_lat'), sr)


# Clean up destination points uids
# The traj_ids are not unique. They repeat for each of the 9 sections, plus they
# get duplicated when they settle on mpas that overlap, so I can't just do 
# groupby/agg. I would be grouping at least 9 different points.
arcpy.env.workspace = dest_fcs
fcs_orig = arcpy.ListFeatureClasses('destpts*')
for fc in fcs_orig:
    print("Copying feature "+fc)
    arcpy.CopyFeatures_management(fc, os.path.join(out_gdb, fc))
# go row by row and delete a row if it is the same as the previous one
arcpy.env.workspace = out_gdb
fcs_out = arcpy.ListFeatureClasses('destpts*')
for fc in fcs_out:
    print("Removing duplicates and adding new ID in "+fc)
    previous_value = None
    with arcpy.da.UpdateCursor(fc, ['traj_id']) as cursor:
        for row in cursor:
            current_value = row[0]
            if current_value == previous_value:
                cursor.deleteRow()
            previous_value = current_value
    # assign a new universal traj_id
    arcpy.AddField_management(fc, 'traj_id_u', 'LONG')
    ind = 1
    with arcpy.da.UpdateCursor(fc, ['traj_id_u']) as cursor:
        for row in cursor:
            row[0] = ind
            ind += 1
            cursor.updateRow(row)


# select particles that are not from the MPAs that we are excluding
to_include = []
with arcpy.da.SearchCursor(mpas_exclude, ['uID_20201124', 'exclude']) as cursor:
    for row in cursor:
        if not row[1] == 1:
            to_include.append(row[0])
arcpy.env.workspace = out_gdb
fcs_out = arcpy.ListFeatureClasses('destpts*')
# select out particles that experience mortality
for fc in fcs_out:
    print("Selecting particles that settle or strand in feature "+fc)
    arcpy.MakeFeatureLayer_management(fc, 'temp_lyr')
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'NEW_SELECTION', '"uID_part" IN {}'.format(str(tuple(to_include))))
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"mortstep" <> -1')
    arcpy.CopyFeatures_management('temp_lyr', 'selattr_'+ fc.split('_', 1)[1])
    arcpy.Delete_management('temp_lyr')

# create features for clipping points that strand/settle
con30 = arcpy.SelectLayerByAttribute_management(contours, 'NEW_SELECTION', '"ContourMax" = 30')
if not arcpy.Exists('landmask_buff30'):
    arcpy.Buffer_analysis(land, 'landmask_buff30', 30)
if not arcpy.Exists('land_mpas_con30'):
    arcpy.Merge_management(['landmask_buff30', mpas, con30], 'land_mpas_con30')
arcpy.Delete_management('landmask_buff30')

# clip ones that settle/strand
fcs = arcpy.ListFeatureClasses(wild_card='selattr_*')
for fc in fcs:
    print('clipping to land_mpas_con30 for {}'.format(fc))
    arcpy.Clip_analysis(fc, 'land_mpas_con30', 'clip_'+ fc.split('_', 1)[1])

# make rasters of points that settle/strand
snapras = 'land_r_500'
arcpy.env.snapRaster = snapras
fcs = arcpy.ListFeatureClasses(wild_card='clip*')
newSnapRas = True
for fc in fcs:
    print(f'Converting {fc} to raster')
    pld = fc.split('_')[-1]
    date = fc.split('_')[-2]
    outname = 'recruit_{}_{}'.format(date, str(pld))
    arcpy.PointToRaster_conversion(fc, '', outname, cellsize=1000)
    if newSnapRas:
        arcpy.env.snapRaster = outname
        newSnapRas = False


# DEALING WITH PARTICLES ON LAND
# There are many scenarios I need to account for here.
# Not all raster cells will overlap with the ocean. Some are completely contained
# by land. In this case there won't be a centroid in the ocean and the distance
# analysis won't work.
# Ther are also some land cells that are only 1 cell wide with a settlement cell inside them, and in this
# case particles on either side may get assigned wrong (see Calvert and the beach)
# Therefore, I am taking just the cells that overlap with the ocean and creating
# centroids for those. Then I am doing a Near analysis to find the closest reference
# point for every destination particle.
# Also, I need to work with a different resolution than the land raster resolution.
# Otherwise, too many cells line up perfectly fully on land, and there are long
# stretches were there are tons of particles but the raster cells are all on land.
# Using a larger cell size account for this (for the most part, see MPA south of Westport).
# However, this doesn't solve everything. See next group of comments...

# make centroids of just ones that overlap the ocean
rasts = arcpy.ListRasters('recruit_1*')
arcpy.env.extent = 'land_r_500' # super important to set this, cells weren't being included initially
outCellStats = CellStatistics(rasts, 'VARIETY', 'DATA')
outCellStats.save('recruit_All')
arcpy.AddField_management('recruit_All', 'uID_rast', 'LONG')
arcpy.CalculateField_management('recruit_All', 'uID_rast', '!OBJECTID!') # values need to be unique or else RasterToPolygon merges non-unique values
arcpy.RasterToPolygon_conversion('recruit_All', 'recruit_All_poly', 'NO_SIMPLIFY', 'uID_rast')
arcpy.Erase_analysis('recruit_All_poly', 'land_r_500_poly', 'recruit_All_erase')
arcpy.FeatureToPoint_management('recruit_All_erase', 'recruit_All_centroid1')
# Even with a near analysis, there are still some areas where raster cells are all
# within land and there is no nearby ocean cell. So...
# also make centroids out of buffered versions of cells within land.
arcpy.Clip_analysis('recruit_All_poly', 'land_r_500_poly', 'recruit_All_polyclip')
arcpy.Buffer_analysis('recruit_All_polyclip', 'recruit_All_polyclipbuff', 100)
arcpy.Erase_analysis('recruit_All_polyclipbuff', 'land_r_500_poly', 'recruit_All_polyclipbufferase1')
arcpy.Buffer_analysis('recruit_All_erase', 'recruit_All_erasebuff', 110)
arcpy.Erase_analysis('recruit_All_polyclipbufferase1', 'recruit_All_erasebuff', 'recruit_All_polyclipbufferase2')
arcpy.FeatureToPoint_management('recruit_All_polyclipbufferase2', 'recruit_All_centroid2')
# merge the two versions and manage attributes
arcpy.Merge_management(['recruit_All_centroid1', 'recruit_All_centroid2'], 'CENTROIDS_dest')
arcpy.AddField_management('CENTROIDS_dest', 'uID_ref_dest', 'SHORT')
arcpy.CalculateField_management('CENTROIDS_dest', 'uID_ref_dest', '!OBJECTID!')

# assign nearest distance reference point to each particle
fcs = arcpy.ListFeatureClasses(wild_card='clip*')
for fc in fcs:
    out_name = 'near_'+ fc.split('_', 1)[1]
    if not arcpy.Exists(out_name):
        print('Copying '+fc)
        arcpy.CopyFeatures_management(fc, out_name)
        print('Near analysis on '+fc)
        arcpy.Near_analysis(out_name, 'CENTROIDS_dest', 1500, 'NO_LOCATION', 'NO_ANGLE', 'PLANAR')
        # the resulting NEAR_FID relates the the OBJECTID in CENTROIDS_dest, 
        # which matches the 'uID_ref' field I created, so just add a new field
        # in the near fc and calculate.
        print('Adding and calculating new near field id for '+fc)
        arcpy.AddField_management(out_name, 'uID_ref_dest', 'SHORT')
        arcpy.CalculateField_management(out_name, 'uID_ref_dest', '!NEAR_FID!') 


####################
# Do the same thing for the ORIGIN points

# Make raster
snapras = 'land_r_500'
arcpy.env.snapRaster = snapras
arcpy.env.extent = snapras
arcpy.env.outputCoordinateSystem = 'land_r_500'
arcpy.PointToRaster_conversion('origin_coordinates', '', 'origin_coordinates_r', cellsize=1000)
# Centroids
arcpy.RasterToPolygon_conversion('origin_coordinates_r', 'origin_coordinates_rp', 'NO_SIMPLIFY', 'Value')
arcpy.Erase_analysis('origin_coordinates_rp', 'land_r_500_poly', 'origin_coordinates_rp_erase')
arcpy.FeatureToPoint_management('origin_coordinates_rp_erase', 'origin_centroid1')
# also make centroids out of buffered versions of cells within land
arcpy.Clip_analysis('origin_coordinates_rp', 'land_r_500_poly', 'origin_coordinates_rp_clip')
arcpy.Buffer_analysis('origin_coordinates_rp_clip', 'origin_coordinates_rp_clipbuff', 100)
arcpy.Erase_analysis('origin_coordinates_rp_clipbuff', 'land_r_500_poly', 'origin_coordinates_rp_clipbufferase1')
arcpy.Buffer_analysis('origin_coordinates_rp_erase', 'origin_coordinates_rp_erasebuff', 110)
arcpy.Erase_analysis('origin_coordinates_rp_clipbufferase1', 'origin_coordinates_rp_erasebuff', 'origin_coordinates_rp_clipbufferase2')
arcpy.FeatureToPoint_management('origin_coordinates_rp_clipbufferase2', 'origin_centroid2')
# merge the two versions and manage attributes
arcpy.Merge_management(['origin_centroid1', 'origin_centroid2'], 'CENTROIDS_origin')
arcpy.AddField_management('CENTROIDS_origin', 'uID_ref_origin', 'SHORT')
arcpy.CalculateField_management('CENTROIDS_origin', 'uID_ref_origin', '!OBJECTID!')
# Assign back to origin points
print('Copying origin points')
outname = 'near_ORIGIN'
arcpy.CopyFeatures_management('origin_coordinates', outname)
print('Near analysis on '+outname)
arcpy.Near_analysis(outname, 'CENTROIDS_origin', 5000, 'NO_LOCATION', 'NO_ANGLE', 'PLANAR')
print('Adding and calculating new near field id for '+outname)
arcpy.AddField_management(outname, 'uID_ref_origin', 'LONG')
arcpy.CalculateField_management(outname, 'uID_ref_origin', '!NEAR_FID!')


#########################################

# ONLY RUN THIS ONCE I KNOW EVERYTHING WORKED

# # delete intermmediate fcs, this gdb is massive
# arcpy.Delete_management('land_fc')
# arcpy.Delete_management('land_mpas_con30')
# fcs = arcpy.ListFeatureClasses(wild_card='origin_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)
# fcs = arcpy.ListRasters(wild_card='origin_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)
# fcs = arcpy.ListFeatureClasses(wild_card='destpts_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)
# fcs = arcpy.ListFeatureClasses(wild_card='selattr_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)
# fcs = arcpy.ListFeatureClasses(wild_card='clip_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)
# fcs = arcpy.ListFeatureClasses(wild_card='recruit_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)
# fcs = arcpy.ListRasters(wild_card='recruit_*')
# for fc in fcs:
#     arcpy.Delete_management(fc)


#########################################

# I now have the two centroid datasets to calculate distances between:
# "CENTROIDS_origin(_exclude)" and "CENTROIDS_dest"
# In the next script, I will step through each point in the origins centroid 
# dataset, and then I can do a selection to get the centroids that are its destination

# I also have the original origin and destination points with an attribute that
# relates them to the centroid datasets (near_ORIGIN and near_{date}_pld{})
# So, once I calculate the distances, I can join them back to one of these datasets.
