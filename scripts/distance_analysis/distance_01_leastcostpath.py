# calculate the overwater distance between MPAs

# There is now an updated way to do this, since I last did it with the Cost Path tools:
# https://www.esri.com/arcgis-blog/products/spatial-analyst/analytics/doing-more-with-euclidean-distance-barriers-and-paths/
# https://community.esri.com/t5/arcgis-pro-questions/measure-distance-between-points-across-polygon-least-cost-path/td-p/66607

# See notes at the bottom for using PathCost instead of Shape_length when plotting


import arcpy
import os


default_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis\distance_analysis_mapping\Default.gdb'
arcpy.env.workspace = default_gdb
mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'
land = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_FINAL'
euc_lines_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis\distance_analysis_mapping\euc_lines.gdb'
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'


# select ones that are not deep in inlets
to_include = []
with arcpy.da.SearchCursor(mpas_exclude, ['uID_20201124', 'exclude']) as cursor:
    for row in cursor:
        if not row[1] == 1:
            to_include.append(row[0])
mpas_sel = arcpy.SelectLayerByAttribute_management(mpas, 'NEW_SELECTION', '"uID_20201124" IN {}'.format(str(tuple(to_include))))


# get centroids - make sure points fall within polys
arcpy.FeatureToPoint_management(mpas_sel, 'mpa_centroid_fc', 'INSIDE')


# copy land, clip, add attribute, convert to raster
arcpy.CopyFeatures_management(land, 'land_fc')
arcpy.AddField_management('land_fc', 'land', 'SHORT')
arcpy.CalculateField_management('land_fc', 'land', 1)
arcpy.FeatureToRaster_conversion('land_fc', 'land', 'land_r_500', '500')

# check if there are any points that overlap with the land raster
# raster to poly
arcpy.RasterToPolygon_conversion('land_r_500', 'land_r_500_fc', 'NO_SIMPLIFY', 'Value', 'MULTIPLE_OUTER_PART')
# copy centroids
arcpy.CopyFeatures_management('mpa_centroid_fc', 'mpa_centroid_fc_manualedit')
# !!!MANUALLY!!!: select by location where they intersect, manually move these few points

# rasterize with Feature to Raster tool to make sure I can control cell size and nothing gets lost
# snap to land raster
arcpy.env.snapRaster = 'land_r_500'
arcpy.FeatureToRaster_conversion('mpa_centroid_fc_manualedit', 'uID_20201124', 'mpa_centroid_r', '500')

# !!!MANUALLY!!!!:
# Some points overlap with each other (e.g. 2 points per 1 raster cell and
# therefore it won't treat these seperately). You can use the Find Identical tool
# on the Shape field with a search radius to find where these are. Move one of the points
# out.
# Also, if certain points fall within the same raster cell but are sufficiently apart
# then find identical won't necessarily identify these. I just searched for these
# manually.
# When I finally create the lines below, there is a line fc created for each point.
# The number of lines should match the number of points (i.e. 351)


# go through each point and find distances
# the tools are set up to find all the distances between ONE starting point and multiple destination points
# First, in the Euclidean distance tool, find the distances of every cell back to the one starting point
# Then, in Coast Path as Polyline, input the Destination points to get all the distances from start to dest.

with arcpy.da.SearchCursor('mpa_centroid_fc_manualedit', ['uID_20201124']) as cursor:
    for row in cursor:
        uID = row[0]
        # if uID < 248:  # I think I had this in for when I needed to run things in bits.
        #     continue
        print('processing uID {}'.format(uID))
        centroid_sel = arcpy.SelectLayerByAttribute_management('mpa_centroid_fc_manualedit', 'NEW_SELECTION', 'uID_20201124={}'.format(row[0]))
        arcpy.env.snapRaster = 'land_r_500'
        origin = 'mpa_centroid_r_{}'.format(uID)
        arcpy.FeatureToRaster_conversion(centroid_sel, 'uID_20201124', origin, '500')

        # Euclidean distance with barriers
        in_rast = origin
        max_dist = ''
        cellsize = 'land_r_500'
        outDirectionRaster = '' # you would use this one if you don't have barriers
        distance_method = 'PLANAR'
        in_barrier_data = 'land_r_500'
        outBackDirectionRaster = 'eucbackdirect' # this direction raster is for when you have barriers

        outEucDistance = arcpy.sa.EucDistance(in_rast, max_dist, cellsize, outDirectionRaster, distance_method, in_barrier_data, outBackDirectionRaster)
        # Save the output 
        outEucDistance.save('eucdistance')

        # Cost Path as Polyline
        inputDestinationLayer = 'mpa_centroid_r'
        inputCostLayer = 'eucdistance'
        inputDirectionLayer = 'eucbackdirect'
        outLines = os.path.join(euc_lines_gdb, 'euc_lines_{}'.format(uID))
        pathType = 'EACH_CELL'
        destfield = ''
        arcpy.sa.CostPathAsPolyline(inputDestinationLayer, inputCostLayer, inputDirectionLayer, outLines, pathType, destfield)

        # add origin_uid attribute
        arcpy.AddField_management(outLines, 'origin_id', 'SHORT')
        arcpy.CalculateField_management(outLines, 'origin_id', uID)

        arcpy.Delete_management(origin)
        arcpy.Delete_management('eucbackdirect')
        arcpy.Delete_management('eucdistance')



# merge FCs
arcpy.env.workspace = euc_lines_gdb
fcs = arcpy.ListFeatureClasses()
arcpy.Merge_management(fcs, os.path.join(default_gdb, 'euc_lines_ALL'))


# CHECK the shape length of the lines
# There are many random line segments that would just end
# in the middle of nowhere. I would move the point where it was supposed to
# end by one grid cell and this would fix it. There doesn't seem to be any rhyme
# or reason to it. You just need to check that everything makes sense.
# What is interesting is that the path cost value is close to what the actual
# distance should be, but for some reason it does not draw the line correctly.

# However, but now I see that it is not the same 2 for each dataset. I tested 
# euc_lines_2 and euc_lines_399 and there were only issues in euc_lines_2.
# Which makes me think that it might be different issues in each datasets.

# CHECK THE ENTIRE DATASET:
# I did this manually: in the ecu_lines_ALL dataset, I added a line_difference
# field and did PathCost - Shape_length. Any difference over ~6000 was usually a
# line that was cut short and did not draw all the way between each point.
# There were a lot of these actually.

# THEREFORE....!!!
# it is not worth moving all of these just to get a proper shape_length.
# I would also have to figure out which point is the problem. It also doesn't
# seem like a point doesn't have consistent problems where there are no lines
# going in/out, so I move it, it might create new problems. So...
# JUST USE PATHCOST for my distance.
# It is close to the line distance and it is all relative anyways. The biggest
# differences are those that have a lot of turns in them or those that are very
# long.