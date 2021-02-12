# create raster of the count of points that strand or are over the 30m bathymetry
# at the end of a PLD
# This will be recruitment.

# This script is currently set up to include all particles that strand before their
# chosen mortality time. It is a cumulative look at all potential positions while
# still applying some mortality.
# I will use that version for the prioritization chapter.
# However, for the analysis where I compare recruitment success by PLD, I may want
# to apply mortality differently.
# IT MIGHT BE EASIER TO JUST MAKE A SEPARATE SCRIPT AND GDB FOR THIS.


import arcpy
from arcpy.sa import *
import os


root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results'
destpts_dir = 'sim{}/outputs/shp/dest_pts'
destpts = 'dest_biology_pts_{}_pld{}.shp'
output_gdb = 'DEST_RAST.gdb'
plds = [60]  # only have 60 day pld data for now
sims = [
    '1101','1105','1108',
    '1401','1405','1408',
    '1701','1705','1708'
]


# get paths of dest_pts
arcpy.env.workspace = os.path.join(root, 'COMBINED.gdb')
files = arcpy.ListFeatureClasses(wild_card='destpts_*')
fcs = []
for fc in files:
    fcs.append(os.path.join(root, 'COMBINED.gdb', fc))

# create raster gdb and set as workspace
if not arcpy.Exists(os.path.join(root, output_gdb)):
    arcpy.CreateFileGDB_management(root, output_gdb)
arcpy.env.workspace = os.path.join(root, output_gdb)

# select points by location
landmask = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_FINAL'
mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'
contours = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\B07_bathy_merge'
con30 = arcpy.SelectLayerByAttribute_management(contours, 'NEW_SELECTION', '"ContourMax" = 30')
if not arcpy.Exists('contour30'):
    arcpy.CopyFeatures_management(con30, 'contour30')


for fc in fcs:
    print('selecting by location and attribute for {}'.format(os.path.basename(fc)))
    arcpy.MakeFeatureLayer_management(fc, 'temp_lyr')
    arcpy.SelectLayerByLocation_management('temp_lyr', 'INTERSECT', landmask, search_distance=50, selection_type='NEW_SELECTION')
    arcpy.SelectLayerByLocation_management('temp_lyr', 'INTERSECT', mpas, selection_type='ADD_TO_SELECTION')
    arcpy.SelectLayerByLocation_management('temp_lyr', 'INTERSECT', 'contour30', selection_type='ADD_TO_SELECTION')
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"uID_part" = "dest_id"')
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"mortstep" < "time_int" and "mortstep" <> -1') # include ones that settle before their mortality time
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"uID_part" = 170') # the particles for this feature are all messed up. Not sure what happened. Best to just remove it.
    arcpy.CopyFeatures_management('temp_lyr', 'selection_'+ os.path.basename(fc))
    arcpy.Delete_management('temp_lyr')
# so for a pld of 60, once I select for ones that didn't die, this pretty much removes everything
# that didn't strand, which means the spatial select by mpas and contour are not that meaningful.
# For a PLD of 60-75, with 1.6million particles, only 8 particles would be remaining at the end.
# This sounds extreme, but it makes sense. Therefore, settlement is super important.
# I will see different patterns with shorte PLDs. Some in open water may still be living.
# This is also a good reason to perhaps have a mortality rate that changes through time,
# then I am not losing so many right away.


# convert to raster
snapras = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_NN'
fcs = arcpy.ListFeatureClasses(wild_card='selection_*')
for fc in fcs:
    arcpy.env.snapRaster = snapras
    pld = fc.split('_')[-1]
    date = fc.split('_')[-2]
    outname = 'recruit_count_{}_{}'.format(date, str(pld))
    arcpy.PointToRaster_conversion(fc, '', outname, 'COUNT', cellsize=1000)


# add rasters by PLD to create one overall raster by PLD
rasters = arcpy.ListRasters(wild_card='recruit_count*')
for pld in plds:
    ras_plds = []
    for ras in rasters:
        base_pld = ras.split('_')[-1]
        if base_pld == 'pld{}'.format(str(pld)):
            ras_plds.append(ras)
    outras = CellStatistics(ras_plds, 'SUM', 'DATA')
    outras.save('recruit_count_ALL_{}'.format(str(pld)))