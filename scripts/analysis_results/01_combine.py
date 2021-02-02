# combine shapefile outputs for biology script
# env: arcpro clone

# for each time period, and for each of 9 groups:
# there will be connection line shapefiles for each of 8 PLDs
# there will be 1 destination point shapefiles for the 60 day PLD


import arcpy
import os


root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results'
shps = os.path.join(root, 'sim{}/outputs/shp')
sims = [
    '1101','1105','1108',
    '1401','1405','1408',
    '1701','1705','1708'
]
plds = [1,3,7,10,21,30,40,60]
gdb = r'COMBINED.gdb'
arcpy.env.workspace = os.path.join(root, gdb)
# create gdb if it doesn't exist
if not arcpy.Exists(arcpy.env.workspace):
    arcpy.CreateFileGDB_management(root, gdb)


# copy in centroids once
centroids = os.path.join(shps.format(sims[0]), 'patch_centroids.shp')
arcpy.CopyFeatures_management(centroids, 'patch_centroids')


# for each time period
# combine the line shapefiles for each group
# output to fgdb so that they can all be in 1 place
for sim in sims:
    folder = shps.format(sim)
    fl = os.listdir(folder)
    for pld in plds:
        files = []
        for file in fl:
            base_pld = file.split('_')[-1]
            if file.startswith('connectivity') and base_pld == 'pld{}.shp'.format(str(pld)):
                files.append(os.path.join(folder, file))
        arcpy.Merge_management(files, 'connectivity_{}_pld{}'.format(sim, str(pld)))


# merge destination points per time period
# this is going to be 3 million points per fc
for sim in sims:
    folder = shps.format(sim)
    fl = os.listdir()
    for pld in plds:
        files = []
        for file in fl:
            base_pld = file.split('_')[-1]
            if file.startswith('dest') and base_pld == 'pld{}.shp'.format(str(pld)):
                files.append(os.path.join(folder, file))
        arcpy.Merge_management(files, 'destpts_{}_pld{}'.format(sim, str(pld)))
