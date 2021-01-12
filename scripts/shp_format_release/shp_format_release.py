# script to format the mpa shapefiles that are used for particle release in 
# Opendrift

# env: arcgispro-py3-clone_MPACONNECTIVITY

import arcpy

mpas_fc = (
    r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA'
    '\mpas.gdb\M08_mpa_20201124_FINALCOPY')
mpas_shp_folder = (
    r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA'
    '\mpas_shp_release')
gdb = 'mpas.gdb'
fc = os.path.basename(os.path.normpath(mpas_fc))
arcpy.env.workspace = os.path.join(mpas_shp_folder, gdb)



if not arcpy.Exists(arcpy.env.workspace):
    arcpy.CreateFileGDB_management(mpas_shp_folder, gdb)
if not arcpy.Exists(fc):
    arcpy.CopyFeatures_management(mpas_fc, fc)
    arcpy.AddField_management(fc, 'part_num', 'LONG')

############################################
# equation for calculating particles to release

# get area min and max
with arcpy.da.SearchCursor(fc, ['Shape_Area']) as cursor:
    area_max = max(cursor)
    cursor.reset()
    area_min = min(cursor)

# min: 9
# max: 1,224,332,320 (Olympic Coast)
# mean: 15,814,497
# median: 722,623
# the largest Canadian MPA is Gwaii Haanas (well, one part of it): 310,867,762

# similar to my seagrass chapter, I will have 84 releases (release every 4 hours
# for two weeks)
# I want my base amount of particles per release to be 10
# If I do things propotionally, that equals 1.3 billion per release for the
# largest MPA. Obviously way too much.

####### Notes for how I calculated particle amount to release:

"""
I can't do a linear relationship between particle count and area. There is a 
massive range in order of magnitude of patch sizes. I wanted to make sure I 
released enough for the smaller patches without having to release a massive 
amount for larger patches.

I explored some simple clean equations, but none of them worked that well 
(log, ln, square root, etc). To get the relationship that I want, I need to 
define a few points. I the end this will be the easiest to justify.

I set the amounts I want to release for the min, max, and median patch size. 
For min, I wanted to release at least 10 points per release. This would be 840
total, which means my min connection strength would be 0.0012. For max, I 
examined the few largest patches. They are massive. I set this at 3,000 (I can 
then detect probabilities as low as 0.0000039). This would be for Olympic 
Coast, which is an order of magnitude larger than the others. I'm simply not 
as concerned with this one anyways. The biggest Canadian mpa is Gwaii Haanas, 
which would have 2300 particles. For the median, which is many order of 
magnitude less than the max, I set this at 750 particles per release. Based 
on a visual assessment of these patches, I think this is adequate.

I put these 3 data points into mycurvefit.com (quick and dirty), and fit a power
 curve. This allows for quick growth at the beginning and then growth slows.
Equation is y = 55.67552x^0.1905742

Compared to my seagrass chapter, it is proportionately more particles released. 
I'll be releasing a total of 27 million particles (over the 84 releases), 
whereas if I used the seagrass proportion I would only release 20 million.
And, since I did some sensitivity tests for chapter 1, I can say that I will at 
least know that it takes those into account indirectly.

When compared to Masson 2015, I am releasing way more particles than them per 
unit area. They also didn't offer any justification for their particle count.
I can say something like "it exceeds the particles released from past studies on
 the BC coast (e.g. ...; ...)"

However, I should also eventually do sensitivity tests on the smallest, median, 
and largest patches (or patches near these sizes that are spaced evenly 
throughout). Start with the particle counts I am using and go up from there. 
Even though I am now interested in all places that particles end up, the primary
 goal is still connections between MPAs. However, I guess since they are more 
 spaced out then this may not tell me much for a sensitivity analysis.

I looked into how other people have done this and it is never explained well. 
See Snauffer et al 2014 and Simons et al 2013.
I think my justifications are fine. It's clear that I am releasing a lot of 
particles.

I can also just skate over this in the manuscript and deal with it if a reviewer
 says something.

equation from mycurvefit.com (yes, quick and dirty)
y = 55.67552 * x^0.1905742
"""

############################################

# calculate particles based on area
with arcpy.da.UpdateCursor(fc, ['Shape_Area', 'part_num']) as cursor:
    for row in cursor:
        particles = 55.67552 * (row[0] ** 0.1905742)
        row[1] = particles
        cursor.updateRow(row)

############################################

# split out mpas into individual shapefiles

# export to shapefile
shp = os.path.join(mpas_shp_folder, 'mpa_.shp')
arcpy.CopyFeatures_management(fc, shp)

# create folder if it doesn't exists
mpas_shp = os.path.join(mpas_shp_folder, 'mpas_shp')
if not os.path.exists(mpas_shp):
    os.makedirs(mpas_shp)

# split
arcpy.SplitByAttributes_analysis(shp, mpas_shp, 'uID_202011')