# Calculate fetch(exposure) for each reference point
# This is mostly Ed Gregr's code. I simplified a few things and had to change
# a few things for Python 3.
# Ed Gregr REFERENCE:
# 

import arcpy
import os
import sys
import math
import time
import datetime

####################################################
# Inputs

ptFC = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch1_dataprep.gdb\refcoords_5remove'
landFC = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch1_dataprep.gdb\land_2merge'
geoDBLocation = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch2_measure.gdb'
workingDir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\fetch\fetch_csv'
bearingInterval = 12  # make evenly divisble by 360
fetchDistance = 500000
# It's not about the absolute max distance to other land - otherwise, you would
# have to consider all the land around the entire Pacific, it's more so about
# the distance of what we would consider that there is no longer any coastal
# effect "protecting" the exposure of a site. Certainly distances like 2000 or
# 20,000km are pretty much equal at that point (aside from something like a
# tsunami from Japan).
# Therefore, I'll use a value of 500km. This would allow some effect from
# Alaska, but not all the way out to the pandhandle. If I make this number too 
# large, then I am heavily skewing the data to mean that exposure at these upper 
# levels matter. I think the only long distance exposure that may have influence
# is north-south movement through Hecate Straight. So for instance, while the 
# points on northern vancouver island would still receive significant waves from
# those areas, but maybe not as routinely as they would from directly West.
# Also, remember, I am looking at a 2D projection. On a sphere, Alaska would 
# provide no protection from circulating ocean currents.

#####################################################
# Set up

startTime = datetime.datetime.now()
print(startTime)
timeStamp = startTime.strftime("%m%d_%H%M")

arcpy.env.workspace = geoDBLocation
arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = False
arcpy.env.XYResolution = "0.00000001 Meters"

startBearing = 0 #Define the start and end bearings
endBearing = 359.999
features = [] # Create a list that will hold each of the Polyline objects (fetch lines) created around each point
g = arcpy.Geometry() # Create an empty Geometry object. This is used to convert features into Geometry.

#Get the projection of the land feature class, and use it as the base coordinate system.
coord_sys = arcpy.Describe(landFC).spatialReference
arcpy.env.outputCoordinateSystem = coord_sys

#Create a copy of the input points to make sure OID starts at 1 instead of 0. 
pointFile = arcpy.FeatureClassToFeatureClass_conversion(ptFC,geoDBLocation,"InputPoints_" + timeStamp)
ptFC = pointFile
countInputPoints = arcpy.GetCount_management(pointFile)
print(str(countInputPoints) + " points read.")

#Prepare the land polygons for processing. Dissolve. It needs to be all one feature.
landDissolve = "land_dissolved" + timeStamp
arcpy.Dissolve_management(landFC,landDissolve,"","","MULTI_PART")
landGeomList = arcpy.CopyFeatures_management(landDissolve, g)
arcpy.Delete_management(landDissolve)
#This is the land geometry object that will be used for all "difference" operations to erase fetch line geometry. 
landGeom = landGeomList[0]
print( "There are " + str(landGeom.partCount) + " land polygons in " + str(len(landGeomList)) + " multi-part feature.")

#Create a temp feature class for the line features (in memory objects are used to reduce processing time). This will be used to assign the attributes for each fetch line, before copying the fetch lines for each point to the output feature class.
outFC = "FetchLines_" + timeStamp
geomErrName = "GeometryErrors_" + timeStamp
geomErrInMemoryName = "in_memory/GeomErrors_" + timeStamp
linesTemp = arcpy.CreateFeatureclass_management("in_memory", "linesTemp" + timeStamp, "POLYLINE", "", "", "", coord_sys)
arcpy.CreateFeatureclass_management(geoDBLocation, outFC, "POLYLINE", "", "", "", coord_sys)
geomErrors = arcpy.CreateFeatureclass_management(geoDBLocation, geomErrName, "POINT", "", "", "", coord_sys)
arcpy.AddField_management(linesTemp,"Point_ID","LONG")
arcpy.AddField_management(outFC,"Point_ID","LONG")

#Empty arcpy Point object
point = arcpy.Point()
#Empty arcpy array objects (to hold coordinate pairs for the line ends and multi-part lines)
array = arcpy.Array()
MultiPart = arcpy.Array()

#Prepare the output csv file header for calculated fetch lengths. 
outFileName = os.path.join(workingDir, "Fetch_" + timeStamp + ".csv")
outFile = open(outFileName, 'a+')
outData = "Point_ID,X,Y"
bearing = startBearing
while bearing <= endBearing:
    outData += ",bearing" + str(int(bearing))
    bearing = bearing + bearingInterval
outData += "\n"
#print outData
outFile.write(outData)
outFile.close()	
outData = ""
#Open the output file and keep it open for the duration.
outFile = open(outFileName, 'a+')

#Initialize the dictionary to store the fetch distances, and the bearing list to store the list of bearings. Initialize the fetch distances to a default value of zero. (For edge conditions, where the point is exactly on the coastline, some fetch bearings may have a valid calculated distance of zero. Also for points that are inside the coastline polygons--i.e. on land--the fetch distances for all bearings will be zero.)
d = {}
bearingList = []
countFetchLines = 0

bearing = startBearing
while bearing <= endBearing:
    bearingList.append(bearing)
    d[bearing] = 0
    bearing = bearing + bearingInterval



#########################################################
# Processing

#Start processing each point, one at a time.
startTime = datetime.datetime.now()
lineTimer = datetime.datetime.now()
print("\n*****************************************")
print("Starting to process input points")
processedPoints = 0

for row in arcpy.da.SearchCursor(ptFC, ["OID@", "SHAPE@XY"]):

    #Get the x,y coordinates for the current point.
    oid = row[0]
    x, y = row[1]	
    #Create a point geometry object.
    ptGeom = arcpy.PointGeometry(arcpy.Point(x,y))
    
    #Write the OID, x and y coordinates to the output string
    outData = ""
    outData += str(oid) + "," + str(x) + "," + str(y)
    
    #Initialize the bearing
    bearing = startBearing
    #Iterate through the bearings and create a fetch line for each bearing. Store the fetch lines as parts of a multi-part line geometry. 
    while bearing <= endBearing:
        
        #Calculate the end points of the current fetch line, based on the bearing and source x,y coordinates. 
        a = math.radians(bearing)
        deltaY = fetchDistance*(math.cos(a))
        deltaX = fetchDistance*(math.sin(a))
        x2 = x + deltaX
        y2 = y + deltaY
        
        #Get the START point of the bearing line
        array.add(arcpy.Point(x,y))
        
        #Get the END point of the bearing line
        point.X = x2
        point.Y = y2
        array.add(point)
        
        #Add a line segment to the multi-part array. This will eventually contain one line part for each bearing. 
        MultiPart.add(array)

        array.removeAll()
        #Increment the bearing, move on to the next fetch line.
        bearing = bearing + bearingInterval
        countFetchLines = countFetchLines + 1
            
    #Create a polyline geometry using the array of calculated line endpoints. Assign a projection, "coord_sys" which has already been defined.  
    currentLine = arcpy.Polyline(MultiPart,coord_sys)
    #reset the multi-part line array
    MultiPart.removeAll()
    
    #Erase the current set of fetch lines using the difference operator on the land Geometry.
    t = datetime.datetime.now()
    try:
        lineErased = currentLine.difference(landGeom)
    except:
        print("\nGeometry difference error on point " + str(oid))
        rejectedPoint = arcpy.PointGeometry(arcpy.Point(x, y))
        geomErrorsInMemory = arcpy.CopyFeatures_management(rejectedPoint, geomErrInMemoryName)
        #Save the rejected point to a file GDB feature class
        arcpy.Append_management(geomErrorsInMemory, geomErrors)
        raw_input("Press Enter to continue...")
        continue
        
    el = ( datetime.datetime.now() - t).total_seconds()
    #print("\n" + str(el) + " sec geometry difference. " + str(lineErased.partCount) + " parts.")

    #Assume that there are multiple parts to this erased line geometry. Iterate through all parts and find the ones that touch the original point geometry. Often, the highest part number contains the line segment touching the original point, so we'll start iterating from there and go backwards.

    t = datetime.datetime.now()
    for i in reversed(range(lineErased.partCount)): 

        #get the current part of the erased line, convert it to a Polyline geometry
        currentLineSegment = arcpy.Polyline(lineErased.getPart(i))
        currentLineArray = lineErased.getPart(i)

        if currentLineSegment.touches(ptGeom):
            #We found a fetch line segment touching the input point!
            #Append this line segment to the list of features that will be written to the output feature class (fetch lines around each point).
            features.append(currentLineSegment)
            
            #Disassemble the part, calculate its bearing. 
            deltaX = currentLineArray[1].X - currentLineArray[0].X
            deltaY = currentLineArray[1].Y - currentLineArray[0].Y
            a = math.atan2(deltaX, deltaY)
            
            #Round off the bearing to the nearest integer so it matches the exiting list of bearings
            bearing = round(math.degrees(a))
            if (bearing < 0):
                #correct negative angles for quadrants where deltaX < 0
                bearing = 360 + bearing
            
            #Record the fetch distance for the current bearing in the fetch distance dictionary
            d[int(bearing)] = str(int(round(currentLineSegment.length)))


    el = ( datetime.datetime.now() - t).total_seconds()
    #print(str(el) + " sec select touching fetch segments, calculate bearing/fetch dist.")

    #print("Multi-part selection/erease done  ... ")
    #We are now done with fetch distance calulcations for this point!
    #Write the fetch distances for this point to the output csv string. Reset the fetch distance dictionary at the same time.
    for b in bearingList:
        #get the fetch distance stored in the dictionary under bearing b, convert to string, and append to the output csv string.
        outData += "," + str(d[b])
        #reset the dictionary
        d[b] = 0

    #Add an EOL line break
    outData += "\n"
    
    #Write fetch distances for this point to disk (append to csv file).
    #outFile = open(outFileName, 'a+')
    outFile.write(outData)
    outFile.close()
    outData = ""
    #immediately re-open the file again, to establish a lock.
    outFile = open(outFileName, 'a+')
    
    #If there are fetch lines in the features buffer to be written to a feature class, then add the attribute Point_ID to each one and append to the in-memory workspace (in case of a point that is entirely on land, the features buffer will be empty, and there will be no fetch lines)
    if (len(features) > 0):
        t = datetime.datetime.now()
        #Copy the features[] polyline array to a new feature class. This step is necessary because you cannot append an array of line features using Append_management. It has to go through a copy first. 
        fetchTmp = arcpy.CopyFeatures_management(features, "in_memory/fetch" + timeStamp)
        arcpy.AddField_management(fetchTmp,"Point_ID","LONG")
        #Update the attributes for Point_ID
        with arcpy.da.UpdateCursor(fetchTmp,["Point_ID"]) as cursor:
            for row in cursor:
                row[0] = int(oid)
                cursor.updateRow(row)

        el = ( datetime.datetime.now() - t).total_seconds()
        #print(str(el) + " sec updated " + str(len(features)) + " fetch lines with Point_ID")
        
        t = datetime.datetime.now()
        arcpy.Append_management(fetchTmp, linesTemp)
        
        #Cleanup
        features = []
        arcpy.Delete_management(fetchTmp)
        el = ( datetime.datetime.now() - t).total_seconds()
        #print(str(el) + " sec append fetch lines to in-memory workspace")

    #Finished a point. Calculate the timing.
    pointTime = (datetime.datetime.now() - lineTimer).total_seconds()
    linesPerSec = (360 / bearingInterval)/pointTime
    #Reset the timer.
    lineTimer = datetime.datetime.now()	
    print("Point " + str(oid) + " completed at " + str("%.2f" % linesPerSec) + " lines/sec.")
    
    #Increment the number of processed points
    processedPoints = processedPoints + 1
    
    #If there are 50 or more points completed, then dump the existing fetch lines to disk. 
    if (processedPoints > 49):
        print("\nBatch of 10 points completed. Saving " + str(countFetchLines) + " fetch lines to disk")
        #Save the temp fetch lines from in-memory to disk
        arcpy.Append_management(linesTemp, outFC, "NO_TEST")
        #Delete all features from the temp feature class
        arcpy.DeleteFeatures_management(linesTemp)
        #Reset counters
        countFetchLines = 0
        processedPoints = 0
    


# outside for loop, Finished processing points
endTime = datetime.datetime.now()
timeDelta = (endTime - startTime).total_seconds()
print("\nFinished processing points. " + str(timeDelta) + " seconds")
outFile.close() #Close the output file, finally.

#Save any remaining temp fetch lines from in-memory to disk
countFetchLines = arcpy.GetCount_management(linesTemp)
if (int(countFetchLines[0]) > 0):
    print("Saving remaining fetch lines")
    arcpy.Append_management(linesTemp, outFC, "NO_TEST")




#########################################################
# Post processing

t = datetime.datetime.now()
#Join the output csv file with the input points, then export to a final shapefile.
print("Joining fetch table to final output...")
featureLayerName =  "ptsToJoin" + timeStamp
arcpy.MakeFeatureLayer_management(ptFC, featureLayerName)
oidField = arcpy.Describe(ptFC).OIDFieldName
arcpy.AddJoin_management(featureLayerName, oidField, outFileName, "Point_ID")
arcpy.CopyFeatures_management(featureLayerName, os.path.join(geoDBLocation, "FetchPoints_" + timeStamp))
el = ( datetime.datetime.now() - t).total_seconds()
print(str(el) + " sec final join")

#If there are no Geometry Errors, delete the error file.
countGeomErrors = str(arcpy.GetCount_management(geomErrors))
if countGeomErrors == "0":
    #If there were no geometry errors, delete the GeometryErrors and InputPoints feature classes. Otherwise, keep them both to help track down the source of error. 
    arcpy.AddMessage("No geometry errors found.")
    arcpy.Delete_management(geomErrors)
    arcpy.Delete_management(ptFC)
else:
    arcpy.AddMessage("There were " + str(countGeomErrors) + " points with geometry errors.")

#Cleanup stray data
arcpy.Delete_management(linesTemp)
finish = datetime.datetime.now()
print(finish)
