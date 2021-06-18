# ABANDONING THIS APPROACH
# it isn't go to tell me what I initially thought it was


# analyzing particle count sensitivity

# selecting MPAs of different size from different locations
# select min, median, max from outer coast and from Salish Sea
# median size = 722,624

# Salish Sea:
# small uID_202011: 353
# median uID_202011: 89
# large: uID_202011: 110
# Outer:
# small uID_202011: 318
# median uID_202011: 178
# large uID_202011: 232

import os
import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


out_ncs = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particle_count\outputs\nc'
particle_factor = [1] # which gets multiplied by 84 releases
#part_uID = [353,89,110,318,178,232] # MPA unique part IDs
#part_uID = [353,89] # MPA unique part IDs
areas = {
    353: 8.964972,
    89: 722855.592613,
    110: 112931048.647291,
    318: 105499.319932,
    178: 775941.62311,
    232: 1224332319.996438
} # just manually entered from feature class

df_columns = ['particle_factor', 'part_uID', 'partperarea', 'total_strand', 'ratio_strandtototal']

df_strand = pd.DataFrame(columns=df_columns)
for mpa in part_uID:
    # get nc files for uID
    nc_files = os.listdir(out_ncs)
    mpa_ncs = []
    for ncf in nc_files:
        id = int(ncf.split('_')[1])
        if id == mpa:
            mpa_ncs.append(ncf)
    for n in mpa_ncs:
        pf = int((n.split('_')[2]).split('.')[0])
        ds = nc.Dataset(os.path.join(out_ncs, n))
        traj = ds.variables['trajectory']
        total_particles = len(traj)
        partperarea = total_particles / areas[mpa]
        status = ds.variables['status']
        total_strand = 0
        for i in status:
            if 1 in i:
                total_strand += 1        
        ratio_strandtotal = total_strand / total_particles
        df_strand_temp = pd.DataFrame([[pf, mpa, partperarea, total_strand, ratio_strandtotal]], columns=df_columns)
        df_strand = df_strand.append(df_strand_temp, ignore_index=True)


# plot line for each part_uID of strand/total by particle factor
# once a line flattens, stop opendrift simulations for this part_uID

sns.set()
g = sns.relplot(
    x='particle_factor', 
    y='ratio_strandtototal',
    hue='part_uID', 
    kind='line', 
    data=df_strand,
    markers=True, 
    ci=None)
g.savefig('ratio_strandtotal.png')