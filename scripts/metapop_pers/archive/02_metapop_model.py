# deterministic model
# so basically, with this version, without a growth rate, the entire system
# will eventually empty out

# Based on Patricks metacommunity model:
# https://github.com/plthompson/mcomsimr/blob/master/R/MC_simulate.R

import arcpy
import pandas as pd
import numpy as np
import seaborn as sns


#### Inputs ####

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity'
conn = r'cluster_results\scripts\COMBINED.gdb\conn_avg_pld{}' # connectivity probabilities per PLD averaged over time
plds = [1,3,7,10,21,30,40,60]
patches = os.path.join(root, r'spatial\MPA\mpas.gdb\M09_mpa_joined')
out_gdb = os.path.join(root, r'scripts\metapop_pers\metapop_pers.gdb')

timesteps = 10000

popn_disp_props = [0.05] # Patrick did 10^-5 to 0.5 in log space

# I need to set this, otherwise for isolated patches, they will always be persistent.
# They would either stay the same or they would always grow based on retention rates.
# However, if I make ALL particles leave at every timestep then stuff empties out super quickly.
# That isn't realistic either.

growth_rate = 1.0

# at 1, births = deaths 
# should make this density dependent so a patch can only grow so big?
# Could also just set a max pop size as a proportion of starting value, then 
# just cut it back to that value at each timestep if it grows too big.

remove_retention = False




#### Get starting population based on patch area ####

field_names = [i.name for i in arcpy.ListFields(patches) if i.type != 'OID']
field_names = [i.name for i in arcpy.ListFields(patches) if i.name in ['uID_20201124', 'Shape_Area']]
cursor = arcpy.da.SearchCursor(patches, field_names)
patches = pd.DataFrame(data=[row for row in cursor], columns=field_names)

# assign proportionally with largest patch receiving 100 million individuals
# I'll probably want to scale this differently later
max_area = patches.Shape_Area.max()
patches['popn'] = (patches.Shape_Area * 100000000) / max_area
# for any patches with fewer than 100 individuals, assign them 100
patches.loc[patches.popn < 100, 'popn'] = 1000
patches = patches.round({'popn':0})

# convert to array
patches = patches.sort_values(by=['uID_20201124'])
popn_base = np.array(patches.popn)




#### Simulation for each PLD ####

# df to hold output data
out_patches = patches

for pld in plds:

    #### read in connectivity data as matrix ####

    # feature class to pandas df
    conns = os.path.join(root, conn.format(pld))
    field_names = [i.name for i in arcpy.ListFields(conns) if i.type != 'OID']
    cursor = arcpy.da.SearchCursor(conns, field_names)
    df = pd.DataFrame(data=[row for row in cursor], columns=field_names)

    # add missing meadows
    # there are a few meadows that have nothing coming into them, which means
    # I don't get a square matrix
    i = 1
    while i < 407:
        if i not in df.to_id.unique():
            df = df.append({'to_id':i, 'from_id':i}, ignore_index=True)
        i += 1

    # to pandas matrix df
    df_p = df.pivot(index='to_id', columns='from_id', values='prob_avg')

    # pandas df to numpy matrix
    conn_matrix = df_p.to_numpy()
    conn_matrix = np.nan_to_num(conn_matrix)
    if np.shape(conn_matrix)[0] != np.shape(conn_matrix)[1]:
        raise ValueError('Not a square matrix')

    # testing: remove retention:
    if remove_retention:
        np.fill_diagonal(conn_matrix, 0)

    for proportion in popn_disp_props:

        # reset population to start
        popn = popn_base

        # dataframe for correlations
        df_cc = pd.DataFrame(columns={'timestep', 'cc'})

        for i in list(range(timesteps)):

            #### calc amount to disperse ####
            popn_ts = popn * proportion
            popn_ts = np.around(popn_ts, 0)
            # if less than 1, set as 1
            popn_ts[popn_ts<1] = 1


            #### immigration ####
            immigrate = np.matmul(conn_matrix, popn_ts)


            #### net population after time step ####
            net = popn + immigrate - popn_ts
            net[net<0] = 0


            #### growth rate ####
            popn = net * growth_rate
            popn = np.around(popn, 0)


            #### track equilibrium ####
            # set any patches with individuals to 1
            # compare matrices between timesteps
            # patches will blink in and out, but eventually the matrix will be reach equilibrium
            pers = np.where(popn > 0, 1, 0)
            if (i != 0) and (i % 5 == 0): # not the first timestep and only every 10 time steps
                cc_m = np.corrcoef([pers_prev, pers])
                cc = cc_m[1,0]
                if np.isnan(cc):
                    cc = 1 # if there is no variance (e.g. all values are 1), then numpy outputs nan
                df_cc = df_cc.append({'timestep':i, 'cc': cc}, ignore_index=True)
            pers_prev = pers


        ### plot correlations ###
        df_cc = df_cc.astype({'timestep':'int64', 'cc':'float32'})
        g = sns.lineplot(data=df_cc, x='timestep', y='cc')
        g.set(ylim=(0.98, 1.01))
        str_prop = str(proportion).split('.')[1]
        g.figure
        if remove_retention:
            out_string = f'cor_pld{pld}_{str_prop}_noret.png'
        else:
            out_string = f'cor_pld{pld}_{str_prop}_ret.png'
        g.figure.savefig(os.path.join(root,'scripts/metapop_pers/figs', out_string))
        
        # Check figures if we have reached equilibrium


        #### add persistence to out_patches ####
        out_patches[f'pld{pld}_{str_prop}'] = pers

    break


# output to feature class table
# out_patches = out_patches.drop(columns=['Shape_Area'])
# x = np.array(np.rec.fromrecords(out_patches.values))
# names = out_patches.dtypes.index.tolist()
# x.dtype.names = tuple(names)
# if remove_retention:
#     out_string = 'metapop_pers_noretention.png'
# else:
#     out_string = 'metapop_pers_retention.png'
# arcpy.da.NumPyArrayToTable(x, os.path.join(out_gdb, out_string))


# In next script:
# copy in patch centroids and join table to points
# see how persistent patches cluster and if they line up with cycles (components)
# calc eigenvector centrality for each PLD conns and see how they compare
# A main result: which patches are consistent across PLDs? How can I show the
# variation? Is there some kind of summary plot? Or just multiple maps?

# Once I know my persistent cycles, I could then calc the leading eigenvalue
# just for those connections and then I can compare the levels of persistence
# between them.

# calculating eigenvalues just gives me a bunch of imaginary numbers
# my guess is that it is looking at the network as a whole, which doesn't
# operate as 1 metapopulation
# Perhaps, if I did it on a set of strongly connected components then I would
# get real numbers.