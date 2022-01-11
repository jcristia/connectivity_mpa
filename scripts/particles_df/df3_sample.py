# Sample from the larger data for building and testing
# use env: plotting. Arcpy enviros don't work.

import pandas as pd
import numpy as np

output_dir_all = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particles_df\output_all'
output_dir_ind = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particles_df\output'


##########################################
# Use pandas and sample from individual feather files
# They are separated out by date and pld, so this makes it easy
# This will be just for the DISTANCE MODEL

fraction = 0.0005
output_label = '0005'

files = [os.path.join(output_dir_ind, f) for f in os.listdir(output_dir_ind)]
# Get unique list of MPAs
df_mpa = pd.read_feather(files[0])
mpas = list(df_mpa.mpa_part_id_orig.unique())

# Process each feather file
df_pd = pd.DataFrame()
for file in files:
    df = pd.read_feather(file)
    # Remove rows that don't have a distance value
    df_dist = df[~df.distance.isna()]
    # Remove ones that have exposure values of 0 (there are only ~4000 out of 2.7 million)
    df_exp = df_dist[df_dist.total_exposure != 0]
    # For distance values of 0, assign a value between 1 and 250
    # Distances can be zero when the origin and dest ref points are the same. We make them 1 so that we can use a gamma distribution in the model. Points are spaced up to 500m apart, so they could have traveled up to that distance anyways.
    #df_sel_length = len(df_exp.loc[(df_exp.distance == 0), 'distance'])
    #df_exp.loc[(df_exp.distance == 0), 'distance'] = np.random.randint(1, 353, df_sel_length) # uniform dist.
    for mpa in mpas:
        df_mpa = df_exp[df_exp.mpa_part_id_orig == mpa]
        if len(df_mpa) < 1: # some mpas didn't have any that settled and therefore there is no distance value. Skip these.
            continue
        # sample fraction
        df_smp = df_mpa.sample(frac=fraction)
        # if none are selected, then just select 1 row so that MPA is represented
        if len(df_smp) < 1:
            df_smp = df_mpa.sample(n=1)
        # append
        df_pd = df_pd.append(df_smp, ignore_index=True)

# drop columns I no longer need
df_pd = df_pd.drop(['traj_id_series', 'dest_x', 'dest_y', 'year', 'period_id', 'provstate', 'settle'], axis=1)
# output to feather
df_pd.to_feather(os.path.join(output_dir_all, f'sample_{output_label}_training.feather'))

print(len(df_pd))


# now get holdout data that is equivalent to 30% (as if the above is 70% of all data)
# divide the length of the dataframe above by 2.333

length_new = len(df_pd) / 2.333
fraction_new =  (fraction * length_new) / len(df_pd)
df_ho = pd.DataFrame()
for file in files:
    df = pd.read_feather(file)
    # Remove rows that don't have a distance value
    df_dist = df[~df.distance.isna()]
    # Remove ones that have exposure values of 0 (there are only ~4000 out of 2.7 million)
    df_exp = df_dist[df_dist.total_exposure != 0]
    # For distance values of 0, assign a value between 1 and 353 (length of diaganol of 250m square)
    # Distances can be zero when the origin and dest ref points are the same. We make them >0 so that we can use a gamma distribution in the model. Points are spaced up to 500m apart, so they could have traveled up to that distance anyways.
    #df_sel_length = len(df_exp.loc[(df_exp.distance == 0), 'distance'])
    #df_exp.loc[(df_exp.distance == 0), 'distance'] = np.random.randint(1, 353, df_sel_length) # uniform dist.
    # Remove from the dataset any that I already selected
    df_merge = df_exp.merge(
        df_pd, on='traj_id_overall', 
        indicator=True,
        suffixes=(None,'y'),
        how='left').query('_merge=="left_only"').drop('_merge', axis=1)
    df_merge = df_merge.drop(['traj_id_series', 'dest_x', 'dest_y', 'year', 'period_id', 'provstate', 'settle', 'origin_xy', 'origin_yy', 'mpa_part_id_origy', 'pldy', 'monthy', 'total_exposurey', 'distancey'], axis=1)
    for mpa in mpas:
        df_mpa = df_merge[df_merge.mpa_part_id_orig == mpa]
        if len(df_mpa) < 1: # some mpas didn't have any that settled and therefore there is no distance value. Skip these.
            continue
        # sample fraction
        df_smp = df_mpa.sample(frac=fraction_new)
        # if none are selected, then just select 1 row so that MPA is represented
        if len(df_smp) < 1:
            df_smp = df_mpa.sample(n=1)
        # append
        df_ho = df_ho.append(df_smp, ignore_index=True)

# output to feather
df_ho.to_feather(os.path.join(output_dir_all, f'sample_{output_label}_holdout.feather'))

# AT VERY LOW fractions, the resulting numbers aren't perfect because the length
# of the initial df_mpa is 0 so we add in single rows that weren't initially
# selected. This compounds and we end up with many more rows. 

print(len(df_ho))




##########################################
# # ARCHIVE (Vaex approach)
# # Sample from entire dataframe of 192 million records
# # I don't think I am using Vaex efficiently, this was really slow

# # YOU HAVE TO RUN THIS ALL AT ONCE. I don't know if it is a Visual Studio thing
# # and the way Vaex doesn't keep things in memory, but basically it doesn't work
# # if you run it in pieces.
# # import vaex
# df = vaex.open(os.path.join(output_dir_all, 'particles_all.feather'))

# # Get unique list of MPAs
# mpas = df.unique('mpa_part_id_orig')

# # Go through each mpa id and select a percentage of particles.
# # This will allow me to get a representation from each MPA as opposed to just
# # doing a random sample from the entire dataframe.
# # THIS WAS REALLY SLOW
# df_pd = pd.DataFrame()
# for mpa in mpas:
#     df_sel = df[df.mpa_part_id_orig==mpa]
#     df_smp = df_sel.sample(frac=0.0001)
#     if df_smp.length() < 1:
#         df_smp = df.sel.sample(n=1)
#     df_t = df_smp.to_pandas_df() # this probably isn't ideal, but I whatever
#     df_pd = df_pd.append(df_t, ignore_index=True)    

# df_pd.to_feather(os.path.join(output_dir_all, 'sample_1.feather'))
# # Read in and check:
# df_1 = vaex.open(os.path.join(output_dir_all, 'sample_1.feather'))