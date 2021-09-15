# Make smaller dataframe for building and testing
# YOU HAVE TO RUN THIS ALL AT ONCE. I don't know if it is a Visual Studio thing
# and the way Vaex doesn't keep things in memory, but basically it doesn't work
# if you run it in pieces.

import vaex
import pandas as pd

output_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particles_df\output_all'
df = vaex.open(os.path.join(output_dir, 'particles_all.feather'))

# Get unique list of MPAs
mpas = df.unique('mpa_part_id_orig')

# Go through each mpa id and select a percentage of particles.
# This will allow me to get a representation from each MPA as opposed to just
# doing a random sample from the entire dataframe.
# THIS WAS REALLY SLOW
df_pd = pd.DataFrame()
for mpa in mpas:
    df_sel = df[df.mpa_part_id_orig==mpa]
    df_smp = df_sel.sample(frac=0.0001)
    if df_smp.length() < 1:
        df_smp = df.sel.sample(n=1)
    df_t = df_smp.to_pandas_df() # this probably isn't ideal, but I whatever
    df_pd = df_pd.append(df_t, ignore_index=True)    

df_pd.to_feather(os.path.join(output_dir, 'sample_1.feather'))


# Read in and check:
df_1 = vaex.open(os.path.join(output_dir, 'sample_1.feather'))

