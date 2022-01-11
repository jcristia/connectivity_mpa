# NOT USING THIS ANYMORE

# Combine individual feather files into 1 big file
# env: plotting

import vaex

input_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particles_df\output'
output_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particles_df\output_all'

files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)]
df = vaex.open_many(files)
df.export_feather(os.path.join(output_dir, 'particles_all.feather'))
