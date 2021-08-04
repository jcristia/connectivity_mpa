# a bar/line plot for % of connections that are within 20km that are realized
# do various threshold levels and group by each PLD
# e.g. PLD 1 (then a bar for 0.0001, 0.001, 0.01)

# should be connected in at least 1 direction

# Important citations:
# Shanks 2003:
#   - 2 evolutionary stable dispersal strategies:
#   - disperse less than 1km
#   - disperse > 20km
#   - Therefore, make reserves large enough to capture small movement, and
#   - space them 10-100km apart
# Burt 2014:
#   - compiles a good table of citations (page 8)
#   - page 21: recommends 20-100km spacing in BC
# Martone 2021:
#   - compiles the same info as above, but it is DFO, so a good resource
#   - Does acknowledge that it should be assessed at multiple scales and based
#   on different mpa criteria
#   - So they hesitate to give a number

# SO... I will start by using 20km


import arcpy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

dist_thres = 20  # 20 km distance threshold for mpa spacing
# probability thresholds
prob_thresh = [0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01]
prob_thresh_field = ['e8', 'e7', 'e6', 'e5', 'e4', 'e3', 'e2']
plds = [60,40,30,21,10,7,3,1]

# output to fc???
out_to_fc = False
# only Canadian connections???
can_only = True


# euclidean distances for all combinations of nodes
distances = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis\distance_analysis_mapping\Default.gdb\euc_lines_ALL"
# averaged connections by PLD
conns_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts\COMBINED.gdb'

# out gdb
out_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis\distance_analysis_mapping\rule_of_thumb.gdb'

# all mpa attributes from master dataset
mpa_att = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'


###############################


# read in distances as pandas df
field_names = [i.name for i in arcpy.ListFields(distances)]
cursor = arcpy.da.SearchCursor(distances, field_names)
dist_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
# drop fields
dist_df = dist_df.drop(['Shape'], 1)
# drop self connections from the dataframe
dist_df = dist_df[dist_df.DestID != dist_df.origin_id]
# convert to km
dist_df['distkm'] = (dist_df.PathCost)/1000.0
dist_df = dist_df.drop(['Shape_Length', 'length_difference', 'PathCost'], 1)

# access connection dataset
# drop self connections
arcpy.env.workspace = conns_gdb
conns = arcpy.ListFeatureClasses(wild_card = 'conn_avg_pld*')
conns_all = pd.DataFrame()
for fc in conns:
    field_names = [i.name for i in arcpy.ListFields(fc)]
    cursor = arcpy.da.SearchCursor(fc, field_names)
    conns_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    conns_df = conns_df.drop(['OBJECTID', 'Shape', 'date_start', 'freq', 'totalori', 'totquant', 'Shape_Length'], 1)
    conns_df = conns_df[conns_df.from_id != conns_df.to_id]
    conns_all = conns_all.append(conns_df)

if can_only:
    # access mpa attributes from master dataset
    field_names = [i.name for i in arcpy.ListFields(mpa_att)]
    cursor = arcpy.da.SearchCursor(mpa_att, field_names)
    mpa_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    # get list of uIDs for those in Canada
    mpa_df = mpa_df[mpa_df.provstate == 'British Columbia']
    mpa_uid = mpa_df.uID_20201124.to_list()
    # remove US ones from conns_all
    conns_all = conns_all[
        (conns_all.from_id.isin(mpa_uid)) & (conns_all.to_id.isin(mpa_uid)) 
        ]
    dist_df = dist_df[
        (dist_df.DestID.isin(mpa_uid)) & (dist_df.origin_id.isin(mpa_uid)) 
    ]


###################################



# I am considering it a successful connection if it was made in at least 1 direction.
# A connection might be made in both directions, but I'm not giving extra credit for that.

# If a connection was made in at least 1 direction, code the weaker one as '-1'
# If they are both zero then just code the first as '-1'.

# this is going to have to be done manually with nested for loops.
# I can't think of any other way to do it.

# A more proper way would have been...
# join to each other based on opposite of from/to, then I can compare values in
# the same row and edit their fields in a more pandas-like way.

# Has to be done by PLD, but you have to also include PLDs that are NaN in each selection

# Also, I have to do the join with distances AFTER I select by PLD. A brief description why:
# If I do the join before:
# A certain connection may exist one way for PLD 60, but it doesn't exist in the other direction.
# However, it does exist in the direction for other PLDs, so when I make a separate df
# I am not getting a NaN row for the direction that doesn't exist.
# This specific scenario isn't a problem. But, there are presumably connections
# that don't exist at all in either direction, but they exist in both directions for other PLDs,
# so I would lose that combination completely for a certain PLD.

arcpy.env.workspace = out_gdb
df_plot = pd.DataFrame(columns=['pld', 'threshold', 'perc_established'])

for pld in plds:

    # get connections just for that PLD
    df_pld_all = conns_all[(conns_all.pld == pld)]
    # join distances, keep all
    df_pld_all = df_pld_all.merge(dist_df, how='right', left_on=['from_id', 'to_id'], right_on=['origin_id', 'DestID'])
    # for any distance combinations that do not have a connection strength, fill as 0
    df_pld_all.prob_avg = df_pld_all.prob_avg.fillna(0)
    # remove rows with distances > 20km
    df_pld_sel = df_pld_all[df_pld_all.distkm < dist_thres]

    # since I am using PathCost, there are some instances where the connections differ slightly
    # e.g. 19.9 vs. 20.1, and therefore one is not selected
    # So now, go back and check that the opposite of each connection exists.
    # Add back in ones that do not.
    for row in df_pld_sel.itertuples():
        check = df_pld_sel[(df_pld_sel.DestID == row.origin_id) & (df_pld_sel.origin_id == row.DestID)]
        if  check.empty:
            df_add = df_pld_all[(df_pld_all.DestID == row.origin_id) & (df_pld_all.origin_id == row.DestID)]
            df_pld_sel = df_pld_sel.append(df_add, ignore_index=False)

    # if connection distance is less than 20km and the prob is over the threshold,
    # code as '1'
    # elif, lessn than 20km and not over, code as '0'
    for thresh,field in zip(prob_thresh, prob_thresh_field):
        df_pld_sel[field] = df_pld_sel.prob_avg.apply(lambda x: 1 if x > thresh else 0)

    # look up oppposite connection
    for field in prob_thresh_field:

        print(f"processing pld {str(pld)} and field {field}")

        # I check the index value to see if I have already looked at it in the iteration.
        # The issue is that the itertuples map object isn't getting edited while
        # I loop through it, so as I edit rows I haven't gotten to yet, they
        # don't change in the itertuples instance.
        # Ugh. This is getting messy and not very 'pandas'-vectorization-like.
        checked = []

        for row in df_pld_sel.itertuples():
            # look up the opposite connection and see what its values are for each thresh
            row_look = df_pld_sel[(df_pld_sel.DestID == row.origin_id) & (df_pld_sel.origin_id == row.DestID)]
            
            # Index values
            i_1 = row.Index
            if i_1 in checked:
                continue            
            i_2 = row_look.index[0]
            checked.append(i_2)

            for field in prob_thresh_field:

                # get 1/0 value for each threshold field
                t_1 = getattr(row, field)
                t_2 = row_look[field].values[0]

                # get prob_avg for each direction
                p_1 = row.prob_avg
                p_2 = row_look.prob_avg.values[0]

                if t_1==-1 or t_2==-1: # if one is already -1, skip
                    continue
                elif t_1==0 and t_2==0: # if neither is connected, just code the second as -1
                    df_pld_sel.at[i_2, field] = -1
                elif t_1==1 and t_2==0: # if 1 is connected, code the other as -1
                    df_pld_sel.at[i_2, field] = -1
                elif t_1==0 and t_2==1: # if 1 is connected, code the other as -1
                    df_pld_sel.at[i_1, field] = -1
                elif t_1==1 and t_2==1: # if both are connected, code the weaker one as -1
                    if p_1 >= p_2:
                        df_pld_sel.at[i_2, field] = -1
                    elif p_2 > p_1:
                        df_pld_sel.at[i_1, field] = -1    
    

    # for each threshold level, populate a new df
    # pld, threshold, count of 1s divided by count of 1s and 0s
    for field, prob in zip(prob_thresh_field, prob_thresh):
        counts = df_pld_sel[field].value_counts()
        perc_established = (counts[1] / (counts[0] + counts[1])) * 100
        df_plot = df_plot.append({'pld': pld, 'threshold': prob, 'perc_established': perc_established}, ignore_index=True)
    
    if out_to_fc:
        # output to fc.
        # I want to map the connections within 20km that were NOT established.    
        # 1 fc per combination of pld and threshold
        # All I need are the ones that are zero.
        # From the dist_df dataframe, get the OBJECTID of the to/from combinations that are zero.
        # Select by attribute and save to new fc
        # add fields for pld and threshold
        for field, prob in zip(prob_thresh_field, prob_thresh):
            df_sel_0 = df_pld_sel[df_pld_sel[field] == 0]
            oid_list = df_sel_0.OBJECTID.tolist()
            sel = arcpy.SelectLayerByAttribute_management(
                distances,
                'NEW_SELECTION',
                '"OBJECTID" IN {}'.format(tuple(oid_list))
            )
            fc_name = f'conns_not_est_{pld}_{field}'
            arcpy.CopyFeatures_management(sel, fc_name)
            arcpy.AddField_management(fc_name, 'pld', 'SHORT')
            arcpy.AddField_management(fc_name, 'threshold', 'FLOAT')

            with arcpy.da.UpdateCursor(fc_name, ['pld', 'threshold']) as cursor:
                for row in cursor:
                    row[0] = pld
                    row[1] = prob
                    cursor.updateRow(row)
   

###############################################


# line plot with threshold along the bottom and percent on y axis 1 line per PLD

df_plot['PD'] = df_plot.pld
df_plot['PD'] = df_plot.PD.astype('int')
df_plot['PD'] = df_plot.PD.astype('str') # seaborn hue needs to be a string

sns.set()
sns.set_style('white')
sns.set_context('paper')

f = sns.lineplot(
    data = df_plot,
    x = 'threshold',
    y = 'perc_established',
    hue = 'PD',
    #palette = 'husl'
)
f.set(xscale='log')
f.set(xlabel='Connectivity strength threshold', ylabel=r'% connections within 20km established')
f.figure.savefig(r'conn_established.svg')