# Identify Strongly Connected Components
# (https://en.wikipedia.org/wiki/Strongly_connected_component)
# This will assess the scale of individual networks in the system.



import arcpy
import pandas as pd
import networkx as nx
import numpy as np


input_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts\COMBINED.gdb'
root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis\distance_analysis_mapping'
output_gdb = 'components.gdb'
plds = [1,3,7,10,21,30,40,60]
thresholds = [0.0000001, 0.00001, 0.001]
thres_labels = ['e7', 'e5', 'e3']
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'
remove_usa = True
# all mpa attributes from master dataset
mpa_att = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'



if not arcpy.Exists(os.path.join(root, output_gdb)):
    arcpy.CreateFileGDB_management(root, output_gdb)

# copy to components gdb
arcpy.env.workspace = input_gdb
fcs = arcpy.ListFeatureClasses(wild_card='conn_avg_pld*')
for fc in fcs:
    arcpy.CopyFeatures_management(fc, os.path.join(root, output_gdb, fc))

# add field for conn_id
arcpy.env.workspace = os.path.join(root, output_gdb)
fcs = arcpy.ListFeatureClasses(wild_card='conn_avg_pld*')
for fc in fcs:
    arcpy.AddField_management(fc, 'conn_uid', 'SHORT')
    with arcpy.da.UpdateCursor(fc, ['OBJECTID', 'conn_uid']) as cursor:
        for row in cursor:
            row[1] = row[0]
            cursor.updateRow(row)

# get mpas to exclude as a list
field_names = [i.name for i in arcpy.ListFields(mpas_exclude)]
cursor = arcpy.da.SearchCursor(mpas_exclude, field_names)
mpaex = pd.DataFrame(data=[row for row in cursor], columns=field_names)
mpaex = mpaex[mpaex.exclude==1]
mpaex_list = mpaex.uID_20201124.to_list()

if remove_usa:
    # access mpa attributes from master dataset
    field_names = [i.name for i in arcpy.ListFields(mpa_att)]
    cursor = arcpy.da.SearchCursor(mpa_att, field_names)
    mpa_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    # get list of uIDs for those in Canada
    mpa_df = mpa_df[mpa_df.provstate == 'British Columbia']
    mpa_uid = mpa_df.uID_20201124.to_list()





for fc in fcs:

    pld = fc[12:]
    if int(pld) not in plds:
        continue

    # read in conns as pandas df
    field_names = [i.name for i in arcpy.ListFields(fc)]
    cursor = arcpy.da.SearchCursor(fc, field_names)
    conns_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    conns_df = conns_df.drop(['OBJECTID', 'Shape', 'date_start', 'freq', 'totalori', 'totquant', 'Shape_Length'], 1)
    conns_df = conns_df[conns_df.from_id != conns_df.to_id]

    # exclude MPAs in deep inlets
    conns_df = conns_df[
        ~(conns_df.from_id.isin(mpaex_list)) | (conns_df.to_id.isin(mpaex_list)) 
        ]

    # remove US ones from conns_all
    if remove_usa:
        conns_df = conns_df[
            (conns_df.from_id.isin(mpa_uid)) & (conns_df.to_id.isin(mpa_uid)) 
            ]

    for thresh, label in zip(thresholds, thres_labels):

        # threshold
        conn_thresh = conns_df[conns_df.prob_avg >= thresh]

        # create graph
        G = nx.from_pandas_edgelist(
            conn_thresh, 
            source='from_id', 
            target='to_id', 
            edge_attr='prob_avg', 
            create_using=nx.DiGraph)

        # identify strongly connected components
        comps = nx.strongly_connected_components(G)

        # add field, give all rows -1 to start
        conns_df[f'component_{label}'] = -1

        # each component is listed as SET type, which are unordered, unchangeable and don't allow duplicates
        i = 1
        for comp in comps:
            c = list(comp)
            #print(c)

            if len(c)==1:
                continue

            # add a component id to connection dataframe
            conns_df[f'component_{label}'] = np.where(
                ((conns_df.from_id.isin(c)) & (conns_df.to_id.isin(c))),
                i,  # change it to i
                conns_df[f'component_{label}'] # or keep it the same
            )

            i+=1

    # output to gdb table
    df_out = conns_df.drop(['from_id', 'to_id', 'prob_avg', 'pld'],1)
    x = np.array(np.rec.fromrecords(df_out.values))
    names = df_out.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    arcpy.da.NumPyArrayToTable(x, os.path.join(root, output_gdb, 'temp_tbl')) # for some reason, you always have to give it the full path

    # join to connections fc
    conn_join = arcpy.AddJoin_management(fc, 'conn_uid', 'temp_tbl', 'conn_uid')
    arcpy.CopyFeatures_management(conn_join, fc+'_comps')

    # delete temp fc and table
    arcpy.Delete_management(fc)
    arcpy.Delete_management('temp_tbl')

