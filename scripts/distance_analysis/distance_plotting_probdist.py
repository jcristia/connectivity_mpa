
# plot connection strength vs. distance (with distance being calculated in ArcGIS in a different script)
# for distance, use PATHCOST, not shape_length

import arcpy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
import math


# euclidean distances for all combinations of nodes
distances = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\distance_analysis\distance_analysis_mapping\Default.gdb\euc_lines_ALL"
# averaged connections by PLD
conns_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts\COMBINED.gdb'


# read in distances as pandas df
field_names = [i.name for i in arcpy.ListFields(distances)]
cursor = arcpy.da.SearchCursor(distances, field_names)
dist_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
# drop fields
dist_df = dist_df.drop(['Shape'], 1)
# drop self connections from the dataframe
dist_df = dist_df[dist_df.DestID != dist_df.origin_id]
dist_df = dist_df.drop(['OBJECTID', 'Shape_Length', 'length_difference'], 1)


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

# remove duplicates and take smallest PLD value
# sort by PLD, groupby connection, then take the first value, which should be the smallest PLD
conns_remove = conns_all.sort_values('pld').groupby(['from_id', 'to_id'], as_index=False).first()
# the issue here is...
# I'm just taking the probability of the smallest PLD. I don't think it makes
# sense to average, which means I am missing a lot of info in this plot.
# I'll still see what it looks like, but it probably makes sense to make a line
# for each PLD with all of the data.

# SO FOR HERE...
# I can read in each individual PLD averaged dataset
# Connections in this chapter are a bit different than the seagrass chapter. 
# There can be connections for pld1 that are NOT in 60 because things can keep 
# drifting even if they are over an MPA. They do not settle until the end of the
# pld. In the seagrass case, I only needed 1 point per connection because I had 
# averaged between all PLDs and then all time periods. I took the average 
# connection strength but the first pld that connection was made.
# BUT FOR HERE...
# I will try 2 plots:
# (1) conns_remove df - smallest PLD that a connection is made, no average
# (2) conns_all - all connections, fit a line for each PLD


# Join distances
# pandas merge, keep all records from distance dataframe. In conns_all, there
# are multiple from and to, so it needs to be on the left, but do a right join
# so that you retain all the records from the distance dataframe.
df_merge_all = conns_all.merge(dist_df, how='right', left_on=['from_id', 'to_id'], right_on=['origin_id', 'DestID'])
df_merge_rem = conns_remove.merge(dist_df, how='right', left_on=['from_id', 'to_id'], right_on=['origin_id', 'DestID'])

# for any distance combinations that do not have a connection strength, fill as 0
df_merge_all.prob_avg = df_merge_all.prob_avg.fillna(0)
df_merge_rem.prob_avg = df_merge_rem.prob_avg.fillna(0)
# convert to km
df_merge_all['distkm'] = (df_merge_all.PathCost)/1000.0
df_merge_rem['distkm'] = (df_merge_rem.PathCost)/1000.0

# create dfs without zeros
df_merge_all_noz = df_merge_all[df_merge_all.prob_avg > 0]
df_merge_rem_noz = df_merge_rem[df_merge_rem.prob_avg > 0]




############
# PLOTTING
############

# FITTING WITH MY OWN EQUATION
# Resource:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html

# # I first tried this, thinking that I would define an asymptote:
# def func(x, a, b, c):
#     return a * np.exp(-b * x) + c
# # in this equation:
# # a is the starting y value when x is zero
# # -b is how quickly the slope changes
# # c is the asymptote
# but I can't actually say that there is an asymptote present. That doesn't really make sense, especially if I log transform my data.
# it is more of a power curve (exponential curve)
# This is what Treml 2012 did. "Data fit to a negative exponential curve, with best fit of y = 0.664x**(-1.032), Rsquared = 0.22."
# However, I found that I can't get a good fit at all unless I first log transform my data. It's simply not possible without it.
# At first I did not want to do this because I wanted to compare when I do and do not include connections with a strength of zero (which can't be log transformed).
# However, for the sake of getting a good fit, I will drop that. So...

#########
# PLOT 1:
# df with one connection per combination of from/to
# fit one line

df_merge_rem_noz['probperc'] = df_merge_rem_noz.prob_avg * 100
df_merge_rem_noz['probperclog'] = np.log10(df_merge_rem_noz.probperc)
def func(x, a, b):
    return a * x**(b)
popt, pcov = curve_fit(func, df_merge_rem_noz.distkm, df_merge_rem_noz.probperclog)
# "Use non-linear least squares to fit a function, f, to data."
print(popt) # to see a,b
#print(pcov)

# get 95% confidence interval
# this ended up not being as straight forward as I thought it would be.
# This is the most straighforward approach and is an answer from 2020:
# https://prodevsblog.com/questions/118243/confidence-interval-for-exponential-curve-fit/
a, b = unc.correlated_values(popt, pcov)
px = np.linspace(0,600,300)
py = a * px**(b)
nom = unp.nominal_values(py)
std = unp.std_devs(py)

#### plot data (run this all together so that the line plots on top) ###
sns.set()
sns.set_style('white')
sns.set_context('paper')
f = sns.lmplot(
    x="distkm", 
    y="probperclog", 
    data=df_merge_rem_noz, 
    hue='pld',
    hue_order=[60,40,30,21,10,7,3,1], 
    scatter=True, 
    fit_reg=False, 
    scatter_kws={"s": 1, 'alpha':1},
    legend=True,
    legend_out=False,
    ) # plot points
plt.plot(px, nom, 'dimgray') # plot fitted curve
#plt.plot(px, nom - 2 * std) # if you want to plot just the bounding lines of the CI
#plt.plot(px, nom + 2 * std)
## or plot it as a fill:
#plt.fill_between(px, nom - 2 * std, nom + 2 * std, color='gray', alpha=0.2) # plot CI
# However, I won't plot the CI. it doesn't really show up, and is it relevant on a log scale?

# LEGEND...
# I need to use lmplot instead of regplot because it allows me to use hue, but
# it is a facetgrid type of plot (higher order?), so the access to the legend is
# different (there is not ax= attribute in lmplot like in regplot).
# I need to order my points so that 60 draws on bottom, which means putting it
# first in the list, and then the legend orders this way.
# There is no easy way to reorder the legend (there is if this was a regplot
# though). So, I need to manually create a legend. Oh well.
pal = sns.color_palette() # I'm just using the default
pal_hex = pal.as_hex()[:8]
pal_hex.reverse()
handles = []
labels = ['1', '3', '7', '10', '21', '30', '40', '60']
import matplotlib.lines as mlines
for h, l in zip(pal_hex, labels):
    blue_line = mlines.Line2D([], [], color=h, linestyle='None', marker='o', markersize=2, label=l)
    handles.append(blue_line)
plt.legend(title='PD (days)', frameon=False, handles=handles)

f.set(xlim=(0,600))
f.set(xlabel='Distance (km)', ylabel=r'log$_{10}$ Connection probability (%)')
f.savefig(r'probvdist_nodupes.svg')
###########

# get r squared
residuals = df_merge_rem_noz.probperclog- func(df_merge_rem_noz.distkm, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((df_merge_rem_noz.probperclog-np.mean(df_merge_rem_noz.probperclog))**2)
r_squared = 1 - (ss_res / ss_tot)
print(r_squared)

# EQUATION
# y = -0.73(x^0.26)
# r2 = 0.30


###########
# PLOT 2
# for reference, plot non transformed points:
sns.set()
sns.set_style('white')
sns.set_context('paper')
f = sns.lmplot(x="distkm", y="prob_avg", data=df_merge_rem_noz, scatter=True, scatter_kws={"s": 1, 'alpha':0.3}, fit_reg=False)
f.set(xlim=(0,600), ylim=(0,1))
f.set(xlabel='Distance (km)', ylabel=r'Connection Strength')
sns.despine(top=False, right=False, left=False, bottom=False)
#f.savefig(r'probvdist_nodupes_nontransformed.svg')




###########
# PLOT 3
# df with all connections
# fit a line for each PLD

df_merge_all_noz['probperc'] = df_merge_all_noz.prob_avg * 100
df_merge_all_noz['probperclog'] = np.log10(df_merge_all_noz.probperc)

# equations for each PLD line
def func(x, a, b):
    return a * x**(b)
plds = [60,40,30,21,10,7,3,1]
pxs = []
noms = []
for pld in plds:
    df_pld = df_merge_all_noz[df_merge_all_noz.pld == pld]
    popt, pcov = curve_fit(func, df_pld.distkm, df_pld.probperclog)
    a, b = unc.correlated_values(popt, pcov)
    print(popt)
    max_dist = df_pld.distkm.max()
    #max_dist = 50 * math.ceil(max_dist/50) # round up to nearest 50
    px = np.linspace(0,round(max_dist),200)
    py = a * px**(b)
    nom = unp.nominal_values(py)
    pxs.append(px)
    noms.append(nom)

sns.set()
sns.set_style('white')
sns.set_context('paper')
f = sns.lmplot(
    x="distkm", 
    y="probperclog", 
    data=df_merge_rem_noz, 
    hue='pld',
    hue_order=[60,40,30,21,10,7,3,1], 
    scatter=True,
    fit_reg=False, 
    scatter_kws={"s": 1, 'alpha':1},
    legend=True,
    legend_out=False,
    )

pal = sns.color_palette() # I'm just using the default
pal_hex = pal.as_hex()[:8]

for p,n,h in zip(pxs, noms, pal_hex):
    plt.plot(p, n, h) # plot fitted curve

# legend
pal_hex.reverse()
handles = []
labels = ['1', '3', '7', '10', '21', '30', '40', '60']
import matplotlib.lines as mlines
for h, l in zip(pal_hex, labels):
    blue_line = mlines.Line2D([], [], color=h, linestyle='None', marker='o', markersize=2, label=l)
    handles.append(blue_line)
plt.legend(title='PD (days)', frameon=False, handles=handles)

f.set(xlim=(0,600))
f.set(xlabel='Distance (km)', ylabel=r'log$_{10}$ Connection probability (%)')
f.savefig(r'probvdist_all.svg')