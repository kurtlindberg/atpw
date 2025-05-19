## Arctic Terrestrial Plant Waxes Figure X
# Student's T-tests (or Mann-Whitney U test) between plant growth forms

# Ecological and environmental controls on modern plant wax production and stable isotope fractionation alogn a latitudinal transect of the Eastern Canadian Arctic

# Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds, Helga Bultmann, Jonathan H. Raberg

# DOI: pending

# Author: Kurt R. Lindberg
# Last edited: 02/25/2025

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from copy import copy
import scipy.stats as stats
import atpw_functions as atpw_fun

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 6
plt.rcParams['font.family'] = "Liberation Sans"

# Import master plant wax datafile and map climate info onto it
arc = pd.read_excel(
  'krl_arctic_terrestrial_plantwax_20241022.xlsx',
  sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc)
arc_data = atpw_fun.map_oipc(arc_data)
arc_data = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).mean(numeric_only=True).reset_index()

temp_maf = atpw_fun.maf(arc_data, data_type="temp", maf_name="T")
precip_maf = atpw_fun.maf(arc_data, data_type="precip", maf_name="P")
rh_maf = atpw_fun.maf(arc_data, data_type="rh", maf_name="RH")
pd2h_maf = atpw_fun.maf(arc_data, data_type="pd2h", maf_name="pd2H")

ftconc = atpw_fun.tconc(arc_data, wax_type="f", tconc_name="ftconc", log="yes")
atconc = atpw_fun.tconc(arc_data, wax_type="a", tconc_name="atconc", log="yes")

facl = atpw_fun.wax_acl(arc_data, data_type="f", start_chain=20, end_chain=32, acl_name="facl")
aacl = atpw_fun.wax_acl(arc_data, data_type="a", start_chain=21, end_chain=33, acl_name="aacl")

fd13c = atpw_fun.iso_avg(arc_data, data_type="f", start_chain=22, end_chain=28, iso_type="d13c", iso_name="fd13c")
ad13c = atpw_fun.iso_avg(arc_data, data_type="a", start_chain=23, end_chain=29, iso_type="d13c", iso_name="ad13c")

feapp = atpw_fun.wax_precip_eapp_maf(arc_data, data_type="f", start_chain=22, end_chain=28, eapp_name="feapp")
aeapp = atpw_fun.wax_precip_eapp_maf(arc_data, data_type="a", start_chain=23, end_chain=29, eapp_name="aeapp")

arc_data = pd.concat(
  [
    arc_data,
    temp_maf,
    precip_maf,
    rh_maf,
    pd2h_maf,
    ftconc,
    atconc,
    facl,
    aacl,
    fd13c,
    ad13c,
    feapp,
    aeapp
  ],
  axis=1
)
# arc_data.to_csv('arc_data_mapped.csv')

# arc_data_conc = arc_data[arc_data['conc_ugg_plant'] == 1].reset_index(drop=True)
# arc_fames_conc = arc_data_conc.filter(regex="fconc")
# arc_fames_tconc = arc_fames_conc.sum(axis=1)
# arc_fames_tconc = pd.Series(data=arc_fames_tconc[arc_fames_tconc > 0], name="fconc")
# arc_alks_conc = arc_data_conc.filter(regex="aconc")
# arc_alks_tconc = arc_alks_conc.sum(axis=1)
# arc_alks_tconc = pd.Series(data=arc_alks_tconc[arc_alks_tconc > 0], name="aconc")

# arc_data_conc = pd.concat(
#   [
#     arc_data_conc,
#     arc_fames_tconc,
#     arc_alks_tconc
#   ],
#   axis=1
# )

fconc_ttestp = np.zeros(shape=(8,8))
aconc_ttestp = np.zeros(shape=(8,8))
facl_ttestp = np.zeros(shape=(8,8))
aacl_ttestp = np.zeros(shape=(8,8))
fd13c_ttestp = np.zeros(shape=(8,8))
ad13c_ttestp = np.zeros(shape=(8,8))
feapp_ttestp = np.zeros(shape=(8,8))
aeapp_ttestp = np.zeros(shape=(8,8))

# fames concentration
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.log(np.array(arc_data[arc_data['cat'] == i].ftconc)),
      np.log(np.array(arc_data[arc_data['cat'] == j].ftconc)),
      # equal_var=False,
      nan_policy='omit'
    )
    fconc_ttestp[i-1,j-1] = ttest.pvalue
fconc_ttestp = pd.DataFrame(
  data=fconc_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes concentration
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.log(np.array(arc_data[arc_data['cat'] == i].atconc)),
      np.log(np.array(arc_data[arc_data['cat'] == j].atconc)),
      # equal_var=False,
      nan_policy='omit'
    )
    aconc_ttestp[i-1,j-1] = ttest.pvalue
aconc_ttestp = pd.DataFrame(
  data=aconc_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# fames ACL
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].facl),
      np.array(arc_data[arc_data['cat'] == j].facl),
      # equal_var=False,
      nan_policy='omit'
    )
    facl_ttestp[i-1,j-1] = ttest.pvalue
facl_ttestp = pd.DataFrame(
  data=facl_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes ACL
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].aacl),
      np.array(arc_data[arc_data['cat'] == j].aacl),
      # equal_var=False,
      nan_policy='omit'
    )
    aacl_ttestp[i-1,j-1] = ttest.pvalue
aacl_ttestp = pd.DataFrame(
  data=aacl_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# fames d13c
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].fd13c),
      np.array(arc_data[arc_data['cat'] == j].fd13c),
      # equal_var=False,
      nan_policy='omit'
    )
    fd13c_ttestp[i-1,j-1] = ttest.pvalue
fd13c_ttestp = pd.DataFrame(
  data=fd13c_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes d13c
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].ad13c),
      np.array(arc_data[arc_data['cat'] == j].ad13c),
      # equal_var=False,
      nan_policy='omit'
    )
    ad13c_ttestp[i-1,j-1] = ttest.pvalue
ad13c_ttestp = pd.DataFrame(
  data=ad13c_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# fames eapp
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].feapp),
      np.array(arc_data[arc_data['cat'] == j].feapp),
      # equal_var=False,
      nan_policy='omit'
    )
    feapp_ttestp[i-1,j-1] = ttest.pvalue
feapp_ttestp = pd.DataFrame(
  data=feapp_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes eapp
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].aeapp),
      np.array(arc_data[arc_data['cat'] == j].aeapp),
      # equal_var=False,
      nan_policy='omit'
    )
    aeapp_ttestp[i-1,j-1] = ttest.pvalue
aeapp_ttestp = pd.DataFrame(
  data=aeapp_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)


## Figure Script

growth_mask = np.triu(np.ones_like(facl_ttestp, dtype=bool))
# my_cmap = copy(plt.cm.Blues)
# my_cmap.set_over("white")
# my_cmap.set_under("white")

fig, axs = plt.subplots(4,3,layout='constrained')

ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
  ax=ax, data=fconc_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask
)

ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
  ax=ax, data=aconc_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask
)

ax = axs[1,0]
sns.heatmap(
  ax=ax, data=facl_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[1,1]
sns.heatmap(
  ax=ax, data=aacl_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[2,0]
sns.heatmap(
  ax=ax, data=fd13c_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[2,1]
sns.heatmap(
  ax=ax, data=ad13c_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[3,0]
sns.heatmap(
  ax=ax, data=feapp_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[3,1]
sns.heatmap(
  ax=ax, data=aeapp_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

fig.delaxes(axs[0,2])
fig.delaxes(axs[1,2])
fig.delaxes(axs[2,2])
fig.delaxes(axs[3,2])

figure_ttest = plt.gcf()
# figure_ttest.savefig('figures/atpw_figure4_mwup_log.svg')



## ECA Mann-Whitney U tests

eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index(drop=True)

fconc_ttestp = np.zeros(shape=(8,8))
aconc_ttestp = np.zeros(shape=(8,8))
facl_ttestp = np.zeros(shape=(8,8))
aacl_ttestp = np.zeros(shape=(8,8))
fd13c_ttestp = np.zeros(shape=(8,8))
ad13c_ttestp = np.zeros(shape=(8,8))
feapp_ttestp = np.zeros(shape=(8,8))
aeapp_ttestp = np.zeros(shape=(8,8))

# fames concentration
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.log(np.array(eca_data[eca_data['cat'] == i].ftconc)),
      np.log(np.array(eca_data[eca_data['cat'] == j].ftconc)),
      # equal_var=False,
      nan_policy='omit'
    )
    fconc_ttestp[i-1,j-1] = ttest.pvalue
fconc_ttestp = pd.DataFrame(
  data=fconc_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes concentration
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.log(np.array(eca_data[eca_data['cat'] == i].atconc)),
      np.log(np.array(eca_data[eca_data['cat'] == j].atconc)),
      # equal_var=False,
      nan_policy='omit'
    )
    aconc_ttestp[i-1,j-1] = ttest.pvalue
aconc_ttestp = pd.DataFrame(
  data=aconc_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# fames ACL
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(eca_data[eca_data['cat'] == i].facl),
      np.array(eca_data[eca_data['cat'] == j].facl),
      # equal_var=False,
      nan_policy='omit'
    )
    facl_ttestp[i-1,j-1] = ttest.pvalue
facl_ttestp = pd.DataFrame(
  data=facl_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes ACL
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(eca_data[eca_data['cat'] == i].aacl),
      np.array(eca_data[eca_data['cat'] == j].aacl),
      # equal_var=False,
      nan_policy='omit'
    )
    aacl_ttestp[i-1,j-1] = ttest.pvalue
aacl_ttestp = pd.DataFrame(
  data=aacl_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# fames d13c
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(eca_data[eca_data['cat'] == i].fd13c),
      np.array(eca_data[eca_data['cat'] == j].fd13c),
      # equal_var=False,
      nan_policy='omit'
    )
    fd13c_ttestp[i-1,j-1] = ttest.pvalue
fd13c_ttestp = pd.DataFrame(
  data=fd13c_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes d13c
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(eca_data[eca_data['cat'] == i].ad13c),
      np.array(eca_data[eca_data['cat'] == j].ad13c),
      # equal_var=False,
      nan_policy='omit'
    )
    ad13c_ttestp[i-1,j-1] = ttest.pvalue
ad13c_ttestp = pd.DataFrame(
  data=ad13c_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# fames eapp
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(eca_data[eca_data['cat'] == i].feapp),
      np.array(eca_data[eca_data['cat'] == j].feapp),
      # equal_var=False,
      nan_policy='omit'
    )
    feapp_ttestp[i-1,j-1] = ttest.pvalue
feapp_ttestp = pd.DataFrame(
  data=feapp_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)
# alkanes eapp
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(eca_data[eca_data['cat'] == i].aeapp),
      np.array(eca_data[eca_data['cat'] == j].aeapp),
      # equal_var=False,
      nan_policy='omit'
    )
    aeapp_ttestp[i-1,j-1] = ttest.pvalue
aeapp_ttestp = pd.DataFrame(
  data=aeapp_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)


## Figure Script

growth_mask = np.triu(np.ones_like(facl_ttestp, dtype=bool))
# my_cmap = copy(plt.cm.Blues)
# my_cmap.set_over("white")
# my_cmap.set_under("white")

fig, axs = plt.subplots(4,3,layout='constrained')

ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
  ax=ax, data=fconc_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask
)

ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
  ax=ax, data=aconc_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask
)

ax = axs[1,0]
sns.heatmap(
  ax=ax, data=facl_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[1,1]
sns.heatmap(
  ax=ax, data=aacl_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[2,0]
sns.heatmap(
  ax=ax, data=fd13c_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[2,1]
sns.heatmap(
  ax=ax, data=ad13c_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[3,0]
sns.heatmap(
  ax=ax, data=feapp_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[3,1]
sns.heatmap(
  ax=ax, data=aeapp_ttestp,
  cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

fig.delaxes(axs[2,1])
fig.delaxes(axs[3,1])
fig.delaxes(axs[0,2])
fig.delaxes(axs[1,2])
fig.delaxes(axs[2,2])
fig.delaxes(axs[3,2])

sfigure_ttest = plt.gcf()
# sfigure_ttest.savefig('figures/atpw_sfigure_ecamwup.svg')