## Arctic Terrestrial Plant Waxes Figure X
# Plant wax chain-length + stable isotope results "scatterbox"

# Ecological and environmental controls on modern plant wax production and stable isotope fractionation alogn a latitudinal transect of the Eastern Canadian Arctic

# Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds, Helga Bultmann, Jonathan H. Raberg

# DOI: pending

# Author: Kurt R. Lindberg
# Last edited: 11/12/2024

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import atpw_functions as atpw_fun

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 9
plt.rcParams['font.family'] = "Liberation Sans"


## Define function for calculating Shapiro-Wilk tests on each plant growth form
# with a chosen data type
def sw_test(df, df_data, tconc="no"):

  df = pd.Series(data=df, name='data')
  gf = [
    'tree',
    'shrub',
    'forb',
    'fern',
    'graminoid',
    'moss',
    'liverwort',
    'lichen'
  ]
  sw = np.zeros(shape=len(gf))

  for i in range(1,9):

    df_test = pd.concat([df, df_data.cat], axis=1)
    df_cat = df_test[df_test['cat'] == i].reset_index(drop=True)

    match tconc:
      case "no":
        sw[i-1] = scipy.stats.shapiro(np.array(df_cat.data), nan_policy='omit').pvalue
      case "yes":
        sw[i-1] = scipy.stats.shapiro(np.log(np.array(df_cat.data)), nan_policy='omit').pvalue

    df_test = pd.DataFrame()

  sw = pd.Series(data=sw, index=gf)

  return sw


arc = pd.read_excel(
  'krl_arctic_terrestrial_plantwax_20241022.xlsx',
  sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc, use_sample_year="yes")
arc_data = atpw_fun.map_oipc(arc_data)
# arc_data = arc_data[arc_data['growth_form'] == 'shrub'].reset_index(drop=True)
# arc_data = arc_data[arc_data['study'] != 'Balascio_et_al_2018'].reset_index(drop=True)
eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index()
# eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|QPT|3LN")].reset_index(drop=True)

# arc_data_conc = arc_data[arc_data['conc_ugg_plant'] == 1]
# # arc_data_conc = arc_data_conc[arc_data_conc['growth_form'] == 'liverwort']
# arc_fames_conc = arc_data_conc.filter(regex="fconc")
# arc_fames_tconc = arc_fames_conc.sum(axis=1)
# arc_fames_tconc = arc_fames_tconc[arc_fames_tconc > 0]
# arc_alks_conc = arc_data_conc.filter(regex="aconc")
# arc_alks_tconc = arc_alks_conc.sum(axis=1)
# arc_alks_tconc = arc_alks_tconc[arc_alks_tconc > 0]

arc_fames_tconc = atpw_fun.tconc(arc_data, wax_type="f", tconc_name='arc_fames_tconc')
arc_alks_tconc = atpw_fun.tconc(arc_data, wax_type="a", tconc_name='arc_alks_tconc')

# eca_data_conc = eca_data[eca_data['conc_ugg_plant'] == 1]
# # eca_data_conc = eca_data_conc[eca_data_conc['growth_form'] == 'lichen']
# eca_fames_conc = eca_data_conc.filter(regex="fconc")
# eca_fames_tconc = eca_fames_conc.sum(axis=1)
# eca_fames_tconc = eca_fames_tconc[eca_fames_tconc > 0]
# eca_alks_conc = eca_data_conc.filter(regex="aconc")
# eca_alks_tconc = eca_alks_conc.sum(axis=1)
# eca_alks_tconc = eca_alks_tconc[eca_alks_tconc > 0]

eca_fames_tconc = atpw_fun.tconc(eca_data, wax_type="f", tconc_name='eca_fames_tconc')
eca_alks_tconc = atpw_fun.tconc(eca_data, wax_type="a", tconc_name='eca_alks_tconc')

arc_fames_acl = atpw_fun.wax_acl(
  arc_data,
  data_type="f",
  start_chain=20,
  end_chain=32,
  acl_name="facl_c20c32"
)
eca_fames_acl = atpw_fun.wax_acl(
  eca_data,
  data_type="f",
  start_chain=20,
  end_chain=32,
  acl_name="facl_c20c32"
)
arc_alks_acl = atpw_fun.wax_acl(
  arc_data,
  data_type="a",
  start_chain=21,
  end_chain=33,
  acl_name="aacl_c21c33"
)
eca_alks_acl = atpw_fun.wax_acl(
  eca_data,
  data_type="a",
  start_chain=21,
  end_chain=33,
  acl_name="aacl_c21c33"
)

arc_fames_eapp_maf = atpw_fun.wax_precip_eapp_maf(
  arc_data,
  data_type="f",
  start_chain=22,
  end_chain=28,
  eapp_name='feapp_c22c28_maf'
)
eca_fames_eapp_maf = atpw_fun.wax_precip_eapp_maf(
  eca_data,
  data_type="f",
  start_chain=22,
  end_chain=28,
  eapp_name="feapp_c22c28_maf"
)
arc_alks_eapp_maf = atpw_fun.wax_precip_eapp_maf(
  arc_data,
  data_type="a",
  start_chain=23,
  end_chain=29,
  eapp_name='aeapp_c23c29_maf'
)
eca_alks_eapp_maf = atpw_fun.wax_precip_eapp_maf(
  eca_data,
  data_type="a",
  start_chain=23,
  end_chain=29,
  eapp_name="aeapp_c23c29_maf"
)

arc_fames_d13c = atpw_fun.iso_avg(
  arc_data,
  data_type="f",
  iso_type="d13c",
  start_chain=22,
  end_chain=28,
  iso_name='fd13c_c22c28'
)
eca_fames_d13c = atpw_fun.iso_avg(
  eca_data,
  data_type="f",
  iso_type="d13c",
  start_chain=22,
  end_chain=28,
  iso_name='fd13c_c22c28'
)
arc_alks_d13c = atpw_fun.iso_avg(
  arc_data,
  data_type="a",
  iso_type="d13c",
  start_chain=23,
  end_chain=29,
  iso_name='ad13c_c23c29'
)
eca_alks_d13c = atpw_fun.iso_avg(
  eca_data,
  data_type="a",
  iso_type="d13c",
  start_chain=23,
  end_chain=29,
  iso_name='ad13c_c23c29'
)

arc_fames_d2h = atpw_fun.iso_avg(
  arc_data,
  data_type="f",
  iso_type="d2h",
  start_chain=22,
  end_chain=28,
  iso_name='fd2h_c22c28'
)
eca_fames_d2h = atpw_fun.iso_avg(
  eca_data,
  data_type="f",
  iso_type="d2h",
  start_chain=22,
  end_chain=28,
  iso_name='fd2h_c22c28'
)
arc_alks_d2h = atpw_fun.iso_avg(
  arc_data,
  data_type="a",
  iso_type="d2h",
  start_chain=23,
  end_chain=29,
  iso_name='ad2h_c23c29'
)
eca_alks_d2h = atpw_fun.iso_avg(
  eca_data,
  data_type="a",
  iso_type="d2h",
  start_chain=23,
  end_chain=29,
  iso_name='ad2h_c23c29'
)

# eca_feapp = pd.concat(
#   [
#     eca_data,
#     eca_fames_eapp_maf
#   ],
#   axis=1
# ).to_csv('eca_feapp.csv')

arc_fames_tconc_sw = sw_test(arc_fames_tconc, arc_data, tconc="yes")
arc_alks_tconc_sw = sw_test(arc_alks_tconc, arc_data, tconc="yes")
arc_fames_acl_sw = sw_test(arc_fames_acl, arc_data)
arc_alks_acl_sw = sw_test(arc_alks_acl, arc_data)
arc_fames_d13c_sw = sw_test(arc_fames_d13c, arc_data)
arc_alks_d13c_sw = sw_test(arc_alks_d13c, arc_data)
arc_fames_eapp_maf_sw = sw_test(arc_fames_eapp_maf, arc_data)
arc_alks_eapp_maf_sw = sw_test(arc_alks_eapp_maf, arc_data)

eca_fames_tconc_sw = sw_test(eca_fames_tconc, eca_data, tconc="yes")
eca_alks_tconc_sw = sw_test(eca_alks_tconc, eca_data, tconc="yes")
eca_fames_acl_sw = sw_test(eca_fames_acl, eca_data)
eca_alks_acl_sw = sw_test(eca_alks_acl, eca_data)
eca_fames_d13c_sw = sw_test(eca_fames_d13c, eca_data)
eca_alks_d13c_sw = sw_test(eca_alks_d13c, eca_data)
eca_fames_eapp_maf_sw = sw_test(eca_fames_eapp_maf, eca_data)
eca_alks_eapp_maf_sw = sw_test(eca_alks_eapp_maf, eca_data)

# print(np.nanmean(eca_fames_eapp_maf))
# print(np.nanstd(eca_fames_eapp_maf))
# print(np.nanmean(eca_alks_eapp_maf))
# print(np.nanstd(eca_alks_eapp_maf))
# print(np.nanmean(arc_fames_eapp_maf))
# print(np.nanstd(arc_fames_eapp_maf))
# print(np.mean(arc_alks_eapp_maf))
# print(np.nanstd(arc_alks_eapp_maf))


## Create full results plot
fig, axs = plt.subplots(4,2,layout='constrained')

# Total chain-length concentration
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_fames_tconc,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_tconc,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10 
)
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_xticklabels("")
# ax.set_xticks(
#   ticks=[1,2,3,4,5,6,7,8],
#   labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
#   rotation=45
# )
ax.set_xlim([0,9])
ax.set_ylim([1, 50000])
ax.set_yscale('log')
ax.set_xlabel("")
ax.set_ylabel("Conc.")
ax.legend('', frameon=False)

ax = axs[0,1]
ax.set_title("n-alkanes")
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_alks_tconc,
  hue=arc_data.lat, s=25, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_tconc,
  hue=eca_data.lat, s=25, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_xticklabels("")
# ax.set_xticks(
#   ticks=[1,2,3,4,5,6,7,8],
#   labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
#   rotation=45
# )
ax.set_xlim([0,9])
ax.set_ylim([1, 50000])
ax.set_yscale('log')
ax.set_xlabel("")
ax.set_ylabel("Conc.")
ax.legend('', frameon=False)

# ACL
ax = axs[1,0]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_fames_acl,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_acl,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([20,32])
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[20,22,24,26,28,30,32],
  labels=[20,"",24,"",28,"",32]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("ACL")
ax.legend('', frameon=False)

ax = axs[1,1]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_alks_acl,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_acl,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([21,33])
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[21,23,25,27,29,31,33],
  labels=[21,"",25,"",29,"",33]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("ACL")
ax.legend('', frameon=False)

# d13C
ax = axs[2,0]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_fames_d13c,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_d13c,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[-45,-40,-35,-30,-25],
  labels=[-45,"",-35,"",-25]
)
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([-45,-25])
ax.set_xlabel("")
ax.set_ylabel("d13C")
ax.legend('', frameon=False)

ax = axs[2,1]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_alks_d13c,
  hue=arc_data.lat, s=25, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_d13c,
  hue=eca_data.lat, s=25, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[-45,-40,-35,-30,-25],
  labels=[-45,"",-35,"",-25]
)
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([-45,-25])
ax.set_xlabel("")
ax.set_ylabel("d13C")
ax.legend('', frameon=False)

# Eapp (wax-maf precip)
ax = axs[3,0]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_fames_eapp_maf,
  hue=arc_data.lat, s=25, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_eapp_maf,
  hue=eca_data.lat, s=25, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([-250,0])
ax.set_xticks(
  ticks=[1,2,3,4,5,6,7,8],
  labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
  rotation=45
)
# ax.set_yticks(
#   ticks=[-200,-150,-100,-50,0],
#   labels=[-200,"",-100,"",0]
# )
ax.set_xlabel("")
ax.set_ylabel("Eapp")
ax.legend('', frameon=False)

ax = axs[3,1]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_alks_eapp_maf,
  hue=arc_data.lat, s=25, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_eapp_maf,
  hue=eca_data.lat, s=25, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([0,9])
ax.set_ylim([-250,0])
# ax.set_yticks(
#   ticks=[-250,-200,-150,-100,-50,0,50],
#   labels=[-250,"",-150,"",-50,"",50]
# )
ax.set_xticks(
  ticks=[1,2,3,4,5,6,7,8],
  labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
  rotation=45
)
ax.set_xlabel("")
ax.set_ylabel("Eapp")
ax.legend('', frameon=False)

norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='Blues', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[:,1], location='right', shrink=0.6)


figure3_waxscatter = plt.gcf()
# figure3_waxscatter.savefig('figures/atpw_figure3_waxscatter.svg')


## Create ECA-only results plot

fig, axs = plt.subplots(4,2,layout='constrained')

# ACL
ax = axs[0,0]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_acl,
  hue=eca_data.lat, s=20, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([20,32])
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[20,22,24,26,28,30,32],
  labels=[20,"",24,"",28,"",32]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("ACL")
ax.legend('', frameon=False)

ax = axs[0,1]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_acl,
  hue=eca_data.lat, s=20, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([21,33])
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[21,23,25,27,29,31,33],
  labels=[21,"",25,"",29,"",33]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("ACL")
ax.legend('', frameon=False)

# Total chain-length concentration
ax = axs[1,0]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_tconc,
  hue=eca_data.lat, s=20, marker='s', palette='YlOrBr', edgecolors='black', zorder=10 
)
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_xticklabels("")
# ax.set_xticks(
#   ticks=[1,2,3,4,5,6,7,8],
#   labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
#   rotation=45
# )
ax.set_xlim([0,9])
ax.set_ylim([1, 50000])
ax.set_yscale('log')
ax.set_xlabel("")
ax.set_ylabel("Conc.")
ax.legend('', frameon=False)

ax = axs[1,1]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_tconc,
  hue=eca_data.lat, s=25, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_xticklabels("")
# ax.set_xticks(
#   ticks=[1,2,3,4,5,6,7,8],
#   labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
#   rotation=45
# )
ax.set_xlim([0,9])
ax.set_ylim([1, 50000])
ax.set_yscale('log')
ax.set_xlabel("")
ax.set_ylabel("Conc.")
ax.legend('', frameon=False)

# d13C
ax = axs[2,0]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_d13c,
  hue=eca_data.lat, s=20, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[-45,-40,-35,-30,-25],
  labels=[-45,"",-35,"",-25]
)
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([-45,-25])
ax.set_xlabel("")
ax.set_ylabel("d13C")
ax.legend('', frameon=False)

ax = axs[2,1]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_d13c,
  hue=eca_data.lat, s=25, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
  ticks=[-45,-40,-35,-30,-25],
  labels=[-45,"",-35,"",-25]
)
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([-45,-25])
ax.set_xlabel("")
ax.set_ylabel("d13C")
ax.legend('', frameon=False)

# Eapp (wax-maf precip)
ax = axs[3,0]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_eapp_maf,
  hue=eca_data.lat, s=25, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([-250,0])
ax.set_xticks(
  ticks=[1,2,3,4,5,6,7,8],
  labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
  rotation=45
)
ax.set_yticks(
  ticks=[-200,-150,-100,-50,0],
  labels=[-200,"",-100,"",0]
)
ax.set_xlabel("")
ax.set_ylabel("Eapp")
ax.legend('', frameon=False)

ax = axs[3,1]
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_eapp_maf,
  hue=eca_data.lat, s=25, marker='s', palette='YlOrBr', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([0,9])
ax.set_ylim([-250,0])
ax.set_xticks(
  ticks=[1,2,3,4,5,6,7,8],
  labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
  rotation=45
)
ax.set_xlabel("")
ax.set_ylabel("Eapp")
ax.legend('', frameon=False)

norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='YlOrBr', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[:,1], location='right', shrink=0.6)


figure3_waxscatter = plt.gcf()
# figure3_waxscatter.savefig('figures/atpw_figure3_waxscatter_eca.svg')