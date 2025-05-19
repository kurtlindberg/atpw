## Arctic Terrestrial Plant Waxes Figure X
# Scatter plots of environmental parameters vs. latitude

# Ecological and environmental controls on modern plant wax production and stable isotope fractionation alogn a latitudinal transect of the Eastern Canadian Arctic

# Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds, Helga Bultmann, Jonathan H. Raberg

# DOI: pending

# Author: Kurt R. Lindberg
# Last edited: 11/12/2024

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import scipy.stats
import atpw_functions as atpw_fun

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = "Liberation Sans"

## Define function that calculates p-values for each Pearson
# correlation in 'pandas.corr()'
def corr_pvalues(df):

  df = df._get_numeric_data()
  dfcols = pd.DataFrame(columns=df.columns)
  pvalues = dfcols.transpose().join(dfcols, how='outer')

  for r in df.columns:
    for c in df.columns:

      if c == r:
        df_corr = df[[r]].dropna()
      else:
        df_corr = df[[r,c]].dropna()
        
      pvalues[r][c] = scipy.stats.pearsonr(df_corr[r], df_corr[c])[1]

  pvalues = pd.DataFrame(data=pvalues, index=df.columns, columns=df.columns)

  return pvalues

arc = pd.read_excel(r'~/Documents/PACEMAP/CH2_Modern_Terrestrial_Plant_Waxes/figure_codedata/krl_arctic_terrestrial_plantwax_20241022.xlsx',
                     sheet_name = 'plantwax')
arc_data = atpw_fun.map_era5clim(arc, use_sample_year="yes")
arc_data = atpw_fun.map_oipc(arc_data)

# eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index()
eca_data = arc_data[arc_data['study'] == 'Hollister_et_al_2022'].reset_index(drop=True)

eca_temp = eca_data.filter(regex="temp")
eca_temp['ann'] = atpw_fun.env_avg(eca_data, data_type="temp", env_avg_name="ann")
eca_temp['jja'] = atpw_fun.env_avg(eca_data, data_type="temp", months=["jun","jul","aug"], env_avg_name="jja")
eca_temp['maf'] = atpw_fun.maf(eca_data, data_type="temp", maf_name="maf")
arc_temp = arc_data.filter(regex="temp")
arc_temp['ann'] = atpw_fun.env_avg(arc_data, data_type="temp", env_avg_name="ann")
arc_temp['jja'] = atpw_fun.env_avg(arc_data, data_type="temp", months=["jun","jul","aug"], env_avg_name="jja")
arc_temp['maf'] = atpw_fun.maf(arc_data, data_type="temp", maf_name="maf")

eca_precip = eca_data.filter(regex="precip")
eca_precip['ann'] = atpw_fun.env_avg(eca_data, data_type="precip", env_avg_name="ann")
eca_precip['jja'] = atpw_fun.env_avg(eca_data, data_type="precip", months=["jun","jul","aug"], env_avg_name="jja")
eca_precip['maf'] = atpw_fun.maf(eca_data, data_type="precip", maf_name="jja")
arc_precip = arc_data.filter(regex="precip")
arc_precip['ann'] = atpw_fun.env_avg(arc_data, data_type="precip", env_avg_name="ann")
arc_precip['jja'] = atpw_fun.env_avg(arc_data, data_type="precip", months=["jun","jul","aug"], env_avg_name="jja")
arc_precip['maf'] = atpw_fun.maf(arc_data, data_type="precip", maf_name="maf")

eca_rh = eca_data.filter(regex="rh")
eca_rh['ann'] = atpw_fun.env_avg(eca_data, data_type="rh", env_avg_name="ann")
eca_rh['jja'] = atpw_fun.env_avg(eca_data, data_type="rh", months=["jun","jul","aug"], env_avg_name="jja")
eca_rh['maf'] = atpw_fun.maf(eca_data, data_type="rh", maf_name="maf")
arc_rh = arc_data.filter(regex="rh")
arc_rh['ann'] = atpw_fun.env_avg(arc_data, data_type="rh", env_avg_name="ann")
arc_rh['jja'] = atpw_fun.env_avg(arc_data, data_type="rh", months=["jun","jul","aug"], env_avg_name="jja")
arc_rh['maf'] = atpw_fun.maf(arc_data, data_type="rh", maf_name="maf")

eca_oipc = eca_data.filter(regex="pd2h")
eca_oipc['ann'] = atpw_fun.env_avg(eca_data, data_type="pd2h", env_avg_name="ann")
eca_oipc['jja'] = atpw_fun.env_avg(eca_data, data_type="pd2h", months=["jun","jul","aug"], env_avg_name="jja")
eca_oipc['maf'] = atpw_fun.maf(eca_data, data_type="pd2h", maf_name="maf")
arc_oipc = arc_data.filter(regex="pd2h")
arc_oipc['ann'] = atpw_fun.env_avg(arc_data, data_type="pd2h", env_avg_name="ann")
arc_oipc['jja'] = atpw_fun.env_avg(arc_data, data_type="pd2h", months=["jun","jul","aug"], env_avg_name="jja")
arc_oipc['maf'] = atpw_fun.maf(arc_data, data_type="pd2h", maf_name="maf")

envcorr_df = pd.concat(
  [
    arc_data.lat,
    arc_temp.maf,
    arc_precip.maf,
    arc_rh.maf,
    arc_data.elevation,
    arc_oipc.maf
  ], axis=1
)
envcorr = envcorr_df.corr(min_periods=3, numeric_only=True)
corr_mask = np.triu(np.ones_like(envcorr, dtype=bool))
# envcorr_p = corr_pvalues(envcorr_df)
# envcorr_p.to_csv('envcorr_out/envcorr_p.csv')


## ECA environmental parameters vs. latitude

fig, axs = plt.subplots(3,2)

ax = axs[0,0]
sns.scatterplot(ax=ax, x=arc_data.lat, y=arc_temp.maf,
                marker='o', s=50, facecolors='#cb181d', zorder=5)
sns.scatterplot(ax=ax, x=eca_data.lat, y=eca_temp.maf,
                marker='s', s=50, facecolors='#cb181d', edgecolors='black', zorder=10)
ax.set_xlim([55, 75])
ax.set_xticklabels("")
ax.set_ylim([3, 13])
ax.set_yticks(
  ticks=[4,6,8,10,12],
  labels=[4,"",8,"",12]
)
ax.set_xlabel("")
ax.set_ylabel("Temperature C")

ax = axs[0,1]
sns.scatterplot(ax=ax, x=arc_data.lat, y=arc_precip.maf*1000,
                marker='o', s=50, facecolors='#2171b5', zorder=10)
sns.scatterplot(ax=ax, x=eca_data.lat, y=eca_precip.maf*1000,
                marker='s', s=50, facecolors='#2171b5', edgecolors='black', zorder=10)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([55, 75])
ax.set_xticklabels("")
ax.set_ylim([100, 900])
ax.set_yticks(
  ticks=[200,400,600,800],
  labels=[200,400,600,800]
)
ax.set_xlabel("")
ax.set_ylabel("Total Precipitation (mm)")

ax = axs[1,0]
sns.scatterplot(ax=ax, x=arc_data.lat, y=arc_rh.maf,
                marker='o', s=50, facecolors='#6a51a3', zorder=5)
sns.scatterplot(ax=ax, x=eca_data.lat, y=eca_rh.maf,
                marker='s', s=50, facecolors='#6a51a3', edgecolors='black', zorder=10)
ax.set_xlim([55, 75])
ax.set_ylim([55, 85])
ax.set_yticks(
  ticks=[55,65,75,85],
  labels=[55,65,75,85]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("Relative Humidity (%)")

ax = axs[1,1]
sns.scatterplot(
  ax=ax, x=arc_data.lat, y=arc_data.elevation,
  marker='o', s=50, facecolors="grey", zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.lat, y=eca_data.elevation,
  marker='s', s=50, facecolors='grey', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([55,75])
ax.set_ylim([-50,1050])
ax.set_yticks(
  ticks=[0,250,500,750,1000],
  labels=[0,250,500,750,1000]
)
ax.set_xlabel("Latitude")
ax.set_ylabel("Elevation (m)")

ax = axs[2,0]
sns.scatterplot(ax=ax, x=arc_data.lat, y=arc_oipc.maf,
                marker='o', s=50, facecolors='#238b45', zorder=5)
sns.scatterplot(ax=ax, x=eca_data.lat, y=eca_oipc.maf,
                marker='s', s=50, facecolors='#238b45', edgecolors='black', zorder=10)
ax.set_xlim([55, 75])
ax.set_ylim([-160, -65])
ax.set_yticks(
  ticks=[-150,-125,-100,-75],
  labels=[-150,-125,-100,-75]
)
ax.set_xlabel("Latitude")
ax.set_ylabel("Precipitation d2H")

fig.delaxes(axs[2,1])

figure2_env = plt.gcf()
# figure2_env.savefig('figures/atpw_figure2_env.svg')


## Pearson correlation matrix of environmental parameters

fig, axs = plt.subplots(1,1)

ax = axs
sns.heatmap(
  ax=ax, data=envcorr,
  cmap='coolwarm', cbar=True, vmin=-1.0, vmax=1.0, mask=corr_mask,
  annot=True, annot_kws={"fontsize":10}, fmt='.2f',
  xticklabels=['Lat','T','P','RH','E','d2H'],
  yticklabels=['Lat','T','P','RH','E','d2H']
)
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

sfigure_envcorr = plt.gcf()
# sfigure_envcorr.savefig('figures/atpw_sfigure_envcorr.svg')