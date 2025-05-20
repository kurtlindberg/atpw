## Arctic Terrestrial Plant Waxes Supplemental Figure X
# Plant wax vs. MAF precipitation d2H correlation

# Ecological and environmental controls on modern plant wax production and stable isotope fractionation alogn a latitudinal transect of the Eastern Canadian Arctic

# Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds, Helga Bultmann, Jonathan H. Raberg

# DOI: pending

# Author: Kurt R. Lindberg
# Last edited: 02/04/2025

import scipy.stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
import atpw_functions as atpw_fun

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Liberation Sans'

arc = pd.read_excel(
  'krl_arctic_terrestrial_plantwax_20241022.xlsx',
  sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc, use_sample_year="yes")
arc_data = atpw_fun.map_oipc(arc_data)
# arc_data = arc_data[arc_data['study'] != 'Balascio_et_al_2018'].reset_index(drop=True)
arc_data = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).mean(numeric_only=True).reset_index()

# arc_stdev = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).std(ddof=0, numeric_only=True).reset_index()
# arc_stdev.to_csv('arc_stdev_test.csv')

eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index()

arc_pd2h_maf = atpw_fun.maf(
  arc_data,
  data_type="pd2h",
  maf_name='pd2h'
)
eca_pd2h_maf = atpw_fun.maf(
  eca_data,
  data_type="pd2h",
  maf_name='pd2h'
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

fames_df = pd.concat([arc_pd2h_maf, arc_fames_d2h], axis=1)
fames_df = fames_df[fames_df['fd2h_c22c28'].notna()].reset_index(drop=True)
alks_df = pd.concat([arc_pd2h_maf, arc_alks_d2h], axis=1)
alks_df = alks_df[alks_df['ad2h_c23c29'].notna()].reset_index(drop=True)

psn_fames = np.array(scipy.stats.pearsonr(fames_df.pd2h, fames_df.fd2h_c22c28))
psn_alks = np.array(scipy.stats.pearsonr(alks_df.pd2h, alks_df.ad2h_c23c29))

print(psn_fames)
print(psn_alks)

model_fames = LinearRegression().fit(np.array(fames_df.pd2h).reshape((-1,1)), np.array(fames_df.fd2h_c22c28))
model_alks = LinearRegression().fit(np.array(alks_df.pd2h).reshape((-1,1)), np.array(alks_df.ad2h_c23c29))

## Create correlation figure
fig, axs = plt.subplots(2,2,layout='constrained')

ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.scatterplot(
  ax=ax, x=arc_pd2h_maf, y=arc_fames_d2h,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_pd2h_maf, y=eca_fames_d2h,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.plot(
  fames_df.pd2h, (fames_df.pd2h*model_fames.coef_ + model_fames.intercept_),
  color='black', zorder=1
)
ax.set_xlim([-145,-75])
ax.set_ylim([-310,-90])
ax.set_xlabel("MAF Precip. d2H")
ax.set_ylabel("d2H22-28")
ax.legend('', frameon=False)

ax = axs[0,1]
ax.set_title("n-alkanes")
sns.scatterplot(
  ax=ax, x=arc_pd2h_maf, y=arc_alks_d2h,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_pd2h_maf, y=eca_alks_d2h,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.plot(
  alks_df.pd2h, (alks_df.pd2h*model_alks.coef_ + model_alks.intercept_),
  color='black', zorder=1
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([-145,-75])
ax.set_ylim([-310,-90])
ax.set_xlabel("MAF Precip. d2H")
ax.set_ylabel("d2H23-29")
ax.legend('', frameon=False)

norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='Blues', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[0,1], location='right', shrink=0.8)

fig.delaxes(axs[1,0])
fig.delaxes(axs[1,1])


sfigure_d2h = plt.gcf()
# sfigure_d2h.savefig('figures/atpw_sfigure_d2hcorr.svg')