## Arctic Terrestrial Plant Waxes Supplemental Figure X
# Plant wax Carbon Preference Index (CPI) "scatterbox"

# Ecological and environmental controls on modern plant wax production and stable isotope fractionation alogn a latitudinal transect of the Eastern Canadian Arctic

# Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds, Helga Bultmann, Jonathan H. Raberg

# DOI: pending

# Author: Kurt R. Lindberg
# Last edited: 02/04/2025

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
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
# arc_data = arc_data[arc_data['growth_form'] == 'lichen']
eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index()

arc_fames_cpi = atpw_fun.wax_cpi(
  arc_data,
  data_type="f",
  start_chain=20,
  end_chain=32,
  cpi_name="fcpi_c20c32"
)
arc_fames_cpi.replace([np.inf, -np.inf], np.nan, inplace=True)
eca_fames_cpi = atpw_fun.wax_cpi(
  eca_data,
  data_type="f",
  start_chain=20,
  end_chain=32,
  cpi_name="fcpi_c20c32"
)
eca_fames_cpi.replace([np.inf, -np.inf], np.nan, inplace=True)

arc_alks_cpi = atpw_fun.wax_cpi(
  arc_data,
  data_type="a",
  start_chain=21,
  end_chain=33,
  cpi_name="acpi_c21c33"
)
arc_alks_cpi.replace([np.inf, -np.inf], np.nan, inplace=True)
eca_alks_cpi = atpw_fun.wax_cpi(
  eca_data,
  data_type="a",
  start_chain=21,
  end_chain=33,
  cpi_name="acpi_c21c33"
)
eca_alks_cpi.replace([np.inf, -np.inf], np.nan, inplace=True)

# print(np.nanmean(eca_fames_cpi), np.nanstd(eca_fames_cpi))
# print(np.nanmean(eca_alks_cpi), np.nanstd(eca_alks_cpi))
# print(np.nanmean(arc_fames_cpi), np.nanstd(arc_fames_cpi))
# print(np.nanmean(arc_alks_cpi), np.nanstd(arc_alks_cpi))


## Create CPI results plot

fig, axs = plt.subplots(2,2,layout='constrained')

ax = axs[0,0]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_fames_cpi,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_fames_cpi,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.hlines(y=1, xmin=0, xmax=9, linestyle='-', color='black', linewidth=1, zorder=1)
ax.hlines(y=2, xmin=0, xmax=9, linestyle='--', color='black', linewidth=1, zorder=1)
ax.set_xlim([0,9])
ax.set_ylim([0.5,150])
ax.set_yscale('log')
ax.set_xticks(
  ticks=[1,2,3,4,5,6,7,8],
  labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
  rotation=45
)
# ax.set_yticks()
ax.set_xlabel("")
ax.set_ylabel("CPI (C20-C32)")
ax.legend('', frameon=False)

ax = axs[0,1]
sns.scatterplot(
  ax=ax, x=arc_data.cat, y=arc_alks_cpi,
  hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
  ax=ax, x=eca_data.cat, y=eca_alks_cpi,
  hue=eca_data.lat, s=20, marker='s', palette='Blues', edgecolors='black', zorder=10
)
ax.hlines(y=1, xmin=0, xmax=9, linestyle='-', color='black', linewidth=1, zorder=1)
ax.hlines(y=2, xmin=0, xmax=9, linestyle='--', color='black', linewidth=1, zorder=1)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([0,9])
ax.set_ylim([0.5,150])
ax.set_yscale('log')
ax.set_xticks(
  ticks=[1,2,3,4,5,6,7,8],
  labels=["Tree","Shrub","Forb","Fern","Grass","Moss","Liverwort","Lichen"],
  rotation=45
)
# ax.set_yticks()
ax.set_xlabel("")
ax.set_ylabel("CPI (C21-C33)")
ax.legend('', frameon=False)

norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='Blues', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[0,1], location='right', shrink=0.8)

fig.delaxes(axs[1,0])
fig.delaxes(axs[1,1])

sfigure_cpi = plt.gcf()
sfigure_cpi.savefig('figures/atpw_sfigure_cpi.svg')