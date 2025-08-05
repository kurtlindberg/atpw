'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure S1: Plant wax CPI results
(a) n-alkanoic acid CPI (C20-C32)
(b) n-alkane CPI (C21-C23)
'''

# See atpw_conda_env.yml
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import atpw_functions as atpw_fun

# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Liberation Sans'

# Import plant wax data and map climate info onto it
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc, use_sample_year="yes")
arc_data = atpw_fun.map_oipc(arc_data)

arc_fames_cpi = atpw_fun.wax_cpi(arc_data, data_type="f", start_chain=20, end_chain=32, cpi_name="fcpi_c20c32")
arc_fames_cpi.replace([np.inf, -np.inf], np.nan, inplace=True)

arc_alks_cpi = atpw_fun.wax_cpi(arc_data, data_type="a", start_chain=21, end_chain=33, cpi_name="acpi_c21c33")
arc_alks_cpi.replace([np.inf, -np.inf], np.nan, inplace=True)

arc_data = pd.concat(
    [
        arc_data,
        arc_fames_cpi,
        arc_alks_cpi
    ], axis=1
)
eca_data = arc_data[arc_data['study'] == 'Lindberg_et_al'].reset_index(drop=True)
arc_data = arc_data[arc_data['study'] != 'Lindberg_et_al'].reset_index(drop=True)


# Figure script
fig, axs = plt.subplots(2,2,layout='constrained')

# n-alkanoic acid CPI (C20-C32)
ax = axs[0,0]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.fcpi_c20c32,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.fcpi_c20c32,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
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

# n-alkane CPI (C21-C33)
ax = axs[0,1]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.acpi_c21c33,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.acpi_c21c33,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
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

# Create color scale legend
norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='Blues', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[0,1], location='right', shrink=0.8)

fig.delaxes(axs[1,0])
fig.delaxes(axs[1,1])


figs1_cpi = plt.gcf()
# figs1_cpi.savefig('figures/atpw_figs1_cpi.svg')
