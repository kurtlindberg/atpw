'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 2
(a) mean temperature of the months above freezing
(b) mean total precipiation amount of the months above freezing
(c) mean relative humidity of the months above freezing
(d) site elevation
(e) mean amount-weighted precipitation isotope d2H
'''

# See atpw_conda_env.yml
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import scipy.stats
import atpw_functions as atpw_fun

# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = "Liberation Sans"

# Import data from main .xlsx file
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc, use_sample_year="yes")
arc_data = atpw_fun.map_oipc(arc_data)

eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index()
# eca_data = arc_data[arc_data['study'] == 'Hollister_et_al_2022'].reset_index(drop=True)

# Calculate temperature averages
eca_temp = eca_data.filter(regex="temp")
eca_temp['ann'] = atpw_fun.env_avg(eca_data, data_type="temp", env_avg_name="ann")
eca_temp['jja'] = atpw_fun.env_avg(eca_data, data_type="temp", months=["jun","jul","aug"], env_avg_name="jja")
eca_temp['maf'] = atpw_fun.maf(eca_data, data_type="temp", maf_name="maf")
arc_temp = arc_data.filter(regex="temp")
arc_temp['ann'] = atpw_fun.env_avg(arc_data, data_type="temp", env_avg_name="ann")
arc_temp['jja'] = atpw_fun.env_avg(arc_data, data_type="temp", months=["jun","jul","aug"], env_avg_name="jja")
arc_temp['maf'] = atpw_fun.maf(arc_data, data_type="temp", maf_name="maf")

# Calculate total precipitation
eca_precip = eca_data.filter(regex="precip")
eca_precip['ann'] = atpw_fun.env_avg(eca_data, data_type="precip", env_avg_name="ann")
eca_precip['jja'] = atpw_fun.env_avg(eca_data, data_type="precip", months=["jun","jul","aug"], env_avg_name="jja")
eca_precip['maf'] = atpw_fun.maf(eca_data, data_type="precip", maf_name="jja")
arc_precip = arc_data.filter(regex="precip")
arc_precip['ann'] = atpw_fun.env_avg(arc_data, data_type="precip", env_avg_name="ann")
arc_precip['jja'] = atpw_fun.env_avg(arc_data, data_type="precip", months=["jun","jul","aug"], env_avg_name="jja")
arc_precip['maf'] = atpw_fun.maf(arc_data, data_type="precip", maf_name="maf")

# Calculate relative humidity averages
eca_rh = eca_data.filter(regex="rh")
eca_rh['ann'] = atpw_fun.env_avg(eca_data, data_type="rh", env_avg_name="ann")
eca_rh['jja'] = atpw_fun.env_avg(eca_data, data_type="rh", months=["jun","jul","aug"], env_avg_name="jja")
eca_rh['maf'] = atpw_fun.maf(eca_data, data_type="rh", maf_name="maf")
arc_rh = arc_data.filter(regex="rh")
arc_rh['ann'] = atpw_fun.env_avg(arc_data, data_type="rh", env_avg_name="ann")
arc_rh['jja'] = atpw_fun.env_avg(arc_data, data_type="rh", months=["jun","jul","aug"], env_avg_name="jja")
arc_rh['maf'] = atpw_fun.maf(arc_data, data_type="rh", maf_name="maf")

# Calculate precipitation d2H averages
eca_oipc = eca_data.filter(regex="pd2h")
eca_oipc['ann'] = atpw_fun.env_avg(eca_data, data_type="pd2h", env_avg_name="ann")
eca_oipc['jja'] = atpw_fun.env_avg(eca_data, data_type="pd2h", months=["jun","jul","aug"], env_avg_name="jja")
eca_oipc['maf'] = atpw_fun.maf(eca_data, data_type="pd2h", maf_name="maf")
arc_oipc = arc_data.filter(regex="pd2h")
arc_oipc['ann'] = atpw_fun.env_avg(arc_data, data_type="pd2h", env_avg_name="ann")
arc_oipc['jja'] = atpw_fun.env_avg(arc_data, data_type="pd2h", months=["jun","jul","aug"], env_avg_name="jja")
arc_oipc['maf'] = atpw_fun.maf(arc_data, data_type="pd2h", maf_name="maf")


# Figure script
fig, axs = plt.subplots(3,2)

# MAF temperature
ax = axs[0,0]
sns.scatterplot(
    ax=ax, x=arc_data.lat, y=arc_temp.maf,
    marker='o', s=50, facecolors='#cb181d', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.lat, y=eca_temp.maf,
    marker='s', s=50, facecolors='#cb181d', edgecolors='black', zorder=10
)
ax.set_xlim([55,75])
ax.set_xticklabels("")
ax.set_ylim([3,13])
ax.set_yticks(
    ticks=[4,6,8,10,12],
    labels=[4,"",8,"",12]
)
ax.set_xlabel("")
ax.set_ylabel("MAF Temperature C")

# MAF total precipitation
ax = axs[0,1]
sns.scatterplot(
    ax=ax, x=arc_data.lat, y=arc_precip.maf*1000,
    marker='o', s=50, facecolors='#2171b5', zorder=10
)
sns.scatterplot(
    ax=ax, x=eca_data.lat, y=eca_precip.maf*1000,
    marker='s', s=50, facecolors='#2171b5', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([55,75])
ax.set_xticklabels("")
ax.set_ylim([100,900])
ax.set_yticks(
    ticks=[200,400,600,800],
    labels=[200,400,600,800]
)
ax.set_xlabel("")
ax.set_ylabel("MAF Total Precipitation (mm)")

# MAF relative humidity
ax = axs[1,0]
sns.scatterplot(
    ax=ax, x=arc_data.lat, y=arc_rh.maf,
    marker='o', s=50, facecolors='#6a51a3', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.lat, y=eca_rh.maf,
    marker='s', s=50, facecolors='#6a51a3', edgecolors='black', zorder=10
)
ax.set_xlim([55,75])
ax.set_ylim([55,85])
ax.set_yticks(
    ticks=[55,65,75,85],
    labels=[55,65,75,85]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("MAF Relative Humidity (%)")

# Sample site elevation
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

# MAF precipitation d2H
ax = axs[2,0]
sns.scatterplot(
    ax=ax, x=arc_data.lat, y=arc_oipc.maf,
    marker='o', s=50, facecolors='#238b45', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.lat, y=eca_oipc.maf,
    marker='s', s=50, facecolors='#238b45', edgecolors='black', zorder=10
)
ax.set_xlim([55,75])
ax.set_ylim([-160,-65])
ax.set_yticks(
    ticks=[-150,-125,-100,-75],
    labels=[-150,-125,-100,-75]
)
ax.set_xlabel("Latitude")
ax.set_ylabel("MAF Precipitation d2H")

fig.delaxes(axs[2,1])

figure2_env = plt.gcf()
# figure2_env.savefig('figures/atpw_figure2_env.svg')
