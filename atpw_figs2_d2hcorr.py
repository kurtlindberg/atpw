'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure S2: Plant wax CPI results
(a) n-alkanoic acid CPI (C20-C32)
(b) n-alkane CPI (C21-C23)
'''

# See atpw_conda_env.yml
import scipy.stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
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
# arc_data = arc_data[arc_data['study'] != 'Balascio_et_al_2018'].reset_index(drop=True)
arc_data = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).mean(numeric_only=True).reset_index()

# arc_stdev = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).std(ddof=0, numeric_only=True).reset_index()
# arc_stdev.to_csv('arc_stdev_test.csv')

# Calculate climate and plant wax indices
arc_pd2h_maf = atpw_fun.maf(arc_data, data_type="pd2h", maf_name='pd2h_maf')
arc_fames_d2h = atpw_fun.iso_avg(arc_data, data_type="f", iso_type="d2h", start_chain=22, end_chain=28, iso_name='fd2h_c22c28')
arc_alks_d2h = atpw_fun.iso_avg(arc_data, data_type="a", iso_type="d2h", start_chain=23, end_chain=29, iso_name='ad2h_c23c29')

arc_data=pd.concat(
    [
        arc_data,
        arc_pd2h_maf,
        arc_fames_d2h,
        arc_alks_d2h
    ], axis=1
)


# Pearson correlations with pan-Arctic data
fames_df = pd.concat([arc_data.pd2h_maf, arc_data.fd2h_c22c28], axis=1)
fames_df = fames_df[fames_df['fd2h_c22c28'].notna()].reset_index(drop=True)
alks_df = pd.concat([arc_data.pd2h_maf, arc_data.ad2h_c23c29], axis=1)
alks_df = alks_df[alks_df['ad2h_c23c29'].notna()].reset_index(drop=True)

psn_fames = np.array(scipy.stats.pearsonr(fames_df.pd2h_maf, fames_df.fd2h_c22c28))
psn_alks = np.array(scipy.stats.pearsonr(alks_df.pd2h_maf, alks_df.ad2h_c23c29))
print(psn_fames)
print(psn_alks)

model_fames = LinearRegression().fit(np.array(fames_df.pd2h_maf).reshape((-1,1)), np.array(fames_df.fd2h_c22c28))
model_alks = LinearRegression().fit(np.array(alks_df.pd2h_maf).reshape((-1,1)), np.array(alks_df.ad2h_c23c29))


# Pearson correlations with no Hollabattjonnen data
arc_data_nohol = arc_data[arc_data['site_name'] != "Hollabattjonnen"].reset_index(drop=True)
alks_df_nohol = pd.concat([arc_data_nohol.pd2h_maf, arc_data_nohol.ad2h_c23c29], axis=1)
alks_df_nohol = alks_df_nohol[alks_df_nohol['ad2h_c23c29'].notna()].reset_index(drop=True)
psn_alks_nohol = np.array(scipy.stats.pearsonr(alks_df_nohol.pd2h_maf, alks_df_nohol.ad2h_c23c29))
print(psn_alks_nohol)
model_alks_nohol = LinearRegression().fit(np.array(alks_df_nohol.pd2h_maf).reshape((-1,1)), np.array(alks_df_nohol.ad2h_c23c29))


# Isolate new ECA data and remove it from pan-Arctic dataset
eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN")].reset_index(drop=True)
arc_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|3LN") == False].reset_index(drop=True)


# Figure script
fig, axs = plt.subplots(2,2,layout='constrained')

# n-alkanoic acid d2H (C22-C28)
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.scatterplot(
    ax=ax, x=arc_data.pd2h_maf, y=arc_data.fd2h_c22c28,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.pd2h_maf, y=eca_data.fd2h_c22c28,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.plot(
    fames_df.pd2h_maf, (fames_df.pd2h_maf*model_fames.coef_ + model_fames.intercept_),
    color='black', zorder=10
)
ax.set_xlim([-145,-75])
ax.set_ylim([-310,-90])
ax.set_xlabel("MAF Precip. d2H")
ax.set_ylabel("d2H22-28")
ax.legend('', frameon=False)

# n-alkane d2H (C23-C29)
ax = axs[0,1]
ax.set_title("n-alkanes")
sns.scatterplot(
    ax=ax, x=arc_data.pd2h_maf, y=arc_data.ad2h_c23c29,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.pd2h_maf, y=eca_data.ad2h_c23c29,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.plot(
    alks_df.pd2h_maf, (alks_df.pd2h_maf*model_alks.coef_ + model_alks.intercept_),
    color='black', zorder=10
)
ax.plot(
    alks_df_nohol.pd2h_maf, (alks_df_nohol.pd2h_maf*model_alks_nohol.coef_ + model_alks_nohol.intercept_),
    color='grey', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([-145,-75])
ax.set_ylim([-310,-90])
ax.set_xlabel("MAF Precip. d2H")
ax.set_ylabel("d2H23-29")
ax.legend('', frameon=False)

# Create color scale legend
norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='Blues', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[0,1], location='right', shrink=0.8)

fig.delaxes(axs[1,0])
fig.delaxes(axs[1,1])


figs2_d2hcorr = plt.gcf()
# figs2_d2hcorr.savefig('figures/atpw_figs2_d2hcorr.svg')
