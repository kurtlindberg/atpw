'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 7: Pearson correlation matrices of vascular and non-vascular plant wax
data vs. environmental parameters
(a) n-alkanoic acid vascular plants
(b) n-alkane vascular plants
(c) n-alkanoic acid non-vascular plants
(d) n-alkane non-vascular plants
'''


# See atpw_conda_env.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import atpw_functions as atpw_fun

# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = "Liberation Sans"


# Import master plant wax datafile and map climate info onto it
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc)
arc_data = atpw_fun.map_oipc(arc_data)
arc_data = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).mean(numeric_only=True).reset_index()

vas_data = arc_data[arc_data['growth_form'].str.contains("tree|shrub|forb|fern|graminoid")].reset_index(drop=True)
nvas_data = arc_data[arc_data['growth_form'].str.contains("moss|liverwort|lichen")].reset_index(drop=True)

# Create Pearson correlation matrices
vas_fames_corr_df, vas_fames_corr, vas_fames_corr_p = atpw_fun.envcorr(vas_data, data_type="f")
vas_alks_corr_df, vas_alks_corr, vas_alks_corr_p = atpw_fun.envcorr(vas_data, data_type="a")

nvas_fames_corr_df, nvas_fames_corr, nvas_fames_corr_p = atpw_fun.envcorr(nvas_data, data_type="f")
nvas_alks_corr_df, nvas_alks_corr, nvas_alks_corr_p = atpw_fun.envcorr(nvas_data, data_type="a")

# vas_fames_corr_p.to_csv('envcorr_out/vas_fames_corr_p.csv')
# vas_alks_corr_p.to_csv('envcorr_out/vas_alks_corr_p.csv')
# nvas_fames_corr_p.to_csv('envcorr_out/nvas_fames_corr_p.csv')
# nvas_alks_corr_p.to_csv('envcorr_out/nvas_alks_corr_p.csv')


# Figure script
growth_mask = np.triu(np.ones_like(vas_fames_corr, dtype=bool))

fig, axs = plt.subplots(3,2,layout='constrained')

# n-alkanoic acid vascular plants
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
    ax=ax, data=vas_fames_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane vascular plants
ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
    ax=ax, data=vas_alks_corr,
    cmap='coolwarm', cbar=True, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkanoic vascular plants
ax = axs[1,0]
sns.heatmap(
    ax=ax, data=nvas_fames_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane vascular plants
ax = axs[1,1]
sns.heatmap(
    ax=ax, data=nvas_alks_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

fig.delaxes(axs[2,0])
fig.delaxes(axs[2,1])


fig7_vascorr = plt.gcf()
# fig7_vascorr.savefig('figures/atpw_fig7_vascorr.svg')
