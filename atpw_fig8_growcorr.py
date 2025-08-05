'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 8: Pearson correlation matrices of shrub, graminoid, and moss plant wax
data vs. environmental parameters
(a) n-alkanoic acid shrubs
(b) n-alkane shrubs
(c) n-alkanoic acid graminoids
(d) n-alkane graminoids
(e) n-alkanoic acid mosses
(f) n-alkane mosses
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

shrub_data = arc_data[arc_data['growth_form'] == 'shrub'].reset_index(drop=True)
gram_data = arc_data[arc_data['growth_form'] == "graminoid"].reset_index(drop=True)
moss_data = arc_data[arc_data['growth_form'] == "moss"].reset_index(drop=True)

# Create Pearson correlation matrices
shrub_fames_corr_df, shrub_fames_corr, shrub_fames_corr_p = atpw_fun.envcorr(shrub_data, data_type="f")
shrub_alks_corr_df, shrub_alks_corr, shrub_alks_corr_p = atpw_fun.envcorr(shrub_data, data_type="a")

gram_fames_corr_df, gram_fames_corr, gram_fames_corr_p = atpw_fun.envcorr(gram_data, data_type="f")
gram_alks_corr_df, gram_alks_corr, gram_alks_corr_p = atpw_fun.envcorr(gram_data, data_type="a")

moss_fames_corr_df, moss_fames_corr, moss_fames_corr_p = atpw_fun.envcorr(moss_data, data_type="f")
moss_alks_corr_df, moss_alks_corr, moss_alks_corr_p = atpw_fun.envcorr(moss_data, data_type="a")


# Figure script
growth_mask = np.triu(np.ones_like(shrub_fames_corr, dtype=bool))

fig, axs = plt.subplots(3,2,layout='constrained')

# n-alkanoic acid shrubs
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
    ax=ax, data=shrub_fames_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane shrubs
ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
    ax=ax, data=shrub_alks_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkanoic acid graminoids
ax = axs[1,0]
sns.heatmap(
    ax=ax, data=gram_fames_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane graminoids
ax = axs[1,1]
sns.heatmap(
    ax=ax, data=gram_alks_corr,
    cmap='coolwarm', cbar=True, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkanoic acid mosses
ax = axs[2,0]
sns.heatmap(
    ax=ax, data=moss_fames_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane mosses
ax = axs[2,1]
sns.heatmap(
    ax=ax, data=moss_alks_corr,
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)


fig8_growcorr = plt.gcf()
# fig8_growcorr.savefig('figures/atpw_fig8_growcorr.svg')
