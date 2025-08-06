'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 9: Pearson correlation matrices of Salix sp. and Betula sp. plant wax
data vs. environmental parameters
(a) n-alkanoic acid Salix sp.
(b) n-alkane Salix sp.
(c) n-alkanoic acid Betula sp.
(d) n-alkane Betula sp.
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

salix_data = arc_data[arc_data['genus'] == "Salix"].reset_index(drop=True)
betula_data = arc_data[arc_data['genus'] == "Betula"].reset_index(drop=True)
betula_data = betula_data[betula_data['growth_form'] != "tree"].reset_index(drop=True)

# Create Pearson correlation matrices
salix_fames_corr_df, salix_fames_corr, salix_fames_corr_p = atpw_fun.envcorr(salix_data, data_type="f")
salix_alks_corr_df, salix_alks_corr, salix_alks_corr_p = atpw_fun.envcorr(salix_data, data_type="a")

betula_fames_corr_df, betula_fames_corr, betula_fames_corr_p = atpw_fun.envcorr(betula_data, data_type="f")
betula_alks_corr_df, betula_alks_corr, betula_alks_corr_p = atpw_fun.envcorr(betula_data, data_type="a")

# salix_fames_corr_p.to_csv('envcorr_out/salix_fames_corr_p.csv')
# salix_alks_corr_p.to_csv('envcorr_out/salix_alks_corr_p.csv')
# betula_fames_corr_p.to_csv('envcorr_out/betula_fames_corr_p.csv')
# betula_alks_corr_p.to_csv('envcorr_out/betula_alks_corr_p.csv')


# Figure script
growth_mask = np.triu(np.ones_like(salix_fames_corr, dtype=bool))

fig, axs = plt.subplots(3,2,layout='constrained')

# n-alkanoic acid Salix sp.
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
  ax=ax, data=salix_fames_corr,
  cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
  annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane Salix sp.
ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
  ax=ax, data=salix_alks_corr,
  cmap='coolwarm', cbar=True, vmin=-1.0, vmax=1.0, mask=growth_mask,
  annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkanoic acid Betula sp.
ax = axs[1,0]
sns.heatmap(
  ax=ax, data=betula_fames_corr,
  cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
  annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane Betula sp.
ax = axs[1,1]
sns.heatmap(
  ax=ax, data=betula_alks_corr,
  cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=growth_mask,
  annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

fig.delaxes(axs[2,0])
fig.delaxes(axs[2,1])


fig9_gencorr = plt.gcf()
fig9_gencorr.savefig('figures/atpw_fig9_gencorr.svg')
