'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure S3: Principal component analyses scree plots of plant wax chain-length
distributions
(a) ECA n-alkanoic acids (C20-C32)
(b) ECA n-alkanes (C21-C33)
(c) Pan-Arctic n-alkanoic acids (C20-C32)
(d) Pan-Arctic n-alkanes (C21-C33)
'''

# See atpw_conda_env.yml
import pandas as pd
import matplotlib.pyplot as plt
import atpw_functions as atpw_fun


# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 9
plt.rcParams['font.family'] = "Liberation Sans"


# Load data and perform PCA
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
# arc_data = atpw_fun.map_era5clim(arc)
arc_data = arc.fillna(0)
# arc_data = arc_data.sort_values(by=['cat'])

arc_fames, arc_fames_pca = atpw_fun.wax_pca(arc_data, data_type="f", groupby="species")
arc_fames = arc_fames.sort_values(by=['cat']).reset_index(drop=True)
# arc_fames.to_csv('arc_fames_pca.csv')
arc_alks, arc_alks_pca = atpw_fun.wax_pca(arc_data, data_type="a", groupby="species")
arc_alks = arc_alks.sort_values(by=['cat']).reset_index(drop=True)
# arc_alks.to_csv('arc_alks_pca.csv')

eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|QPT|3LN")].reset_index()
eca_fames, eca_fames_pca = atpw_fun.wax_pca(eca_data, data_type="f", groupby="species")
eca_fames = eca_fames.sort_values(by=['cat']).reset_index(drop=True)
# eca_fames.to_csv('eca_fames_pca.csv')
eca_alks, eca_alks_pca = atpw_fun.wax_pca(eca_data, data_type="a", groupby="species")
eca_alks = eca_alks.sort_values(by=['cat']).reset_index(drop=True)
# eca_alks.to_csv('eca_alks_pca.csv')


# Figure script
fig, axs = plt.subplots(2,2,layout='constrained')

# ECA n-alkanoic acids (C20-C32)
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
ax.bar(
    x=eca_fames_pca["pc_values"], height=eca_fames_pca["pca"].explained_variance_ratio_*100,
    color='#a6cee3', edgecolor='k'
)
ax.plot(
    eca_fames_pca["pc_values"], eca_fames_pca["pca"].explained_variance_ratio_*100,
    'o-', linewidth=1.5, color='k'
)
ax.set_ylim([0,60])
ax.set_xticks([1,2,3,4,5,6,7])
ax.set_xticklabels("")
ax.set_ylabel('% Variance Explained')

ax = axs[0,1]
ax.set_title("n-alkanes")
ax.bar(
    x=eca_alks_pca["pc_values"], height=eca_alks_pca["pca"].explained_variance_ratio_*100,
    color='#a6cee3', edgecolor='k'
)
ax.plot(
    eca_alks_pca["pc_values"], eca_alks_pca["pca"].explained_variance_ratio_*100,
    'o-', linewidth=1.5, color='k'
)
ax.set_ylim([0,60])
ax.set_xticks([1,2,3,4,5,6,7])
ax.set_xticklabels("")
ax.set_ylabel('% Variance Explained')
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")

ax = axs[1,0]
ax.bar(
    x=arc_fames_pca["pc_values"], height=arc_fames_pca["pca"].explained_variance_ratio_*100,
    color='#a6cee3', edgecolor='k'
)
ax.plot(
    arc_fames_pca["pc_values"], arc_fames_pca["pca"].explained_variance_ratio_*100,
    'o-', linewidth=1.5, color='k'
)
ax.set_ylim([0,60])
ax.set_xticks([1,2,3,4,5,6,7])
ax.set_xlabel('Principal Component #')
ax.set_ylabel('% Variance Explained')

ax = axs[1,1]
ax.bar(
    x=arc_alks_pca["pc_values"], height=arc_alks_pca["pca"].explained_variance_ratio_*100,
    color='#a6cee3', edgecolor='k'
)
ax.plot(
    arc_alks_pca["pc_values"], arc_alks_pca["pca"].explained_variance_ratio_*100,
    'o-', linewidth=1.5, color='k'
)
ax.set_ylim([0,60])
ax.set_xticks([1,2,3,4,5,6,7])
ax.set_xlabel('Principal Component #')
ax.set_ylabel('% Variance Explained')
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")


figs3_pcascree = plt.gcf()
# figs3_pcascree.savefig('figures/atpw_figs3_pcascree.svg')
