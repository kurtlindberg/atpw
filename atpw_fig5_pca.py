'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 5: Principal component analyses biplots of plant wax chain-length distributions
(a) ECA n-alkanoic acids (C20-C32)
(b) ECA n-alkanes (C21-C33)
(c) Pan-Arctic n-alkanoic acids (C20-C32)
(d) Pan-Arctic n-alkanes (C21-C33)
'''


# See atpw_conda_env.yml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns
import atpw_functions as atpw_fun


# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 9
plt.rcParams['font.family'] = "Liberation Sans"


# Create function for plotting PCA confidence ellipses
def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
    '''
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The Axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to '~matplotlib.patches.Ellipse'

    Returns
    -------
    matplotlib.patches.Ellipse
    '''

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0,1]/np.sqrt(cov[0,0] * cov[1,1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse(
        (0,0), width=ell_radius_x * 2, height=ell_radius_y * 2,
        facecolor=facecolor, **kwargs
    )

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0,0]) * n_std
    mean_x = np.mean(x)

    # Calculating the standard deviation of y
    scale_y = np.sqrt(cov[1,1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    return ax.add_patch(ellipse)


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
eca_site_colors = [
    '#4575b4',
    '#74add1',
    '#f46d43',
    '#d73027'
]

eca_growth_colors = [
    '#00441b',
    '#1b7837',
    '#5aae61',
    '#d9f0d3',
    '#c2a5cf',
    '#9970ab',
    '#762a83'
]

arc_growth_colors = [
    '#00441b',
    '#1b7837',
    '#5aae61',
    '#a6dba0',
    '#d9f0d3',
    '#c2a5cf',
    '#9970ab',
    '#762a83'
]

tree = mpatches.Patch(color=arc_growth_colors[0], label='Tree')
shrub = mpatches.Patch(color=arc_growth_colors[1], label='Shrub')
forb = mpatches.Patch(color=arc_growth_colors[2], label='Forb')
fern = mpatches.Patch(color=arc_growth_colors[3], label='Fern')
graminoid = mpatches.Patch(color=arc_growth_colors[4], label='Graminoid')
moss = mpatches.Patch(color=arc_growth_colors[5], label='Moss')
liverwort = mpatches.Patch(color=arc_growth_colors[6], label='Liverwort')
lichen = mpatches.Patch(color=arc_growth_colors[7], label='Lichen')

eca_fames_gf = eca_fames.growth_form.unique()
eca_alks_gf = eca_alks.growth_form.unique()
arc_fames_gf = arc_fames.growth_form.unique()
arc_alks_gf = arc_fames.growth_form.unique()
# print(arc_fames.cat)
# print(arc_fames_gf)

fig, axs = plt.subplots(2,2, layout='constrained')

# ECA n-alkanoic acids (C20-C32)
ax = axs[0,0]
for i, feature in enumerate(eca_fames_pca["features"]):
    ax.arrow(
        0, 0, eca_fames_pca["ldings"][0,i], eca_fames_pca["ldings"][1,i],
        head_width=0.02, head_length=0.02
    )
    ax.text(
        eca_fames_pca["ldings"][0,i] * 1.15, eca_fames_pca["ldings"][1,i] * 1.15,
        str(feature[0:3]), fontsize=10
    )
sns.scatterplot(
    ax=ax, x=eca_fames.pc1, y=eca_fames.pc2,
    hue=eca_fames.cat, style=eca_fames.site_name,
    style_order=['AFR','CF8','QPT','3LN'], markers=['^','X','P','v'],
    s=25, palette=eca_growth_colors, edgecolors='k', zorder=10
)
# for j in range(0, len(eca_fames_gf)):
#     ell_data = eca_fames[eca_fames['growth_form'] == eca_fames_gf[j]]
#     # print(eca_fames_gf[j])
#     # print(len(ell_data.growth_form))
#     if len(ell_data.growth_form) > 2:
#         confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=eca_growth_colors[j])
ax.legend('', frameon=False)
ax.axhline(y=0, color='k', linestyle='--', linewidth=0.75, zorder=5)
ax.axvline(x=0, color='k', linestyle='--', linewidth=0.75, zorder=5)
ax.set_xlabel('PC1 ' + str(np.round(eca_fames_pca["pca"].explained_variance_ratio_[0]*100, decimals=0)) + '%')
ax.set_ylabel('PC2 ' + str(np.round(eca_fames_pca["pca"].explained_variance_ratio_[1]*100, decimals=0)) + '%')
ax.set_xlim([-0.75, 0.75])
ax.set_ylim([-0.75, 0.75])
ax.set_xticks([-0.5, 0, 0.5])
ax.set_yticks([-0.5, 0, 0.5])
ax.xaxis.set_label_position("top")
ax.xaxis.set_ticks_position("top")
# ax.legend(loc='center left', bbox_to_anchor=(1,0.5))

# ECA n-alkanes (C21-C33)
ax = axs[0,1]
for i, feature in enumerate(eca_alks_pca["features"]):
    ax.arrow(
        0, 0, eca_alks_pca["ldings"][0,i], eca_alks_pca["ldings"][1,i],
        head_width=0.02, head_length=0.02
    )
    ax.text(
        eca_alks_pca["ldings"][0,i] * 1.15, eca_alks_pca["ldings"][1,i] * 1.15,
        str(feature[0:3]), fontsize=10
    )
sns.scatterplot(
    ax=ax, x=eca_alks.pc1, y=eca_alks.pc2,
    hue=eca_alks.cat, style=eca_alks.site_name,
    style_order=['AFR','CF8','QPT','3LN'], markers=['^','X','P','v'],
    s=25, palette=eca_growth_colors, edgecolors='k', zorder=10
)
# for j in range(0, len(eca_alks_gf)):
#     ell_data = eca_alks[eca_alks['growth_form'] == eca_alks_gf[j]]
#     # print(arc_alks_gf[j])
#     # print(len(ell_data.growth_form))
#     if len(ell_data.growth_form) > 2:
#         confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=eca_growth_colors[j])
ax.legend('', frameon=False)
ax.axhline(y=0, color='k', linestyle='--', linewidth=0.75, zorder=5)
ax.axvline(x=0, color='k', linestyle='--', linewidth=0.75, zorder=5)
ax.set_xlabel('PC1 ' + str(np.round(eca_alks_pca["pca"].explained_variance_ratio_[0]*100, decimals=0)) + '%')
ax.set_ylabel('PC2 ' + str(np.round(eca_alks_pca["pca"].explained_variance_ratio_[1]*100, decimals=0)) + '%')
ax.set_xlim([-0.75, 0.75])
ax.set_ylim([-0.75, 0.75])
ax.set_xticks([-0.5, 0, 0.5])
ax.set_yticks([-0.5, 0, 0.5])
ax.xaxis.set_label_position("top")
ax.xaxis.set_ticks_position("top")
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
# ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
ax.legend(
  handles=[tree,shrub,forb,fern,graminoid,moss,liverwort,lichen],
  loc='center left',
  bbox_to_anchor=(1,0.5)
)

# Pan-Arctic n-alkanoic acids (C20-C32)
ax = axs[1,0]
for i, feature in enumerate(arc_fames_pca["features"]):
    ax.arrow(
        0, 0, arc_fames_pca["ldings"][0,i], arc_fames_pca["ldings"][1,i],
        head_width=0.02, head_length=0.02
    )
    ax.text(
        arc_fames_pca["ldings"][0,i] * 1.15, arc_fames_pca["ldings"][1,i] * 1.15,
        str(feature[0:3]), fontsize=10
    )
sns.scatterplot(
    ax=ax, x=arc_fames.pc1, y=arc_fames.pc2,
    hue=arc_fames.cat,
    s=15, palette=arc_growth_colors, edgecolors='k', zorder=10)
# for j in range(0, len(arc_fames_gf)):
#     ell_data = arc_fames[arc_fames['growth_form'] == arc_fames_gf[j]]
#     # print(arc_fames_gf[j])
#     # print(len(ell_data.growth_form))
#     if len(ell_data.growth_form) > 2:
#         confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=arc_growth_colors[j])
ax.legend('', frameon=False)
ax.axhline(y=0, color='k', linewidth=0.75, linestyle='--')
ax.axvline(x=0, color='k', linewidth=0.75, linestyle='--')
ax.set_xlabel('PC1 ' + str(np.round(arc_fames_pca["pca"].explained_variance_ratio_[0]*100, decimals=0)) + '%')
ax.set_ylabel('PC2 ' + str(np.round(arc_fames_pca["pca"].explained_variance_ratio_[1]*100, decimals=0)) + '%')
ax.set_xlim([-0.75, 0.75])
ax.set_ylim([-0.75, 0.75])
ax.set_xticks([-0.5, 0, 0.5])
ax.set_yticks([-0.5, 0, 0.5])

# Pan-Arctic n-alkanes (C21-C33)
ax = axs[1,1]
for i, feature in enumerate(arc_alks_pca["features"]):
    ax.arrow(
        0, 0, arc_alks_pca["ldings"][0,i], arc_alks_pca["ldings"][1,i],
        head_width=0.02, head_length=0.02
    )
    ax.text(
        arc_alks_pca["ldings"][0,i] * 1.15, arc_alks_pca["ldings"][1,i] * 1.15,
        str(feature[0:3]), fontsize=10
    )
sns.scatterplot(
    ax=ax, x=arc_alks.pc1, y=arc_alks.pc2,
    hue=arc_alks.cat,
    s=15, palette=arc_growth_colors, edgecolors='k', zorder=10
)
# for j in range(0, len(arc_alks_gf)):
#     ell_data = arc_alks[arc_alks['growth_form'] == arc_alks_gf[j]]
#     # print(arc_alks_gf[j])
#     # print(len(ell_data.growth_form))
#     if len(ell_data.growth_form) > 2:
#         confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=arc_growth_colors[j])
ax.legend('', frameon=False)
ax.axhline(y=0, color='k', linewidth=0.75, linestyle='--')
ax.axvline(x=0, color='k', linewidth=0.75, linestyle='--')
ax.set_xlabel('PC1 ' + str(np.round(arc_alks_pca["pca"].explained_variance_ratio_[0]*100, decimals=0)) + '%')
ax.set_ylabel('PC2 ' + str(np.round(arc_alks_pca["pca"].explained_variance_ratio_[1]*100, decimals=0)) + '%')
ax.set_xlim([-0.75, 0.75])
ax.set_ylim([-0.75, 0.75])
ax.set_xticks([-0.5, 0, 0.5])
ax.set_yticks([-0.5, 0, 0.5])
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")


fig5_pcabi = plt.gcf()
# fig5_pcabi.savefig('figures/atpw_fig5_pcabi.svg')
