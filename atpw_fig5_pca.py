## Arctic Terrestrial Plant Waxes Figure X
# Principal Component Analysis (PCA) of plant wax chain-length data

# Ecological and environmental controls on modern plant wax production and stable isotope fractionation alogn a latitudinal transect of the Eastern Canadian Arctic

# Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds, Helga Bultmann, Jonathan H. Raberg

# DOI: pending

# Author: Kurt R. Lindberg
# Last edited: 11/11/2024

# import necessary packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from composition_stats import clr, closure, multiplicative_replacement
import atpw_functions as atpw_fun

# set graphical parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 9
plt.rcParams['font.family'] = "Liberation Sans"

# Create function for PCA analysis and storing results for plotting
def wax_pca(df, data_type, groupby="species"):

  df_data = df.copy()

  if data_type == "fconc":
    for i in range(0, len(df.c20_fconc)):
      # if np.sum([df['c20_fconc'].iloc[i], df['c22_fconc'].iloc[i], df['c24_fconc'].iloc[i], df['c26_fconc'].iloc[i], df['c28_fconc'].iloc[i], df['c30_fconc'].iloc[i], df['c32_fconc'].iloc[i]]) == 0:
      if np.sum([df_data.c20_fconc[i], df_data.c22_fconc[i], df_data.c24_fconc[i], df_data.c26_fconc[i], df_data.c28_fconc[i], df_data.c30_fconc[i], df_data.c32_fconc[i]]) == 0:
        # wax_df = df.drop([i])
        df_data.drop([i], inplace=True)
  elif data_type == "aconc":
    for i in range(0, len(df.c21_aconc)):
      # if np.sum([df['c21_aconc'].iloc[i], df['c23_aconc'].iloc[i], df['c25_aconc'].iloc[i], df['c27_aconc'].iloc[i], df['c29_aconc'].iloc[i], df['c31_aconc'].iloc[i], df['c33_aconc'].iloc[i]]) == 0:
      if np.sum([df_data.c21_aconc[i], df_data.c23_aconc[i], df_data.c25_aconc[i], df_data.c27_aconc[i], df_data.c29_aconc[i], df_data.c31_aconc[i], df_data.c33_aconc[i]]) == 0:
        # wax_df = df.drop([i])
        df_data.drop([i], inplace=True)

  wax_df = df_data.reset_index(drop=True)

  match groupby:
    case "species":
      wax_df = wax_df.groupby(
        [
        'genus',
        'species',
        'growth_form',
        'habitat',
        'site_name',
        'sample_year'
        ]
      ).mean(numeric_only=True).reset_index()

    case "genus":
      wax_df = wax_df.groupby(
        [
        'genus',
        'growth_form',
        'habitat',
        'site_name',
        'sample_year'
        ]
      ).mean(numeric_only=True).reset_index()

    case "none":
      wax_df = wax_df

  match data_type:
    case "fconc":
      wax_dat = wax_df[
        [
        'c20_fconc',
        'c22_fconc',
        'c24_fconc',
        'c26_fconc',
        'c28_fconc',
        'c30_fconc',
        'c32_fconc'
        ]
      ]
    case "aconc":
      wax_dat = wax_df[
        [
        'c21_aconc',
        'c23_aconc',
        'c25_aconc',
        'c27_aconc',
        'c29_aconc',
        'c31_aconc',
        'c33_aconc'
        ]
      ]

  wax_frac = closure(multiplicative_replacement(np.array(wax_dat)))
  wax_clr = clr(wax_frac)
  wax_data = pd.DataFrame(data=wax_clr, columns=wax_dat.columns)

  wax_scaler = StandardScaler()
  wax_scaler.fit(wax_data)
  wax_data_scaled = wax_scaler.transform(wax_data)

  wax_pca = PCA(n_components=7)
  wax_PC_scores = pd.DataFrame(wax_pca.fit_transform(wax_data_scaled),
                                columns=['PC1','PC2','PC3','PC4','PC5','PC6','PC7'])

  wax_loadings = pd.DataFrame(wax_pca.components_.T,
                                columns=['PC1','PC2','PC3','PC4','PC5','PC6','PC7'],
                                index=wax_data.columns)

  wax_pc1 = wax_pca.fit_transform(wax_data_scaled)[:,0]
  wax_pc2 = wax_pca.fit_transform(wax_data_scaled)[:,1]
  wax_ldings = wax_pca.components_

  wax_scale_pc1 = 1.0/(wax_pc1.max() - wax_pc1.min())
  wax_scale_pc2 = 1.0/(wax_pc2.max() - wax_pc2.min())
  wax_features = wax_data.columns

  wax_pc_values = np.arange(wax_pca.n_components_) + 1

  wax_pc1_scores = pd.DataFrame(data=(wax_pc1 * wax_scale_pc1), columns=['pc1'])
  wax_pc2_scores = pd.DataFrame(data=(wax_pc2 * wax_scale_pc2), columns=['pc2'])
  wax_df = pd.concat([wax_df, wax_pc1_scores, wax_pc2_scores], axis=1)

  pca_dict = {
    "pca": wax_pca,
    "pc_values": wax_pc_values,
    "features": wax_features,
    "ldings": wax_ldings,
    "pc1": wax_pc1,
    "scale_pc1": wax_scale_pc1,
    "pc2": wax_pc2,
    "scale_pc2": wax_scale_pc2,
  }

  return wax_df, pca_dict

# Create function for plotting PCA confidence ellipses
def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
  """
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
  """

  if x.size != y.size:
    raise ValueError("x and y must be the same size")

  cov = np.cov(x, y)
  pearson = cov[0,1]/np.sqrt(cov[0,0] * cov[1,1])
  # Using a special case to obtain the eigenvalues of this
  # two-dimensional dataset.
  ell_radius_x = np.sqrt(1 + pearson)
  ell_radius_y = np.sqrt(1 - pearson)
  ellipse = Ellipse((0,0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                    facecolor=facecolor, **kwargs)

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

## Load data and perform PCA
arc = pd.read_excel(r'~/Documents/PACEMAP/CH2_Modern_Terrestrial_Plant_Waxes/figure_codedata/krl_arctic_terrestrial_plantwax_20241022.xlsx',
                     sheet_name = 'plantwax')
# arc_data = atpw_fun.map_era5clim(arc)
arc_data = arc.fillna(0)
# arc_data = arc_data.sort_values(by=['cat'])

arc_fames, arc_fames_pca = wax_pca(arc_data, data_type="fconc", groupby="species")
arc_fames = arc_fames.sort_values(by=['cat']).reset_index(drop=True)
arc_fames.to_csv('arc_fames_pca.csv')
arc_alks, arc_alks_pca = wax_pca(arc_data, data_type="aconc", groupby="species")
arc_alks = arc_alks.sort_values(by=['cat']).reset_index(drop=True)
arc_alks.to_csv('arc_alks_pca.csv')

eca_data = arc_data[arc_data['site_name'].str.contains("AFR|CF8|QPT|3LN")].reset_index()
eca_fames, eca_fames_pca = wax_pca(eca_data, data_type="fconc", groupby="species")
eca_fames = eca_fames.sort_values(by=['cat']).reset_index(drop=True)
# eca_fames.to_csv('eca_fames_pca.csv')
eca_alks, eca_alks_pca = wax_pca(eca_data, data_type="aconc", groupby="species")
eca_alks = eca_alks.sort_values(by=['cat']).reset_index(drop=True)
# eca_alks.to_csv('eca_alks_pca.csv')


## Create PCA biplots

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
sns.scatterplot(ax=ax, x=eca_fames.pc1, y=eca_fames.pc2,
                hue=eca_fames.cat, style=eca_fames.site_name,
                style_order=['AFR','CF8','QPT','3LN'], markers=['^','X','P','v'],
                s=25, palette=eca_growth_colors, edgecolors='k', zorder=10)
# for j in range(0, len(eca_fames_gf)):
#   ell_data = eca_fames[eca_fames['growth_form'] == eca_fames_gf[j]]
#   # print(eca_fames_gf[j])
#   # print(len(ell_data.growth_form))
#   if len(ell_data.growth_form) > 2:
#     confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=eca_growth_colors[j])
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
sns.scatterplot(ax=ax, x=eca_alks.pc1, y=eca_alks.pc2,
                hue=eca_alks.cat, style=eca_alks.site_name,
                style_order=['AFR','CF8','QPT','3LN'], markers=['^','X','P','v'],
                s=25, palette=eca_growth_colors, edgecolors='k', zorder=10)
# for j in range(0, len(eca_alks_gf)):
#   ell_data = eca_alks[eca_alks['growth_form'] == eca_alks_gf[j]]
#   # print(arc_alks_gf[j])
#   # print(len(ell_data.growth_form))
#   if len(ell_data.growth_form) > 2:
#     confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=eca_growth_colors[j])
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
sns.scatterplot(ax=ax, x=arc_fames.pc1, y=arc_fames.pc2,
                hue=arc_fames.cat, 
                s=15, palette=arc_growth_colors, edgecolors='k', zorder=10)
# for j in range(0, len(arc_fames_gf)):
#   ell_data = arc_fames[arc_fames['growth_form'] == arc_fames_gf[j]]
#   # print(arc_fames_gf[j])
#   # print(len(ell_data.growth_form))
#   if len(ell_data.growth_form) > 2:
#     confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=arc_growth_colors[j])
ax.legend('', frameon=False)
ax.axhline(y=0, color='k', linewidth=0.75, linestyle='--')
ax.axvline(x=0, color='k', linewidth=0.75, linestyle='--')
ax.set_xlabel('PC1 ' + str(np.round(arc_fames_pca["pca"].explained_variance_ratio_[0]*100, decimals=0)) + '%')
ax.set_ylabel('PC2 ' + str(np.round(arc_fames_pca["pca"].explained_variance_ratio_[1]*100, decimals=0)) + '%')
ax.set_xlim([-0.75, 0.75])
ax.set_ylim([-0.75, 0.75])
ax.set_xticks([-0.5, 0, 0.5])
ax.set_yticks([-0.5, 0, 0.5])

ax = axs[1,1]
for i, feature in enumerate(arc_alks_pca["features"]):
  ax.arrow(0, 0, arc_alks_pca["ldings"][0,i], arc_alks_pca["ldings"][1,i],
            head_width=0.02, head_length=0.02)
  ax.text(arc_alks_pca["ldings"][0,i] * 1.15, arc_alks_pca["ldings"][1,i] * 1.15,
          str(feature[0:3]), fontsize=10)
sns.scatterplot(ax=ax, x=arc_alks.pc1, y=arc_alks.pc2,
                hue=arc_alks.cat,
                s=15, palette=arc_growth_colors, edgecolors='k', zorder=10)
# for j in range(0, len(arc_alks_gf)):
#   ell_data = arc_alks[arc_alks['growth_form'] == arc_alks_gf[j]]
#   # print(arc_alks_gf[j])
#   # print(len(ell_data.growth_form))
#   if len(ell_data.growth_form) > 2:
#     confidence_ellipse(ax=ax, x=ell_data.pc1, y=ell_data.pc2, n_std=2.0, edgecolor=arc_growth_colors[j])
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

figure4_pcabi = plt.gcf()
# figure4_pcabi.savefig('figures/atpw_figure5_pcabi.svg')


## Create PCA scree plots

fig, axs = plt.subplots(2,2,layout='constrained')

ax = axs[0,0]
ax.set_title("n-alkanoic acids")
ax.bar(x=eca_fames_pca["pc_values"], height=eca_fames_pca["pca"].explained_variance_ratio_*100, color='#a6cee3', edgecolor='k')
ax.plot(eca_fames_pca["pc_values"], eca_fames_pca["pca"].explained_variance_ratio_*100, 'o-', linewidth=1.5, color='k')
ax.set_ylim([0,60])
ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
ax.set_xticklabels("")
ax.set_ylabel('% Variance Explained')

ax = axs[0,1]
ax.set_title("n-alkanes")
ax.bar(x=eca_alks_pca["pc_values"], height=eca_alks_pca["pca"].explained_variance_ratio_*100, color='#a6cee3', edgecolor='k')
ax.plot(eca_alks_pca["pc_values"], eca_alks_pca["pca"].explained_variance_ratio_*100, 'o-', linewidth=1.5, color='k')
ax.set_ylim([0,60])
ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
ax.set_xticklabels("")
ax.set_ylabel('% Variance Explained')
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")

ax = axs[1,0]
ax.bar(x=arc_fames_pca["pc_values"], height=arc_fames_pca["pca"].explained_variance_ratio_*100, color='#a6cee3', edgecolor='k')
ax.plot(arc_fames_pca["pc_values"], arc_fames_pca["pca"].explained_variance_ratio_*100, 'o-', linewidth=1.5, color='k')
ax.set_ylim([0,60])
ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
ax.set_xlabel('Principal Component #')
ax.set_ylabel('% Variance Explained')

ax = axs[1,1]
ax.bar(x=arc_alks_pca["pc_values"], height=arc_alks_pca["pca"].explained_variance_ratio_*100, color='#a6cee3', edgecolor='k')
ax.plot(arc_alks_pca["pc_values"], arc_alks_pca["pca"].explained_variance_ratio_*100, 'o-', linewidth=1.5, color='k')
ax.set_ylim([0,60])
ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
ax.set_xlabel('Principal Component #')
ax.set_ylabel('% Variance Explained')
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")


figures1_scree = plt.gcf()
# figures1_scree.savefig('figures/atpw_sfigure_scree.svg')