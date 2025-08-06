'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 4: Mann-Whitney U tests between plant growth forms
(a) n-alkanoic acid log total concentration (ug/g plant)
(b) n-alkane acid log total concentration (ug/g plant)
(c) n-alkanoic acid Average Chain-Length (ACL; C20-C32)
(d) n-alkane acid Average Chain-Length (ACL; C21-C33)
(e) n-alkanoic acid d13C (C22-C28)
(f) n-alkane acid d13C (C23-C29)
(g) n-alkanoic acid apparent fractionation (eapp; C22-C28)
(h) n-alkane acid apparent fractionation (eapp; C23-C29)
'''

# See atpw_conda_env.yml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from copy import copy
import scipy.stats as stats
import atpw_functions as atpw_fun

# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 6
plt.rcParams['font.family'] = "Liberation Sans"


# Import plant wax data and map climate info onto it
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc)
arc_data = atpw_fun.map_oipc(arc_data)
arc_data = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).mean(numeric_only=True).reset_index()

# Calculate plant wax indices
ftconc = atpw_fun.tconc(arc_data, data_type="f", tconc_name="ftconc", log="yes")
atconc = atpw_fun.tconc(arc_data, data_type="a", tconc_name="atconc", log="yes")

facl = atpw_fun.wax_acl(arc_data, data_type="f", start_chain=20, end_chain=32, acl_name="facl")
aacl = atpw_fun.wax_acl(arc_data, data_type="a", start_chain=21, end_chain=33, acl_name="aacl")

fd13c = atpw_fun.iso_avg(arc_data, data_type="f", start_chain=22, end_chain=28, iso_type="d13c", iso_name="fd13c")
ad13c = atpw_fun.iso_avg(arc_data, data_type="a", start_chain=23, end_chain=29, iso_type="d13c", iso_name="ad13c")

feapp = atpw_fun.wax_precip_eapp_maf(arc_data, data_type="f", start_chain=22, end_chain=28, eapp_name="feapp")
aeapp = atpw_fun.wax_precip_eapp_maf(arc_data, data_type="a", start_chain=23, end_chain=29, eapp_name="aeapp")

arc_data = pd.concat(
    [
        arc_data,
        ftconc,
        atconc,
        facl,
        aacl,
        fd13c,
        ad13c,
        feapp,
        aeapp
    ], axis=1
)


# Set up MWU test output p-value arrays
fconc_ttestp = np.zeros(shape=(8,8))
aconc_ttestp = np.zeros(shape=(8,8))
facl_ttestp = np.zeros(shape=(8,8))
aacl_ttestp = np.zeros(shape=(8,8))
fd13c_ttestp = np.zeros(shape=(8,8))
ad13c_ttestp = np.zeros(shape=(8,8))
feapp_ttestp = np.zeros(shape=(8,8))
aeapp_ttestp = np.zeros(shape=(8,8))

# MWU tests: n-alkanoic acid log total concentration
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.log(np.array(arc_data[arc_data['cat'] == i].ftconc)),
            np.log(np.array(arc_data[arc_data['cat'] == j].ftconc)),
            # equal_var=False,
            nan_policy='omit'
        )
        fconc_ttestp[i-1,j-1] = ttest.pvalue

fconc_ttestp = pd.DataFrame(
    data=fconc_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU test: n-alkane log total concentration
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.log(np.array(arc_data[arc_data['cat'] == i].atconc)),
            np.log(np.array(arc_data[arc_data['cat'] == j].atconc)),
            # equal_var=False,
            nan_policy='omit'
        )
        aconc_ttestp[i-1,j-1] = ttest.pvalue

aconc_ttestp = pd.DataFrame(
    data=aconc_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU tests: n-alkanoic acid ACL (C20-C32)
for i in range(1,9):
  for j in range(1,9):
    ttest = stats.mannwhitneyu(
      np.array(arc_data[arc_data['cat'] == i].facl),
      np.array(arc_data[arc_data['cat'] == j].facl),
      # equal_var=False,
      nan_policy='omit'
    )
    facl_ttestp[i-1,j-1] = ttest.pvalue

facl_ttestp = pd.DataFrame(
  data=facl_ttestp,
  columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
  index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU tests: n-alkane ACL (C21-C33)
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.array(arc_data[arc_data['cat'] == i].aacl),
            np.array(arc_data[arc_data['cat'] == j].aacl),
            # equal_var=False,
            nan_policy='omit'
        )
        aacl_ttestp[i-1,j-1] = ttest.pvalue

aacl_ttestp = pd.DataFrame(
    data=aacl_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU tests: n-alkanoic acid d13C (C22-C28)
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.array(arc_data[arc_data['cat'] == i].fd13c),
            np.array(arc_data[arc_data['cat'] == j].fd13c),
            # equal_var=False,
            nan_policy='omit'
        )
        fd13c_ttestp[i-1,j-1] = ttest.pvalue

fd13c_ttestp = pd.DataFrame(
    data=fd13c_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU tests: n-alkane d13C (C23-C29)
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.array(arc_data[arc_data['cat'] == i].ad13c),
            np.array(arc_data[arc_data['cat'] == j].ad13c),
            # equal_var=False,
            nan_policy='omit'
        )
        ad13c_ttestp[i-1,j-1] = ttest.pvalue

ad13c_ttestp = pd.DataFrame(
    data=ad13c_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU tests: n-alkanoic acid eapp (C22-C28)
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.array(arc_data[arc_data['cat'] == i].feapp),
            np.array(arc_data[arc_data['cat'] == j].feapp),
            # equal_var=False,
            nan_policy='omit'
        )
        feapp_ttestp[i-1,j-1] = ttest.pvalue

feapp_ttestp = pd.DataFrame(
    data=feapp_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)

# MWU tests: n-alkane eapp (C23-C29)
for i in range(1,9):
    for j in range(1,9):
        ttest = stats.mannwhitneyu(
            np.array(arc_data[arc_data['cat'] == i].aeapp),
            np.array(arc_data[arc_data['cat'] == j].aeapp),
            # equal_var=False,
            nan_policy='omit'
        )
        aeapp_ttestp[i-1,j-1] = ttest.pvalue

aeapp_ttestp = pd.DataFrame(
    data=aeapp_ttestp,
    columns=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic'],
    index=['Tre','Shr','For','Fer','Gra','Mos','Liv','Lic']
)


# Figure script

growth_mask = np.triu(np.ones_like(facl_ttestp, dtype=bool))

fig, axs = plt.subplots(4,3,layout='constrained')

ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
    ax=ax, data=fconc_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask
)

ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
    ax=ax, data=aconc_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask
)

ax = axs[1,0]
sns.heatmap(
    ax=ax, data=facl_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[1,1]
sns.heatmap(
    ax=ax, data=aacl_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[2,0]
sns.heatmap(
    ax=ax, data=fd13c_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[2,1]
sns.heatmap(
    ax=ax, data=ad13c_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[3,0]
sns.heatmap(
    ax=ax, data=feapp_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

ax = axs[3,1]
sns.heatmap(
    ax=ax, data=aeapp_ttestp,
    cmap=["green","lightgrey"], center=0.05, cbar=False, mask=growth_mask,
)

fig.delaxes(axs[0,2])
fig.delaxes(axs[1,2])
fig.delaxes(axs[2,2])
fig.delaxes(axs[3,2])

figure_ttest = plt.gcf()
# figure_ttest.savefig('figures/atpw_figure4_mwup.svg')
