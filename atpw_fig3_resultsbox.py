'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 3: plant wax results from new ECA data and pan-Arctic compilation
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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import atpw_functions as atpw_fun

# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 9
plt.rcParams['font.family'] = "Liberation Sans"


# Define function for calculating Shapiro-Wilk tests on each plant growth form
# with a chosen data type
def sw_test(df, df_data, tconc="no"):

    df = pd.Series(data=df, name='data')
    gf = [
        'tree',
        'shrub',
        'forb',
        'fern',
        'graminoid',
        'moss',
        'liverwort',
        'lichen'
    ]
    sw = np.zeros(shape=len(gf))

    for i in range(1,9):
        df_test = pd.concat([df, df_data.cat], axis=1)
        df_cat = df_test[df_test['cat'] == i].reset_index(drop=True)

        match tconc:
            case "no":
                sw[i-1] = scipy.stats.shapiro(np.array(df_cat.data), nan_policy='omit').pvalue
            case "yes":
                sw[i-1] = scipy.stats.shapiro(np.log(np.array(df_cat.data)), nan_policy='omit').pvalue

        df_test = pd.DataFrame()

    sw = pd.Series(data=sw, index=gf)

    return sw


# Import plant wax data and map climate info onto it
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc, use_sample_year="yes")
arc_data = atpw_fun.map_oipc(arc_data)


# Calculate plant wax indices
arc_fames_tconc = atpw_fun.tconc(arc_data, data_type="f", tconc_name='ftconc_c20c32')
arc_alks_tconc = atpw_fun.tconc(arc_data, data_type="a", tconc_name='atconc_c21c33')

arc_fames_acl = atpw_fun.wax_acl(arc_data, data_type="f", start_chain=20, end_chain=32, acl_name="facl_c20c32")
arc_alks_acl = atpw_fun.wax_acl(arc_data, data_type="a", start_chain=21, end_chain=33, acl_name="aacl_c21c33")

arc_fames_d13c = atpw_fun.iso_avg(arc_data, data_type="f", iso_type="d13c", start_chain=22, end_chain=28, iso_name='fd13c_c22c28')
arc_alks_d13c = atpw_fun.iso_avg(arc_data, data_type="a", iso_type="d13c", start_chain=23, end_chain=29, iso_name='ad13c_c23c29')

arc_fames_eapp_maf = atpw_fun.wax_precip_eapp_maf(arc_data, data_type="f", start_chain=22, end_chain=28, eapp_name='feapp_c22c28_maf')
arc_alks_eapp_maf = atpw_fun.wax_precip_eapp_maf(arc_data, data_type="a", start_chain=23, end_chain=29, eapp_name='aeapp_c23c29_maf')

arc_fames_d2h = atpw_fun.iso_avg(arc_data, data_type="f", iso_type="d2h", start_chain=22, end_chain=28, iso_name='fd2h_c22c28')
arc_alks_d2h = atpw_fun.iso_avg(arc_data, data_type="a", iso_type="d2h", start_chain=23, end_chain=29, iso_name='ad2h_c23c29')

arc_data = pd.concat(
    [
        arc_data,
        arc_fames_tconc,
        arc_alks_tconc,
        arc_fames_acl,
        arc_alks_acl,
        arc_fames_d13c,
        arc_alks_d13c,
        arc_fames_d2h,
        arc_alks_d2h,
        arc_fames_eapp_maf,
        arc_alks_eapp_maf
    ], axis=1
)

# Shapiro-Wilk tests on plant growth forms
arc_fames_tconc_sw = sw_test(arc_data.ftconc_c20c32, arc_data, tconc="yes")
arc_alks_tconc_sw = sw_test(arc_data.atconc_c21c33, arc_data, tconc="yes")
arc_fames_acl_sw = sw_test(arc_data.facl_c20c32, arc_data)
arc_alks_acl_sw = sw_test(arc_data.aacl_c21c33, arc_data)
arc_fames_d13c_sw = sw_test(arc_data.fd13c_c22c28, arc_data)
arc_alks_d13c_sw = sw_test(arc_data.ad13c_c23c29, arc_data)
arc_fames_eapp_maf_sw = sw_test(arc_data.feapp_c22c28_maf, arc_data)
arc_alks_eapp_maf_sw = sw_test(arc_data.aeapp_c23c29_maf, arc_data)

# Isolate new ECA data and remove it from pan-Arctic dataset
eca_data = arc_data[arc_data['study'] == 'Lindberg_et_al'].reset_index(drop=True)
arc_data = arc_data[arc_data['study'] != 'Lindberg_et_al'].reset_index(drop=True)


# Figure script
fig, axs = plt.subplots(4,2,layout='constrained')

# n-alkanoic acid log total concentration
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.ftconc_c20c32,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.ftconc_c20c32,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([1,50000])
ax.set_yscale('log')
ax.set_xlabel("")
ax.set_ylabel("Conc.")
ax.legend('', frameon=False)

# n-alkane log total concentration
ax = axs[0,1]
ax.set_title("n-alkanes")
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.atconc_c21c33,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.atconc_c21c33,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([1, 50000])
ax.set_yscale('log')
ax.set_xlabel("")
ax.set_ylabel("Conc.")
ax.legend('', frameon=False)

# n-alkanoic acid ACL (C20-C32)
ax = axs[1,0]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.facl_c20c32,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.facl_c20c32,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([20,32])
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
    ticks=[20,22,24,26,28,30,32],
    labels=[20,"",24,"",28,"",32]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("ACL")
ax.legend('', frameon=False)

# n-alkane ACL (C21-C33)
ax = axs[1,1]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.aacl_c21c33,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.aacl_c21c33,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([21,33])
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
    ticks=[21,23,25,27,29,31,33],
    labels=[21,"",25,"",29,"",33]
)
ax.set_xticklabels("")
ax.set_xlabel("")
ax.set_ylabel("ACL")
ax.legend('', frameon=False)

# n-alkanoic acid d13C (C22-C28)
ax = axs[2,0]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.fd13c_c22c28,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.fd13c_c22c28,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
    ticks=[-45,-40,-35,-30,-25],
    labels=[-45,"",-35,"",-25]
)
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([-45,-25])
ax.set_xlabel("")
ax.set_ylabel("d13C")
ax.legend('', frameon=False)

# n-alkane d13C (C23-C29)
ax = axs[2,1]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.ad13c_c23c29,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.ad13c_c23c29,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xticks(ticks=[1,2,3,4,5,6,7,8])
ax.set_yticks(
    ticks=[-45,-40,-35,-30,-25],
    labels=[-45,"",-35,"",-25]
)
ax.set_xticklabels("")
ax.set_xlim([0,9])
ax.set_ylim([-45,-25])
ax.set_xlabel("")
ax.set_ylabel("d13C")
ax.legend('', frameon=False)

# n-alkanoic acid eapp (C22-C28)
ax = axs[3,0]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.feapp_c22c28_maf,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.feapp_c22c28_maf,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.set_xlim([0,9])
ax.set_ylim([-250,0])
ax.set_xticks(
    ticks=[1,2,3,4,5,6,7,8],
    labels=["Tree","Shrub","Forb","Fern","Graminoid","Moss","Liverwort","Lichen"],
    rotation=45
)
ax.set_yticks(
    ticks=[-250,-200,-150,-100,-50,0],
    labels=["",-200,"",-100,"",0]
)
ax.set_xlabel("")
ax.set_ylabel("Eapp")
ax.legend('', frameon=False)

# n-alkane eapp (C23-C29)
ax = axs[3,1]
sns.scatterplot(
    ax=ax, x=arc_data.cat+0.2, y=arc_data.aeapp_c23c29_maf,
    hue=arc_data.lat, s=20, marker='o', palette='Blues', edgecolors='black', zorder=5
)
sns.scatterplot(
    ax=ax, x=eca_data.cat-0.2, y=eca_data.aeapp_c23c29_maf,
    hue=eca_data.lat, s=15, marker='D', palette='Blues', edgecolors='black', zorder=10
)
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position("right")
ax.set_xlim([0,9])
ax.set_ylim([-250,0])
ax.set_yticks(
    ticks=[-250,-200,-150,-100,-50,0],
    labels=["",-200,"",-100,"",0]
)
ax.set_xticks(
    ticks=[1,2,3,4,5,6,7,8],
    labels=["Tree","Shrub","Forb","Fern","Graminoid","Moss","Liverwort","Lichen"],
    rotation=45
)
ax.set_xlabel("")
ax.set_ylabel("Eapp")
ax.legend('', frameon=False)

# Create color scale legend
norm = plt.Normalize(arc_data['lat'].min(), arc_data['lat'].max())
sm = plt.cm.ScalarMappable(cmap='Blues', norm=norm)
sm.set_array([])
ax.legend('', frameon=False)
ax.figure.colorbar(sm, ax=axs[:,1], location='right', shrink=0.6)


fig3_waxscatter = plt.gcf()
# fig3_waxscatter.savefig('figures/atpw_fig3_resultsbox.svg')
