'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Figure 6: Pearson correlation and "leave one out (loo)" matrices on full
pan-Arctic datasets vs. environmental parameters
(a) n-alkanoic acid Pearson correlation matrix
(b) n-alkane Pearson correlation matrix
(c) n-alkanoic acid loo max r-value change
(d) n-alkane loo max r-value change
(e) n-alkanoic acid loo max r-value growth form
(f) n-alkane loo max r-value growth form
'''


# See atpw_conda_env.yml
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import atpw_functions as atpw_fun

# Figure parameters for editing in Inkscape
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = "Liberation Sans"


# Define function that calculates correlation matrices for each
# 'leave one out' test
def envcorr_loo(df, wax_type, data_type):

    match wax_type:
        case "f":
            chain_start = 20
            chain_end = 32
            iso_start = 22
            iso_end = 28

        case "a":
            chain_start = 21
            chain_end = 33
            iso_start = 23
            iso_end = 29

    loo_test = [
        'all',
        'tree',
        'shrub',
        'forb',
        'fern',
        'graminoid',
        'moss',
        'liverwort',
        'lichen'
    ]

    loo_corr_dict = {}

    for test in range(0, len(loo_test)):

        df_test = df[df['cat'] != test]
        df_test = df_test.reset_index(drop=True)

        temp_maf = atpw_fun.maf(df_test, data_type="temp", maf_name='T')
        precip_maf = atpw_fun.maf(df_test, data_type="precip", maf_name='P')
        rh_maf = atpw_fun.maf(df_test, data_type="rh", maf_name='RH')
        # pd2h_maf = atpw_fun.maf(df_test, data_type="pd2h", maf_name='pd2H')

        # Build correlation matrix based on selected plant wax data type
        # fractional abundance = "frac"; d13c = "d13c"; apparent fractionation = "eapp"
        match data_type:
            case "frac":
                arc_frac = atpw_fun.frac_abd(df_test, data_type=wax_type, start_chain=chain_start, end_chain=chain_end)
                arc_acl = atpw_fun.wax_acl(df_test, data_type=wax_type, start_chain=chain_start, end_chain=chain_end, acl_name='ACL')

                corr_df = pd.concat(
                    [
                        temp_maf,
                        precip_maf,
                        rh_maf,
                        # pd2h_maf,
                        arc_frac.C20,
                        arc_frac.C22,
                        arc_frac.C24,
                        arc_frac.C26,
                        arc_frac.C28,
                        arc_frac.C30,
                        arc_frac.C32,
                        arc_acl
                    ], axis=1
                )

            case "d13c":
                c1_d13c = atpw_fun.iso_avg(df_test, data_type=wax_type, start_chain=iso_start, end_chain=iso_start, iso_type="d13c", iso_name='C'+str(iso_start))
                c2_d13c = atpw_fun.iso_avg(df_test, data_type=wax_type, start_chain=iso_start+2, end_chain=iso_start+2, iso_type="d13c", iso_name='C'+str(iso_start+2))
                c3_d13c = atpw_fun.iso_avg(df_test, data_type=wax_type, start_chain=iso_start+4, end_chain=iso_start+4, iso_type="d13c", iso_name='C'+str(iso_start+4))
                c4_d13c = atpw_fun.iso_avg(df_test, data_type=wax_type, start_chain=iso_end, end_chain=iso_end, iso_type="d13c", iso_name='C'+str(iso_end))
                c1c4_d13c = atpw_fun.iso_avg(df_test, data_type=wax_type, start_chain=iso_start, end_chain=iso_end, iso_type="d13c", iso_name='Avg')

                corr_df = pd.concat(
                    [
                        temp_maf,
                        precip_maf,
                        rh_maf,
                        # pd2h_maf,
                        c1_d13c,
                        c2_d13c,
                        c3_d13c,
                        c4_d13c,
                        c1c4_d13c
                    ], axis=1
                )

            case "eapp":
                c1_eapp_maf = atpw_fun.wax_precip_eapp_maf(df_test, data_type=wax_type, start_chain=iso_start, end_chain=iso_start, eapp_name='C'+str(iso_start))
                c2_eapp_maf = atpw_fun.wax_precip_eapp_maf(df_test, data_type=wax_type, start_chain=iso_start+2, end_chain=iso_start+2, eapp_name='C'+str(iso_start+2))
                c3_eapp_maf = atpw_fun.wax_precip_eapp_maf(df_test, data_type=wax_type, start_chain=iso_start+4, end_chain=iso_start+4, eapp_name='C'+str(iso_start+4))
                c4_eapp_maf = atpw_fun.wax_precip_eapp_maf(df_test, data_type=wax_type, start_chain=iso_end, end_chain=iso_end, eapp_name='C'+str(iso_end))
                c1c4_eapp_maf = atpw_fun.wax_precip_eapp_maf(df_test, data_type=wax_type, start_chain=iso_start, end_chain=iso_end, eapp_name='Avg')

                corr_df = pd.concat(
                    [
                        temp_maf,
                        precip_maf,
                        rh_maf,
                        # pd2h_maf,
                        c1_eapp_maf,
                        c2_eapp_maf,
                        c3_eapp_maf,
                        c4_eapp_maf,
                        c1c4_eapp_maf
                    ], axis=1
                )

        corr = corr_df.corr(numeric_only=True)
        loo_corr_dict.update({f"{loo_test[test]}": corr})

    loo_maxr = np.zeros(shape=np.shape(loo_corr_dict[loo_test[0]]))
    loo_maxr_test = np.empty(shape=np.shape(loo_corr_dict[loo_test[0]]), dtype=object)

    for i in range(0, len(loo_corr_dict[loo_test[0]].iloc[:,0])):
        for j in range(0, len(loo_corr_dict[loo_test[0]].iloc[:,0])):

            test_r = np.zeros(shape=len(loo_test))
            test_r_diff = np.zeros(shape=len(loo_test))
            test_r_absdiff = np.zeros(shape=len(loo_test))

            for k in range(0, len(loo_test)):
                test_r[k] = loo_corr_dict[loo_test[k]].iloc[i,j]

            for l in range(0, len(test_r)):
                test_r_diff[l] = test_r[l] - test_r[0]
                test_r_absdiff[l] = np.abs(test_r[l] - test_r[0])

            # loo_maxr[i,j] = np.max(test_r_absdiff)
            loo_maxr[i,j] = test_r_diff[np.argmax(test_r_absdiff)]
            loo_maxr_test[i,j] = str(loo_test[np.argmax(test_r_absdiff)])

    loo_maxr = pd.DataFrame(
        data=loo_maxr,
        columns=loo_corr_dict[loo_test[0]].columns,
        index=loo_corr_dict[loo_test[0]].columns
    )
    loo_maxr_test = pd.DataFrame(
        data=loo_maxr_test,
        columns=loo_corr_dict[loo_test[0]].columns,
        index=loo_corr_dict[loo_test[0]].columns
    )

    return loo_corr_dict, loo_maxr, loo_maxr_test


# Correlation matrices with all datatypes
def envcorr_loosum(df, wax_type):

    match wax_type:
        case "f":
            chain_start = 20
            chain_end = 32
            iso_start = 22
            iso_end = 28

        case "a":
            chain_start = 21
            chain_end = 33
            iso_start = 23
            iso_end = 29

    loo_test = [
        'all',
        'tree',
        'shrub',
        'forb',
        'fern',
        'graminoid',
        'moss',
        'liverwort',
        'lichen'
    ]

    loo_corr_dict = {}

    for test in range(0, len(loo_test)):

        df_test = df[df['cat'] != test].reset_index(drop=True)
        # df_test = df_test.reset_index(drop=True)

        elev = pd.Series(data=df_test.elevation, name='E')
        temp_maf = atpw_fun.maf(df_test, data_type="temp", maf_name='T')
        # swi_maf = atpw_fun.swi(df_test, swi_name="SWI")
        precip_maf = atpw_fun.maf(df_test, data_type="precip", maf_name='P')
        rh_maf = atpw_fun.maf(df_test, data_type="rh", maf_name='RH')
        # pd2h_maf = atpw_fun.maf(df_test, data_type="pd2h", maf_name='d2H')

        tconc = atpw_fun.tconc(df_test, data_type=wax_type, tconc_name='logC', log="yes")
        acl = atpw_fun.wax_acl(df_test, data_type=wax_type, start_chain=chain_start, end_chain=chain_end, acl_name='ACL')
        c1c4_d13c = atpw_fun.iso_avg(df_test, data_type=wax_type, start_chain=iso_start, end_chain=iso_end, iso_type="d13c", iso_name='d13C')
        c1c4_eapp_maf = atpw_fun.wax_precip_eapp_maf(df_test, data_type=wax_type, start_chain=iso_start, end_chain=iso_end, eapp_name='eapp')

        corr_df = pd.concat(
            [
                temp_maf,
                precip_maf,
                rh_maf,
                elev,
                # pd2h_maf,
                tconc,
                acl,
                c1c4_d13c,
                c1c4_eapp_maf
            ], axis=1
        )

        corr = corr_df.corr(min_periods=3, numeric_only=True)
        loo_corr_dict.update({f"{loo_test[test]}": corr})

        if test == 0:
            pval_all = atpw_fun.corr_pvalues(corr_df)
            continue

    loo_maxr = np.zeros(shape=np.shape(loo_corr_dict[loo_test[0]]))
    loo_maxr_test = np.zeros(shape=np.shape(loo_corr_dict[loo_test[0]]))
    # loo_maxr_test = np.empty(shape=np.shape(loo_corr_dict[loo_test[0]]), dtype=object)

    for i in range(0, len(loo_corr_dict[loo_test[0]].iloc[:,0])):
        for j in range(0, len(loo_corr_dict[loo_test[0]].iloc[:,0])):

            test_r = np.zeros(shape=len(loo_test))
            test_r_diff = np.zeros(shape=len(loo_test))
            test_r_absdiff = np.zeros(shape=len(loo_test))

            for k in range(0, len(loo_test)):
                test_r[k] = loo_corr_dict[loo_test[k]].iloc[i,j]

            for l in range(0, len(test_r)):
                test_r_diff[l] = test_r[l] - test_r[0]
                test_r_absdiff[l] = np.abs(test_r[l] - test_r[0])

            # loo_maxr[i,j] = np.max(test_r_absdiff)
            loo_maxr[i,j] = test_r_diff[np.argmax(test_r_absdiff)]
            loo_maxr_test[i,j] = np.argmax(test_r_absdiff)

    loo_maxr = pd.DataFrame(
        data=loo_maxr,
        columns=loo_corr_dict[loo_test[0]].columns,
        index=loo_corr_dict[loo_test[0]].columns
    )
    loo_maxr_test = pd.DataFrame(
        data=loo_maxr_test,
        columns=loo_corr_dict[loo_test[0]].columns,
        index=loo_corr_dict[loo_test[0]].columns
    )

    return loo_corr_dict, loo_maxr, loo_maxr_test, pval_all


# Import master plant wax datafile and map climate info onto it
arc = pd.read_excel(
    'Lindberg_Arctic_terrestrial_plantwax.xlsx',
    sheet_name = 'plantwax'
)
arc_data = atpw_fun.map_era5clim(arc)
arc_data = atpw_fun.map_oipc(arc_data)
arc_data = arc_data.groupby(['genus','species','growth_form','habitat','site_name','sample_year']).mean(numeric_only=True).reset_index()


# Run summary function for n-alkanoic acids and n-alkanes
fames_sum_loo_corr, fames_sum_loo_maxr, fames_sum_loo_tests, fames_sum_loo_p = envcorr_loosum(arc_data, wax_type="f")
# fames_sum_loo_p.to_csv('envcorr_out/fames_sum_loo_p.csv')
# fames_sum_loo_tests.to_csv('envcorr_out/fames_sum_loo_tests.csv')
alks_sum_loo_corr, alks_sum_loo_maxr, alks_sum_loo_tests, alks_sum_loo_p = envcorr_loosum(arc_data, wax_type="a")
# alks_sum_loo_p.to_csv('envcorr_out/alks_sum_loo_p.csv')
# alks_sum_loo_tests.to_csv('envcorr_out/alks_sum_loo_tests.csv')


# Run function for desired outputs
# fames_frac_loo_corr, fames_frac_loo_maxr, fames_frac_loo_tests = envcorr_loo(arc_data, wax_type="f", data_type="frac")
# fames_frac_loo_tests.to_csv('fames_frac_loo_tests.csv')
#
# fames_d13c_loo_corr, fames_d13c_loo_maxr, fames_d13c_loo_tests = envcorr_loo(arc_data, wax_type="f", data_type="d13c")
# fames_d13c_loo_tests.to_csv('fames_d13c_loo_tests.csv')
#
# fames_eapp_loo_corr, fames_eapp_loo_maxr, fames_eapp_loo_tests = envcorr_loo(arc_data, wax_type="f", data_type="eapp")
# fames_eapp_loo_tests.to_csv('fames_eapp_loo_tests.csv')
#
# frac_mask = np.triu(np.ones_like(fames_frac_loo_corr['all'], dtype=bool))
# iso_mask = np.triu(np.ones_like(fames_d13c_loo_corr['all'], dtype=bool))


# Figure script
sum_mask  = np.triu(np.ones_like(fames_sum_loo_corr['all'], dtype=bool))

arc_growth_colors_dict = {
    0: '#FFFFFF',
    1: '#00441b',
    2: '#1b7837',
    3: '#5aae61',
    4: '#a6dba0',
    5: '#d9f0d3',
    6: '#c2a5cf',
    7: '#9970ab',
    8: '#762a83'
}
arc_growth_colors_cm = mpl.colors.ListedColormap([arc_growth_colors_dict[x] for x in arc_growth_colors_dict.keys()])


fig, axs = plt.subplots(3,2,layout='constrained')

# n-alkanoic acid Pearson correlation matrix
ax = axs[0,0]
ax.set_title("n-alkanoic acids")
sns.heatmap(
    ax=ax, data=fames_sum_loo_corr['all'],
    cmap='coolwarm', cbar=False, vmin=-1.0, vmax=1.0, mask=sum_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane Pearson correlation matrix
ax = axs[0,1]
ax.set_title("n-alkanes")
sns.heatmap(
    ax=ax, data=alks_sum_loo_corr['all'],
    cmap='coolwarm', cbar=True, vmin=-1.0, vmax=1.0, mask=sum_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkanoic acid loo max r-value change
ax = axs[1,0]
sns.heatmap(
    ax=ax, data=fames_sum_loo_maxr,
    cmap='PuOr', cbar=False, vmin=-0.6, vmax=0.6, mask=sum_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkane loo max r-value change
ax = axs[1,1]
sns.heatmap(
    ax=ax, data=alks_sum_loo_maxr,
    cmap='PuOr', cbar=True, vmin=-0.6, vmax=0.6, mask=sum_mask,
    annot=True, annot_kws={"fontsize":6}, fmt='.2f'
)

# n-alkanoic acid loo max r-value growth form
ax = axs[2,0]
sns.heatmap(
    ax=ax, data=fames_sum_loo_tests,
    cmap=arc_growth_colors_cm, cbar=False, vmin=0, vmax=8, mask=sum_mask,
    annot=False
)

# n-alkane loo max r-value growth form
ax = axs[2,1]
sns.heatmap(
    ax=ax, data=alks_sum_loo_tests,
    cmap=arc_growth_colors_cm, cbar=False, vmin=0, vmax=8, mask=sum_mask,
    annot=False
)


fig6_loo = plt.gcf()
# fig6_loo.savefig('figures/atpw_figure6_loo.svg')
