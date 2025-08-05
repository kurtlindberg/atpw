'''
Manuscript title: Ecological and environmental controls on plant wax production
and stable isotope fractionation in modern terrestrial Arctic vegetation

Manuscript authors: Kurt R. Lindberg, Elizabeth K. Thomas, Martha K. Raynolds,
Helga Bultmann, Jonathan H. Raberg

DOI: pending

Code author: Kurt R. Lindberg

Supporting functions for scripts in github.com/kurtlindberg/atpw
'''


# Import necessary packages
import numpy as np
import pandas as pd
import math
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from composition_stats import clr, closure, multiplicative_replacement


# Calculate total plant wax concentration
def tconc(df, data_type, tconc_name, log="no"):

    conc_data_type = data_type + "conc"
    conc_df = df.filter(regex=conc_data_type)
    conc_df_arr = np.array(conc_df)
    tconc = np.zeros(shape=len(df.cat))

    for i in range(0, (len(df.cat))):
        if df.conc_ugg_plant[i] == 0:
            tconc[i] = np.nan
        else:
            tconc[i] = np.nansum(conc_df_arr[i,:])

        if tconc[i] == 0:
            tconc[i] = np.nan
            continue

    match log:
        case "yes":
            tconc = np.log(tconc)
        case "no":
            tconc = tconc

    tconc = pd.Series(data=tconc, name=tconc_name)

    return tconc


# Calculate plant wax Average Chain-Length (ACL) for specified chain-length range
def wax_acl(df, data_type, start_chain, end_chain, acl_name):
    '''
    Input arguments:

    df: input dataframe of plant wax data
    rows = individual samples, columns = data categories
    regex aruments based on column headers for the spreadsheet included in this study

    data_type: regex argument for selecting plant wax compound class to filter by
    n-alkanoic acids/FAMEs = 'f', n-alkanes = 'a'

    start_chain: shortest chosen carbon chain-length number

    end_chain: longest chosen carbon chain-length number

    acl_name: column header name (str) for the returned series
    this makes it easier to concatenate the return to an existing dataframe

    Returns:

    acl: pandas Series datatype with ACL calculated for each input dataframe row (sample)
    '''

    conc_data_type = data_type + "conc"
    wax_conc_all = df.filter(regex=conc_data_type).fillna(0)
    # wax_conc_all = wax_conc_all.fillna(0)
    chain_lengths = list(range(start_chain, end_chain+1))
    wax_conc = pd.DataFrame()

    # filter for chain-length concentration data within start-end range
    for n in chain_lengths:
        chain = pd.DataFrame(data=np.array(wax_conc_all.filter(regex=(str(n)))), columns=[str(n)])
        wax_conc = pd.concat([wax_conc, chain], axis=1)

    wax_conc = np.array(wax_conc)
    acl_numer = np.zeros(len(wax_conc[:,0]))
    acl = np.zeros(len(wax_conc[:,0]))

    # calculate ACL for each row (sample)
    for i in range(0, len(wax_conc[:,0])):
        for j in range(0, len(wax_conc[0,:])):
            acl_numer[i] += wax_conc[i,j] * chain_lengths[j]
        acl[i] = acl_numer[i]/np.sum(wax_conc[i,:])

    acl = pd.Series(acl, name=acl_name)

    return acl


# Calculate plant wax Carbon Preference Index (CPI) for specified chain-length range
def wax_cpi(df, data_type, start_chain, end_chain, cpi_name):
    '''
    Input arguments:

    df: input dataframe of plant wax data
        rows = individual samples, columns = data categories
        regex aruments based on column headers for the spreadsheet included in this study

    data_type: regex argument for selecting plant wax compound class to filter by
        n-alkanoic acids/FAMEs = 'f', n-alkanes = 'a'

    start_chain: shortest chosen carbon chain-length number

    end_chain: longest chosen carbon chain-length number

    cpi_name: column header name (str) for the returned series
        this makes it easier to concatenate the return to an existing dataframe

    Returns:

    cpi: pandas Series datatype with CPI calculated for each input dataframe row (sample)
    '''

    conc_data_type = data_type + "conc"
    wax_conc_all = df.filter(regex=conc_data_type).fillna(0)
    chain_lengths = list(range(start_chain, end_chain+1))
    chain_lengths_even = [num for num in chain_lengths if num % 2 == 0]
    chain_lengths_odd = [num for num in chain_lengths if num % 2 == 1]
    wax_conc = pd.DataFrame()

    # filter for chain-length concentration data within start-end range
    for n in chain_lengths:
        chain = pd.DataFrame(data=np.array(wax_conc_all.filter(regex=(str(n)))), columns=[str(n)])
        wax_conc = pd.concat([wax_conc, chain], axis=1)

    wax_conc_even = np.array(wax_conc.filter(items=(map(str, chain_lengths_even))))
    wax_conc_odd = np.array(wax_conc.filter(items=(map(str, chain_lengths_odd))))
    wax_conc = np.array(wax_conc)
    cpi = np.zeros(len(wax_conc[:,0]))

    match conc_data_type:
        case 'fconc':
            for i in range(0, len(wax_conc[:,0])):
                cpi[i] = (np.sum(wax_conc_even[i,0:-1]) + np.sum(wax_conc_even[i,1:])) / (2 * np.sum(wax_conc_odd[i,:]))

            cpi = pd.Series(cpi, name=cpi_name)

            return cpi

        case 'aconc':
            for i in range(0, len(wax_conc[:,0])):
                cpi[i] = (np.sum(wax_conc_odd[i,0:-1]) + np.sum(wax_conc_odd[i:1:])) / (2 * np.sum(wax_conc_even[i,:]))

            cpi = pd.Series(cpi, name=cpi_name)

            return cpi


# Calculate plant wax chain-length fractional abundance
def frac_abd(df, data_type, start_chain, end_chain):

    conc_data_type = data_type + "conc"
    wax_conc_all = df.filter(regex=conc_data_type).fillna(0)
    chain_lengths = list(range(start_chain, end_chain+1))
    wax_conc = pd.DataFrame()

    match data_type:
        case "f":
            chain_lengths = [num for num in chain_lengths if num % 2 == 0]
        case "a":
            chain_lengths = [num for num in chain_lengths if num % 2 == 1]

    for n in chain_lengths:
        wax_chain = pd.DataFrame(data=np.array(wax_conc_all.filter(items=["c"+str(n)+"_"+conc_data_type])), columns=["C"+str(n)])
        wax_conc = pd.concat([wax_conc, wax_chain], axis=1)

    wax_conc_arr = np.array(wax_conc)
    wax_frac = np.zeros(np.shape(wax_conc_arr))

    for i in range(0, len(wax_conc_arr[:,0])):
        for j in range(0, len(wax_conc_arr[0,:])):
            wax_frac[i,j] = wax_conc_arr[i,j]/np.sum(wax_conc_arr[i,:])

    wax_frac = pd.DataFrame(data=wax_frac, columns=wax_conc.columns)

    return wax_frac


# Convert dew point temperature to relative humidity
def dew_to_rh(dew, temp):

    rh = 100*(math.exp((17.625*dew)/(243.04+dew))/math.exp((17.625*temp)/(243.04+temp)))

    return rh


# Map environmental data to plant wax data
def map_era5clim(df, use_sample_year="yes", start_year=1991, end_year=2020):

    era5_columns = [
        "jan_temp",
        "feb_temp",
        "mar_temp",
        "apr_temp",
        "may_temp",
        "jun_temp",
        "jul_temp",
        "aug_temp",
        "sep_temp",
        "oct_temp",
        "nov_temp",
        "dec_temp",
        "jan_precip",
        "feb_precip",
        "mar_precip",
        "apr_precip",
        "may_precip",
        "jun_precip",
        "jul_precip",
        "aug_precip",
        "sep_precip",
        "oct_precip",
        "nov_precip",
        "dec_precip",
        "jan_rh",
        "feb_rh",
        "mar_rh",
        "apr_rh",
        "may_rh",
        "jun_rh",
        "jul_rh",
        "aug_rh",
        "sep_rh",
        "oct_rh",
        "nov_rh",
        "dec_rh"
    ]
    df_clim = df.reindex(df.columns.tolist() + era5_columns, axis=1)
    sites = df_clim.site_name.unique()

    for i in range(0, len(sites)):
        era5_temp = pd.read_csv("climate_reanalyzer/"+sites[i]+"_temp.csv", header=8, index_col="Year")
        era5_precip = pd.read_csv("climate_reanalyzer/"+sites[i]+"_precip.csv", header=8, index_col="Year")
        era5_dew = pd.read_csv("climate_reanalyzer/"+sites[i]+"_dew.csv", header=8, index_col="Year")

        era5_temp_range = era5_temp.query(f"{start_year-1} < index < {end_year+1}")
        era5_precip_range = era5_precip.query(f"{start_year-1} < index < {end_year+1}")
        era5_dew_range = era5_dew.query(f"{start_year-1} < index < {end_year+1}")


        for j in range(0, len(df_clim.site_name)):
            match use_sample_year:
                case "yes":
                    sample_year = df_clim.sample_year[j]

                    if df_clim.site_name[j] == sites[i]:
                        df_clim.loc[j, "jan_temp"] = np.average(np.array(era5_temp.Jan[sample_year]))
                        df_clim.loc[j, "feb_temp"] = np.average(np.array(era5_temp.Feb[sample_year]))
                        df_clim.loc[j, "mar_temp"] = np.average(np.array(era5_temp.Mar[sample_year]))
                        df_clim.loc[j, "apr_temp"] = np.average(np.array(era5_temp.Apr[sample_year]))
                        df_clim.loc[j, "may_temp"] = np.average(np.array(era5_temp.May[sample_year]))
                        df_clim.loc[j, "jun_temp"] = np.average(np.array(era5_temp.Jun[sample_year]))
                        df_clim.loc[j, "jul_temp"] = np.average(np.array(era5_temp.Jul[sample_year]))
                        df_clim.loc[j, "aug_temp"] = np.average(np.array(era5_temp.Aug[sample_year]))
                        df_clim.loc[j, "sep_temp"] = np.average(np.array(era5_temp.Sep[sample_year]))
                        df_clim.loc[j, "oct_temp"] = np.average(np.array(era5_temp.Oct[sample_year]))
                        df_clim.loc[j, "nov_temp"] = np.average(np.array(era5_temp.Nov[sample_year]))
                        df_clim.loc[j, "dec_temp"] = np.average(np.array(era5_temp.Dec[sample_year]))

                        df_clim.loc[j, "jan_precip"] = np.average(np.array(era5_precip.Jan[sample_year]))
                        df_clim.loc[j, "feb_precip"] = np.average(np.array(era5_precip.Feb[sample_year]))
                        df_clim.loc[j, "mar_precip"] = np.average(np.array(era5_precip.Mar[sample_year]))
                        df_clim.loc[j, "apr_precip"] = np.average(np.array(era5_precip.Apr[sample_year]))
                        df_clim.loc[j, "may_precip"] = np.average(np.array(era5_precip.May[sample_year]))
                        df_clim.loc[j, "jun_precip"] = np.average(np.array(era5_precip.Jun[sample_year]))
                        df_clim.loc[j, "jul_precip"] = np.average(np.array(era5_precip.Jul[sample_year]))
                        df_clim.loc[j, "aug_precip"] = np.average(np.array(era5_precip.Aug[sample_year]))
                        df_clim.loc[j, "sep_precip"] = np.average(np.array(era5_precip.Sep[sample_year]))
                        df_clim.loc[j, "oct_precip"] = np.average(np.array(era5_precip.Oct[sample_year]))
                        df_clim.loc[j, "nov_precip"] = np.average(np.array(era5_precip.Nov[sample_year]))
                        df_clim.loc[j, "dec_precip"] = np.average(np.array(era5_precip.Dec[sample_year]))

                        if df_clim.jan_temp[j] != np.nan:
                            df_clim.loc[j, "jan_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Feb[sample_year])), temp=df_clim.jan_temp[j])
                            df_clim.loc[j, "feb_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Jan[sample_year])), temp=df_clim.feb_temp[j])
                            df_clim.loc[j, "mar_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Mar[sample_year])), temp=df_clim.mar_temp[j])
                            df_clim.loc[j, "apr_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Apr[sample_year])), temp=df_clim.apr_temp[j])
                            df_clim.loc[j, "may_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.May[sample_year])), temp=df_clim.may_temp[j])
                            df_clim.loc[j, "jun_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Jun[sample_year])), temp=df_clim.jun_temp[j])
                            df_clim.loc[j, "jul_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Jul[sample_year])), temp=df_clim.jul_temp[j])
                            df_clim.loc[j, "aug_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Aug[sample_year])), temp=df_clim.aug_temp[j])
                            df_clim.loc[j, "sep_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Sep[sample_year])), temp=df_clim.sep_temp[j])
                            df_clim.loc[j, "oct_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Oct[sample_year])), temp=df_clim.oct_temp[j])
                            df_clim.loc[j, "nov_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Nov[sample_year])), temp=df_clim.nov_temp[j])
                            df_clim.loc[j, "dec_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew.Dec[sample_year])), temp=df_clim.dec_temp[j])


                case "no":
                    if df_clim.site_name[j] == sites[i]:
                        df_clim.loc[j, "jan_temp"] = np.average(np.array(era5_temp_range.Jan))
                        df_clim.loc[j, "feb_temp"] = np.average(np.array(era5_temp_range.Feb))
                        df_clim.loc[j, "mar_temp"] = np.average(np.array(era5_temp_range.Mar))
                        df_clim.loc[j, "apr_temp"] = np.average(np.array(era5_temp_range.Apr))
                        df_clim.loc[j, "may_temp"] = np.average(np.array(era5_temp_range.May))
                        df_clim.loc[j, "jun_temp"] = np.average(np.array(era5_temp_range.Jun))
                        df_clim.loc[j, "jul_temp"] = np.average(np.array(era5_temp_range.Jul))
                        df_clim.loc[j, "aug_temp"] = np.average(np.array(era5_temp_range.Aug))
                        df_clim.loc[j, "sep_temp"] = np.average(np.array(era5_temp_range.Sep))
                        df_clim.loc[j, "oct_temp"] = np.average(np.array(era5_temp_range.Oct))
                        df_clim.loc[j, "nov_temp"] = np.average(np.array(era5_temp_range.Nov))
                        df_clim.loc[j, "dec_temp"] = np.average(np.array(era5_temp_range.Dec))

                        df_clim.loc[j, "jan_precip"] = np.average(np.array(era5_precip_range.Jan))
                        df_clim.loc[j, "feb_precip"] = np.average(np.array(era5_precip_range.Feb))
                        df_clim.loc[j, "mar_precip"] = np.average(np.array(era5_precip_range.Mar))
                        df_clim.loc[j, "apr_precip"] = np.average(np.array(era5_precip_range.Apr))
                        df_clim.loc[j, "may_precip"] = np.average(np.array(era5_precip_range.May))
                        df_clim.loc[j, "jun_precip"] = np.average(np.array(era5_precip_range.Jun))
                        df_clim.loc[j, "jul_precip"] = np.average(np.array(era5_precip_range.Jul))
                        df_clim.loc[j, "aug_precip"] = np.average(np.array(era5_precip_range.Aug))
                        df_clim.loc[j, "sep_precip"] = np.average(np.array(era5_precip_range.Sep))
                        df_clim.loc[j, "oct_precip"] = np.average(np.array(era5_precip_range.Oct))
                        df_clim.loc[j, "nov_precip"] = np.average(np.array(era5_precip_range.Nov))
                        df_clim.loc[j, "dec_precip"] = np.average(np.array(era5_precip_range.Dec))

                        if df_clim.jan_temp[j] != np.nan:
                            df_clim.loc[j, "jan_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Jan)), temp=df_clim.jan_temp[j])
                            df_clim.loc[j, "feb_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Feb)), temp=df_clim.feb_temp[j])
                            df_clim.loc[j, "mar_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Mar)), temp=df_clim.mar_temp[j])
                            df_clim.loc[j, "apr_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Apr)), temp=df_clim.apr_temp[j])
                            df_clim.loc[j, "may_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.May)), temp=df_clim.may_temp[j])
                            df_clim.loc[j, "jun_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Jun)), temp=df_clim.jun_temp[j])
                            df_clim.loc[j, "jul_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Jul)), temp=df_clim.jul_temp[j])
                            df_clim.loc[j, "aug_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Aug)), temp=df_clim.aug_temp[j])
                            df_clim.loc[j, "sep_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Sep)), temp=df_clim.sep_temp[j])
                            df_clim.loc[j, "oct_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Oct)), temp=df_clim.oct_temp[j])
                            df_clim.loc[j, "nov_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Nov)), temp=df_clim.nov_temp[j])
                            df_clim.loc[j, "dec_rh"] = dew_to_rh(dew=np.average(np.array(era5_dew_range.Dec)), temp=df_clim.dec_temp[j])

    return df_clim


# Map OIPC precipitation isotope data to plant wax data
def map_oipc(df):

    oipc_columns = [
        "jan_pd2h",
        "feb_pd2h",
        "mar_pd2h",
        "apr_pd2h",
        "may_pd2h",
        "jun_pd2h",
        "jul_pd2h",
        "aug_pd2h",
        "sep_pd2h",
        "oct_pd2h",
        "nov_pd2h",
        "dec_pd2h"
    ]
    df_iso = df.reindex(df.columns.tolist() + oipc_columns, axis=1)

    water_iso = pd.read_excel(
        'Lindberg_Arctic_terrestrial_plantwax.xlsx',
        sheet_name = 'wateriso'
    )

    for i in range(0, len(water_iso.site_name)):
        for j in range(0, len(df_iso.site_name)):
            if df_iso.site_name[j] == water_iso.site_name[i]:
                df_iso.loc[j, "jan_pd2h"] = water_iso.jan_d2h[i]
                df_iso.loc[j, "feb_pd2h"] = water_iso.feb_d2h[i]
                df_iso.loc[j, "mar_pd2h"] = water_iso.mar_d2h[i]
                df_iso.loc[j, "apr_pd2h"] = water_iso.apr_d2h[i]
                df_iso.loc[j, "may_pd2h"] = water_iso.may_d2h[i]
                df_iso.loc[j, "jun_pd2h"] = water_iso.jun_d2h[i]
                df_iso.loc[j, "jul_pd2h"] = water_iso.jul_d2h[i]
                df_iso.loc[j, "aug_pd2h"] = water_iso.aug_d2h[i]
                df_iso.loc[j, "sep_pd2h"] = water_iso.sep_d2h[i]
                df_iso.loc[j, "oct_pd2h"] = water_iso.oct_d2h[i]
                df_iso.loc[j, "nov_pd2h"] = water_iso.nov_d2h[i]
                df_iso.loc[j, "dec_pd2h"] = water_iso.dec_d2h[i]

    return df_iso


# Calculate mean environmental parameter value for a specified range of months
# includes temperature, total precipitation, relative humidity, amount-weight precipitation d2H
def env_avg(df, data_type, env_avg_name, months=["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]):

    df_temp_all = df.filter(regex="temp")
    df_precip_all = df.filter(regex="precip")
    df_rh_all = df.filter(regex="rh")
    df_oipc_all = df.filter(regex="pd2h")

    df_temp = pd.DataFrame()
    df_precip = pd.DataFrame()
    df_rh = pd.DataFrame()
    df_oipc = pd.DataFrame()

    for n in months:
        temp_month = pd.DataFrame(data=np.array(df_temp_all.filter(regex=n)), columns=[n])
        df_temp = pd.concat([df_temp, temp_month], axis=1)

        precip_month = pd.DataFrame(data=np.array(df_precip_all.filter(regex=n)), columns=[n])
        df_precip = pd.concat([df_precip, precip_month], axis=1)

        rh_month = pd.DataFrame(data=np.array(df_rh_all.filter(regex=n)), columns=[n])
        df_rh = pd.concat([df_rh, rh_month], axis=1)

        oipc_month = pd.DataFrame(data=np.array(df_oipc_all.filter(regex=n)), columns=[n])
        df_oipc = pd.concat([df_oipc, oipc_month], axis=1)

    df_temp = np.array(df_temp)
    df_precip = np.array(df_precip)
    df_rh = np.array(df_rh)
    df_oipc = np.array(df_oipc)
    env_avg = np.zeros(len(df_temp[:,0]))

    match data_type:
        case "temp":
            for i in range(0, len(df_temp[:,0])):
                env_avg[i] = np.average(df_temp[i,:])

            env_avg = pd.Series(data=np.array(env_avg), name=env_avg_name)

            return env_avg

        case "precip":
            for i in range(0, len(df_temp[:,0])):
                env_avg[i] = np.sum(df_precip[i,:])

            env_avg = pd.Series(data=np.array(env_avg), name=env_avg_name)

            return env_avg

        case "rh":
            for i in range(0, len(df_temp[:,0])):
                env_avg[i] = np.average(df_rh[i,:])

            env_avg = pd.Series(data=np.array(env_avg), name=env_avg_name)

            return env_avg

        case "pd2h":
            for i in range(0, len(df_temp[:,0])):
                env_avg_numer = 0

                for j in range(0, len(df_temp[0,:])):
                    env_avg_numer += df_oipc[i,j] * df_precip[i,j]

                env_avg[i] = env_avg_numer/np.sum(df_precip[i,:])

            env_avg = pd.Series(data=np.array(env_avg), name=env_avg_name)

            return env_avg


# Calculate mean environmental parameter value for the months above freezing (MAF)
# includes temperature, total precipitation, relative humidity, amount-weight precipitation d2H
def maf(df, data_type, maf_name):

    df_temp_all = np.array(df.filter(regex="temp"))
    df_precip_all = np.array(df.filter(regex="precip"))
    df_rh_all = np.array(df.filter(regex="rh"))
    df_oipc_all = np.array(df.filter(regex="pd2h"))
    maf_data = np.zeros(len(df_temp_all[:,0]))

    match data_type:
        case "temp":
            for i in range(0, len(df_temp_all[:,0])):
                maf_data_numer = 0
                maf_data_denom = 0

                for j in range(0, len(df_temp_all[0,:])):
                    if df_temp_all[i,j] > 0:
                        maf_data_numer += df_temp_all[i,j]
                        maf_data_denom += 1

                maf_data[i] = maf_data_numer/maf_data_denom

            maf_data = pd.Series(data=maf_data, name=maf_name)

            return maf_data

        case "precip":
            for i in range(0, len(df_temp_all[:,0])):
                for j in range(0, len(df_temp_all[0,:])):
                    if df_temp_all[i,j] > 0:
                        maf_data[i] += df_precip_all[i,j]

            maf_data = pd.Series(data=maf_data, name=maf_name)

            return maf_data

        case "rh":
            for i in range(0, len(df_temp_all[:,0])):
                maf_data_numer = 0
                maf_data_denom = 0

                for j in range(0, len(df_temp_all[0,:])):
                    if df_temp_all[i,j] > 0:
                        maf_data_numer += df_rh_all[i,j]
                        maf_data_denom += 1

                maf_data[i] = maf_data_numer/maf_data_denom

            maf_data = pd.Series(data=maf_data, name=maf_name)

            return maf_data

        case "pd2h":
            for i in range(0, len(df_temp_all[:,0])):
                precip_total = 0

                for j in range(0, len(df_temp_all[0,:])):
                    if df_temp_all[i,j] > 0:
                        maf_data[i] += df_oipc_all[i,j] * df_precip_all[i,j]
                        precip_total += df_precip_all[i,j]

                maf_data[i] = maf_data[i]/precip_total

            maf_data = pd.Series(data=maf_data, name=maf_name)

            return maf_data


# Calculate Summer Warmth Index (SWI)
def swi(df, swi_name):

    df_temp = np.array(df.filter(regex="temp"))
    swi = np.zeros(len(df_temp[:,0]))

    for i in range(0, len(df_temp[:,0])):
        for j in range(0, len(df_temp[0,:])):
            if df_temp[i,j] > 0:
                swi[i] += df_temp[i,j]

    swi = pd.Series(data=swi, name=swi_name)

    return swi


# Calculate Plant Wax Apparent Fractionation for range of months
def wax_precip_eapp(df, data_type, start_chain, end_chain, eapp_name, months=["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]):

    conc_data_type = data_type + "conc"
    d2h_data_type = data_type + "d2h"

    wax_conc_all = df.filter(regex=conc_data_type).fillna(0)
    wax_d2h_all = df.filter(regex=d2h_data_type).fillna(0)
    chain_lengths = list(range(start_chain, end_chain+1))
    wax_conc = pd.DataFrame()
    wax_d2h = pd.DataFrame()

    match data_type:
        case "f":
            chain_lengths = [num for num in chain_lengths if num % 2 == 0]
        case "a":
            chain_lengths = [num for num in chain_lengths if num % 2 == 1]

    for n in chain_lengths:
        wax_chain = pd.DataFrame(data=np.array(wax_conc_all.filter(items=["c"+str(n)+"_"+conc_data_type])), columns=[str(n)])
        wax_conc = pd.concat([wax_conc, wax_chain], axis=1)

        d2h_chain = pd.DataFrame(data=np.array(wax_d2h_all.filter(items=["c"+str(n)+"_"+d2h_data_type])), columns=[str(n)])
        wax_d2h = pd.concat([wax_d2h, d2h_chain], axis=1)

    wax_conc = np.array(wax_conc)
    wax_d2h = np.array(wax_d2h)
    eapp_numer = np.zeros(len(wax_d2h[:,0]))

    for i in range(0, len(wax_d2h[:,0])):
        for j in range(0, len(wax_d2h[0,:])):
            if wax_d2h[i,j] == 0:
                wax_conc[i,j] = 0

    for i in range(0, len(wax_d2h[:,0])):
        for j in range(0, len(wax_d2h[0,:])):
            eapp_numer[i] += wax_conc[i,j] * wax_d2h[i,j]

        eapp_numer[i] = eapp_numer[i]/np.sum(wax_conc[i,:])

    era5_precip_all = df.filter(regex="precip")
    oipc_iso_all = df.filter(regex="pd2h")
    era5_precip = pd.DataFrame()
    oipc_iso = pd.DataFrame()

    for n in months:
        # print(n+"_precip")
        precip_month = pd.DataFrame(data=np.array(era5_precip_all.filter(regex=n)), columns=[n])
        era5_precip = pd.concat([era5_precip, precip_month], axis=1)

        iso_month = pd.DataFrame(data=np.array(oipc_iso_all.filter(regex=n)), columns=[n])
        oipc_iso = pd.concat([oipc_iso, iso_month], axis=1)

    era5_precip = np.array(era5_precip)
    oipc_iso = np.array(oipc_iso)
    eapp_denom = np.zeros(len(era5_precip[:,0]))

    for i in range(0, len(era5_precip[:,0])):
        for j in range(0, len(era5_precip[0,:])):
            eapp_denom[i] += era5_precip[i,j] * oipc_iso[i,j]

        eapp_denom[i] = eapp_denom[i]/np.sum(era5_precip[i,:])

    eapp = np.zeros(len(eapp_numer))

    for i in range(0, len(eapp_numer)):
        eapp[i] = (((1000+eapp_numer[i])/(1000+eapp_denom[i]))-1)*1000

    eapp = pd.Series(eapp, name=eapp_name)

    return eapp


# Calculate plant wax apparent fractionation (eapp) for months above freezing (maf)
def wax_precip_eapp_maf(df, data_type, start_chain, end_chain, eapp_name):

    conc_data_type = data_type + "conc"
    d2h_data_type = data_type + "d2h"

    wax_conc_all = df.filter(regex=conc_data_type, axis=1).fillna(0)
    wax_d2h_all = df.filter(regex=d2h_data_type, axis=1).fillna(0)
    chain_lengths = list(range(start_chain, end_chain+1))
    wax_conc = pd.DataFrame()
    wax_d2h = pd.DataFrame()

    match data_type:
        case "f":
            chain_lengths = [num for num in chain_lengths if num % 2 == 0]
        case "a":
            chain_lengths = [num for num in chain_lengths if num % 2 == 1]

    for n in chain_lengths:
        wax_chain = pd.DataFrame(data=np.array(wax_conc_all.filter(items=["c"+str(n)+"_"+conc_data_type], axis=1)), columns=[str(n)])
        wax_conc = pd.concat([wax_conc, wax_chain], axis=1)
        d2h_chain = pd.DataFrame(data=np.array(wax_d2h_all.filter(items=["c"+str(n)+"_"+d2h_data_type], axis=1)), columns=[str(n)])
        wax_d2h = pd.concat([wax_d2h, d2h_chain], axis=1)

    wax_conc = np.array(wax_conc)
    wax_d2h = np.array(wax_d2h)
    eapp_numer = np.zeros(len(wax_d2h[:,0]))

    for i in range(0, len(wax_d2h[:,0])):
        for j in range(0, len(wax_d2h[0,:])):
            if wax_d2h[i,j] == 0:
                wax_conc[i,j] = 0

    for i in range(0, len(wax_d2h[:,0])):
        for j in range(0, len(wax_d2h[0,:])):
            eapp_numer[i] += wax_conc[i,j] * wax_d2h[i,j]

        eapp_numer[i] = eapp_numer[i]/np.sum(wax_conc[i,:])

    era5_temp = np.array(df.filter(regex="temp", axis=1))
    era5_precip = np.array(df.filter(regex="precip", axis=1))
    oipc_iso = np.array(df.filter(regex="pd2h", axis=1))

    eapp_denom = np.zeros(len(eapp_numer))

    for i in range(0, len(era5_temp[:,0])):
        precip_total = 0

        for j in range(0, len(era5_temp[0,:])):
            if era5_temp[i,j] > 0:
                eapp_denom[i] += oipc_iso[i,j] * era5_precip[i,j]
                precip_total += era5_precip[i,j]

        eapp_denom[i] = eapp_denom[i]/precip_total

    eapp_maf = np.zeros(len(eapp_numer))

    for i in range(0, len(eapp_numer)):
        eapp_maf[i] = (((1000+eapp_numer[i])/(1000+eapp_denom[i]))-1)*1000

    eapp_maf = pd.Series(eapp_maf, name=eapp_name)

    return eapp_maf


# Calculate chain-length concentration-weighted isotope value
def iso_avg(df, data_type, iso_type, start_chain, end_chain, iso_name):

    conc_data_type = data_type + "conc"
    iso_data_type = data_type + iso_type

    wax_conc_all = df.filter(regex=conc_data_type).fillna(0)
    wax_iso_all = df.filter(regex=iso_data_type).fillna(0)
    chain_lengths = list(range(start_chain, end_chain+1))

    wax_conc = pd.DataFrame()
    wax_iso = pd.DataFrame()

    match data_type:
        case "f":
            chain_lengths = [num for num in chain_lengths if num % 2 == 0]
        case "a":
            chain_lengths = [num for num in chain_lengths if num % 2 == 1]

    for n in chain_lengths:
        wax_chain = pd.DataFrame(data=np.array(wax_conc_all.filter(items=["c"+str(n)+"_"+conc_data_type])), columns=[str(n)])
        wax_conc = pd.concat([wax_conc, wax_chain], axis=1)
        iso_chain = pd.DataFrame(data=np.array(wax_iso_all.filter(items=["c"+str(n)+"_"+iso_data_type])), columns=[str(n)])
        wax_iso = pd.concat([wax_iso, iso_chain], axis=1)

    wax_conc = np.array(wax_conc)
    wax_iso = np.array(wax_iso)
    iso_avg = np.zeros(len(wax_iso[:,0]))

    for i in range(0, len(wax_iso[:,0])):
        for j in range(0, len(wax_iso[0,:])):
            if wax_iso[i,j] == 0:
                wax_conc[i,j] = 0

    for i in range(0, len(wax_iso[:,0])):
        for j in range(0, len(wax_iso[0,:])):
            iso_avg[i] += wax_iso[i,j] * wax_conc[i,j]

        iso_avg[i] = iso_avg[i]/np.sum(wax_conc[i,:])

    iso_avg = pd.Series(data=iso_avg, name=iso_name)

    return iso_avg


# Calculate the range of chain-length isotope values per sample
def iso_diff(df, data_type, iso_type, start_chain, end_chain, diff_name):

    # conc_data_type = data_type + "conc"
    iso_data_type = data_type + iso_type
    wax_iso_all = df.filter(regex=iso_data_type)
    wax_iso = pd.DataFrame()
    chain_lengths = list(range(start_chain, end_chain+1))

    match data_type:
        case "f":
            chain_lengths = [num for num in chain_lengths if num % 2 == 0]
        case "a":
            chain_lengths = [num for num in chain_lengths if num % 2 == 1]

    for n in chain_lengths:
        iso_chain = pd.DataFrame(data=np.array(wax_iso_all.filter(items=["c"+str(n)+"_"+iso_data_type])), columns=[str(n)])
        wax_iso = pd.concat([wax_iso, iso_chain], axis=1)

    wax_iso = np.array(wax_iso)
    iso_diff = np.zeros(len(wax_iso[:,0]))

    for i in range(0, len(wax_iso[:,0])):
        # for j in range(0, len(wax_iso[0,:])):
        iso_diff[i] = np.max(wax_iso[i,:]) - np.min(wax_iso[i,:])

    iso_diff = pd.Series(data=iso_avg, name=diff_name)

    return iso_diff


# PCA of plant wax chain-length distributions and storing results for plotting
def wax_pca(df, data_type, groupby="species"):

    df_data = df.copy()
    conc_data_type = data_type + "conc"

    match data_type:
        case "f":
            for i in range(0, len(df.c20_fconc)):
                # if np.sum([df['c20_fconc'].iloc[i], df['c22_fconc'].iloc[i], df['c24_fconc'].iloc[i], df['c26_fconc'].iloc[i], df['c28_fconc'].iloc[i], df['c30_fconc'].iloc[i], df['c32_fconc'].iloc[i]]) == 0:
                if np.sum([df_data.c20_fconc[i], df_data.c22_fconc[i], df_data.c24_fconc[i], df_data.c26_fconc[i], df_data.c28_fconc[i], df_data.c30_fconc[i], df_data.c32_fconc[i]]) == 0:
                    # wax_df = df.drop([i])
                    df_data.drop([i], inplace=True)

        case "a":
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
        case "f":
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

        case "a":
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
    wax_PC_scores = pd.DataFrame(
        wax_pca.fit_transform(wax_data_scaled),
        columns=['PC1','PC2','PC3','PC4','PC5','PC6','PC7']
    )

    wax_loadings = pd.DataFrame(
        wax_pca.components_.T,
        columns=['PC1','PC2','PC3','PC4','PC5','PC6','PC7'],
        index=wax_data.columns
    )

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
