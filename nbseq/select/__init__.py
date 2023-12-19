"""
Functions to evaluate the process of phage display selection: enrichment, enrichment specificity, etc.
"""
import pickle
from scipy.interpolate import interp2d
import itertools
import sys
from pathlib import *
import subprocess
import tempfile

import pandas as pd
import numpy as np
import anndata

from ..ft import dataframe_to_anndata, fortify, query as ft_query, get_identifier
from ..utils import *


def calculate_bias_df(
        ft, group_by=['name'],
        comparison=('pre', 'post'),
        comparison_col='amp',
        bias_col='bias',
        round_col='s',
        add_rank=False,
        fillna=None, dropna=False, 
        add_log=False
        ):
    """ Calculate the change in abundance each nanobody in a feature table at each round of selection

    This is intended for calculating amplification bias, e.g. change in abundance between R(n)o and R(n+1)i, but can be abused to calculate the per-round enrichment
    by setting round_col='r', comparison_col='io', comparison=('i','o')

    Parameters
    ----------
    ft : anndata.AnnData
            feature table
    group_by : list of str, optional
            which column(s) of `ft.obs` should be used to identify sample conditions
    comparison : tuple of str, optional
            which values of ``comparison_col`` should be compared. defaults to ('pre', 'post')
    comparison_col : str
            column of `ft.obs` that is used to make comparisons. defaults to 'amp', which is 'pre' for io == 'o' and 'post' for io == 'i'
    round_col : str
            column of `ft.obs` that is used to identify different rounds of selection. 
            defaults to 's' which is calculated as (r-1) for samples where io == 'i' and r for samples where io == 'o'. 
    bias_col : str
            what should the calculated column be called in the output DataFrame
    input_library_query : str, optional
            if given, identifies rows of the feature table corresponding to the raw
            (input) library. These rows will be used to fill missing (`NA`) values for
            ``comparison_col == 'R1i'``; features with missing abundances for this
            round will be filled with the mean relative abundance in the input library.
    input_comparison_name : str, default = 'R1i'
            What value should samples from the input library have `comparison_col` set to?
    fillna : str or None
            how to fill values for nanobodies where ``comparison_col == comparison[0]`` is NaN
            'min': fills with the minimum relative abundance for that sample
            None : give NaN for the bias of such nanobodies
            other: (e.g. 'R3i'): give the value where ``comparison_col == fillna``. For example, if a VHH is not seen in R2i, give the value in R3i
    dropna : bool
        True to drop enrichments that are na
    add_log : bool
        True to add a column `log_enrichment` as `np.log10(enrichment)`

    Returns
    -------
    pd.DataFrame
            One row for each feature for each selection for each round. 
            ``bias_col`` contains the calculated bias value; ``round_col`` indicates the round where the calculation was performed;
            columns from ``group_by`` identify the selection and feature.

    """

    identifier = get_identifier(ft, axis='var')


    # fortify to a long DataFrame, one row per sample/nanobody pair
    ft_df = fortify(ft, relative=True)

    # attach minimal metadata about grouping variable and panning round
    # df = pd.merge(ft_df, ft.obs[group_by + ['round']], left_on='ID', right_index=True)
    r_col = 'r'
    io_col = 'io'

    df = ft_df.join(ft.obs[group_by + [r_col, io_col]], on='ID')

    # in some expts, I use io == 'k' for blanks
    df = df.loc[df[io_col].isin(['i','o']),:]


    # r | io | s | amp
    # 1 | i  | 0 | post
    # 1 | o  | 1 | pre
    # 2 | i  | 1 | post
    # 2 | o  | 2 | pre
    # 3 | i  | 2 | post
    df['s'] =   df[r_col].astype(int) + df[io_col].map({'i':-1, 'o':0}).astype(int)
    df['amp'] = df[io_col].map({'i':'post', 'o':'pre'})

    # construct table, one row per (selection name)/feature, with columns for abundance in R1i, R1o, R2i, R2o, etc.
    dfp = pd.pivot_table(df, index=(
        group_by + [identifier, round_col]), columns=comparison_col, values='abundance')

    # check that all values of `comparison` are represented:
    missing_comparison_values = set(comparison) - set(dfp.columns)
    if len(missing_comparison_values) > 0:
        raise ValueError(
            f"Error computing enrichment: dataset contains no values where {comparison_col} = {missing_comparison_values}. "
            f"Unique values of {comparison_col}: {list(dfp.columns)}. Check that you have correctly filtered the dataset."
        )

    # optionally, assume un-observed nanobodies in the comparison round have same relative abundance as
    # the least abundant nanobody that was observed in that sample.
    if fillna == 'min':
        # find the minimum relative abundance for each sample at round == comparison[0]
        # note: must convert .to_frame() to make .combine_first work sensibly below
        c0_min = dfp.reset_index().groupby(
            group_by)[comparison[0]].min().to_frame()

        # for any nanobodies not observed at round == comparison[0], set their abundances
        # to the minimum relative abundance for that sample
        # important: must keep this a pd.DataFrame (hence dfp[[comparison[0]]] instead of dfp[comparison[0]])
        # for .combine_first to work sensibly, but then grab first column so
        # bias calculation by division works
        c0 = dfp[[comparison[0]]].combine_first(c0_min).iloc[:, 0]
    elif fillna is None:
        c0 = dfp[comparison[0]]
    else:
        # for any nanobodies not observed at round == comparison[0], set their abundances
        # to round == fillna
        c0 = dfp[[comparison[0]]].combine_first(dfp[[fillna]]).iloc[:, 0]

    c1 = dfp[comparison[1]]

    # calculate bias
    dfp[bias_col] = c1 / c0

    # add log10 bias
    if add_log:
        dfp[f'log_{bias_col}'] = np.log10(dfp[bias_col])

    if add_rank:
        dfp = bias_rank(dfp, inplace=True)

    # optionally drop na biases
    if dropna:
        dfp = dfp.dropna(subset=[bias_col])

    out = dfp.reset_index()
    out.columns.name = None
    return out




def calculate_enrichment_df(ft, group_by=['name'],
                            comparison=('R2i', 'R5i'), comparison_col='round',
                            enrichment_col='enrichment',
                            add_start_end=False,
                            add_rank=False,
                            input_library_query=None,
                            input_comparison_name='R1i',
                            fillna='min', dropna=False, add_log=False):
    """ Calculate the enrichment of each nanobody in a feature table over rounds of selection.

    Generally, enrichment = relative abundance at last round / relative abundance at first round. Enrichment is calculated for each nanobody, for each selection condition.
    The feature table has one row for each selection condition-round.

    Parameters
    ----------
    ft : anndata.AnnData
            feature table
    group_by : list of str, optional
            which column(s) of `ft.obs` should be used to identify sample conditions
    comparison : tuple of str, optional
            which values of ``comparison_col`` should be compared.
    comparison_col : str
            column of `ft.obs` that is used to identify different rounds of selection
    input_library_query : str, optional
            if given, identifies rows of the feature table corresponding to the raw
            (input) library. These rows will be used to fill missing (`NA`) values for
            ``comparison_col == 'R1i'``; features with missing abundances for this
            round will be filled with the mean relative abundance in the input library.
    input_comparison_name : str, default = 'R1i'
            What value should samples from the input library have `comparison_col` set to?
    fillna : str or None
            how to fill values for nanobodies where ``comparison_col == comparison[0]`` is NaN
            'min': fills with the minimum relative abundance for that sample
            None : give NaN for the enrichment of such nanobodies
            other: (e.g. 'R3i'): give the value where ``comparison_col == fillna``. For example, if a VHH is not seen in R2i, give the value in R3i
    dropna : bool
        True to drop enrichments that are na
    add_log : bool
        True to add a column `log_enrichment` as `np.log10(enrichment)`

    Returns
    -------
    pd.DataFrame
            Columns are named by values of ``comparison_col``, typically ``R2i``, ``R3i``, etc. ``enrichment_col`` contains the enrichment value.

    """

    identifier = get_identifier(ft, axis='var')
    df_r1 = None

    # the raw library may have been sequenced as a set of separate samples. If so, it could be considered
    # to contain the "R1 input" abundances for all samples. input_library_query describes which samples
    if input_library_query is not None:
        # if given, identify subset of observations that correspond to the input (un-panned) library.
        # calculate the mean
        df_r1 = (
            fortify(
                ft_query(ft, input_library_query, axis='obs'),
                relative=True)
            .groupby(identifier)['abundance'].mean()
            .rename(input_comparison_name).to_frame()
        )

    # fortify to a long DataFrame, one row per sample/nanobody pair
    ft_df = fortify(ft, relative=True)

    # attach minimal metadata about grouping variable and panning round
    # df = pd.merge(ft_df, ft.obs[group_by + ['round']], left_on='ID', right_index=True)
    df = ft_df.join(ft.obs[group_by + [comparison_col]], on='ID')

    # construct table, one row per (selection name)/feature, with columns for abundance in R1i, R1o, R2i, R2o, etc.
    dfp = pd.pivot_table(df, index=(
        group_by + [identifier]), columns=comparison_col, values='abundance')

    # if given, fill (missing) values for R1i with those of the raw library. This allows for some samples
    # to have R1i sequenced separately, if they had a different input library for some reason
    if df_r1 is not None:
        # alternative:
        # dfp.drop(['R1i'], axis='columns').join(df_r1)
        dfp = dfp.combine_first(df_r1)

    if add_start_end:
        dfp = enrichment_start_end(dfp, inplace=True)

    # check that all values of `comparison` are represented:
    missing_comparison_values = set(comparison) - set(dfp.columns)
    if len(missing_comparison_values) > 0:
        raise ValueError(
            f"Error computing enrichment: dataset contains no values where {comparison_col} = {missing_comparison_values}. "
            f"Unique values of {comparison_col}: {list(dfp.columns)}. Check that you have correctly filtered the dataset."
        )

    # optionally, assume un-observed nanobodies in the comparison round have same relative abundance as
    # the least abundant nanobody that was observed in that sample.
    if fillna == 'min':
        # find the minimum relative abundance for each sample at round == comparison[0]
        # note: must convert .to_frame() to make .combine_first work sensibly below
        c0_min = dfp.reset_index().groupby(
            group_by)[comparison[0]].min().to_frame()

        # for any nanobodies not observed at round == comparison[0], set their abundances
        # to the minimum relative abundance for that sample
        # important: must keep this a pd.DataFrame (hence dfp[[comparison[0]]] instead of dfp[comparison[0]])
        # for .combine_first to work sensibly, but then grab first column so
        # enrichment calculation by division works
        c0 = dfp[[comparison[0]]].combine_first(c0_min).iloc[:, 0]
    elif fillna is None:
        c0 = dfp[comparison[0]]
    else:
        # for any nanobodies not observed at round == comparison[0], set their abundances
        # to round == fillna
        c0 = dfp[[comparison[0]]].combine_first(dfp[[fillna]]).iloc[:, 0]

    c1 = dfp[comparison[1]]

    # calculate enrichment
    dfp[enrichment_col] = c1 / c0

    # add log10 enrichment
    if add_log:
        dfp[f'log_{enrichment_col}'] = np.log10(dfp[enrichment_col])

    if add_rank:
        dfp = enrichment_rank(dfp, inplace=True)

    # optionally drop na enrichments
    if dropna:
        dfp = dfp.dropna(subset=[enrichment_col])

    return dfp.reset_index()


def calculate_enrichment_ad(ft, name_col=None, group_by=['name'], feature_col=None, dropna=True, log=False, **kwargs):
    """ Calculate enrichment matrix as an AnnData table

    See ``calculate_enrichment_df``
    """
    if feature_col is None:
        feature_col = get_identifier(ft, axis='var')
        if feature_col == '':
            feature_col = 'feature'

    if name_col is not None and group_by is not None:
        group_by = [name_col]
        print(
            "Warning: both `name_col` and `group_by` specified. `name_col` is deprecated; use `group_by` instead. Using group_by={group_by}")
    if len(group_by) > 1:
        raise NotImplementedError(
            "len(group_by) > 1 not supported for calculate_enrichment_ad; use calculate_enrichment_df instead")

    obs = ft.obs.groupby(group_by).head(1).set_index(group_by)
    var = ft.var
    df = calculate_enrichment_df(ft, group_by=group_by, **kwargs)

    if dropna:
        df = df.dropna(subset=['enrichment'])

    enrichment = dataframe_to_anndata(
        df, obs_col=group_by[0], var_col=feature_col, value_col='enrichment', obs=obs, var=var)

    if log:
        enrichment.X.data = np.log10(enrichment.X.data)
    enrichment.X = sparse_drop_na(enrichment.X)

    return enrichment


def calculate_enrichment_lme4(ft=None,
                              df_path=None, ag_col=None, feature_col='CDR3ID', output_path=None,
                              conda=None, verbose=True, threads=1, return_output=None):

    with tempfile.TemporaryDirectory() as tmpdir:

        # TODO: accept feature table?
        if df_path is None:
            df_path = Path(tmpdir) / 'df.csv'
            # TODO: save df to path
            raise NotImplementedError()

        if output_path is None:
            tmp_output_path = Path(tmpdir) / 'models.csv'
        else:
            tmp_output_path = output_path
        mkdirp_file(tmp_output_path)

        # call lme4 from R
        cmd = ['Rscript', Path(__file__).parent / 'lme4.R',
               str(df_path), str(tmp_output_path), str(
                   feature_col), str(ag_col), str(int(threads))
               ]

        if verbose:
            print("Running R...")
        run_cmd(cmd, conda=conda, verbose=verbose)

        if output_path is None:
            return pd.read_csv(tmp_output_path)


enrichment_lme4 = calculate_enrichment_lme4


def calculate_enrichment_deseq2(ft, design, coefficients,
                                output_prefix=None,
                                conda=None, verbose=True, threads=1):

    with tempfile.TemporaryDirectory() as tmpdir:

        # write (dense) feature table matrix to TSV
        df_path = Path(tmpdir) / 'df.tsv'
        ft.to_df().to_csv(df_path, sep="\t")

        # write metadata to TSV
        metadata_path = Path(tmpdir) / 'obs.tsv'
        obs = ft.obs
        if obs.index.name is None or obs.index.name == '':
            index_label = 'SAMPLE__'
        else:
            index_label = obs.index.name
        obs.to_csv(metadata_path, header=True, index=True,
                   index_label=index_label, sep="\t")

        # output path
        if output_prefix is None:
            output_prefix = Path(tmpdir) / 'out'

        # call DESeq2 in R
        cmd = ['Rscript', Path(__file__).parent / 'deseq2.R',

               # col_data_path
               str(metadata_path),

               # count_data_path
               str(df_path),

               # (sample) identifier_name
               index_label,

               # design_formula
               design,

               # coef_names
               coefficients,

               # out_data_path
               str(output_prefix)
               #    str(df_path), str(tmp_output_prefix), str(
               #        feature_col), str(ag_col), str(int(threads))
               #    ]
               ]

        if verbose:
            print("Running R...")
        run_cmd(cmd, conda=conda, verbose=verbose)

        out = {}
        for coeff in coefficients:
            out[coeff] = pd.read_csv(f"{output_prefix}-{coeff}.csv")
    return out


enr_methods = {
    'df': calculate_enrichment_df,
    'ad': calculate_enrichment_ad,
    'lme4': calculate_enrichment_lme4,
    'deseq2': calculate_enrichment_deseq2
}


def enr(ft, method, *args, **kwargs):
    if method not in enr_methods:
        raise NotImplementedError(
            f"Unknown method '{method}' for calculating enrichment. Options: {list(enr_methods.keys())}")
    func = enr_methods[method]
    return func(ft, *args, **kwargs)


def summarize_enrichment(enrichment):
    print(f"Calculated enrichment for {enrichment.shape[0]} samples across {enrichment.shape[1]} features; {enrichment.X.nnz} "
          f"({enrichment.X.nnz/(enrichment.shape[0]*enrichment.shape[1]):.1%}) non-empty enrichment values")


def calculate_specificity(ft, q1, q2):
    from ..ft import query as ft_query
    e1 = ft_query(ft, q1, axis='sample')
    e2 = ft_query(ft, q2, axis='sample')

    e1.X = sparse_drop_na(e1.X)
    e2.X = sparse_drop_na(e2.X)

    assert (e1.var_names.values == e2.var_names.values).all()

    df = pd.DataFrame({'enrichment_1': e1.X.mean(axis=0).A1,
                       'enrichment_2': e2.X.mean(axis=0).A1},
                      index=e1.var_names)

    # df = (pd.Series(e1.X.mean(axis=0).A1,
    #                 index=e1.var_names, name='e1').join(
    #       pd.Series(e2.X.mean(axis=0).A1,
    #                 index=e2.var_names, name='e2'))
    #      )
    df['specificity'] = df['e1']/df['e2']
    return df


# def calculate_specificities(ft, enrichment, q1, q2, mean='geometric'):
#     """ Calculate the enrichment and abundance specificity in two different groups
#     """

#     from .ft import query as ft_query

#     # if isinstance(partitions, dict):
#     #     names = list(partitions.keys())
#     #     queries = list(partitions.values())
#     # else:
#     #     queries = partitions
#     #     names = list(range(0, len(queries)))

#     a1 = ft_query(ft, q1, axis='sample')
#     a2 = ft_query(ft, q2, axis='sample')

#     enrichment = enrichment[:, a1.var_names.values]
#     e1 = ft_query(enrichment, q1, axis='sample')
#     e2 = ft_query(enrichment, q2, axis='sample')

#     e1X = sparse_drop_na(e1.X)
#     e2X = sparse_drop_na(e2.X)
#     a1X = sparse_drop_na(a1.X)
#     a2X = sparse_drop_na(a2.X)

#     if mean == 'geometric':
#         _mean = sparse_gmean_nonzero
#         _std = sparse_gstd
#     else:
#         _mean = sparse_mean_nonzero
#         _std = sparse_std

#     df = pd.DataFrame({'enrichment_1':           _mean(e1X, axis=0).A1,
#                        'enrichment_1_std':        _std(e1X, axis=0).A1,
#                        'enrichment_2':           _mean(e2X, axis=0).A1,
#                        'enrichment_2_std':        _std(e2X, axis=0).A1,
#                        'abundance_1':            _mean(a1X, axis=0).A1,
#                        'abundance_1_std':         _std(a1X, axis=0).A1,
#                        'abundance_2':            _mean(a2X, axis=0).A1,
#                        'abundance_2_std':         _std(a2X, axis=0).A1,
#                        },
#                       index=a1.var_names)
#     df['abundance_specificity'] = df[f'abundance_1'] / df[f'abundance_2']
#     df['enrichment_specificity'] = df[f'enrichment_1'] / df[f'enrichment_2']
#     return df

def calculate_specificities_mean(ft, enrichment, q1, q2, 
enrichment_model=None, std=False,
mean='geometric', suffixes=('_1','_2')):
    """Calculate the enrichment and abundance specificity in two different groups

    Parameters
    ----------
    ft : anndata.AnnData
        feature table of abundances
    enrichment : anndata.AnnData
        feature table of enrichments
    q1 : str
        query identifying samples in group 1
    q2 : str
        query identifying samples in group 2
    enrichment_model : callable(abundance : float, enrichment : float), optional
        calculate P(enrichment >= y | abundance = x) for each enrichment
    mean : str, optional
        which mean to use, 'geometric' or 'arithmetic'; in either case, mean will be taken over *non-zero* values; by default 'geometric'

    Returns
    -------
    pd.DataFrame
        with columns:
            - f'enrichment{suffixes[0]}'
            - f'enrichment{suffixes[1]}'
            - f'abundance{suffixes[0]}'
            - f'abundance{suffixes[1]}'
            - f'enrichment{suffixes[0]}_pvalue'
            - f'enrichment{suffixes[1]}_pvalue'
            - 'abundance_specificity'
            - 'enrichment_specificity'
    """

    from ..ft import query as ft_query
    from ..utils import intersection_ordered, sparse_mean_nonzero, sparse_std, sparse_gmean_nonzero, sparse_gstd

    # if isinstance(partitions, dict):
    #     names = list(partitions.keys())
    #     queries = list(partitions.values())
    # else:
    #     queries = partitions
    #     names = list(range(0, len(queries)))

    features = intersection_ordered(
        ft.var_names.values, enrichment.var_names.values)

    a1 = ft_query(ft, q1, axis='sample')[:, features]
    a2 = ft_query(ft, q2, axis='sample')[:, features]

    enrichment = enrichment[:, features]
    e1 = ft_query(enrichment, q1, axis='sample')
    e2 = ft_query(enrichment, q2, axis='sample')

    e1X = sparse_drop_na(e1.X)
    e2X = sparse_drop_na(e2.X)
    a1X = sparse_drop_na(a1.X)
    a2X = sparse_drop_na(a2.X)

    if mean == 'geometric':
        _mean = sparse_gmean_nonzero
        _std = sparse_gstd
    else:
        _mean = sparse_mean_nonzero
    _std = sparse_std

    m_e1X = _mean(e1X, axis=0)
    m_e2X = _mean(e2X, axis=0)
    m_a1X = sparse_mean_nonzero(a1X, axis=0)
    m_a2X = sparse_mean_nonzero(a2X, axis=0)
    df = pd.DataFrame({f'mean_enrichment{suffixes[0]}': m_e1X,
                       f'mean_enrichment{suffixes[1]}': m_e2X,
                       f'mean_abundance{suffixes[0]}':  m_a1X,
                       f'mean_abundance{suffixes[1]}':  m_a2X,
                       },
                      index=a1.var_names)

    if std:
        df[f'std_enrichment{suffixes[0]}'] = _std(e1X, axis=0)
        df[f'std_enrichment{suffixes[1]}'] = _std(e2X, axis=0)

    if enrichment_model is not None:
        df[f'enrichment{suffixes[0]}_pvalue'] =  enrichment_model(m_a1X, m_e1X)
        df[f'enrichment{suffixes[1]}_pvalue'] =  enrichment_model(m_a2X, m_e2X)
        
    df['abundance_specificity']  = df[f'mean_abundance{suffixes[0]}']  / df[f'mean_abundance{suffixes[1]}']
    df['enrichment_specificity'] = df[f'mean_enrichment{suffixes[0]}'] / df[f'mean_enrichment{suffixes[1]}']
    return df


calculate_specificities = calculate_specificities_mean


def enrichment_rank_products(df):
    pass


def calculate_enrichment_rank(enrichment: anndata.AnnData):
    """calculates the rank of each non-zero feature within each sample; ranks 
    are ascending (lowest rank is the smallest value). returns a feature 
    table of ranks
    """
    enrichment_rnk = enrichment.copy()
    enrichment_rnk.X = sparse_rank(enrichment_rnk.X, pct=False, asc=True)
    return enrichment_rnk


def calculate_enrichment_pct(enrichment: anndata.AnnData):
    """calculates the percentile of each non-zero feature within each sample; 
    percentile(i,j) = rank(i,j) / nnz(i). sum(percentile(i,j) for j in nz(i)) = 1
    percentiles are ascending (lowest percentile is the smallest value). 
    returns a feature table of percentiles
    """
    enrichment_pct = enrichment.copy()
    enrichment_pct.X = sparse_rank(enrichment_pct.X, pct=True, asc=True)
    return enrichment_pct

# ***


def calculate_pct_product(enrichment, query, enrichment_pct=None, p_value=True, plot_simulation=False, **kwargs):
    """calculate the percentile product and optionally p-value for all features among a subset of samples

    Parameters
    ----------
    enrichment : anndata.AnnData
        feature table of enrichment values
    query : str
        query string to identify relevant samples
    enrichment_pct : anndata.AnnData, optional
        feature table of enrichment percentiles; if omitted, will be calculated by `calculate_enrichment_pct`
    p_value : bool, optional
        if True, perform simulations to calculate a p-value for each enrichment percentile, by default True

    Returns
    -------
    pd.DataFrame
        with columns:
            - feature: name of the feature; determined by ``get_identifier(enrichment)``
            - 'n_samples': how many samples does the feature have a finite enrichment
            - 'pct_product': percentile product of feature
            - 'p_value': if ``p_value == True``, p-value for H0: enrichments are randomly ranked within this group of samples
    """

    from ..ft import query_ids, query as ft_query, add_feature_data, fortify_features, get_identifier

    if enrichment_pct is None:
        enrichment_pct = calculate_enrichment_pct(enrichment)

    # query enrichment and enrichment percentile_product; drop samples not seen in query
    en1 = add_feature_data(
        ft_query(enrichment, query, axis='sample'), 'nsamples')

    if len(en1) == 0:
        raise ValueError(
            f"Unable to calculate enrichment percentile product: \n"
            f"no selections in `enrichment` for query `{query}`.\n"
            f"`enrichment` = {repr(enrichment)}"
        )

    # determine how many samples each feature appeared in _after applying the filtering query_
    en1 = add_feature_data(en1, 'nsamples')

    # find enrichment percentiles
    ep1 = ft_query(enrichment_pct, query, axis='sample')
    
    # calculate percentile product for each group and add column for that metric
    en1.var['pct_product'] = sparse_product(ep1.X, log=False)
    en1.var['pct_sum'] = ep1.X.sum(axis=0).A1

    # convert to dataframe, one row per feature per selection
    identifier = get_identifier(enrichment)
    dfe1 = fortify_features(
        en1, abundance_col='enrichment').reset_index().set_index(identifier)

    # make a contingency table of `n_features`` and `n_samples`,
    # i.e. for each value of `n_samples`, determine how many features
    # have finite enrichments in that many samples
    nse1 = dfe1.reset_index().groupby('nsamples')[
        identifier].count().to_frame().reset_index()

    # simulate rank product, conditioned on number of samples, to generate ecdf and p-value
    c_ecdf1 = simulate_rank_product_conditional(
        n_samples=nse1['nsamples'], n_genes=nse1[identifier], plot=plot_simulation, **kwargs)

    # calculate p-value(s) from ecdf(s); add -log10(p-value)s
    if p_value:
        dfe1['p_value'] = dfe1.apply(
            lambda x: 1-c_ecdf1(x['nsamples'], x['pct_product']), axis='columns')
    else:
        dfe1['p_value'] = float('nan')

    dfe1['nlogp'] = -np.log10(dfe1['p_value'])

    return dfe1


def calculate_specificities_pct_product(enrichment, enrichment_pct,
                                        q1, q2, 
                                        pvalue_1=True, pvalue_2=True, 
                                        suffixes=('_1', '_2'), **kwargs):
    """Calculate the enrichment percentile product in two different groups

    See `calculate_pct_product`; this function just calls this for two groups (e.g. ag+ and ag-) and joins the results

    Parameters
    ----------
    enrichment : anndata.AnnData
        
    enrichment_pct : anndata.AnnData
        
    q1 : str
        query for group 1
    q2 : str
        query for group 1
    pvalue_1 : bool, optional
        True to calculate p-value for percentile products in group 1, by default True
    pvalue_2 : bool, optional
        True to calculate p-value for percentile products in group 2, by default True
    suffixes : tuple of str, optional
        suffix to append to columns for each group, by default ('_1', '_2')

    Returns
    -------
    pd.DataFrame
        see `calculate_pct_product`; the result of joining the two tables, with columns from each group renamed by `suffixes`
    """



    from ..ft import query_ids, query as ft_query, add_feature_data, fortify_features

    # # query enrichment and enrichment percentile_product; drop samples not seen in query
    # en1 = add_feature_data(
    #     ft_query(enrichment, q1, axis='sample'), 'nsamples')
    # en2 = add_feature_data(
    #     ft_query(enrichment, q2, axis='sample'), 'nsamples')
    # # en1 = add_feature_data(drop_empty(ft_query(enrichment, q1, axis='sample'), axis='var'),'nsamples')
    # # en2 = add_feature_data(drop_empty(ft_query(enrichment, q2, axis='sample'), axis='var'),'nsamples')

    # # determine how many samples each feature appeared in _after applying the filtering query_
    # en1 = add_feature_data(en1, 'nsamples')
    # en2 = add_feature_data(en2, 'nsamples')

    # features = list(set(en1.var_names.values).union(en2.var_names.values))

    # # find enrichment percentiles
    # ep1 = ft_query(enrichment_pct, q1, axis='sample')  # [:,features]
    # ep2 = ft_query(enrichment_pct, q2, axis='sample')  # [:,features]

    # # calculate percentile product for each group and add column for that metric
    # en1.var['pct_product'] = sparse_product(ep1.X, log=False)
    # en2.var['pct_product'] = sparse_product(ep2.X, log=False)

    # # convert to dataframe
    # identifier = enrichment.var_names.name
    # dfe1 = fortify_features(
    #     en1, abundance_col='enrichment').reset_index().set_index(identifier)
    # dfe2 = fortify_features(
    #     en2, abundance_col='enrichment').reset_index().set_index(identifier)

    # # for each feature, count number of samples where the enrichment was finite
    # nse1 = dfe1.reset_index().groupby('nsamples')[
    #     identifier].count().to_frame().reset_index()
    # nse2 = dfe2.reset_index().groupby('nsamples')[
    #     identifier].count().to_frame().reset_index()

    # # simulate rank product, conditioned on number of samples, to generate ecdf and p-value
    # c_ecdf1 = simulate_rank_product_conditional(
    #     nse1['nsamples'], nse1[identifier], plot=False, **kwargs)
    # c_ecdf2 = simulate_rank_product_conditional(
    #     nse2['nsamples'], nse2[identifier], plot=False, **kwargs)

    # # calculate p-value(s) from ecdf(s); add -log10(p-value)s
    # if pvalue_1:
    #     dfe1['p_value'] = dfe1.apply(
    #         lambda x: 1-c_ecdf1(x['nsamples'], x['pct_product']), axis='columns')
    # else:
    #     dfe1['p_value'] = float('nan')
    # if pvalue_2:
    #     dfe2['p_value'] = dfe2.apply(
    #         lambda x: 1-c_ecdf2(x['nsamples'], x['pct_product']), axis='columns')
    # else:
    #     dfe2['p_value'] = float('nan')

    # dfe1['nlogp'] = -np.log10(dfe1['p_value'])
    # dfe2['nlogp'] = -np.log10(dfe2['p_value'])

    dfe1 = calculate_pct_product(
        enrichment=enrichment, enrichment_pct=enrichment_pct, query=q1, p_value=pvalue_1, **kwargs)
    dfe2 = calculate_pct_product(
        enrichment=enrichment, enrichment_pct=enrichment_pct, query=q2, p_value=pvalue_2, **kwargs)

    # join data for q1 and q2
    dfe1 = dfe1[['pct_product', 'pct_sum', 'nsamples', 'p_value', 'nlogp']].rename(
        columns=lambda c: c+suffixes[0])
    dfe2 = dfe2[['pct_product', 'pct_sum', 'nsamples', 'p_value', 'nlogp']].rename(
        columns=lambda c: c+suffixes[1])

    df = dfe1.join(dfe2)

    # df['abundance_specificity']  = df[f'abundance_1']  / df[f'abundance_2']
    # df['enrichment_specificity'] = df[f'enrichment_1'] / df[f'enrichment_2']
    return df


def simulate_rank_product_conditional(n_samples, n_genes, n=1000, n_boots=5, plot=False, alpha=0.01):
    """generate an ECDF for P(percentile product > x | n_samples = y)

    for each pair of (n_samples, n_genes), perform `n/n_genes*n_boots` calls 
    to `simulate_percentile_product(n_samples, n_genes)`. Use the results to
    fit an ECDF for each value of `n_samples`

    Parameters
    ----------
    n_samples : iterable of int
        an ECDF for  will be fit for each value of `n_samples`, assuming that there are `n_samples[i]` samples which have `n_genes[i]` each. 
    n_genes : iterable of int
        number of genes to use for simulation, for each value of `n_samples`
    n : int, optional
        number of simulation iterations, by default 1000
    n_boots : int, optional
        multiplier for simulation iterations, by default 5
    plot : bool, optional
        show a boxplot of simulation results, with samples below; by default True
    alpha : float, default = 0.01
        if plotting, show points with p-values smaller than alpha

    Returns
    -------
    callable
        conditional ECDF; accepts parameters ``(y, x)`` and returns P(percentile product <= ``x`` | ``n_samples`` == ``y``).
        Also has attribute `inverse` which a callable giving the inverse ECDF, i.e. accepts ``(y, p)`` and returns ``x`` such that ``p`` == P(percentile product <= ``x`` | ``n_samples`` == ``y``).
    """

    from statsmodels.distributions.empirical_distribution import ECDF

    dfs = []
    for _n_samples, _n_genes in zip(n_samples, n_genes):
        n_iters = int(round((n/_n_genes)*n_boots, 0))
        #print(_n_samples, _n_genes, n_iters)
        for i in range(n_iters):
            spp = simulate_percentile_product(_n_samples, _n_genes)
            dfs.append(pd.DataFrame(
                {'pct_product': spp, 'n_samples': _n_samples, 'n_genes': _n_genes}))
    dff = pd.concat(dfs)
    dfg = dff.groupby(['n_samples', 'n_genes']).sample(n, replace=True)

    conditional_ecdfs = (dfg
        .groupby(['n_samples'])
        .apply(lambda r: ECDF(r['pct_product']))
        .to_dict()
    )

    if plot:
        import seaborn as sns
        import matplotlib.pyplot as plt

        dfg_x_from_p = dfg.groupby(['n_samples', 'n_genes']).apply(lambda r: ecdf_x_from_p(
            1-alpha, r['pct_product'])).rename('pct_product').to_frame().reset_index()

        g = sns.catplot(data=dfg, y='pct_product', x='n_samples', kind='boxen')
        g.ax.scatter(dfg_x_from_p['n_samples'],
                     dfg_x_from_p['pct_product'], c='red')
        g.set_xlabels('# samples where feature appears')
        g.set_ylabels('Percentile product')
        plt.title('Simulated percentile product \n(conditioned on # samples)')
        plt.show()

    conditional_ecdfs[0] = lambda x: float('nan')

    def conditional_ecdf(n_sample, x):
        ecdf = conditional_ecdfs[n_sample]
        return ecdf(x)

    def conditional_ecdf_inverse(n_sample, p):
        ecdf = conditional_ecdfs[n_sample]
        return ecdf_x_from_p(p, ecdf=ecdf)

    setattr(conditional_ecdf, 'inverse', conditional_ecdf_inverse)
    return conditional_ecdf

# ***


def simulate_percentile_product(k, n):
    """randomly permute values 1...n, k times; calculate the percentile product for each gene in this simulated experiment

    Parameters
    ----------
    k : int
        number replicates (samples)
    n : int
        number of genes

    Returns
    -------
    np.array
        with shape ``(n,)``, giving the percentile products for each gene
    """
    from scipy.stats import rankdata
    from numpy.random import default_rng
    x = np.arange(0, n)
    rng = default_rng()

    # perms.shape = (k, n)
    perms = rng.permuted(np.tile(x, k).reshape(k, n), axis=1)

    percentiles = perms / n
    # percentiles = rankdata(perms, axis=0) / n

    # return 10**(np.log10(percentiles).sum(axis=1))
    # rp.shape = (n,)
    rp = 10**(np.log10(percentiles).sum(axis=0))
    return rp

# ***


def enrichment_rank(df_enr, groupby='name', inplace=True, **kwargs):
    if inplace:
        out = df_enr
    else:
        out = pd.DataFrame(
            index=df_enr.index,
            columns=['enrichment_rank'])

    out['enrichment_rank'] = df_enr.groupby(groupby)['enrichment'].rank(**kwargs)
    return out



def enrichment_start_end(df_enr, rounds=None, inplace=True):
    """add column indicating the starting and ending abundance for this feature

    starting abundance is the abundance from the earliest round where the abundance is not NA;
    ending abundance is abundance from the last round where the abundance is not NA

    Parameters
    ----------
    df_enr : pd.DataFrame
        enrichment dataframe from `calculate_enrichment_df`

    Returns
    -------
    pd.DataFrame
        `df_enr` with columns 'start' and 'end' appended
    """
    if rounds is None:
        rounds = df_enr.columns[df_enr.columns.str.startswith('R')]

    if inplace:
        out = df_enr
    else:
        out = pd.DataFrame(
            index=df_enr.index,
            columns=['start', 'end'])

    # rowwise backfill NAs with the first non-NA value; take the first column
    out['start'] = df_enr[rounds].fillna(method='bfill', axis=1).iloc[:, 0]
    # rowwise forward fill NAs with the last non-NA value; take the last column
    out['end']   = df_enr[rounds].fillna(method='ffill', axis=1).iloc[:, -1]
    return out

    # return pd.concat([
    #     df_enr,
    #     df_enr[_rounds].fillna(
    #         method='bfill', axis=1).iloc[:, 0].rename('start'),
    #     df_enr[_rounds].fillna(
    #         method='bfill', axis=1).iloc[:, -1].rename('end'),
    # ], axis=1, verify_integrity=False, copy=False)
    # # return (
    # #     df_enr
    # #      .join(df_enr[_rounds].fillna(method='bfill', axis=1).iloc[:, 0].rename('start'))
    # #      .join(df_enr[_rounds].fillna(method='bfill', axis=1).iloc[:, -1].rename('end'))
    # # )

def enrichment_pvalues(enrichment_df, cond_ecdf=None, abundance_col='abundance', space='cdr3', alpha=0.01, inplace=True):
    """use empirical condition CDF of enrichment and abundance to assign p-values to observed enrichment/abundances

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        with columns 'enrichment', `abundance_col`
    cond_ecdf : callable, optional
        function that takes arguments `abundance`, `enrichment` and yields `P(enrichment >= x | abundance = y)`; 
        alternatively, a str or Path to a pickled version of the above; if not given, will be loaded from a default file path , by default None
    abundance_col : str, optional
        name of the column that contains the relevant abundance measure, by default 'abundance'
    space : str, optional
        if cond_ecdf is None, this is used to locate a pickled version, by default 'cdr3'
    alpha : float, optional
        significance threshold, by default 0.01
    inplace : bool, optional
        True to append columns to `enrichment_df`, False to return a new DataFrame with the same index

    Returns
    -------
    pd.DataFrame
        with columns:
            - 'p_value'
            - 'nlogp': -log10(p_value)
            - 'sig': p_value > alpha
    """
    # load ECDF from input library simulations
    cond_ecdf = load_cond_ecdf(cond_ecdf, space)

    # calculate p-value
    if inplace:
        out = enrichment_df
    else:
        out = pd.DataFrame(index=enrichment_df.index, columns=[
                           'p_value', 'nlogp', 'sig'])

    if type(cond_ecdf).__name__ == 'interp2':
        # slow code path; older scipy.interpolate.interp2d can't do vectorized
        # evaluation of arbitrary points; must evaluate a dense grid of 
        # x/y coordinates
        def cond_pvalue(abundance, enrichment):
            return 1 - cond_ecdf(np.log10(abundance), np.log10(enrichment))

        out['p_value'] = np.vectorize(cond_pvalue)(
            enrichment_df[abundance_col], enrichment_df['enrichment'])
    else:
        # newer scipy.interpolate.RectBivariateSpline.__call__(xs, ys, grid=False) 
        # allows evaluating zs[i] = xs[i], ys[i]
        out['p_value'] = (1 - cond_ecdf(
            np.log10(enrichment_df[abundance_col]), 
            np.log10(enrichment_df['enrichment'])
        ))

    out['nlogp'] = -np.log10(out['p_value'])
    out['sig'] = out['p_value'] < alpha
    return out

def load_cond_ecdf(cond_ecdf=None, space='cdr3'):
    """load a pickled conditional ECDF, previously trained with :func:`enrichment_abundance_ecdf`

    Parameters
    ----------
    cond_ecdf : str or Path
        path to pickled conditional ECDF
    space : str, optional
        if given and cond_ecdf is omitted, will pick a default path.

    Returns
    -------
    Callable[[float,float], float]
        function mapping log(abundance), log(enrichment) to probability
    """
    if cond_ecdf is None:
        cond_ecdf = f'results/tables/{space.lower()}/enrichment/null/ecdf.pickle'

    if not callable(cond_ecdf): #type(cond_ecdf) == str:
        import pickle
        with open(cond_ecdf, 'rb') as f:
            cond_ecdf = pickle.load(f)
    return cond_ecdf


# ***

def compare_binary_phenotypes(ft, phenotypes,
                              n_boots=100, pos_query=None, neg_query=None, pos_value=1, neg_value=None, identifier=None,
                              plot=False, 
                              plot_sparkline = False,
                              n_features_to_plot=10, 
                              alpha=0.05,
                              subplot_size = (3, 1.2)
                              ):
    """calculate probability of observing feature in N phenotype+ samples given distribution in phenotype- samples

    For each phenotype, partition samples into phenotype-positive and phenotype-
    negative samples. For each feature, determine how many times the feature appears
    in the phenotype-positive partition. Then, repeatedly (`n_boots` times) sample
    the phenotype-negative partition, and count how many times the feature appears there.
    Calculate an empirical CDF

    Parameters
    ----------
    ft : anndata.AnnData
        binary feature table
    phenotypes : list of str
        iterable of phenotypes for which to compare positive/negative values. 
        - if `pos_query` and `neg_query` are None: `phenotypes` is an iterable 
            naming columns in `ft.obs`; each column should have values 0 or 1. 
            Phenotype-positive samples are rows where `ft.obs[phenotype] == pos_value`; 
        - if neg_value is None, then phenotype-negative values are all non-positive 
            rows. If neg_value is given, phenotype-negative rows are those 
            were `ft.obs[phenotype] == neg_value`.
        - Alternatively, is `pos_query` and/or `neg_query` is given, then positive
            rows are those selected by those queries, formatted with `phenotype`
    n_boots : int, optional
        number of simulations to run; minimum p-value reported will be `1/n_boots`, by default 100
    pos_value : int, optional
        if pos_query is not given, phenotype-positive samples are rows where `ft.obs[phenotype] == pos_value`, by default 1
    neg_value : bool, optional
        if given, phenotype-negative values will be those where `ft.obs[phenotype] == neg_value`.
        If None, phenotype-negative values will be `~(ft.obs[phenotype] == 1)`
    pos_query : str, optional
        if given, use to choose positive samples by `ft.obs.query(pos_query.format(phenotype=phenotype))`; by default None
    neg_query : str, optional
        if given, use to choose negative samples by `ft.obs.query(neg_query.format(phenotype=phenotype))`; by default None
    identifier : str, optional
        name of the column to use to identify features; by default, chosen from `ft`
    plot : bool, optional
        if True, produce a plot summarizing the simulations for the top features; by default False
    plot_sparkline : bool, optional
        if True, render a plot summarizing the simulations for each feature as a sparkline; by default False
    n_features_to_plot : int, optional
        if plotting, produce subplots for (at most) this many features; features with p_value < alpha are sorted by number of positive samples in descending order; by default 10
    alpha : float, optional
        if plotting, plot features with p_value < alpha, by default 0.05

    Returns
    -------
    pd.DataFrame
        with one row per feature per phenotye; with columns:
            - identifier: identifies the feature
            - 'phenotype': identifies the phenotype
            - 'p_value': P(n Ph+ samples | m Ph- samples)
            - 'n_samples': number of phenotype+ samples that have this feature
            - 'f_samples': fraction of phenotype+ samples that have this feature
    """

    from numpy.random import default_rng
    from ..ft import sum as ft_sum, get_ids
    rng = default_rng()

    # store output
    p_values = np.zeros(shape=(len(phenotypes), ft.shape[1]))
    n_samples = np.zeros(shape=(len(phenotypes), ft.shape[1]))
    f_samples = np.zeros(shape=(len(phenotypes), ft.shape[1]))

    features = get_ids(ft, axis='var')
    if identifier is None:
        identifier = get_identifier(ft)

    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(ncols=len(phenotypes), nrows=n_features_to_plot, 
                                figsize=(subplot_size[0] * len(phenotypes), subplot_size[1] * n_features_to_plot), squeeze=False)

    if plot_sparkline:
        sparklines = np.empty(shape=(len(phenotypes), ft.shape[1]), dtype='object')
        sparklines.fill('')

    p_features = ft.shape[1]

    for i, phenotype in enumerate(phenotypes):
        if pos_query is not None:
            ph_pos = ft.obs.eval(pos_query.format(
                phenotype=phenotype), inplace=False).values
        elif pos_query is None:
            ph_pos = ft.obs[phenotype] == pos_value

        if neg_query is not None:
            ph_neg = ft.obs.eval(neg_query.format(
                phenotype=phenotype), inplace=False).values
        elif neg_query is None:
            if neg_value is None:
                ph_neg = ~ph_pos
            else:
                ph_neg = ft.obs[phenotype] == neg_value

        # partition feature table into phenotype+ and phenotype- samples
        ftf_pos = ft[ph_pos, :]
        ftf_neg = ft[ph_neg, :]

        n_pos = sum(ph_pos)
        n_neg = sum(ph_neg)

        # choose `n_pos` random samples from ftf_neg (e.g. the antigen- samples).
        # do this `n_boots` times.
        idx_neg = rng.integers(low=0, high=n_neg, size=(n_boots, n_pos))

        Xn_sums = np.zeros(shape=(n_boots, p_features), dtype='float64')

        # for each bootstrap iteration:
        for j in range(n_boots):
            Xn = ftf_neg[idx_neg[j, :], :]

            # calculate how many of the randomly-chosen Ag- samples meet the criteria
            Xn_sums[j, :] = Xn.X.sum(axis=0)

        # for each feature, calculate how many ag+ samples meet the critera
        # extract numpy array because multidimensional indexing on pd.Series is not supported
        # "Multi-dimensional indexing (e.g. `obj[:, None]`) is no longer supported. Convert to a numpy array before indexing instead."
        Xp_sum = ft_sum(ftf_pos, axis='var').values

        # for each feature, calculate the probability that the observed sum exceeds
        # those in the simulation
        p_values[i, :] = np.maximum(
            # Xp_sum.shape = (1, ft.shape[1])
            # Xn_sum.shape = (n_boots, ft.shape[1])
            # np.greater(Xp_sum[np.newaxis, :], Xn_sums)
            #   .shape = (n_boots, ft.shape[1])
            #   [i, j] = True if, for feature j, in simulation i # positive samples > # negative samples
            #   Xp_sum is broadcasted over Xn_sums axis 0
            # .mean(axis=0) : average over all simulations
            # 1-x : get fraction of observations where Xn_sum > Xp_sum (have to place Xp_sum first for broadcasting to work)
            1 - np.greater(Xp_sum[np.newaxis, :], Xn_sums).mean(axis=0),

            # even if Xn_sum>Xp_sum in 0 simulatiuons, do not report p-value
            # less than 1/n_boots, since we only did that many simulations
            1.0/float(n_boots))
        n_samples[i, :] = Xp_sum
        f_samples[i, :] = Xp_sum/n_pos

        if plot:
            # find features where p-value is significant, then sort by number of samples
            # np.argpartition(p_values[i,:], n_features_to_plot)[:n_features_to_plot]
            features_to_plot = np.argwhere(p_values[i, :] < alpha).reshape(-1)
            features_to_plot = features_to_plot[np.argsort(
                n_samples[i, features_to_plot])[::-1]]
            features_to_plot = features_to_plot[:min(
                n_features_to_plot, len(features_to_plot))]

            for k, feature in enumerate(features_to_plot):
                feature_hash = features[feature][:6]
                axs[k, i].hist(Xn_sums[:, feature], bins=np.arange(
                    0, n_pos + 1.5) - 0.5, color=f"#{feature_hash}")
                axs[k, i].axvline(Xp_sum[feature], color='r')

                axs[k, i].text(0.99, 0.95, f"{feature_hash} \n {Xp_sum[feature]:.0f} sample(s) \n p ≤ {p_values[i,feature]:0.2f}",
                               horizontalalignment='right', verticalalignment='top',
                               transform=axs[k, i].transAxes)

                axs[k, i].set(ylabel=None, yticklabels=[])

            axs[0, i].set_title(phenotype)
            axs[-1, i].set(xlabel='# samples')

        if plot_sparkline:
            from ..viz.utils import mini_histogram, sparkline

            # only bother making these for features where p_values < alpha
            features_to_plot = np.argwhere(p_values[i, :] < alpha).reshape(-1)

            for feature in features_to_plot:
                feature_hash = features[feature][:6]

                sparklines[i, feature] = mini_histogram(
                    Xn_sums[:, feature], 
                    vmax=Xp_sum[feature],
                    bins=(np.arange(0, n_pos + 1.5) - 0.5), 
                    color=f"#{feature_hash}",
                    formatter="{:.0f}",
                    marks = [Xp_sum[feature]], 
                    mark_kws = dict(color='r')
                )

        
    if plot:
        for k in range(n_features_to_plot):
            axs[k, 0].set(ylabel='% simulations')
        fig.set_tight_layout(True)

    out = pd.DataFrame({
        # C = row-major order = last index (columns = features) changes fastest
        #   phenotype_0, feature_0, p_values[0,0],
        #   phenotype_0, feature_1, p_values[0,1],
        #   ...
        #   phenotype_0, feature_P, p_values[0,P],
        #   phenotype_1, feature_0, p_values[1,0],
        #   ...
        #   phenotype_N, feature_P, p_values[N,P],
        'p_value':   p_values.reshape(-1, order='C'),
        'n_samples': n_samples.reshape(-1, order='C'),
        'f_samples': f_samples.reshape(-1, order='C'),
        # tile: repeats the entire array
        identifier:  np.tile(features, len(phenotypes)),
        'phenotype': np.repeat(phenotypes, len(features))
    })
    if plot_sparkline:
        out['plot'] = sparklines.reshape(-1, order='C')
    return out

# ***


def partition_library(ex, library_query, pan_query=None, show_summary=False):
    if pan_query is None:
        pan_query = f"~({library_query})"

    fts = dict(
        ft_pan_cdr3=ex.query(library_query, space="cdr3"),
        ft_pan_aa=ex.query(library_query, space="aa"),

        ft_lib_cdr3=ex.query(pan_query,     space="cdr3"),
        ft_lib_aa=ex.query(pan_query,     space="aa")
    )

    if show_summary:
        from textwrap import dedent
        print(dedent("""\
        Total # AAs observed:   {ex.fts.aa.shape[1]:>12}
            - in raw library:   {(ft_lib_aa.X.sum(axis=0) > 0).sum():>12}
        
        Total # CDR3s observed: {ex.fts.cdr3.shape[1]:>12}
            - in raw library: {(ft_lib_cdr3.X.sum(axis=0) > 0).sum():>12}
        """.format(fts)))
    return fts


def features_not_in_raw_library(ft, library_query, pan_query=None, last_query="round == 'R5i'"):
    """show a venn diagram and some statistics about how many features appeared in selections but not in the input library (and vice-versa)

    Parameters
    ----------
    ft : anndata.AnnData
        feature table
    library_query : str
        query to identify which samples in `ft` are in the input library
    pan_query : str, optional
        query to identify which samples in `ft` are from panned selections; by default, f"~({library_query})"
    last_query : str, optional
        query to identify which samples in `ft` are in the last round of selection, by default "round == 'R5i'"

    Returns
    -------
    _type_
        _description_
    """
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3
    import upsetplot
    from ..ft import query as ft_query, get_ids

    if pan_query is None:
        pan_query = f"~({library_query})"

    ft_pan = ft_query(ft, pan_query, axis='obs')
    ft_lib = ft_query(ft, library_query, axis='obs')
    ft_lst = ft_query(ft, last_query, axis='obs')

    feat_in_lib = (ft_lib.X.sum(axis=0) > 0).A1
    feat_in_pan = (ft_pan.X.sum(axis=0) > 0).A1
    feat_in_lst = (ft_lst.X.sum(axis=0) > 0).A1

    features_in_lib = get_ids(ft_lib, 'var')[feat_in_lib]
    features_not_lib = get_ids(ft_pan, 'var')[~feat_in_lib]

    features_in_pan = get_ids(ft_pan, 'var')[feat_in_pan]
    features_in_lst = get_ids(ft_lst, 'var')[feat_in_lst]

    print(f"Total features:      {ft.shape[1]:>12.0f}")
    print(f"- Seen in lib:       {len(features_in_lib):>12.0f}")
    print(
        f"   - Seen in last:   {(ft_lst[:,features_in_lib].X.sum(axis=0) > 0).sum():>12.0f}")
    print(
        f"- Not in lib:        {len(features_not_lib):>12.0f} = ({len(features_not_lib) / ft.shape[1]:.2%})")
    print(
        f"   - Seen in last:   {(ft_lst[:,features_not_lib].X.sum(axis=0) > 0).sum():>12.0f}")
    print()

    ft_pan_reads = ft_pan.X.sum()
    print(f"Total reads:         {ft.X.sum():>12.0f}")
    print(f"- Total reads (pan): {ft_pan_reads:>12.0f}")
    print(f"  - Seen in lib:     {ft_pan[:,features_in_lib].X.sum():>12.0f}")
    print(f"    - Seen in last:  {ft_lst[:,features_in_lib].X.sum():>12.0f}")
    last_round_reads = ft_lst.X.sum()
    ft_pan_not_lib_reads = ft_pan[:, features_not_lib].X.sum()
    print(
        f"  - Not in lib:      {ft_pan_not_lib_reads:>12.0f} = ({ft_pan_not_lib_reads / ft_pan_reads:.2%}) ")
    print(
        f"    - Seen in last:  {ft_lst[:,features_not_lib].X.sum():>12.0f} = ({ft_lst[:,features_not_lib].X.sum() / last_round_reads:.2%})")
    print(f"- Total reads (lib): {ft_lib.X.sum():>12.0f}")

    features = set(get_ids(ft, 'var'))
    features_in_lib = set(features_in_lib)
    features_in_pan = set(features_in_pan)
    features_in_lst = set(features_in_lst)

    fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
    venn3((features_in_lib,
           features_in_pan,
           features_in_lst),
          set_labels=('Input', 'Experiment', 'Round 5'), ax=axs[0])
    axs[0].set_title("Features")

    subsets = {
        '100': int(ft_lib[:, list(features_in_lib - features_in_pan - features_in_lst)].X.sum()),
        '010': int(ft_pan[:, list(features_in_pan - features_in_lib - features_in_lst)].X.sum()),
        # int(ft_lst[:, list(features_in_lst - features_in_lib - features_in_pan)].X.sum()),
        '001': 0,
        '110': int(ft_lib[:, list(features_in_lib - features_in_lst)].X.sum()),
        '011': int(ft_pan[:, list(features_in_pan - features_in_lib)].X.sum()),
        # int(ft_lib[:, list(features_in_lib - features_in_pan)].X.sum()),
        '101': 0,
        '111': (int(ft_pan[:, list(features_in_lib & features_in_pan & features_in_lst)].X.sum())),
    }

    venn3(
        subsets=subsets,
        set_labels=('Input', 'Experiment', 'Round 5'),
        ax=axs[1],
        subset_label_formatter=lambda x: f"{x:.2g}"
    )

    return upsetplot.plot(
        upsetplot.from_indicators(
            pd.DataFrame([[bool(int(v)) for v in k] for k in subsets.keys()], columns=(
                'Input', 'Experiment', 'Round 5')),
            data=list(subsets.values())
        )
    )

    # venn3(
    #     subsets={
    #         '100':np.log10(int(ft_lib[:, list(features_in_lib - features_in_pan - features_in_lst)].X.sum())),
    #         '010':np.log10(int(ft_pan[:, list(features_in_pan - features_in_lib - features_in_lst)].X.sum())),
    #         '001':0,#int(ft_lst[:, list(features_in_lst - features_in_lib - features_in_pan)].X.sum()),
    #         '110':np.log10(int(ft_lib[:, list(features_in_lib - features_in_lst)].X.sum())),
    #         '011':np.log10(int(ft_pan[:, list(features_in_pan - features_in_lib)].X.sum())),
    #         '101':0,#int(ft_lib[:, list(features_in_lib - features_in_pan)].X.sum()),
    #         '111':np.log10((int(ft_pan[:, list(features_in_lib & features_in_pan & features_in_lst)].X.sum()))),
    #     },
    #     set_labels=('Input','Experiment', 'Round 5'),
    #     ax=axs[1],
    #     subset_label_formatter= lambda x: f"10^{x:.1f}"
    # )
    axs[1].set_title('Reads')


def enrichment_selection_summary(df_enr, phenotypes, min_enrichment = 2, alpha=0.01, selection_metadata=None, identifier='CDR3ID'):

    def summarize_ph(row):
        
        n_sig = row['p_value']
        n_enr = row['enrichment']
        phenotype = row['phenotype']
        
        _ph_out = []
        if n_sig > 0:
            _ph_out.append(f"{n_sig:.0f}*")
        if n_enr > 0:
            _ph_out.append("{:.0f}".format(n_enr))
        if n_sig + n_enr > 0:
            return("+".join(_ph_out) + " " + phenotype + "+")
        else: return ''

    if selection_metadata is not None:
        df_enr = df_enr.join(selection_metadata, on='name')
    
    # df = df_enr[[identifier, 'enrichment', 'p_value'] + list(phenotypes)]
    # return df.groupby(identifier).apply(_summarize_selections).rename('selection_summary')

    df = df_enr[[identifier, 'p_value','enrichment'] + list(phenotypes)].copy()

    df['p_value'] = df['p_value'] < alpha
    df['enrichment'] = df['enrichment'] > min_enrichment

    dfm = pd.melt(df, id_vars=[identifier, 'p_value','enrichment'], var_name='phenotype')
    # dfm: [identifier,'p_value','enrichment','phenotype','value]
    dfm = dfm[dfm.value > 0]
    # nsels: index = [identifier,'phenotype']; columns=['p_value','enrichment'], values = # of sig/enriched selections
    nsels = dfm.groupby([identifier,'phenotype'])[['p_value','enrichment']].sum().sort_values(['p_value','enrichment'], ascending=False)

    def summarize_nsels_per_feature(nsels_feature):
        nsels_per_ph = nsels_feature.reset_index().apply(summarize_ph, axis=1)
        return ", ".join(nsels_per_ph[nsels_per_ph.astype('bool')])

    return nsels.groupby(identifier).apply(summarize_nsels_per_feature)


# ***

def simulate_enrichments(ft,ft2=None,
                         r1_label='abundance', r2_label='end_abundance', feature_label='feature',
                         zeros='drop',
                         log=True,
                         inv=True):
    """simulate fake enrichment values by calculating enrichment for all nanobodies in all pairs of samples in this feature table

    Parameters
    ----------
    ft : anndata.AnnData
        feature table; all samples should be similar, e.g. input library
    r1_label : str, optional
        what to call the column containing "starting" abundance, by default 'abundance'
    r2_label : str, optional
        what to call the column containing "ending" abundance, by default 'end_abundance'
    feature_label : str, optional
        what to call the column containing feature, by default 'feature'
    zeros : str, optional
        how to handle circumstances where one of the abundances is zero:
            - 'drop': drop points where one abundance is zero
            - 'min':  if one abundance is zero, fill it with the minimum value 

    inv : bool, optional
        invert the starting and ending abundance, by default True

    Returns
    -------
    pd.DataFrame
        # rows ≈ (# features) * (# pairwise combinations of samples), except infinite or NA enrichments are dropped
    """

    from ..ft import to_relative, get_ids
    from ..utils import sparse_min_nz

    ft = to_relative(ft)
    ft1 = ft
    if ft2 is None:
        ft2 = ft
        features = get_ids(ft, 'var')
        samples = get_ids(ft, 'obs')
        sample_pairs = itertools.combinations(samples, 2)
        mins1 = mins2 = dict(zip(samples, sparse_min_nz(ft2.X, axis=1)))
    else:
        ft2 = to_relative(ft2)
        features1 = get_ids(ft1,  'var')
        features2 = get_ids(ft2, 'var')
        features = np.intersect1d(features1, features2)
        ft1 = ft1[:,  features]
        ft2 = ft2[:, features]

        samples1 = get_ids(ft1, 'obs')
        samples2 = get_ids(ft2, 'obs')
        
        sample_pairs = itertools.product(samples1, samples2)

        mins1 = dict(zip(samples1, sparse_min_nz(ft1.X, axis=1)))
        mins2 = dict(zip(samples2, sparse_min_nz(ft2.X, axis=1)))
    del ft

    dfs = []
    

    for (s1, s2) in sample_pairs:
        if zeros == 'min':
            x1 = ft1[s1, :].X.todense().A1
            x2 = ft2[s2, :].X.todense().A1
            x1[x1 == 0] = mins1[s1]
            x2[x2 == 0] = mins2[s2]
            e1 = x1 / x2
            e2 = x2 / x1
        else:
            e1 = (ft1[s1, :].X / ft2[s2, :].X).A1
            e2 = 1/e1
            x1 = ft1[s1, :].X.todense().A1
            x2 = ft2[s2, :].X.todense().A1
        dfs.append(pd.DataFrame({'enrichment': e2 if inv else e1,
                                 'feature':       features,
                                 'abundance':     x1,
                                 'end_abundance': x2
                                 }))
        dfs.append(pd.DataFrame({'enrichment': e1 if inv else e2,
                                 'feature':       features,
                                 'abundance':     x2,
                                 'end_abundance': x1
                                 }))

    df = pd.concat(dfs)
    if zeros == 'drop':
        df = df[df['enrichment'] != 0]
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    if log:
        df['log_enrichment'] = np.log10(df['enrichment']).replace([float('inf'), float('-inf')],0)
    return df


# def simulate_enrichments(ft,
#                          r1_label='abundance', r2_label='end_abundance', feature_label='feature',
#                          zeros='drop',
#                          inv=True):
#     """simulate fake enrichment values by calculating enrichment for all nanobodies in all pairs of samples in this feature table

#     Parameters
#     ----------
#     ft : anndata.AnnData
#         feature table; all samples should be similar, e.g. input library
#     r1_label : str, optional
#         what to call the column containing "starting" abundance, by default 'abundance'
#     r2_label : str, optional
#         what to call the column containing "ending" abundance, by default 'end_abundance'
#     feature_label : str, optional
#         what to call the column containing feature, by default 'feature'
#     zeros : str, optional
#         how to handle circumstances where one of the abundances is zero:
#             - 'drop': drop points where one abundance is zero
#             - 'min':  if one abundance is zero, fill it with the minimum value 

#     inv : bool, optional
#         invert the starting and ending abundance, by default True

#     Returns
#     -------
#     pd.DataFrame
#         # rows ≈ (# features) * (# pairwise combinations of samples), except infinite or NA enrichments are dropped
#     """

#     from ..ft import to_relative, get_ids
#     from ..utils import sparse_min_nz

#     ft = to_relative(ft)

#     features = get_ids(ft, 'var')
#     samples = get_ids(ft, 'obs')

#     dfs = []
#     mins = dict(zip(samples, sparse_min_nz(ft.X, axis=1)))

#     for (s1, s2) in itertools.combinations(samples, 2):
#         if zeros == 'min':
#             x1 = ft[s1, :].X.todense().A1
#             x2 = ft[s2, :].X.todense().A1
#             x1[x1 == 0] = mins[s1]
#             x2[x2 == 0] = mins[s2]
#             e1 = x1 / x2
#             e2 = x2 / x1
#         else:
#             e1 = (ft[s1, :].X / ft[s2, :].X).A1
#             e2 = 1/e1
#             x1 = ft[s1, :].X.todense().A1
#             x2 = ft[s2, :].X.todense().A1
#         dfs.append(pd.DataFrame({'enrichment': e2 if inv else e1,
#                                  'feature':       features,
#                                  'abundance':     x1,
#                                  'end_abundance': x2
#                                  }))
#         dfs.append(pd.DataFrame({'enrichment': e1 if inv else e2,
#                                  'feature':       features,
#                                  'abundance':     x2,
#                                  'end_abundance': x1
#                                  }))

#     df = pd.concat(dfs)
#     if zeros == 'drop':
#         df = df[df['enrichment'] != 0]
#     df = df.replace([np.inf, -np.inf], np.nan).dropna()
#     return df


def ecdf_p_plot(xs, ps=[0.95, 0.99, 0.999]):
    """plot empirical CDF with vertical/horizontal lines showing threshold for particular significance

    for each probability `p` in `ps`, find the value `x` where ``P(xs <= x) = p``; show a vertical
    line from `x` to the curve, then a horizonal line from the curve to `p`

    Parameters
    ----------
    xs : iterable of float
        sequence of values of which to show ECDF
    ps : list of float, optional
        mark position of these p-values, by default [0.95, 0.99, 0.999]
    """

    # df = simulate_enrichment(ft_lib_cdr3)
    # ecdf_p_plot(df.loc[df['enrichment'] > 1,'enrichment'])
    # ecdf_p_plot(df.loc[df['abundance'] > 0,'abundance'])

    from statsmodels.distributions.empirical_distribution import ECDF
    import matplotlib.pyplot as plt

    ecdf = ECDF(xs)

    plt.plot(ecdf.x, ecdf.y, label='ECDF')
    plt.xscale('log')
    plt.xlabel(xs.name)
    plt.ylabel('proportion')
    # plt.axhline(0.95,c='red')

    xmin = min(xs)
    for i, p in enumerate(ps):
        e_p = min(ecdf.x[ecdf.y > p])

        # plt.axhline(p,c='red')
        label = f"P({xs.name} < {e_p:0.2g}) = {p}"
        plt.hlines(y=p, xmin=xmin, xmax=e_p, color=f'C{i}', label=label)
        plt.vlines(x=e_p, ymin=0, ymax=p, color=f'C{i}')
        print(label)


def abundance_enrichment_plot(df, nz=True, n_points=10000, kind='kde', scatter=True, log_scale=True, height=10, scatter_kw=None, **kwargs):
    """Plot joint distribution of enrichment vs abundance

    Parameters
    ----------
    df : pd.DataFrame
            with columns 'enrichment' and 'abundance'
    nz : bool, optional
            True to remove points where enrichment or abundance is 0, by default True
    n_points : int, optional
            if given, downsample df to this many rows, by default 10000
    scatter : bool, optional
            if True, overlay a scatterplot over the histogram/KDE, by default True
    log_scale : bool, optional
            if True, plot the joint histogram/KDE with log scales, by default True
    height : int, optional
            height of the plot, by default 10

    Returns
    -------
    sns.FacetGrid
            plot
    """
    import seaborn as sns

    # dff_nz = df.query("(enrichment >= 1) & (abundance > 0)")

    if nz:
        dff = df.query(
            "(enrichment != 0) & (abundance > 0)").reset_index(drop=True)
    else:
        dff = df

    if n_points is not None:
        dffs = dff.sample(n=n_points)
    else:
        dffs = dff

    # sns.jointplot(data=dffs, x='abundance', y='enrichment',
    #               kind='hist', log_scale=True, height=10)

    g = sns.jointplot(data=dffs, x='abundance', y='enrichment',
                      kind=kind, log_scale=log_scale, height=height, **kwargs)

    if scatter:
        if scatter_kw is None:
            scatter_kw = dict()
        scatter_kw = {'alpha':0.5*(1000/len(dffs)), **scatter_kw}
        dffs.plot.scatter(x='abundance', y='enrichment', ax=g.ax_joint, 
                          **scatter_kw)

    return g


def enrichment_abundance_ecdf(
        df,
        abundance_bins=None, enrichment_bins=None, p_values=[0.95, 0.99, 0.999, 0.9999],
        plot=None,
        plot_joint_pmf=False,
        plot_marginal_abundances=False,
        plot_conditional_pmf=False,
        plot_conditional_cdf=False, 
        fig_kw = None,
        force_interp2d=False,
        verbose=False):
    """Calculate an empirical CDF of log10(enrichment) conditioned on (starting) log10(abundance).

    Given a series of (abundance, enrichment) points: fit a joint Gaussian KDE;
    evaluate it at points on a grid (given by `abundance_bins` and `enrichment_bins`); 
    divide by marginal (abundance) probability at each point to calculate the conditional 
    probability; cumulative sum to calculate the conditional PMF; interpolate

    Parameters
    ----------
    df : pd.DataFrame
        with columns 'enrichment' and 'abundance', e.g. from ``df = simulate_enrichment(ft_input_lib_cdr3)``
    abundance_bins : np.array, optional
        boundaries of bins for abundance values; if omitted, will be chosen automatically
    enrichment_bins : np.array, optional
        boundaries of bins for enrichment values; if omitted, will be chosen automatically
    p_values : list, optional
        if plotting conditional CDF, by default [0.95, 0.99, 0.999]
    plot_joint_pmf : bool, optional
        True to plot the joint probability mass function in a new figure, by default False
    plot_marginal_abundances : bool, optional
        True to plot the distribution of abundances in a new figure, by default False
    plot_conditional_pmf : bool, optional
        True to plot the conditional probability mass function, by default False
    plot_conditional_cdf : bool, optional
        True to plot the conditional cumulative distribution function, by default False
    plot : bool, optional
        True as shorthand for ``plot_joint_pmf = plot_marginal_abundances = plot_conditional_pmf = plot_conditional_cdf = True``
    verbose : bool, optional
        True to print progress information


    Returns
    -------
    Callable[[float,float], float]
        callable which returns P(log10(enrichment) >= y | log10(abundance) = x)
    """

    if plot:
        plot_joint_pmf = True
        plot_marginal_abundances = True
        plot_conditional_pmf = True
        plot_conditional_cdf = True
    if plot_joint_pmf or plot_marginal_abundances or plot_conditional_pmf or plot_conditional_cdf:
        import matplotlib.pyplot as plt
        if fig_kw is None: fig_kw = dict()

    # df = simulate_enrichment(ft_lib_cdr3)
    from scipy.stats import gaussian_kde

    def midpoints(x):
        return (x[1:] + x[:-1]) / 2

    # fit a gaussian kernel density estimator to log10(abundance) vs log10(enrichment)
    # ldf = np.log10(df[['abundance', 'enrichment']]).dropna()

    if verbose:
        print("Fitting joint KDE...")

    ldf = np.log10(df[['abundance', 'enrichment']]).replace(
        [np.inf, -np.inf], 0).dropna()
    joint_kde = gaussian_kde(ldf.values.T)

    # evaluate the KDE at a grid of abundance/enrichment values to calculate the joint probability distribution
    if abundance_bins is None:
        abundance_bins = np.linspace(-6, -1, 6*40)
    if enrichment_bins is None:
        enrichment_bins = np.linspace(-1, 4, 4*40)

    # evaluate at the midpoint of each bin so we can use plt.colormesh to define the boundaries
    abundances = midpoints(abundance_bins)
    enrichments = midpoints(enrichment_bins)

    # estimate the joint probability at each point using the KDE
    if verbose:
        print("Estimating joint probabilities on grid from KDE...")
    xx, yy = np.meshgrid(abundances, enrichments)
    zz = joint_kde([xx.ravel(), yy.ravel()])
    joint_probabilities = zz.reshape(xx.shape) / zz.sum()

    if plot_joint_pmf:
        fig, ax = plt.subplots(**fig_kw)
        im = ax.pcolormesh(abundance_bins, enrichment_bins,
                           joint_probabilities)
        ax.set_xlabel('log10(abundance)')
        ax.set_ylabel('log10(enrichment)')
        ax.set_title("Joint PMF")
        fig.colorbar(im, ax=ax, label='p(abundance = x & enrichment = y)')

    # the marginal probability (e.g. P(abundance = x)) is the sum over enrichments of the joint probability
    marginal_probabilities = joint_probabilities.sum(axis=0)

    if plot_marginal_abundances:
        fig, ax = plt.subplots(**fig_kw)
        ax.plot(abundances, marginal_probabilities)
        ax.set_title("Marginal PMF (log10(abundance))")
        ax.set_xlabel("log10(abundance)")
        ax.set_ylabel("P(log10(abundance) = x)")

    # P(A|B) = P(A & B)/P(B)
    # i.e. conditional probability is joint probability / marginal probability
    conditional_probabilities = joint_probabilities / marginal_probabilities

    if plot_conditional_pmf:
        fig, ax = plt.subplots(**fig_kw)
        im = ax.pcolormesh(abundance_bins, enrichment_bins,
                           conditional_probabilities)
        ax.set_xlabel('log10(abundance)')
        ax.set_ylabel('log10(enrichment)')
        ax.set_title("Conditional distribution")
        fig.colorbar(im, ax=ax, label='p(enrichment = y | abundance = x)')

    conditional_cdf = conditional_probabilities.cumsum(axis=0).cumsum(axis=1)
    conditional_cdf = conditional_cdf / conditional_cdf.sum()

    if plot_conditional_cdf:
        fig, ax = plt.subplots(**fig_kw)
        im = ax.pcolormesh(abundance_bins, enrichment_bins,
                           conditional_probabilities.cumsum(axis=0))
        ax.set_xlabel('log10(abundance)')
        ax.set_ylabel('log10(enrichment)')
        ax.set_title("Conditional CDF")
        fig.colorbar(im, ax=ax, label='p(enrichment <= y | abundance = x)')

        # for each abundance, calculate the minimum enrichment y such that P(enrichment <= y | abundance = x) = p
        for p in p_values:
            ax.plot(abundances, enrichments[(
                conditional_probabilities.cumsum(axis=0) > p).argmax(axis=0)], label=f'p = {p}')
        ax.legend()

    if verbose:
        print("Fitting 2D spline to interpolate ECDF...")


    if force_interp2d:
        try:
            from scipy.interpolate import interp2d
            cond_ecdf = interp2d(abundances, enrichments,
                                conditional_probabilities.cumsum(axis=0))
            return cond_ecdf
        except (ImportError):
            pass

    from scipy.interpolate import RectBivariateSpline
    cond_ecdft = RectBivariateSpline(
        abundances, enrichments, conditional_probabilities.cumsum(axis=0).T)

    # cond_ecdf = lambda xnew, ynew: cond_ecdft(xnew, ynew).T
    cond_ecdf = SplineEvaluator(cond_ecdft, transpose=True)

    return cond_ecdf

class SplineEvaluator():
    def __init__(self, spline, transpose=True):
        self.spline = spline
        self.transpose = transpose
    
    def __call__(self, *args, grid=False, **kwargs):
        if self.transpose:
            return self.spline(*args, grid=grid, **kwargs).T
        else:
            return self.spline(*args, grid=grid, **kwargs)



# def cond_pvalue(abundance, enrichment):
#     return 1 - cond_ecdf(np.log10(abundance),np.log10(enrichment))

# nbseq.utils.mkdirp('results/tables/cdr3/enrichment/null')
# with open('results/tables/cdr3/enrichment/null/ecdf.pickle', 'wb') as f:
#     pickle.dump(cond_ecdf, f)


# sig_df = pd.DataFrame({'abundance': 10**abundances, **{ str(p): 10**enrichments[(conditional_probabilities.cumsum(axis=0) > p).argmax(axis=0)] for p in ps }})
# sig_df
