import pandas as pd
import numpy as np

from .ft import dataframe_to_anndata, fortify, query as ft_query
import anndata
from anndata import AnnData

"""
convert feature tables to design matrices for statistical inference/machine learning

a feature table is an NR x P matrix, where
- N is the number of selection conditions
- R is the number of rounds of selection
- P is the number of VHH features (e.g. NA sequences, AA sequences, CDR3s, etc.)

a design matrix is an N x F matrix, where F is a number of features. F may be equal to P, but not necessarily.

rows of the design matrix can be registered to the antigen matrix for CCA, fitting classifiers, etc.
"""


def concat_rounds_colwise(ft, group_by=['name'], comparison_col='round', fillna=None, identifier=None, format='anndata'):
    """
    accepts a NR x P feature table where rows are selection condition-rounds, columns are VHH features; produces an N x P matrix where rows are sample-conditions, columns are VHH feature-rounds
    """
    if identifier is None:
        identifier = ft.var_names.name

    df = fortify(ft, obs=True).sort_values(by=[comparison_col,identifier])

    if fillna is not None:
        df['abundance'] = df['abundance'].fillna(fillna)

    id_comparison_col = f"{identifier}_{comparison_col}"
    df[id_comparison_col] = df[identifier] + '_' + df[comparison_col]

    # no matter what else I do, pandas insists on sorting the columns in the pivot_table step. apply this to preserve the desired sorting order
    feature_order = df[id_comparison_col].drop_duplicates().values

    if format == 'pandas':
        dfp = pd.pivot_table(df, index=group_by, columns=id_comparison_col, values='abundance', sort=False)
        return dfp[feature_order]
    elif format == 'anndata':
        obs_col = '-'.join(group_by)
        df[obs_col] = df[group_by].agg('-'.join, axis=1)
        dfp = dataframe_to_anndata(df, obs_col=obs_col, var_col=id_comparison_col, value_col='abundance')
        return dfp[:,feature_order]

def last_round(ft, group_by=['name'], comparison_col='round', round='R5i', **kwargs):
    """
    produce a design matrix where rows are selection conditions and columns are VHH features at the last round of selection
    """
    ft = ft_query(ft, f"{comparison_col} == '{round}'", axis='obs')
    return concat_rounds_colwise(ft, group_by=group_by, comparison_col=comparison_col, **kwargs)

def enrichment(ft, **kwargs):
    """
    accepts a NR x P feature table, produces an N x P matrix representing the change in abundance of each nanobody feature for each round
    """
    from .select import calculate_enrichment_ad
    ft_enrichment = calculate_enrichment_ad(ft)
    return enrichment_ft_to_design(ft_enrichment, **kwargs)


def enrichment_ft_to_design(ft_enrichment, fillna=None, group_by=['name'], identifier=None, format='anndata'):
    if identifier is None:
        identifier = ft_enrichment.var_names.name

    df_enrichment = fortify(ft_enrichment, abundance_col='enrichment')
    df_enrichment[identifier] = df_enrichment[identifier] + '_enr'
    if fillna is not None:
        df_enrichment['enrichment'] = df_enrichment['enrichment'].fillna(fillna)

    if format == 'pandas':
        return pd.pivot_table(df_enrichment, index=group_by, columns=identifier, values='enrichment', sort=False)
    elif format == 'anndata':
        obs_col = '-'.join(group_by)
        df_enrichment[obs_col] = df_enrichment[group_by].agg('-'.join, axis=1)
        return dataframe_to_anndata(df_enrichment, obs_col=obs_col, var_col=identifier, value_col='enrichment')


def concat_designs(designs, axis=1, **kwargs):
    if isinstance(designs[0], pd.DataFrame):
        return pd.concat(designs, axis=axis)
    elif isinstance(designs[0], AnnData):
        return anndata.concat(designs, axis=axis, **kwargs)


def design(ft, formula):
    """
    accepts a formula in patsy format
    """
    pass


def subset_design_for_label(design, label_matrix, label, shuffle=False, densify=False):
    """ subset rows of a design matrix according to whether the corresponding row of `label_matrix`[`label`] is NaN or not.
    Used for identifying which rows of `design` can be used to train a model for the metadata variable `label` (as opposed to those rows where `label` is unknown)

    if shuffle=True or shuffle='labels', randomly permute the values of `label_matrix[label]`
    if shuffle='rows', randomly choose rows from among `design`

    returns: (train_idx, X_train, X_unknown, y)
    """
    from scipy.sparse import issparse

    y_values = label_matrix[label]
    train_idx = ~y_values.isna()
    y = (y_values[train_idx] > 0).astype(int)
    if shuffle == True or shuffle == 'labels':
        y = y.sample(frac=1)
    elif shuffle == 'rows':
        #np.random.shuffle(train_idx)
        train_idx = train_idx.sample(frac=1)

    if isinstance(design, pd.DataFrame):
        X = design.values
    elif isinstance(design, AnnData):
        X = design.X

    X_train   = X[ train_idx,:]
    X_unknown = X[~train_idx,:]

    if densify and issparse(X):
        X_train = X_train.toarray()
        X_unknown = X_unknown.toarray()

    return (train_idx, X_train, X_unknown, y)

def get_design_for_ag(design, ag_matrix, ag, shuffle=False, **kwargs):
    return subset_design_for_label(design, ag_matrix, ag, shuffle=shuffle, **kwargs)
