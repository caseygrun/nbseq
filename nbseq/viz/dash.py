# from altair_transform import transform_chart
import numpy as np
import pandas as pd
import altair as alt
from altair import datum
import matplotlib.pyplot as plt
import panel as pn
import scipy.stats.mstats

from ..utils import sample_metadata_to_selection_metadata, clamp_finite, get_rounds, replace_multiple
from ..ft import dataframe_to_anndata, to_relative, fortify, fortify_top_asvs, query as ft_query
from .. import select as nbselect
from ..asvs import get_identifier
from .utils import (
    hash_to_color, pretty_hex, pack_hex, hash_to_mn_short,
    buffer_to_figure, figure_to_buffer, sparkline, figure_to_sparkline,
    mini_histogram, plot_mini_histogram, plot_abundance_sparkline,
    rug_bar,
)

from .sample import alt_scale_features



# # monkey patch panel.pane.vega._get_selections, because now in vega-lite v5/altair v5,
# # selections are replaced with more general parameters. Surprisingly with this hack,
# # panel's handling otherwise still works!
# import panel.pane.vega
# from panel.pane.vega import _get_type, _containers, _isin
# def _get_selections(obj):
#     # import pdb; pdb.set_trace()
#     selections = {}
#     if _isin(obj, 'params'):
#         try:
#             selections.update({
#                 spec.name: _get_type(spec)
#                 for spec in obj['params']
#             })
#         except (AttributeError, TypeError):
#             pass
#     for c in _containers:
#         if _isin(obj, c):
#             for subobj in obj[c]:
#                 selections.update(_get_selections(subobj))
#     return selections
# panel.pane.vega._get_selections = _get_selections




def selection_group_vhh_overview(df, point_selector=None, feature_scale=None, tree=None):

    selector = alt.selection_interval(name='brush', encodings=['x', 'y'], resolve='intersect',
                                      on="[mousedown[event.altKey], mouseup] > mousemove",
                                      translate="[mousedown[event.altKey], mouseup] > mousemove",)

    if point_selector is None:
        point_selector = alt.selection_point(
            fields=['feature'], on="click[!event.altKey]", clear="dblclick[!event.altKey]")

    # selector = alt.param(name='brush', select=dict(type='interval', encodings=['x','y'], fields=['feature']))
    # selector = alt.selection_point(fields=['feature'])

    if feature_scale is None:
        feature_scale = alt_scale_features(df['feature'].unique(), sort=True)
    _opacity = alt.condition(selector & point_selector,
                             alt.value(1), alt.value(0.2))
    _size = alt.condition(selector & point_selector,
                          alt.value(100), alt.value(20))

    df['pct_prod_nlogp_neg'] = clamp_finite(df['pct_prod_nlogp_neg'])
    df['pct_prod_nlogp_pos'] = clamp_finite(df['pct_prod_nlogp_pos'])
    df = df.query("binary_pvalue < 1")


    # binary enrichment
    chart1 = alt.Chart(df).mark_point().encode(
        x=alt.X(
            # 'f_samples_pos_sig:Q',
            'f_samples_pos_sig_jitter:Q',
                title='VHH enriched in %samples', 
                axis=alt.Axis(format='%'),
                scale=alt.Scale(padding=+2)
                ),
        y=alt.Y(
            # 'binary_nlogp:Q', 
            'binary_nlogp_jitter:Q', 
            title='-log10(binary p-value)', 
            scale=alt.Scale(padding=+2),
            ),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('f_samples_pos_sig', format='.2g'),
                 alt.Tooltip('n_samples_pos_sig', format='.0f'),
                 alt.Tooltip('binary_pvalue', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_pos', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_neg', format='.2g'), 
        ]
    ).properties(width=200, height=200)

    # pct product p-value, pos/neg
    chart2 = alt.Chart(df).mark_point().encode(
        # x=alt.X('pct_prod_nlogp_pos:Q', scale=alt.Scale(clamp=True)),
        x=alt.X('pct_sum_pos:Q', title="percentile sum (antigen+ samples)", scale=alt.Scale(clamp=True)),
        # y=alt.Y('pct_prod_nlogp_neg:Q', scale=alt.Scale(clamp=True)),
        y=alt.Y('pct_sum_neg:Q', title="percentile sum (antigen- samples)", scale=alt.Scale(clamp=True)),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('pct_prod_pos', format='.2g'),
                 alt.Tooltip('pct_sum_pos', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_pos', format='.2g'),
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('pct_prod_neg', format='.2g'),
                 alt.Tooltip('pct_sum_neg', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_neg', format='.2g'),
                 alt.Tooltip('n_samples_finite_neg', format='.0f'),
                 ],
    ).properties(width=200, height=200)

    # mean enrichment, pos/neg
    chart3 = alt.Chart(df).mark_point().encode(
        # scale=alt.Scale(clamp=True)),
        x=alt.X('mean_enrichment_pos:Q', title="geometric mean enrichment (antigen+ samples)", scale=alt.Scale(type='log')),
        y=alt.Y('mean_enrichment_neg:Q', title="geometric mean enrichment (antigen- samples)", scale=alt.Scale(type='log')),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('mean_enrichment_pos', format='.2g'),
                 alt.Tooltip('n_samples_finite_neg', format='.0f'),
                 alt.Tooltip('mean_enrichment_neg', format='.2g')],
    ).properties(width=200, height=200)

    # start/end abundance
    _extents = (
        df[['mean_start', 'mean_end']].min().min(), 
        df[['mean_start', 'mean_end']].max().max()
    )
    chart4 = alt.Chart(df).mark_point().encode(
        x=alt.X('mean_start:Q', title="geometric mean starting abundance", 
                scale=alt.Scale(type='log', domain=_extents,
                padding=+2), axis=alt.Axis(format='.0e')),  # scale=alt.Scale(clamp=True)),
        y=alt.Y('mean_end:Q',   title="geometric mean final abundance", 
                scale=alt.Scale(type='log', domain=_extents, padding=+2), 
                axis=alt.Axis(format='.0e')),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('mean_start', format='.2g'),
                 alt.Tooltip('mean_end', format='.2g'),
                 alt.Tooltip('mean_enrichment_pos', format='.2g')],
    ).properties(width=200, height=200)

    charts = [chart1, chart2, chart3, chart4]


    chart = (
        alt.hconcat(*charts)
        .transform_calculate(
            f_samples_pos_sig_jitter="datum.f_samples_pos_sig+0.02*random()",
            binary_nlogp_jitter="datum.binary_nlogp+0.1*random()"
        )
    )


    if tree is not None:
        from .tree import plot_tree_alt, fortify_tree
        nodes_df, links_df = fortify_tree(
            df[['feature','mn', 'n_samples_finite_pos', 'f_samples_pos_sig', 'binary_nlogp', 'seq']].set_index('feature'), 
            tree, identifier='feature', features=df['feature'])

        # for some reason, if NAs are allowed for these fields, it causes chart1 to 
        # have broken scales. probably related to an infinite value introduced in the
        # jitter transform
        nodes_df[['f_samples_pos_sig', 'binary_nlogp']] = nodes_df[['f_samples_pos_sig', 'binary_nlogp']].fillna(0)
        chart5 = plot_tree_alt(nodes_df, links_df, identifier='feature', label='hex',
            feature_scale=feature_scale,
            nodes=dict(
                opacity=_opacity, #alt.condition(point_selector, alt.value(1), alt.value(0.1)),
                # tooltip=['feature'],
                tooltip=['feature','mn', alt.Tooltip('n_samples_finite_pos', format='.0f'), 'seq']
            )
        ).interactive()
        # chart5 = plot_tree_alt(
        #     tree, 
        #     features = df['feature'], 
        #     identifier = 'feature', 
        #     feature_data = df[['feature','mn', 'n_samples_finite_pos', 'f_samples_pos_sig', 'binary_nlogp']].set_index('feature'),
        #     feature_scale=feature_scale,
        #     nodes=dict(
        #         opacity=_opacity, #alt.condition(point_selector, alt.value(1), alt.value(0.1)),
        #         # tooltip=['feature'],
        #         tooltip=['feature','mn', alt.Tooltip('n_samples_finite_pos', format='.0f')]
        #     )
        # )
        # charts.append(chart5)
        chart = alt.hconcat(chart, chart5)

    chart = chart.add_params(selector, point_selector)
    return chart


def selection_group_vhh_overview_simple(df, point_selector=None, feature_scale=None, tree=None):

    selector = alt.selection_interval(name='brush', encodings=['x', 'y'], resolve='intersect',
                                      on="[mousedown[event.altKey], mouseup] > mousemove",
                                      translate="[mousedown[event.altKey], mouseup] > mousemove",)

    if point_selector is None:
        point_selector = alt.selection_point(
            fields=['feature'], on="click[!event.altKey]", clear="dblclick[!event.altKey]")

    # selector = alt.param(name='brush', select=dict(type='interval', encodings=['x','y'], fields=['feature']))
    # selector = alt.selection_point(fields=['feature'])

    if feature_scale is None:
        feature_scale = alt_scale_features(df['feature'].unique(), sort=True)
    _opacity = alt.condition(selector & point_selector,
                             alt.value(1), alt.value(0.2))
    _size = alt.condition(selector & point_selector,
                          alt.value(100), alt.value(20))

    df['pct_prod_nlogp_neg'] = clamp_finite(df['pct_prod_nlogp_neg'])
    df['pct_prod_nlogp_pos'] = clamp_finite(df['pct_prod_nlogp_pos'])
    df = df.query("binary_pvalue < 1")


    # binary enrichment
    chart1 = alt.Chart(df).mark_point().encode(
        x=alt.X(
            # 'f_samples_pos_sig:Q',
            'f_samples_pos_sig_jitter:Q',
                title='VHH enriched in %samples', 
                axis=alt.Axis(format='%'),
                scale=alt.Scale(padding=+2)
                ),
        y=alt.Y(
            # 'binary_nlogp:Q', 
            'binary_nlogp_jitter:Q', 
            title='-log10(binary p-value)', 
            scale=alt.Scale(padding=+2),
            ),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('f_samples_pos_sig', format='.2g'),
                 alt.Tooltip('n_samples_pos_sig', format='.0f'),
                 alt.Tooltip('binary_pvalue', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_pos', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_neg', format='.2g'), 
        ]
    ).properties(width=200, height=200)

    # pct product p-value, pos/neg
    chart2 = alt.Chart(df).mark_point().encode(
        # x=alt.X('pct_prod_nlogp_pos:Q', scale=alt.Scale(clamp=True)),
        x=alt.X('pct_sum_pos:Q', title="percentile sum (antigen+ samples)", scale=alt.Scale(clamp=True)),
        # y=alt.Y('pct_prod_nlogp_neg:Q', scale=alt.Scale(clamp=True)),
        y=alt.Y('pct_sum_neg:Q', title="percentile sum (antigen- samples)", scale=alt.Scale(clamp=True)),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('pct_prod_pos', format='.2g'),
                 alt.Tooltip('pct_sum_pos', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_pos', format='.2g'),
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('pct_prod_neg', format='.2g'),
                 alt.Tooltip('pct_sum_neg', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_neg', format='.2g'),
                 alt.Tooltip('n_samples_finite_neg', format='.0f'),
                 ],
    ).properties(width=200, height=200)

    charts = [chart1, chart2]


    chart = (
        alt.hconcat(*charts)
        .transform_calculate(
            f_samples_pos_sig_jitter="datum.f_samples_pos_sig+0.02*random()",
            binary_nlogp_jitter="datum.binary_nlogp+0.1*random()"
        )
    )


    chart = chart.add_params(selector, point_selector)
    return chart


def feature_traceplot(
        df, detail='feature', tooltip=['feature', 'name', 'description'],
        features=None,
        selector=None, feature_scale=None):
    df
    # feature_scale = alt.Scale(domain=_features, range=[hash_to_color(f) for f in _features])

    _opacity = alt.value(1)
    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.1))

    if features is None:
        features = df['feature'].unique()

    if feature_scale is None:
        feature_scale = alt_scale_features(features)

    traceplot = alt.Chart(df).mark_line().encode(
        x=alt.X('r:O', title=None),
        y=alt.Y('abundance:Q'),
        color=alt.Color('feature:N', scale=feature_scale, sort="ascending"),
        opacity=_opacity,
        detail=detail,
        # tooltip=[
        #     'feature:N',
        #     'mn:N',
        #     'description:N',
        #     alt.Tooltip('enrichment:Q', format=".2g"),
        #     alt.Tooltip('log_enrichment:Q', format=".2f")
        # ],

        facet=alt.Facet('name:O', 
            title=None, 
            # header=alt.Header(labels=False, labelPadding=0, labelFontSize=0)
        )
    ).properties(width=100, height=100).transform_filter(datum.feature != 'others').interactive()

    return traceplot


def feature_barplot(df, x='r:O', selector=None,
                    features=None, feature_scale=None, phylo=False, fancy_sorting=False,  **kwargs):
    """barplot of most abundant feature

    Parameters
    ----------
    df : pd.DataFrame
        fortified from a feature table, with columns:
        - 'feature'
        - 'abundance'
        - 'r' (or whatever is passed as `x`)
    x : str
        field expression for altair to identify which column of `df` to use
    selector : alt.Selection or bool
        if an altair/Vega selection is passed, it will be bound to this chart, and used to adjust selected points
        if True, a selection will be created which will allow `feature`s to be selected

    Returns
    -------
    alt.Chart
        altair chart
    """

    if features is None:
        features = df['feature'].unique()
    if feature_scale is None:
        feature_scale = alt_scale_features(features)
    n = len(features)

    _opacity = alt.value(1)
    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.1))

    barplot = alt.Chart(df).mark_bar().encode(
        x=alt.X('r:O', title='round'),
        y='abundance:Q',
        tooltip=[
            'feature:N',
            'mn:N',
            alt.Tooltip('abundance:Q', format=".2g"),
            alt.Tooltip('log_abundance:Q', format=".2f"),
            alt.Tooltip('enrichment:Q', format=".2g"),
        ],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.1)),
        color=alt.Color('feature:N',
                        legend=alt.Legend(
                            columns=n//20, symbolLimit=0, symbolType="circle",
                            labelExpr="slice(datum.value,0,6)", 
                            # sort=sorted(features),
                            # labelExpr="slice(datum.mn,0,6)"
                        ),
                        # sort="ascending",
                        sort = sorted(features),
                        scale=feature_scale
                        ),
        facet=alt.Facet(
            'name:O', 
            title=None, 
            # header=alt.Header(labels=False, labelPadding=0, labelFontSize=0)
        )
    ).properties(width=100, height=200)

    return barplot

def whittaker_plot(df, selector, n_head=100, n_sample=50, n_tail=100):
    dff = pd.concat([
        df.groupby(['ID', 'round']).head(n_head),
        df.groupby(['ID', 'round']).sample(n_sample),
        df.groupby(['ID', 'round']).tail(n_tail)
    ])
    dff['log_rank'] = np.log10(dff['rank'])

    whittakerplot = alt.Chart(dff).mark_point().encode(
        x='rank',
        y='log_abundance',
        color=alt.Color('round', scale=alt.Scale(scheme='viridis')),
        opacity=alt.condition(selector, alt.value(1), alt.value(0.005)),
        tooltip=['feature', 'round', 'rank', 'abundance'],
        facet=alt.Facet('name:O', header=alt.Header(labels=False), title=None)
    ).properties(width=100, height=100).interactive()
    return whittakerplot

def enrichment_rank_plot(df, selector, feature_scale, features=None):
    _opacity = alt.value(1)
    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.05))

    enrrankplot = alt.Chart(df).mark_point().encode(
        x='enrichment_rank',
        y='log_enrichment',
        tooltip=['feature', 'mn', 
            alt.Tooltip('enrichment_rank', format=".0f"),
            alt.Tooltip('log_enrichment', format=".2g"),
            alt.Tooltip('enrichment', format=".2g"),
            # alt.Tooltip('start', format=".2g"),
            # alt.Tooltip('end', format=".2g"),
        ],
        fill=alt.Color('feature:N', scale=feature_scale),
        color=alt.value(None),
        opacity=_opacity,
        facet=alt.Facet('name:O', header=alt.Header(labels=False), title=None)
    ).properties(width=100, height=100)#.interactive()
    return enrrankplot

def enrichment_abundance_plot(df, selector, feature_scale, features=None):
    _opacity = alt.value(1)
    _size = alt.value(20)
    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.05))
        _size = alt.condition(selector, alt.value(100), alt.value(20))

    _tooltips = ['feature', 'mn',
                alt.Tooltip('start', format=".2g"),
                alt.Tooltip('end', format=".2g"),
                alt.Tooltip('enrichment', format=".2g"),
    ]
    if 'enr_p_value' in df.columns:
        _tooltips.append(alt.Tooltip('enr_p_value', format=".2g"))
    if 'log_enrichment' in df.columns:
        _tooltips.append(alt.Tooltip('log_enrichment', format=".2g"))
    if 'enrichment_rank' in df.columns:
        _tooltips.append(alt.Tooltip('enrichment_rank', format=".0f"))

    enrabdplot = alt.Chart(df).mark_point().encode(
        x=alt.X('end', title='end abundance',
                scale=alt.Scale(type='log')),
        # x='end',
        # y='log_enrichment',
        y=alt.Y('enrichment', scale=alt.Scale(type='log')),
        tooltip=_tooltips,
        fill=alt.Color('feature:N', scale=feature_scale, sort='ascending' ),
        color=alt.condition(selector, alt.value('black'), alt.value(None)),
        opacity=_opacity,
        size=_size,
        facet=alt.Facet('name:O', header=alt.Header(labels=False), title=None)
    ).properties(width=100, height=100).interactive()  # .interactive()
    return enrabdplot


# todo
def SelectionGroupData():
    def __init__(self,ft_pos, df_selections, df_samples, df_samples_top, df_features, df_enr_features, sel_fts_q, features):
        self.ft_pos = ft_pos
        self.df_selections = df_selections
        self.df_samples = df_samples
        self.df_samples_top = df_samples_top
        self.df_features = df_features
        self.df_enr_features = df_enr_features
        self.sel_fts_q = sel_fts_q
        self.features = features

    @staticmethod
    def _calculate(*args):
        return SelectionGroupData(*_calculate_selection_group_dashboard_phenotype(*args))


def _calculate_selection_group_dashboard_phenotype(
    ft, ft_enr, df_enr, sel_fts, pos_query, neg_query,
    feature_col, space, phenotype):
    ft_pos = to_relative(
        ft_query(ft, pos_query, axis='sample')
    )

    # one row per sample per feature; for Whittaker plots
    df_samples = fortify(
        ft_pos, feature_col=feature_col, obs=True)
    df_selections = sample_metadata_to_selection_metadata(
        df_samples)

    # process features to calculate enrichment statistics
    _df_bin = nbselect.compare_binary_phenotypes(
        sel_fts['sig'], phenotypes=[phenotype], 
        pos_query=pos_query,
        neg_query=neg_query,
        plot=False, 
        # plot_sparkline=True, 
        plot_sparkline=False,
        n_boots=1000
    )
    _df_bin['nlogp'] = -np.log10(_df_bin['p_value'])
    _df_mean = nbselect.calculate_specificities_mean(
        ft, ft_enr, q1=pos_query, q2=neg_query, suffixes=('_pos', '_neg'))
    _df_pct = nbselect.calculate_specificities_pct_product(
        ft_enr, sel_fts['pct'], q1=pos_query, q2=neg_query, suffixes=('_pos', '_neg'), plot_simulation=False)
    # calculate mean start/end abundance per feature
    _df_start_end = (df_enr
                    .loc[df_enr['name'].isin(df_selections.index), [feature_col, 'start', 'end']]
                    .groupby(feature_col)
                    .agg(scipy.stats.mstats.gmean)
                    .rename(columns=lambda x: f'mean_{x}'))

    # join various metrics together; one row per feature
    df_features = (
        _df_bin.rename(columns={
            'p_value': 'binary_pvalue',
            'nlogp': 'binary_nlogp',
            'n_samples': 'n_samples_pos_sig',
            'f_samples': 'f_samples_pos_sig',
            'plot': 'plot_binary'
        })
        .join(
            _df_pct.rename(columns={
                'pct_product_pos': 'pct_prod_pos',
                'pct_sum_pos': 'pct_sum_pos',
                'p_value_pos': 'pct_prod_pvalue_pos',
                'nlogp_pos': 'pct_prod_nlogp_pos',
                'nsamples_pos': 'n_samples_finite_pos',
                'pct_product_neg': 'pct_prod_neg',
                'pct_sum_neg': 'pct_sum_neg',
                'p_value_neg': 'pct_prod_pvalue_neg',
                'nlogp_neg': 'pct_prod_nlogp_neg',
                'nsamples_neg': 'n_samples_finite_neg',
            }),
            on=feature_col)
        .join(_df_mean, on=feature_col)
        .join(_df_start_end, on=feature_col)
        .join(ft_pos.var, on=feature_col)
        .rename(columns={space.upper(): 'seq'})
    )
    df_features['mn'] = df_features[feature_col].apply(hash_to_mn_short)

    # find relevant features
    features = df_features.query('binary_pvalue < 1')[feature_col]

    # get fortified dataframes
    # one row per sample per feature; relevant features + top n most abundant features
    df_samples_top = (
        fortify_top_asvs(ft, pos_query, n=100, select_from_round=8,
                                other_features=features)
        .join(
            (df_enr
                .set_index(['name', feature_col])
                [['enrichment', 'log_enrichment','enrichment_rank','p_value','start', 'end']]
                .rename(columns={'p_value': 'enr_p_value'})
            ),
            on=['name', 'feature']
        )
    )
    df_samples_top['log_abundance'] = np.log10(
        df_samples_top['abundance'])
    df_samples_top['mn'] = df_samples_top['feature'].apply(
        hash_to_mn_short)

    # todo: do we need these p-values and such?
    df_samples['log_abundance'] = np.log10(df_samples['abundance'])
    df_samples['rank'] = (df_samples.groupby(
        'ID')['log_abundance'].rank(ascending=False, pct=False))
    df_samples = (
        df_samples
        .join(
            _df_bin
            .set_index(feature_col), on=feature_col)
        .join(
            df_enr.set_index(['name', feature_col])
            [['enrichment', 'log_enrichment', 'p_value', 'start', 'end']]
            .rename(columns={'p_value': 'enr_p_value'}),
            on=['name', feature_col]
        )
    )

    # one row per sample per feature, restricted to features of interest
    df_enr_features = df_enr.loc[df_enr[feature_col].isin(features), :]

    # various feature tables indexed by selection, subset by this query
    sel_fts_q = {
        k: ft_query(sel_fts[k], pos_query, axis='obs') for k in sel_fts
    }

    return (ft_pos, df_selections, df_samples, df_samples_top, df_features, df_enr_features, sel_fts_q, features)


def _make_features_table(
        df_samples, df_features, df_enr_features, features, feature_col,  sel_fts_q, 
        n_jobs=1, 
        rounds = ['R5i', 'R6i', 'R7i', 'R8i'],
        show_histograms=True,
        show_traces=True,
        ):
    from ..utils import reorder_columns
    from ..viz.utils import plot_mini_histogram, figure_to_sparkline

    df_samples_abundance = df_samples.sort_values([feature_col, 'name', 'round']).pivot(
        index=[feature_col, 'name'], columns='round', values='abundance')
    _rounds = sorted(list(df_enr_features.columns[df_enr_features.columns.str.startswith('R')]))

    def make_traces():
        def trace_feature(dfg, feature):
            fig, ax = plt.subplots(figsize=(2, 0.5))
            ax.plot(dfg[_rounds].values.T,
                    color='C0', linewidth=0.5)
            ax.plot(df_samples_abundance.loc[(feature, slice(
                None)), :].values.T, color='red', marker='.', linewidth=1)
            return figure_to_sparkline(fig)

        traces = pd.Series(index=features, dtype='object',
                            name=f'abundance_traces')
        for feature, abundances in df_enr_features.groupby(feature_col)[['name'] + _rounds]:
            traces[feature] = trace_feature(abundances, feature)
        return traces

    def make_feature_hist(data, feature, column='log_enrichment'):
        hist_color = hash_to_color(feature)
        fig = plot_mini_histogram(data, color=hist_color)

        # find abundance of samples in query
        marks = sel_fts_q[column][:, feature].X.data

        color = 'red'
        for m in marks:
            fig.axes[0].axvline(m, c=color)
        fig.axes[0].axvline(0, c='black')
        return figure_to_sparkline(fig, axis_off=False, subplots_adjust=True)

    def make_feature_hists(column):
        # from joblib import Parallel, delayed
        
        # _pool = Parallel(n_jobs=n_jobs)
        # def _do_make_feature_hist(feature, abundances):
        #     return (feature, make_feature_hist(
        #         abundances, feature, column=column))

        # _hists = (
        #     _do_make_feature_hist(feature, abundances)
        #     for feature, abundances in df_enr_features.groupby(feature_col)[column]
        # )
        # _hists = _pool(delayed(_do_make_feature_hist)(feature, abundances)
        #       for feature, abundances in df_enr_features.groupby(feature_col)[column])
        
        # idx, values = zip(*_hists)
        # hists = pd.Series(values, idx, dtype='object', name=f'{column}_hist')

        hists = pd.Series(index=features, dtype='object',
                            name=f'{column}_hist')
        for feature, abundances in df_enr_features.groupby(feature_col)[column]:
            hists[feature] = make_feature_hist(
                abundances, feature, column=column)
        return hists

    table = (df_features.query('binary_pvalue < 1'))
    col_order = ['color', 'mn', 'n_samples_pos_sig', 'seq', 'samples',
                    'log_enrichment_hist', f'{rounds[0]}_hist', 'abundance_traces', f'{rounds[-1]}_hist',
                    'binary_pvalue', 'pct_product', 'enrichment', 'abundance']


    if show_histograms:
        table = (table       
                .join(make_feature_hists('log_enrichment'), on=feature_col)
                .join(make_feature_hists(rounds[0]), on=feature_col)
                .join(make_feature_hists(rounds[-1]), on=feature_col)
        )
    else:
        col_order.remove('log_enrichment_hist')
        col_order.remove(f'{rounds[0]}_hist') 
        col_order.remove(f'{rounds[-1]}_hist')

    if show_traces:
        table = table.join(make_traces(), on=feature_col)
    else:
        col_order.remove('abundance_traces')

    table['samples'] = table.apply(lambda row: (
        "⊕: {n_samples_pos_sig:.0f} / {n_samples_finite_pos:.0f} = {f_samples_pos_sig:.0%}, p = {binary_pvalue:.2g}<br/>"
        "⊖: {n_samples_finite_neg:.0f}").format(**row), axis=1)

    table['pct_product'] = table.apply(lambda row: (
        "⊕: {pct_prod_pos:8.2g} in {n_samples_finite_pos:3.0f} samples, p = {pct_prod_pos:.2g}<br/>"
        "⊖: {pct_prod_neg:8.2g} in {n_samples_finite_neg:3.0f} samples, p = {pct_prod_neg:.2g}").format(**row), axis=1)

    table['enrichment'] = table.apply(lambda row: (
        "⊕: {mean_enrichment_pos:8.2g}<br/>"
        "⊖: {mean_enrichment_neg:8.2g}").format(**row), axis=1)

    table['abundance'] = table.apply(lambda row: (
        "{mean_start:8.2g} → <br/>"
        "{mean_end:8.2g}").format(**row), axis=1)

    table['color'] = table[feature_col].apply(pretty_hex)

    from .syntax import aa_highlighter
    table['seq'] = table['seq'].apply(lambda x: aa_highlighter.highlight(x))

    from bokeh.models.widgets.tables import NumberFormatter
    widget = pn.widgets.Tabulator(
        # table[['log_enrichment',
        #        'n_samples', 'f_samples', 'p_value']],
        # reorder_columns(table, col_order).set_index(feature_col),
        table.set_index(feature_col)[col_order],
        titles={
            'n_samples_pos_sig': '⊕*'
        },
        formatters={
            'n_samples_pos_sig': NumberFormatter(format='0'),
            'seq': {'type': 'html'},
            'log_enrichment_hist': {'type': 'html'},
            f'{rounds[0]}_hist': {'type': 'html'},
            'abundance_traces': {'type': 'html'},
            f'{rounds[-1]}_hist': {'type': 'html'},
            'plot_binary': {'type': 'html'},
            'samples': {'type': 'html'},
            'pct_product': {'type': 'html'},
            'enrichment': {'type': 'html'},
            'abundance': {'type': 'html'},
        },
        sorters=[
            {'field': 'n_samples_pos_sig', 'dir': 'desc'}
            # {'field': 'binary_pvalue',     'dir': 'asc'},
        ],
        disabled=True,
        configuration={
            'columnDefaults': {
                'titleFormatter': 'html',
                'formatter': 'html'
            },
        }
    )
    widget.style.applymap(lambda c: f"background-color: #{pack_hex(c)};", subset=['color'])

    return widget


def _make_overview(df_samples_top, df_features, feature_col, tree=None, initial_selection=None):
    df = df_samples_top.rename(columns={feature_col: 'feature'})

    selector = alt.selection_point(name='feature_selection', value=initial_selection, fields=[
        'feature'], bind='legend', on='click', clear='dblclick')
    
    features = df['feature'].unique()
    feature_scale = alt_scale_features(sorted(features), sort=True)

    chart = alt.vconcat(
        selection_group_vhh_overview(
            df_features.rename(columns={feature_col: 'feature'}),
            feature_scale=feature_scale,
            point_selector=selector,
            tree=tree
        ),
        # whittaker_plot(
        #     df_samples.rename(columns={feature_col: 'feature'}), selector=selector),
        alt.vconcat(
            # enrichment_rank_plot(
            enrichment_abundance_plot(
                df, selector=selector,
                features=features, feature_scale=feature_scale),
            feature_traceplot(
                df, selector=selector,
                features=features, feature_scale=feature_scale).transform_filter(selector),
            feature_barplot(
                df, selector=selector,
                features=features, feature_scale=feature_scale)
        ),
        resolve=alt.Resolve(scale={'color': 'independent'})
    ).add_params(selector)  # .configure_header(labels=False, title=None)

    # print(chart.to_json())
    return chart

def _make_overview_simple(df_samples_top, df_features, feature_col, tree=None, initial_selection=None):
    df = df_samples_top.rename(columns={feature_col: 'feature'})

    selector = alt.selection_point(name='feature_selection', value=initial_selection, fields=[
        'feature'], bind='legend', on='click', clear='dblclick')
    
    features = df['feature'].unique()
    feature_scale = alt_scale_features(sorted(features), sort=True)

    chart = selection_group_vhh_overview_simple(
            df_features.rename(columns={feature_col: 'feature'}),
            feature_scale=feature_scale,
            point_selector=selector,
            tree=tree
        ).add_params(selector)  # .configure_header(labels=False, title=None)

    # print(chart.to_json())
    return chart


def selection_group_dashboard(ex, 
                              starting_phenotype=None, 
                              global_query="io == 'i' & kind == '+'",
                              pos_query = "{phenotype} == 1",neg_query = "({phenotype} == 0 | {phenotype} == -1)", 
                              enr_comparison=('start', 'end'), 
                              rounds=None, #['R5i','R6i','R7i','R8i'],
                              tree=True,
                              initial_selection=[],
                              n_jobs=1):

    def setup_selection_group_dashboard(global_query="kind == '+' & io == 'i'", space='cdr3'):
        """build dashboard for a selection group having selected a feature space and sample query"""

        # setup full dataset
        # --------------------------------------------------------
        # subset feature table according to global_query
        feature_col = get_identifier(space)
        ft = ex.query(global_query, space=space, axis='sample')

        nonlocal rounds
        if rounds is None:
            rounds = get_rounds(ft.obs)

        # calculate enrichment df, add p-values and start/end abundance
        df_enr = nbselect.enr(ft, method='df', comparison=enr_comparison,
                                  dropna=True, 
                                  add_start_end=True, 
                                  add_rank=True,
                                  add_log=True)
        df_enr = nbselect.enrichment_pvalues(
            df_enr, abundance_col=enr_comparison[0], inplace=True)
        # df_enr = nbselect.enrichment_start_end(df_enr, inplace=True)
        # df_enr = nbselect.enrichment_rank(df_enr, inplace=True)

        if (feature_col not in df_enr.columns) and ('feature' in df_enr.columns):
            df_enr[feature_col] = df_enr['feature']

        # make several feature tables grouped by selection
        sel_fts = {
            col: dataframe_to_anndata(df_enr, obs=ex.selection_metadata, obs_col='name', var_col=feature_col, value_col=col) for col in ['enrichment', 'log_enrichment', 'sig', 'p_value', rounds[0], rounds[-1], 'start', 'end']
        }
        ft_enr = sel_fts['enr'] = sel_fts['enrichment']
        sel_fts['pct'] = nbselect.calculate_enrichment_pct(ft_enr)

        def setup_selection_group_dashboard_phenotype(
            phenotype,
            pos_query = "{phenotype} == 1",
            neg_query = "({phenotype} == 0 | {phenotype} == -1)"):
            import scipy.stats

            # setup subset
            pos_query = pos_query.format(phenotype=phenotype)
            neg_query = neg_query.format(phenotype=phenotype)

            # check that dataset contains >0 samples matching both pos_query and neg_query
            if (len(ft.obs.query(pos_query))) < 1 or (len(ft.obs.query(neg_query)) < 1):
                return pn.pane.Alert((
                    "Bad phenotype:"
                    f"{len(ft.obs.query(pos_query))} {phenotype}+ samples (query: `{pos_query}`) "
                    f"{len(ft.obs.query(neg_query))} {phenotype}- samples (query: `{neg_query}`). "
                    "Check global query and phenotype value.")
                )

            (ft_pos, df_selections, df_samples, df_samples_top, df_features, df_enr_features,
             sel_fts_q, features) = _calculate_selection_group_dashboard_phenotype(
                ft, ft_enr, df_enr, sel_fts, pos_query, neg_query, feature_col, space, phenotype)


            def make_samples_table():
                return sample_metadata_to_selection_metadata(df_samples)            

            overview = pn.pane.Vega(
                _make_overview(df_samples_top, df_features, feature_col, 
                               tree=(ex.tree[space] if space in ex.tree and tree else None),
                               initial_selection=initial_selection
                               ), 
                debounce=1000, show_actions=True)
            table = pn.panel(
                _make_features_table(df_samples, df_features, df_enr_features,
                                     features, feature_col,  sel_fts_q, n_jobs, 
                                     rounds=rounds), loading_indicator=True)

            # def filter_table(df, selection):
            #     if selection is None or len(selection) == 0:
            #         return df
            #     query = ' & '.join(
            #         f'{crange[0]:.3f} <= `{col}` <= {crange[1]:.3f}'
            #         for col, crange in selection.items()
            #     )
            #     return df.query(query)

            # table.add_filter(
            #     pn.bind(filter_table, selection=overview.selection.param.brush))
            return pn.Column(
                overview, 
                table, 
                sizing_mode='stretch_width')

        nonlocal starting_phenotype
        if starting_phenotype is None:
           starting_phenotype = ex.pheno_names[0]
        phenotype_selector = pn.widgets.AutocompleteInput(
            name='phenotype', options=list(ex.pheno_names), placeholder='Phenotype to compare', 
            value=starting_phenotype
        )

        pos_query_input = pn.widgets.TextInput(
            value=pos_query,
            name="positive query"
        )

        neg_query_input = pn.widgets.TextInput(
            value=neg_query,
            name="negative query"
        )

        return pn.Column(
            pn.Row(
                phenotype_selector,
                pos_query_input,
                neg_query_input,
                sizing_mode='stretch_width'
            ),
            pn.panel(
                pn.bind(setup_selection_group_dashboard_phenotype,
                    phenotype_selector,
                    pos_query_input,
                    neg_query_input
                ),
                loading_indicator=True
            ),
            sizing_mode='stretch_width'
        )

    global_query_input = pn.widgets.TextInput(
        name='global_query', value=global_query,
        placeholder='Filter entire dataset')

    space_input = pn.widgets.Select(
        name='space', options=['cdr3', 'aa'], value='cdr3'
    )



    # neg_query = pn.widgets.TextInput(
    #     name='neg_query', value="expt == '027i' & io == 'i' & kind == '+'"
    #     placeholder='Select antigen-negative samples')

    return pn.Column(
        pn.Row(
            global_query_input,
            space_input,
            sizing_mode='stretch_width'
        ),
        # pn.bind(summarize_clonotype, autocomplete),
        # pn.bind(make_overview, autocomplete),
        # pn.bind(make_table, autocomplete),
        pn.panel(
            pn.bind(setup_selection_group_dashboard,
                    global_query_input, space_input),
            loading_indicator=True),
        sizing_mode='stretch_width'
    ).servable()

def selection_group_dashboard_simple(ex, 
                              starting_phenotype=None, 
                              global_query="expt == '027j' & io == 'i' & kind == '+'",
                              pos_query = "{phenotype} == 1",neg_query = "({phenotype} == 0 | {phenotype} == -1)", 
                              enr_comparison=('start', 'end'), 
                              tree=True,
                              initial_selection=[],
                              n_jobs=1):

    def setup_selection_group_dashboard(global_query="kind == '+' & io == 'i'", space='cdr3'):
        """build dashboard for a selection group having selected a feature space and sample query"""

        # setup full dataset
        # --------------------------------------------------------
        # subset feature table according to global_query
        feature_col = get_identifier(space)
        ft = ex.query(global_query, space=space, axis='sample')

        # calculate enrichment df, add p-values and start/end abundance
        df_enr = nbselect.enr(ft, method='df', comparison=enr_comparison,
                                  dropna=True, 
                                  add_start_end=True, 
                                  add_rank=True,
                                  add_log=True)
        df_enr = nbselect.enrichment_pvalues(
            df_enr, abundance_col=enr_comparison[0], inplace=True)
        # df_enr = nbselect.enrichment_start_end(df_enr, inplace=True)
        # df_enr = nbselect.enrichment_rank(df_enr, inplace=True)

        # make several feature tables grouped by selection
        sel_fts = {
            col: dataframe_to_anndata(df_enr, obs=ex.selection_metadata, obs_col='name', var_col=feature_col, value_col=col) for col in ['enrichment', 'log_enrichment', 'sig', 'p_value', 'R5i', 'R8i', 'start', 'end']
        }
        ft_enr = sel_fts['enr'] = sel_fts['enrichment']
        sel_fts['pct'] = nbselect.calculate_enrichment_pct(ft_enr)

        def setup_selection_group_dashboard_phenotype(
            phenotype,
            pos_query = "{phenotype} == 1",
            neg_query = "({phenotype} == 0 | {phenotype} == -1)"):
            import scipy.stats

            # setup subset
            pos_query = pos_query.format(phenotype=phenotype)
            neg_query = neg_query.format(phenotype=phenotype)

            # check that dataset contains >0 samples matching both pos_query and neg_query
            if (len(ft.obs.query(pos_query))) < 1 or (len(ft.obs.query(neg_query)) < 1):
                return pn.pane.Alert((
                    "Bad phenotype:"
                    f"{len(ft.obs.query(pos_query))} {phenotype}+ samples (query: `{pos_query}`) "
                    f"{len(ft.obs.query(neg_query))} {phenotype}- samples (query: `{neg_query}`). "
                    "Check global query and phenotype value.")
                )

            (ft_pos, df_selections, df_samples, df_samples_top, df_features, df_enr_features,
             sel_fts_q, features) = _calculate_selection_group_dashboard_phenotype(
                ft, ft_enr, df_enr, sel_fts, pos_query, neg_query, feature_col, space, phenotype)


            def make_samples_table():
                return sample_metadata_to_selection_metadata(df_samples)            

            overview = pn.pane.Vega(
                _make_overview_simple(df_samples_top, df_features, feature_col, 
                               tree=(ex.tree[space] if space in ex.tree and tree else None),
                               initial_selection=initial_selection
                               ), 
                debounce=1000, show_actions=True)

            # def filter_table(df, selection):
            #     if selection is None or len(selection) == 0:
            #         return df
            #     query = ' & '.join(
            #         f'{crange[0]:.3f} <= `{col}` <= {crange[1]:.3f}'
            #         for col, crange in selection.items()
            #     )
            #     return df.query(query)

            # table.add_filter(
            #     pn.bind(filter_table, selection=overview.selection.param.brush))
            return pn.Column(overview, sizing_mode='stretch_width')

        nonlocal starting_phenotype
        if starting_phenotype is None:
           starting_phenotype = ex.pheno_names[0]
        phenotype_selector = pn.widgets.AutocompleteInput(
            name='phenotype', options=list(ex.pheno_names), placeholder='Phenotype to compare', 
            value=starting_phenotype
        )

        pos_query_input = pn.widgets.TextInput(
            value=pos_query,
            name="positive query"
        )

        neg_query_input = pn.widgets.TextInput(
            value=neg_query,
            name="negative query"
        )

        return pn.Column(
            pn.Row(
                phenotype_selector,
                pos_query_input,
                neg_query_input,
                sizing_mode='stretch_width'
            ),
            pn.panel(
                pn.bind(setup_selection_group_dashboard_phenotype,
                    phenotype_selector,
                    pos_query_input,
                    neg_query_input
                ),
                loading_indicator=True
            ),
            sizing_mode='stretch_width'
        )

    global_query_input = pn.widgets.TextInput(
        name='global_query', value=global_query,
        placeholder='Filter entire dataset')

    space_input = pn.widgets.Select(
        name='space', options=['cdr3', 'aa'], value='cdr3'
    )



    # neg_query = pn.widgets.TextInput(
    #     name='neg_query', value="expt == '027i' & io == 'i' & kind == '+'"
    #     placeholder='Select antigen-negative samples')

    return pn.Column(
        pn.Row(
            global_query_input,
            space_input,
            sizing_mode='stretch_width'
        ),
        # pn.bind(summarize_clonotype, autocomplete),
        # pn.bind(make_overview, autocomplete),
        # pn.bind(make_table, autocomplete),
        pn.panel(
            pn.bind(setup_selection_group_dashboard,
                    global_query_input, space_input),
            loading_indicator=True),
        sizing_mode='stretch_width'
    ).servable()




def enrichment_binary_volcano_plot(df, selector=None, feature_scale=None, legend=None):
    _opacity = alt.value(1)
    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.05))
    if feature_scale is None:
        feature_scale = alt_scale_features(df['feature'].unique())
    return alt.Chart(df.query("binary_pvalue < 1")).mark_point().encode(
        x=alt.X(
            'f_samples_pos_sig_jitter:Q',
                title='VHH enriched in %samples', axis=alt.Axis(format='%')),
        y=alt.Y(
            'binary_nlogp_jitter:Q', 
            title='-log10(binary p-value)'
            ),
        fill=alt.Color('feature:N', scale=feature_scale, legend=legend),
        color=alt.value(None),
        opacity=_opacity,
        tooltip=['feature', 'f_samples_pos_sig', 'n_samples_pos_sig',
                 'binary_pvalue', 'pct_prod_nlogp_pos', 'pct_prod_nlogp_neg']
    ).transform_calculate(
        f_samples_pos_sig_jitter="datum.f_samples_pos_sig+0.02*random()",
        binary_nlogp_jitter="datum.binary_nlogp+0.1*random()"
    )


# ************************************************

def plot_dist_buffer(data, fig=None):
    fig = plot_mini_histogram(data, figsize=(2, 0.5))
    return figure_to_buffer(fig)


def mark_enr_hist(row, buffers, color=None):
    # fig = buffer_to_figure(buffers.loc[row.loc['name']])
    fig = buffer_to_figure(buffers.loc[row['name']])
    if color is None:
        color = 'red'
    fig.axes[0].axvline(row.loc['log_enrichment'], color=color)
    return figure_to_sparkline(fig, axis_off=False, subplots_adjust=False)


from .utils import shorten_descriptions


class FeatureData():
    def __init__(self, ft, space, phenotypes, df_enr, df_sample, df_selection):
        self.ft = ft
        self.space = space
        self.phenotypes = phenotypes
        self.df_enr = df_enr
        self.df_sample = df_sample
        self.df_selection = df_selection

    def subset(self, feature):
        identifier = get_identifier(self.space)
        return FeatureData(ft=self.ft, 
                           space=self.space,
                           phenotypes=self.phenotypes, 
                           df_enr=self.df_enr.query(f"{identifier} == '{feature}'"),
                           df_sample=self.df_sample,
                           df_selection=self.df_selection)

    @staticmethod
    def from_experiment(ex, global_query, space='cdr3', enr_comparison = ('start', 'end'), abbreviations={}):
        identifier = get_identifier(space)

        ft = ft_query(ex.rfts[space], global_query, axis='obs')
        df_enr = nbselect.enr(ft,'df',comparison=enr_comparison, add_start_end=True, add_log=True, dropna=True)
        df_enr = (
            nbselect.enrichment_pvalues(df_enr, abundance_col=enr_comparison[0], inplace=True)
            .rename(columns={'p_value': 'enr_p_value'})
        )
        # df_enr = nbselect.enrichment_start_end(df_enr, inplace=True)


        df_sample = (fortify(ft, relative=True, obs=True)
                .join(
                    df_enr.set_index(['name',identifier])[['enrichment','log_enrichment','start','end','enr_p_value']], 
                    on=['name',identifier])
            )

        df_sample = shorten_descriptions(df_sample, replacements=abbreviations)

        df_selection = sample_metadata_to_selection_metadata(df_sample)

        return FeatureData(ft=ft, space=space, phenotypes=ex.pheno_names, df_enr=df_enr, df_sample=df_sample, df_selection=df_selection)

def feature_selection_overview(feature, data, starting_phenotype=None):
    """show some plots summarizing the behavior of a feature across multiple selections

    Parameters
    ----------
    feature : str
        feature identifier
    data : FeatureData
        data bundle describing features

    Returns
    -------
    alt.Chart
        altair chart
    """
    data = data.subset(feature)
    df = data.df_enr
    identifier = get_identifier(data.space)

    # df = df_enr.query(f"CDR3ID == '{feature}'")

    # scale for coloring phenotype marks
    _phenotype_scale = alt.Scale(
        domain=['?', -1, 0, 1], range=['#bab0ac', '#4c78a8', '#72b7b2', '#f58518'])

    # enrichment histograms
    # ---------------------

    # allow cross-filtering by selecting histogram
    brush = alt.selection_interval(name='brush',
                                    on="[mousedown[event.altKey], mouseup] > mousemove",
                                    translate="[mousedown[event.altKey], mouseup] > mousemove")
    dropdown = alt.binding_select(
        options=[identifier] + list(data.phenotypes),
        name='Color by Phenotype',   
    )
    phenotype_selector = alt.param(
        name='phenotype_selector', value=(starting_phenotype if starting_phenotype is not None else identifier), bind=dropdown)

    phenotype_value_selector = alt.selection_point(
        fields=['grouping_color'], bind='legend'
    )

    chart_enr_base = alt.Chart(df).mark_bar().encode(
        x=alt.X("log_enrichment:Q", bin=True),
        y=alt.Y('count()', title='# selections'),
    ).properties(height=200, width=200)

    chart_enr_bg = chart_enr_base.encode(
        color=alt.value('#ddd')
    )  # .add_params(brush)

    # highlight the selected bars
    chart_enr_hl = chart_enr_base.encode(
        color=alt.value(hash_to_color(feature))
    ).transform_filter(brush)

    chart_enr_rug = (alt.Chart(df)
                        .mark_point(
        shape='stroke', angle=90, y=0,
    )
        .encode(
            x=alt.X("log_enrichment:Q"),
            color=alt.condition(
                f"phenotype_selector != '{identifier}'",
                alt.Color('grouping_color:N',
                            title='Phenotype', scale=_phenotype_scale),
                alt.value(hash_to_color(feature)),
                # 'grouping_color:N'
            ),
    )
    )

    chart_enr = chart_enr_bg + chart_enr_hl + chart_enr_rug

    chart_enr = chart_enr.transform_calculate(
        # 'datum[phenotype_selector]'
        grouping_color='datum[phenotype_selector]'
    )

    # abundance traces
    # ----------------
    df_abd = data.df_sample.query(f"{identifier} == '{feature}' & kind == '+' & io == 'i'")[
        [identifier, 'r', 'name', 'desc_short', 'abundance', 'enrichment', 'log_enrichment', 'start', 'end', 'enr_p_value'] + list(data.phenotypes)].fillna('?')

    chart_abd_trace = (
        alt.Chart(df_abd).mark_line(point=True).encode(
            x=alt.X('r:O', title='round'),
            y=alt.Y('abundance:Q'),
            color=alt.condition(
                f"phenotype_selector != '{identifier}'",
                alt.Color('grouping_color:N',
                            title='Phenotype', scale=_phenotype_scale),
                alt.value(hash_to_color(feature)),
                # 'grouping_color:N'
            ),
            detail='name',
            tooltip=['name', 'desc_short', 'grouping_color:N']
        )
        .properties(height=200, width=200)
        .transform_filter(brush & phenotype_value_selector)
    )  # .interactive()
    chart_abd = chart_abd_trace

    _opacity = alt.condition(brush & phenotype_value_selector, alt.value(1), alt.value(0.05))
    _size = alt.condition(brush, alt.value(100), alt.value(20))

    _tooltips = ['name', 'desc_short', 'grouping_color:N',
                    alt.Tooltip('start', format=".2g"),
                    alt.Tooltip('end', format=".2g"),
                    alt.Tooltip('enrichment', format=".2g"),
                    ]
    if 'enr_p_value' in df_abd.columns:
        _tooltips.append(alt.Tooltip('enr_p_value', format=".2g"))
    if 'log_enrichment' in df_abd.columns:
        _tooltips.append(alt.Tooltip('log_enrichment', format=".2g"))
    if 'enrichment_rank' in df_abd.columns:
        _tooltips.append(alt.Tooltip('enrichment_rank', format=".0f"))

    _scales = alt.selection_interval(bind='scales',
                                        on="[mousedown[!event.altKey], mouseup] > mousemove",
                                        translate="[mousedown[!event.altKey], mouseup] > mousemove"
                                        )

    # abundance-enrichment plot
    height = 200
    enrabdplot = alt.Chart(df_abd).mark_point().encode(
        y=alt.Y('end', title='abundance, final round',
                scale=alt.Scale(type='log')),
        # x='end',
        # y='log_enrichment',
        x=alt.X('enrichment', scale=alt.Scale(type='log')),
        tooltip=_tooltips,
        color=alt.condition(
            f"phenotype_selector == '{identifier}'",
            alt.value(hash_to_color(feature)),
            alt.Color('grouping_color:N', title='Phenotype',
                        scale=_phenotype_scale)
            # 'grouping_color:N'
        ),
        opacity=_opacity,
        size=_size,
    ).properties(height=height, width=height).add_params(_scales)

    abd_hist = alt.Chart(df_abd).mark_bar().encode(
        y=alt.Y("end:Q", bin=True, scale=alt.Scale(type='log')),
        x=alt.X('count()', title='# selections'),
    ).properties(height=height, width=height)

    # chart = (
    #     (chart_abd_trace | (enrabdplot.add_params(brush) & chart_enr))
    # ).add_params(phenotype_selector).transform_calculate(grouping_color='datum[phenotype_selector]')

    chart = (
        chart_abd_trace.add_params(_scales) |
        enrabdplot.add_params(brush)
    ).transform_calculate(grouping_color='datum[phenotype_selector]').add_params(phenotype_selector, phenotype_value_selector)

    return chart

def vhh_dashboard(
        ex, global_query="kind == '+' & io == 'i'", space='cdr3', 
        rounds=['R5i', 'R6i', 'R7i', 'R8i'], 
        feature=None, 
        enr_comparison = ('start', 'end'), abbreviations={},
        show_histograms=True,
        show_traces=True,
        show_table=True,
        ):

    # setup
    full_data = FeatureData.from_experiment(ex, global_query=global_query, space='cdr3', enr_comparison=enr_comparison, abbreviations=abbreviations)

    # pre-build cache of histograms; for each sample, hist of feature enrichments
    sample_enr_dist_buffers = full_data.df_enr.groupby(
        'name')['log_enrichment'].apply(plot_dist_buffer)
    
    def make_vhh_dashboard(feature):
        data = full_data.subset(feature)

        if len(data.df_enr) == 0:
            return pn.pane.Alert((
                f"No selections found with finite enrichments for {space.upper()} `{feature}`"
            ))

        def feature_selection_table(feature, data, show_traces=True, show_histograms=True, rounds=None):
            from matplotlib.figure import Figure
            import matplotlib.pyplot as plt
            from .utils import stars

            data = data.subset(feature)
            df = data.df_enr

            if rounds is None:
                rounds = sorted(list(df.columns[df.columns.str.startswith('R')]))

            # _abundances_fig_ax = _abundances_fig.subplots(1,1)
            # _abundances_fig.add_subplot(1,1,1)

            table = df[['name','enrichment', 'log_enrichment', 'start', 'end', 'enr_p_value']]
            
            columns = ['desc_short', 'log_enrichment', 'enr_p_value', 'stars',
                       'sample_enrichments', 'start', 'abundance', 'end']

            if show_traces:
                _abundances_fig = Figure(figsize=(2, 0.25))
                abundances = df[rounds].apply(
                                sparkline, plot=plot_abundance_sparkline, fig=_abundances_fig, marker='.', axis=1).rename('abundance')
                plt.close(_abundances_fig)

                table = table.join(abundances)
            else:
                columns.remove('abundance')

            if show_histograms:
                sample_enr_hists = df.apply(mark_enr_hist, axis=1,
                                            buffers=sample_enr_dist_buffers,
                                            color=hash_to_color(feature)).rename('sample_enrichments')

                table = table.join(sample_enr_hists)
            else:
                columns.remove('sample_enrichments')
                
            table = table.set_index('name').join(data.df_selection['desc_short'])
            table['stars'] = table['enr_p_value'].apply(stars)

            # abundance_vmin = round(
            #     np.log10(df[['start', 'end']].min(axis=1).min(axis=0)), 1)
            # abundance_vmax = round(
            #     np.log10(df[['start', 'end']].max(axis=1).max(axis=0)), 1)

            widget = pn.widgets.Tabulator(
                table[columns],
                titles={
                    'desc_short': 'description',
                    'enr_p_value': 'p(enrichment|start)',
                    'stars': '*',
                    'sample_enrichments': 'sample enrichments',
                    # 'R5i': f"start<br/>{mini_histogram(np.log10(table['R5i']), vmin=abundance_vmin, vmax=abundance_vmax)}",
                    'log_enrichment': f"log10(enrichment)<br/>{mini_histogram(table['log_enrichment'])}",
                    # 'R8i':   f"end<br/>{mini_histogram(np.log10(table['R8i']), vmin=abundance_vmin, vmax=abundance_vmax)}"
                },
                formatters={'abundance': {'type': 'html'},
                            'sample_enrichments': {'type': 'html'}},
                sorters=[{'field': 'log_enrichment', 'dir': 'desc'}],
                disabled=True,
                configuration={
                    'columnDefaults': {
                        'titleFormatter': 'html',
                    },
                }
            )
            widget.style.applymap(rug_bar, subset='log_enrichment', vmin=table['log_enrichment'].min(
            ), vmax=table['log_enrichment'].max())

            # widget.style.applymap(lambda x, **kwargs: _bar(np.log10(x), **kwargs), subset='R5i', vmin=np.log10(table['R5i'].min()), vmax=np.log10(table['R5i'].max()))
            # widget.style.applymap(lambda x, **kwargs: _bar(np.log10(x), **kwargs), subset='R8i', vmin=np.log10(table['R8i'].min()), vmax=np.log10(table['R8i'].max()))
            # widget.style.applymap(_bar, subset='enrichment',     vmin=table['enrichment'].min(),     vmax=table['enrichment'].max())
            return widget


        def make_sequences_card(feature):
            table = ex.viz.summarize_similar_features([feature], space=space)
            return pn.Card(
                table, 
                title=table.caption,
                sizing_mode='stretch_width'
            )

        def make_title(feature):
            from .utils import tag
            return tag(feature, space=space)

        overview = pn.pane.Vega(feature_selection_overview(
            feature, data), debounce=1000, show_actions=True)
        
        nonlocal show_table
        if show_table:
            table = feature_selection_table(
                feature, data,
                show_histograms=show_histograms,
                show_traces=show_traces)

        def filter_table(df, selection):
            if selection is None or len(selection) == 0:
                return df
            query = ' & '.join(
                f'{crange[0]:.3f} <= `{col}` <= {crange[1]:.3f}'
                for col, crange in selection.items()
            )
            return df.query(query)

        # table.add_filter(
        #     pn.bind(filter_table, selection=overview.selection.param.brush))

        items = [
            make_title(feature),
            # make_sequences_card(feature),
            overview, 
        ]

        if show_table:
            items.append(pn.panel(table, loading_indicator=True))

        return pn.Column(*items, 
            sizing_mode='stretch_width')

    features = list(ex.fts[space].var_names.values)
    if feature is None:
        feature = features[0]
    autocomplete = pn.widgets.AutocompleteInput(
        name=f'{space}', options=features, value=feature,
        placeholder=f'Select {space.upper()}')

    return pn.Column(
        autocomplete,
        # pn.bind(summarize_clonotype, autocomplete),
        # pn.bind(make_overview, autocomplete),
        # pn.bind(make_table, autocomplete),
        pn.panel(pn.bind(make_vhh_dashboard, autocomplete), loading_indicator=True),
        sizing_mode='stretch_width'
    ).servable()
