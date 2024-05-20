import plotnine as p9
from plotnine import *

from ..ft import *
from .asv import hash_to_color, pretty_hex, pack_hex
from .utils import contrasting_color


import anndata
from anndata import AnnData
import pandas as pd

def _fortify_samples(ft_ad):
    df = fortify(ft_ad, obs=True, var=True)
    # return pd.merge(
    #     pd.merge(df,
    #              df.groupby('sample')['abundance'].rank(method='dense',ascending=False).rename('rank'),
    #              left_index=True, right_index=True)
    #         .sort_values('abundance',ascending=False),
    #     fda_cdr3,
    #     left_on='feature',right_on='CDR3ID')
    return (df.join(
                 df.groupby('sample')['abundance'].rank(method='dense',ascending=False).rename('rank')
            ).sort_values('abundance',ascending=False))


def collapse_top_shared_asvs(ft_ad, ids_1, ids_2, n=20, top_from_samples=None, relative=False):
    from ..ft import to_relative as ft_to_relative_abundance
    ft_samples = set(ft_ad.obs_names)
    ids_1 = list(ft_samples.intersection(ids_1))
    ids_2 = list(ft_samples.intersection(ids_2))
    samples = ids_1 + ids_2

    s1 = ft_ad[ids_1,:]
    s2 = ft_ad[ids_2,:]

    shared_vars = (s1.X.sum(axis=0) > 0).A1 & (s2.X.sum(axis=0) > 0).A1

    if relative:
        ft_ad = ft_to_relative_abundance(ft_ad)

    ft2 = ft_ad[ids_1 + ids_2, shared_vars]

    asv_counts = ft_sum(ft2, axis='var').sort_values(ascending=False)
    top_asvs = asv_counts.iloc[0:n]
    other_asvs = asv_counts.iloc[n:]

    fts = ft_ad[samples,:]
    ft_top = fts[:,top_asvs.index.values]
    ft_top_sum = ft_top.X.sum(axis=1)

    other_shared_sum = ft2.X.sum(axis=1) - ft_top_sum
    fts_sum = fts.X.sum(axis=1)

    # interestingly, much faster to get sum of non-top-ASVs by summing everything,
    # then summing top ASVs and subtracting, as opposed to summing non-top-ASVs...
    # the below options are ~10x slower
    #     other_sum = fts[:,other_asvs.index.values].X.sum(axis=1)
    #     other_sum = ft_sum(fts[:,other_asvs.index.values])
    unshared_sum = fts_sum - other_shared_sum - ft_top_sum

    ft_other = AnnData(
        np.hstack([
            other_shared_sum,
            unshared_sum]),
       obs=pd.DataFrame(index=fts.obs_names),
       var=pd.DataFrame(index=['other shared', 'unshared']))

    return anndata.concat([ft_top, ft_other], axis=1)

from ..ft import collapse_top_asvs, fortify_top_asvs

def make_md5_colors(hashes, special_values={'others':"#aaaaaa", 'unshared':"#aaaaaa", 'other shared': "#999999"}):
    return [hash_to_color(h, **special_values) for h in hashes]

def _scale_manual_md5(hashes, sort=True, special_values={'others':"#aaaaaa", 'unshared':"#aaaaaa", 'other shared': "#999999"}, **kwargs):
    if sort:
        hashes = sorted(list(hashes[~np.isin(hashes, list(special_values))]) + list(special_values.keys()))
    labels = [h if h in special_values else pretty_hex(h) for h in hashes]
    values = make_md5_colors(hashes, special_values=special_values)
    return {**dict(values=values, limits=hashes, breaks=hashes, labels=labels), **kwargs}

def scale_fill_md5(hashes, **kwargs):
    return p9.scale_fill_manual(**_scale_manual_md5(hashes, **kwargs))

def scale_color_md5(hashes, **kwargs):
    return p9.scale_color_manual(**_scale_manual_md5(hashes, **kwargs))


def alt_scale_features(features, sort=False, **kwargs):
    import altair as alt
    import pandas as pd

    if isinstance(features, pd.Series):
        features = features.values
    if sort:
        features = sorted(features)

    return alt.Scale(domain=features, range=[hash_to_color(f) for f in features])


def top_asv_barplot(df, special_values={'others':"#aaaaaa"}, limits=None, x='r', feature_name="feature", facet=None, title=None, **kwargs):

    # not sure why this is necessary, but get this error if omitted:
    # TypeError: object of type 'NoneType' has no len()
    from pandas.api.types import is_numeric_dtype, is_string_dtype
    if is_numeric_dtype(df[x]):
        x_scale = p9.scale_x_continuous(name='Round')
    else:
        x_scale = p9.scale_x_discrete(name='Round')

    gg = (p9.ggplot(df, p9.aes(x=x,y='abundance',fill='feature'))
            + p9.geom_col()
            + scale_fill_md5(df['feature'].unique(), special_values=special_values, name=feature_name)
            + x_scale
            + p9.scale_y_continuous(name='Abundance')
         )

    if facet is not None:
        if isinstance(facet, dict):
            gg = gg + p9.facet_wrap(**facet)
        else:
            gg = gg + p9.facet_wrap(facet)
    if limits is not None:
        gg = gg + p9.coord_cartesian(ylim=limits)
    if title is not None:
        gg = gg + p9.ggtitle(title)
    return gg

def top_asv_barplot_alt(ft, query, n=30, select_from_round=8, x='r:O', fancy_sorting=False, phylo=False, **kwargs):
    import altair as alt

    df = fortify_top_asvs(ft, query, n=n, select_from_round=select_from_round)
    
    # establish sorting order for features; sort by geometric mean abundance 
    # across all samples in query
    df['log_abundance'] = np.log10(df['abundance'])
    features_by_mean_abundance = (df.groupby('feature')
                                  ['log_abundance'].mean()
                                  .rename('mean_log_abundance')
                                  .sort_values(ascending=False))
    
    # attach to df so we can use this for sorting bars
    df = df.join(features_by_mean_abundance, on='feature')
    
    if phylo:
        raise NotImplementedError()
    else: 
        _features = features_by_mean_abundance.index.values
        selector = alt.selection_point(fields=['feature'], bind='legend', on='click', clear='dblclick')
        feature_scale = alt.Scale(domain=_features, 
                                  range=[hash_to_color(f) for f in _features])
        
        return (alt.Chart(df)
         .mark_bar()
             .encode(x=x, 
                     y='abundance', 
                     tooltip=[
                         'feature',
                         alt.Tooltip('abundance', format=".2g"),
                         alt.Tooltip('log_abundance', format=".2f")],
                     opacity=alt.condition(selector, alt.value(1), alt.value(0.1)),
                     # order=(alt.Order('mean_log_abundance',sort='ascending') if fancy_sorting else None),
                     color=alt.Color('feature:N',
                                     legend=alt.Legend(columns=n//20,symbolLimit=0,labelExpr="slice(datum.value,0,6)"),
                                     scale=feature_scale,
                                     # sort=(_features if fancy_sorting else None)
                                    ))
         # .facet(column='description')
        .add_params(selector)
        )

def top_asv_barplot_hist(ex, query, n=100, space='cdr3', select_from_round=8, x='r:O'):
    import altair as alt
    from altair_transform import transform_chart

    # df2 = fortify_top_asvs(ex,query='True', n=n)
    df = fortify_top_asvs(ex, query, space, n=n, select_from_round=select_from_round)
    
    # establish sorting order for features; sort by geometric mean abundance 
    # across all samples in query
    df['log_abundance'] = np.log10(df['abundance'])
    features_by_mean_abundance = (df.groupby('feature')
                                  ['log_abundance'].mean()
                                  .rename('mean_log_abundance')
                                  .sort_values(ascending=False))
    
    # attach to df so we can use this for sorting bars
    df = df.join(features_by_mean_abundance, on='feature')
    
    df2_features = list(set(ex.fts[space].var_names).intersection(features_by_mean_abundance.index.values))
    
    
    df2 = fortify(ex.fts[space][:,df2_features], feature_col='feature')
    df2['log_abundance'] = np.log10(df2['abundance'])


    # _features = pd.concat([df['feature'], df2['feature']]).unique()
    _features = list(features_by_mean_abundance.index.values)# + ['others']
    feature_scale = alt.Scale(domain=_features, range=[hash_to_color(f) for f in _features])
    
    selector = alt.selection_point(fields=['feature'], bind='legend', on='click', clear='dblclick')
    # selector = alt.selection_point(fields=['feature'])
    
    bar_chart = (alt.Chart(df)
     .mark_bar()
         .encode(x=x, 
                 y='abundance', 
                 tooltip=[
                     # 'feature',
                     # alt.Tooltip('feature', labelExpr="slice(datum.value,0,6) + ' ' + slice(datum.value,6)"),
                     'feature',
                     alt.Tooltip('abundance', format=".2g"),
                     alt.Tooltip('log_abundance', format=".2f", title='log10(abundance)')],
#                  color=alt.condition(
#                     selector,
#                     'feature:N',
#                     alt.value('lightgray'),

#                     legend=alt.Legend(columns=3,symbolLimit=0),
#                     scale=feature_scale
#                  )),
                 color=alt.Color(
                    'feature:N',
                    legend=alt.Legend(columns=n//20,symbolLimit=0, 
                                      labelExpr="slice(datum.value,0,6)"
                                     ),
                    scale=feature_scale, 
                 ),                 
                 opacity=alt.condition(selector, alt.value(1), alt.value(0.1))
                )
     .add_params(selector)
    )

    feature_hist = alt.Chart(df2).mark_bar().encode(
        alt.X('log_abundance', bin=True),
        y=alt.Y('count()', title='# samples', stack=None),
        fill=alt.Fill('feature', scale=feature_scale),
        tooltip=['feature'],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.01))
    )
    feature_hist_t = transform_chart(feature_hist).add_params(selector)
    # .transform_filter(
    #     selector
    # )

    return bar_chart | feature_hist_t    
    

def top_asv_lineplot(
    df, special_values={'others':"#aaaaaa"}, limits=None, x='r', facet=None, title=None,
    figsize=(12,8), **kwargs):

    gg = (p9.ggplot(df, p9.aes(x=x,y='abundance',fill='feature'))
            + p9.geom_line(p9.aes(group='feature'))
            + scale_color_md5(df['feature'].unique(), special_values=special_values)
            + p9.scale_x_discrete(name='Round')
            + p9.scale_y_continuous(name='Abundance')
         )

    if facet is not None:
        gg = gg + p9.facet_wrap(facet)
    if limits is not None:
        gg = gg + p9.coord_cartesian(ylim=limits)
    if title is not None:
        gg = gg + p9.ggtitle(title)
    return gg

def feature_traceplot(
        ft, query, detail='feature', tooltip=['feature', 'name', 'description'],
        features=None,
        n=50, select_from_round=None,
        selector=None, feature_scale=None):
    import altair as alt
    # feature_scale = alt.Scale(domain=_features, range=[hash_to_color(f) for f in _features])

    df = fortify_top_asvs(ft, query, n=n, select_from_round=select_from_round)

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

def top_asv_plot_phylo(
    df, tree, special_values={'others':"#aaaaaa"}, limits=None, x='r', facet=None, title=None,
    figsize=(12,8), tree_width=0.5):

    import patchworklib as pw
    from .tree import draw as draw_tree

    features = frozenset(df['feature'].unique())
    tip_names = [t.name for t in tree.tips()]

    tree = tree.shear(features.intersection(tip_names))

    features_in_tree = list(reversed([f.name for f in tree.tips() if f.name in features]))
    features_not_in_tree = sorted(features.difference(features_in_tree))

    features_sorted = features_in_tree + features_not_in_tree

    df['feature'] = pd.Categorical(df['feature'], categories=features_sorted)
    gg = (p9.ggplot(df, p9.aes(x=x,y='abundance',fill='feature'))
            + p9.geom_col()
            + scale_fill_md5(features_sorted, sort=False, special_values=special_values)
            + p9.scale_x_discrete(name='Round')
            + p9.scale_y_continuous(name='Abundance')
            + p9.theme(legend_position='none')
         )

    if facet is not None:
        gg = gg + p9.facet_wrap(facet)
    if limits is not None:
        gg = gg + p9.coord_cartesian(ylim=limits)
    if title is not None:
        gg = gg + p9.ggtitle(title)

    tree_figsize = (figsize[0] * tree_width, figsize[1])

    ax_tree = pw.Brick(figsize=tree_figsize)
    # fig_tree, ax_tree = plt.subplots(figsize=(10,10));
    # ax_tree = subfigs[1].subplots(1,1)
    draw_tree(tree,
         label_func=lambda c: pretty_hex(c.name) if (c.is_tip() and (c.name is not None)) else None,
         label_colors=lambda c: pack_hex(c) if c is not None else None,
         show_confidence=False,
         do_show=False,
         axes=ax_tree)


    bar_figsize = (figsize[0] * (1-tree_width), figsize[1])
    return (pw.load_ggplot(gg, figsize=bar_figsize) | ax_tree).savefig()




def rank_abundance_plot(ft, log=True, pct=False, obs=True, transform=None,
                        point=False, line=True,
                        n_head=10,n_sample=None,n_tail=10,
                         **kwargs):
    """make plot of VHH rank vs. [log-]abundance (Whittaker plots)

    Parameters
    ----------
    ft : anndata.AnnData
        feature table
    log : bool, optional
        plot log(abundance) or linear abundance, by default True
    pct : bool, optional
        plot percentile on X axis? if false, plot rank on x axis, by default False
    obs : bool, optional
        if True, when feature table is fortified to a DataFrame, join the `obs` metadata to the output so that its columns can be referenced, by default True

    Returns
    -------
    ggplot
    """
    from ..ft import fortify
    df = fortify(ft, relative=True, obs=obs)
    
    if log:
        df['abundance'] = np.log10(df['abundance'])
    df['rank'] = df.groupby('ID')['abundance'].rank(ascending=False, pct=pct)

    if n_sample is not None:
        if not isinstance(n_sample, dict):
            n_sample = {'n': n_sample}

        df = pd.concat([
            df.groupby(['ID', 'round']).head(n_head),
            df.groupby(['ID', 'round']).sample(**n_sample),
            df.groupby(['ID', 'round']).tail(n_tail)
        ])

    if transform is not None:
        df = transform(df)

    gg = (p9.ggplot(df, p9.aes(x='rank',y='abundance',**kwargs)) 
            + p9.xlab('feature percentile' if pct else 'feature rank') 
            + p9.ylab('log10(relative abundance)' if log else 'relative abundance'))
    if point:
        if not isinstance(point, dict):
            point = {}
        gg += p9.geom_point(**point)
    if line:
        if not isinstance(line, dict):
            line = {}
        gg += p9.geom_line(**line)
    return gg


def summarize_features(df, features=None, seq_col='CDR3', id_col='CDR3ID', max_features=100, title=None):
    from .syntax import aa_highlighter
    caption = []
    if title is not None:
        caption = [f"<strong>{title}</strong>"]
    if features is not None:
        dfs = df.loc[df[id_col].isin(features),:]
        selected_reads = dfs.abundance.sum()
        total_reads = df.abundance.sum()
        caption.append(f"Selected {len(features)} / {len(df)} features ({selected_reads:.0f} / {total_reads:.0f} reads, {selected_reads/total_reads:.1%})")
    else:
        dfs = df

    if len(df) > max_features:
        dff = dfs.iloc[0:max_features,:].copy()
        dfm = dfs.iloc[max_features:,:]
        masked_reads = dfm.abundance.sum()
        total_reads = dfs.abundance.sum()
        caption.append(f"Showing top {max_features}; hiding {len(dfm)} more features, {masked_reads:.0f} / {total_reads:.0f} ({masked_reads/total_reads:.1%}) more reads")
    else:
        dff = df
        caption = ''

    dff[seq_col] = dff[seq_col].apply(lambda x: aa_highlighter.highlight(x))
    dff['color'] = dff[id_col].str.slice(0,6)
    table = (
        dff[['color',seq_col,'abundance','rank']].style
            .applymap(lambda c: f"background-color: #{c}; color: {contrasting_color('#'+c)}", subset=['color'])
            .applymap(lambda c: f"font-family: monospace;", subset=[seq_col])
            .format({'rank':'{:.0f}', 'abundance':'{:.0f}'})
            .set_caption("<br />".join(caption))
    )
    return table


def compare_samples(ids_1, ids_2, ft, n=20, names=['1','2'], show_tables=True, **kwargs):
    from ..distance import shared_reads, shared_reads_frac, shared_features
    from .utils import display_accordion, display_table
    from IPython.display import display

    print(f"Comparing samples {ids_1} vs. {ids_2}...")


    if len(ids_1) != len(ids_2):
        raise Exception("ids_1 and ids_2 must have same length")

    n1 = names[0]
    n2 = names[1]

    res = pd.DataFrame(index = pd.MultiIndex.from_tuples(zip(ids_1, ids_2), names=['id1','id2']),
                       columns=[
                           'features',
                           f'features_{n1}',
                           f'features_{n2}',
                           'features_shared',
                           'features_shared_frac',
                           f'features_only_{n1}',
                           f'features_only_{n1}_frac',
                           f'features_only_{n2}',
                           f'features_only_{n2}_frac',
                           'reads',
                           f'reads_{n1}',
                           f'reads_{n2}',
                           'reads_shared',
                           'reads_shared_frac',
                           f'reads_only_{n1}',
                           f'reads_only_{n1}_frac',
                           f'reads_only_{n2}',
                           f'reads_only_{n2}_frac',
                       ])
    for id1, id2 in zip(ids_1, ids_2):
        ft1 = ft[id1,:]
        ft2 = ft[id2,:]

        v1 = ft1.X
        v2 = ft2.X

        reads_v1 = v1.sum()
        reads_v2 = v2.sum()

        intersect, shared_v1, shared_v2 = np.intersect1d(v1.indices, v2.indices,
            assume_unique=True,
            return_indices=True)

        res.loc[(id1, id2), 'reads'] = reads_v1 + reads_v2
        res.loc[(id1, id2), f'reads_{n1}'] = reads_v1
        res.loc[(id1, id2), f'reads_{n2}'] = reads_v2

        res.loc[(id1, id2), 'reads_shared'] = shared_reads(v1, v2)
        res.loc[(id1, id2), f'reads_shared_{n1}'] = v1.data[shared_v1].sum()
        res.loc[(id1, id2), f'reads_shared_{n2}'] = v2.data[shared_v2].sum()
        res.loc[(id1, id2), 'reads_shared_frac'] = shared_reads_frac(v1, v2)

        res.loc[(id1, id2), f'reads_only_{n1}'] = reads_v1 - res.loc[(id1, id2), f'reads_shared_{n1}']
        res.loc[(id1, id2), f'reads_only_{n1}_frac'] = res.loc[(id1, id2), f'reads_only_{n1}'] / reads_v1
        res.loc[(id1, id2), f'reads_only_{n2}'] = reads_v2 - res.loc[(id1, id2), f'reads_shared_{n2}']
        res.loc[(id1, id2), f'reads_only_{n2}_frac'] = res.loc[(id1, id2), f'reads_only_{n2}'] / reads_v2

        res.loc[(id1, id2), 'features'] = (v1 + v2).nnz
        res.loc[(id1, id2), f'features_{n1}'] = v1.nnz
        res.loc[(id1, id2), f'features_{n2}'] = v2.nnz

        res.loc[(id1, id2), 'features_shared'] = shared_features(v1, v2)
        res.loc[(id1, id2), 'features_shared_frac'] = res.loc[(id1, id2), 'features_shared'] / res.loc[(id1, id2), 'features']

        res.loc[(id1, id2), f'features_only_{n1}'] = v1.nnz - shared_v1.size
        res.loc[(id1, id2), f'features_only_{n1}_frac'] = res.loc[(id1, id2), f'features_only_{n1}'] / v1.nnz
        res.loc[(id1, id2), f'features_only_{n2}'] = v2.nnz - shared_v2.size
        res.loc[(id1, id2), f'features_only_{n2}_frac'] = res.loc[(id1, id2), f'features_only_{n2}'] / v2.nnz

        if show_tables:
            display(
                display_accordion(display_table([[
                    summarize_features(_fortify_samples(ft1), features=ft1.var_names.values[np.delete(v1.indices, shared_v1)], max_features=20, title=f"{n1}-only"),
                    summarize_features(_fortify_samples(ft2), features=ft2.var_names.values[np.delete(v2.indices, shared_v2)], max_features=20, title=f"{n2}-only")
    #                 display_accordion(summarize_features(_fortify_samples(ft2[:,np.delete(v2.indices, shared_v2)]), max_features=20), n2)
                ]]), title=f"{id1} vs. {id2}")
            )


    ft_top_shared = collapse_top_shared_asvs(ft, ids_1, ids_2, relative=True, n=n)
    df = fortify(ft_top_shared, obs=metadata, relative=False)
    display(top_asv_plot(df, special_values={'unshared':"#aaaaaa", 'other shared': "#888888"}, facet='expt', **kwargs))
#     return res

    res = res.infer_objects()
    df = res.reset_index()
    df['pair'] = df['id1'] + '/' + df['id2']

    df_features = (df
                .loc[:,['pair','features_shared',f'features_only_{n1}',f'features_only_{n2}']]
                .melt(id_vars=['pair'],var_name='partition',value_name='features'))
    df_features['partition'] = df_features['partition'].astype('category').cat.reorder_categories([f'features_only_{n1}','features_shared',f'features_only_{n2}'])

    df_reads = (df
                .loc[:,['pair',f'reads_shared_{n1}',f'reads_shared_{n2}',f'reads_only_{n1}',f'reads_only_{n2}']]
                .melt(id_vars=['pair'],var_name='partition',value_name='reads'))
    df_reads['partition'] = df_reads['partition'].astype('category').cat.reorder_categories([f'reads_only_{n1}',f'reads_shared_{n1}',f'reads_shared_{n2}',f'reads_only_{n2}'])


    from plotnine import ggplot, aes, geom_col, theme, element_text
#     display_table([[
    display(ggplot(df_features,
               aes(x='pair',y='features',fill='partition')) +
            geom_col(position='stack') +
            theme(figure_size=(4,4), axis_text_x=element_text(rotation=-45))),

    display(ggplot(df_reads,
               aes(x='pair',y='reads',fill='partition')) +
            geom_col(position='stack') +
            theme(figure_size=(4,4), axis_text_x=element_text(rotation=-45)))
#     ]])
    display(res.style.format(
        dict(zip(res.columns,
                 ['{:,.1%}' if '_frac' in c else '{:,.0f}' for c in res.columns]))
    ))

    return ft_top_shared


def design_heatmap(des, log=False, row_cluster=False, col_cluster=False, figsize=(30,10), **kwargs):
    from matplotlib.colors import LogNorm, SymLogNorm, Normalize
    import matplotlib.pyplot as plt
    import seaborn as sns

    kwargs = { 'norm': LogNorm() if log else None, **kwargs }
    if row_cluster or col_cluster:
        sns.clustermap(des, figsize=figsize, **kwargs)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(des, ax=ax, **kwargs)


def plot_design_distributions(designs, n_features=100, n_points = 10000):
    """ Plots the marginal and per-feature distributions for a subset of features, for each of a series of design matrices
    """

    from numpy.linalg import LinAlgError
    import matplotlib.pyplot as plt
    import seaborn as sns

    palette=sns.color_palette(['tab:blue']*n_features)

    fig, axs = plt.subplots(nrows=len(designs), ncols=2, figsize=(12, 1.5*len(designs)))
    for i, (name, mat) in enumerate(designs.items()):
        if isinstance(mat, AnnData):
            X = mat.X.data
        else:
            X = mat.values

        log_scale = 10 if (X >= 0).all() else False

        if log_scale:
            X_points = np.random.choice(X[X>0], size=n_points)
        else:
            X_points = np.random.choice(X.reshape(-1), size=n_points)

        sns.histplot(X_points, log_scale=log_scale, ax=axs[i, 0])

        # Xn = X[X<0]
        # Xp = X[X>0]

        # sns.histplot(np.random.choice(X[X>0], size=n_points), log_scale=True, ax=axs[i, 1])

        # if len(Xn) > 0:
        #     sns.histplot(-np.random.choice(Xn, size=min(len(Xn),n_points)), log_scale=True, ax=axs[i, 0])
        #     axs[i, 0].invert_xaxis()
        # axs[i, 0].set_xscale('symlog')

        for attempt in range(3):
            try:
                if isinstance(mat, AnnData):
                    # randomly pick `n_features` columns from `mat`
                    idxs = np.array([1] * n_features + [1] * (mat.shape[1] - n_features))
                    np.random.shuffle(idxs)
                    data = mat[:, idxs].to_df()
                else:
                    data = mat.sample(n_features, axis=1)
                sns.kdeplot(data=data, alpha=0.1,
                            log_scale=log_scale,
                            legend=False, palette=palette, warn_singular=False, ax=axs[i, 1])
                axs[i, 1].set_xscale('symlog')
            except LinAlgError:
                pass
            break

        axs[i,0].set_title(f"{name}\n(marginal)")
        axs[i,1].set_title(f"{name}\n(features)")

    plt.tight_layout()
