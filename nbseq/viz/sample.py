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

def collapse_top_asvs(ft_ad, samples, n=20, top_from_samples=None, relative=False):

    if top_from_samples is None:
        top_from_samples = samples

    ft2 = ft_ad[top_from_samples,:]

    asv_counts = ft_sum(ft2, axis='var').sort_values(ascending=False)
    top_asvs = asv_counts.iloc[0:n]

    if relative:
        ft_ad = to_relative(ft_ad)

    ft_samples = ft_ad[samples,:]
    ft_top = ft_samples[:,top_asvs.index.values]

    # sum top ASVs for each sample
    ft_top_sum = ft_top.X.sum(axis=1)

    # sum all ASVs for each sample
    ft_samples_sum = ft_samples.X.sum(axis=1)

    # calculate sum of all non-top ASVs for each sample
    # interestingly, much faster to get sum of non-top-ASVs by summing everything,
    # then summing top ASVs and subtracting, as opposed to summing non-top-ASVs...
    # the below options are ~10x slower
    #     other_sum = fts[:,other_asvs.index.values].X.sum(axis=1)
    #     other_sum = ft_sum(fts[:,other_asvs.index.values])
    other_asv_sum = ft_samples_sum - ft_top_sum

    ft_other = AnnData(
       other_asv_sum,
       obs=pd.DataFrame(index=ft_samples.obs_names),
       var=pd.DataFrame(index=['others'])
    )
    return anndata.concat([ft_top, ft_other], axis=1, merge="first")

def make_md5_colors(hashes, special_values={'others':"#aaaaaa", 'unshared':"#aaaaaa", 'other shared': "#999999"}):
    return [hash_to_color(h, **special_values) for h in hashes]

def _scale_manual_md5(hashes, sort=True, special_values={'others':"#aaaaaa", 'unshared':"#aaaaaa", 'other shared': "#999999"}):
    if sort:
        hashes = sorted(list(hashes[~np.isin(hashes, list(special_values))]) + list(special_values.keys()))
    labels = [h if h in special_values else pretty_hex(h) for h in hashes]
    values = make_md5_colors(hashes, special_values=special_values)
    return dict(values=values, limits=hashes, breaks=hashes, labels=labels)

def scale_fill_md5(hashes, **kwargs):
    return p9.scale_fill_manual(**_scale_manual_md5(hashes, **kwargs))

def scale_color_md5(hashes, **kwargs):
    return p9.scale_color_manual(**_scale_manual_md5(hashes, **kwargs))

def top_asv_barplot(df, special_values={'others':"#aaaaaa"}, limits=None, x='r', facet=None, title=None, **kwargs):

    # not sure why this is necessary, but get this error if omitted:
    # TypeError: object of type 'NoneType' has no len()
    from pandas.api.types import is_numeric_dtype, is_string_dtype
    if is_numeric_dtype(df[x]):
        x_scale = p9.scale_x_continuous(name='Round')
    else:
        x_scale = p9.scale_x_discrete(name='Round')

    gg = (p9.ggplot(df, p9.aes(x=x,y='abundance',fill='feature'))
            + p9.geom_col()
            + scale_fill_md5(df['feature'].unique(), special_values=special_values)
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

def top_asv_lineplot(
    df, special_values={'others':"#aaaaaa"}, limits=None, x='r', facet=None,
    figsize=(12,8), **kwargs):

    gg = (p9.ggplot(df, p9.aes(x=x,y='abundance',fill='feature'))
            + p9.geom_line(aes(group='feature'))
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




def summarize_features(df, features=None, seq_col='CDR3', id_col='CDR3ID', max_features=100, title=None):
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
