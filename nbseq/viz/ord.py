from matplotlib.colors import LinearSegmentedColormap
import numpy as np

def ordination_plotly(
    ordn, obs=None, var=None, biplot=False,
    df=None, 
    dims=[0,1,2], color_discrete_sequence=None, 
    **kwargs):
    """diplay ordination using plotly

    Parameters
    ----------
    ordn : skbio.stats.ordination.OrdinationResults
        ordination results
    obs : pd.DataFrame, optional
        metadata about samples (observations/rows), by default None
    df : pd.DataFrame, optional
        if given, short-circuits join to `obs`, by default None
    dims : list, optional
        which axes of ordn to plot, by default [0,1,2]
    color_discrete_sequence : _type_, optional
        color palette, by default px.colors.qualitative.Alphabet


    Examples
    --------

        ft_cdr3_norm = nbseq.norm.scran(ft_path='results/tables/cdr3/feature_table.biom', conda='R-DESeq2', verbose=True)
        ft_cdr3_norm = nbseq.ft.join_metadata(ft_cdr3_norm, ex.fts.cdr3.obs, axis='obs')
        ft_cdr3_norm_r5 = nbseq.ft.query(ft_cdr3_norm, '(expt == "027i.lib") | (r == 5 & kind == "+" & io == "i")', axis='sample')
        ftf = nbseq.ft.query(ft_cdr3_norm_r5, "~(ID in ['08-H12'])", axis='obs')
        ord_skl, ord_skbio = nbseq.ordination.ordinate(ftf, method='TSVD', **{ 'n_components':100 })
        ordination_plotly(ord_skbio, ftf.obs, color='category', hover_name='description', height=600)

    """
    import plotly.express as px
    import plotly.graph_objects as go

    if df is None:
        if obs is not None:
            df_samples = ordn.samples.join(obs)
        else:
            df_samples = ordn.samples
    else:
        df_samples = df_samples
    
    df_samples = df_samples.reset_index()

    
    pe = ordn.proportion_explained[dims]
    aspect = np.sqrt(pe)

    if color_discrete_sequence is None:
        color_discrete_sequence = px.colors.qualitative.Alphabet

    fig = px.scatter_3d(df_samples, x=dims[0],y=dims[1],z=dims[2],
                  color_discrete_sequence=color_discrete_sequence,
                        labels={str(a): f'PC{a+1} [{p:0.1%}]' for a,p in zip(dims,pe)}, **kwargs)

    if biplot:
        if var is None:
            df_features = ordn.features
        else:
            df_features = ordn.features.join(var)
        # https://plotly.com/python-api-reference/generated/plotly.express.line_3d.html?highlight=color_discrete_sequence
        fig2 = px.line_3d(df_features, x=dims[0], y=dims[1], z=dims[2], 
            color='CDR3ID')
        fig.add_trace(fig2.data[0])


        # https://plotly.com/python-api-reference/generated/plotly.graph_objects.Scatter3d.html#plotly.graph_objects.Scatter3d
        # https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html#plotly.graph_objects.Figure.add_scatter3d

    fig.update_layout(
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(zip(['x','y','z'], aspect / max(aspect))),
            xaxis=dict(backgroundcolor='#cccccc',showticklabels=False),
            yaxis=dict(backgroundcolor='#bbbbbb',showticklabels=False),
            zaxis=dict(backgroundcolor='#dddddd',showticklabels=False)
        )
    )
#     fig.show()
    
    return fig


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    import matplotlib.pyplot as plt

    base = plt.colormaps[base_cmap]
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return LinearSegmentedColormap.from_list(cmap_name, color_list, N)


def cmap_to_dict(color_order, cmap):
    import matplotlib.pyplot as plt
    from collections import defaultdict

    if isinstance(cmap, dict):
        if None in cmap:
            _default = cmap[None]
            cmap = defaultdict(lambda: _default) | cmap
        color_dict = {level: cmap[level]
                      for i, level in enumerate(color_order)}
    elif isinstance(cmap, list):
        color_dict = {level: cmap[i] for i, level in enumerate(color_order)}
    else:
        if cmap is None:
            color_dict = {level: f'C{i}' for i,
                          level in enumerate(color_order)}
        else:
            if isinstance(cmap, str):
                cmap = plt.colormaps[cmap]
            color_dict = {level: cmap(i)
                          for i, level in enumerate(color_order)}

    return color_dict


def ordination_mpl(
        ordn,
        obs=None, var=None, biplot=False,
        df=None,
        color=None,
        color_order=None,
        color_kws=None,
        cmap='tab10',
        dims=[0, 1, 2],
        log=False,
        camera=None,
        xlim=None,
        ylim=None,
        zlim=None,
        show_ticks=False,
        biplot_color_order=None,
        biplot_cmap=None,
        zoom=1,
        fig_kw={},
        **kwargs):
    """diplay ordination using plotly

    Parameters
    ----------
    ordn : skbio.stats.ordination.OrdinationResults
        ordination results
    obs : pd.DataFrame, optional
        metadata about samples (observations/rows), by default None
    df : pd.DataFrame, optional
        if given, short-circuits join to `obs`, by default None
    dims : list, optional
        which axes of ordn to plot, by default [0,1,2]
    color_discrete_sequence : _type_, optional
        color palette, by default px.colors.qualitative.Alphabet


    Examples
    --------

        ft_cdr3_norm = nbseq.norm.scran(ft_path='results/tables/cdr3/feature_table.biom', conda='R-DESeq2', verbose=True)
        ft_cdr3_norm = nbseq.ft.join_metadata(ft_cdr3_norm, ex.fts.cdr3.obs, axis='obs')
        ft_cdr3_norm_r5 = nbseq.ft.query(ft_cdr3_norm, '(expt == "027i.lib") | (r == 5 & kind == "+" & io == "i")', axis='sample')
        ftf = nbseq.ft.query(ft_cdr3_norm_r5, "~(ID in ['08-H12'])", axis='obs')
        ord_skl, ord_skbio = nbseq.ordination.ordinate(ftf, method='TSVD', **{ 'n_components':100 })
        ordination_plotly(ord_skbio, ftf.obs, color='category', hover_name='description', height=600)

    """
    import matplotlib.pyplot as plt

    if df is None:
        if obs is not None:
            df_samples = ordn.samples.join(obs)
        else:
            df_samples = ordn.samples
    else:
        df_samples = df_samples

    # if biplot:
    #     if var is not None:
    #         df_features = ordn.biplot_scores.join(var)
    # else:
    #     df_features = ordn.biplot_scores

    df_samples = df_samples.reset_index()

    pe = list(ordn.proportion_explained[dims])
    aspect = np.sqrt(pe)

    fig = plt.figure(**fig_kw)
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect(aspect, zoom=zoom)

    if camera is not None:
        ax.view_init(**camera)

    if color_order is None:
        # sorted(df_samples[color].unique())
        color_order = (df_samples[color].unique())

    color_dict = cmap_to_dict(color_order, cmap)
    if color_kws is None:
        color_kws = dict()

    groups = df_samples.groupby(color)[dims]

    # for name, group in groups:
    for name in color_order:
        if name not in groups.groups:
            continue
        group = groups.get_group(name)

        xs = group[dims[0]]
        ys = group[dims[1]]
        zs = group[dims[2]]

        if log:
            xs = np.log10(xs)
            ys = np.log10(ys)
            zs = np.log10(zs)

        if name in color_kws:
            kws = color_kws[name]
        else:
            kws = dict()
        if name in color_dict:
            kws['color'] = color_dict[name]
        elif None in color_dict:
            kws['color'] = color_dict[None]

        ax.scatter(xs, ys, zs, label=name, **{**kwargs, **kws})

    if biplot:
        if biplot_color_order is None:
            biplot_color_order = ordn.biplot_scores.index

        biplot_color_dict = cmap_to_dict(biplot_color_order, biplot_cmap)
        for feature in biplot_color_order:
            if feature not in ordn.biplot_scores.index:
                continue
            xs = [0, ordn.biplot_scores.loc[feature, dims[0]]]
            ys = [0, ordn.biplot_scores.loc[feature, dims[1]]]
            zs = [0, ordn.biplot_scores.loc[feature, dims[2]]]

            if log:
                xs = [0, np.log10(xs[1])]
                ys = [0, np.log10(ys[1])]
                zs = [0, np.log10(zs[1])]

            ax.plot(
                xs, ys, zs,
                label=feature,
                color=biplot_color_dict[feature]
            )

    ax.set_xlabel('{a} [{p:0.1%}]'.format(
        a=(f"PC{dims[0]+1}" if isinstance(dims[0], int) else dims[0]), p=pe[0]))
    ax.set_ylabel('{a} [{p:0.1%}]'.format(
        a=(f"PC{dims[1]+1}" if isinstance(dims[1], int) else dims[1]), p=pe[1]))
    ax.set_zlabel('{a} [{p:0.1%}]'.format(
        a=(f"PC{dims[2]+1}" if isinstance(dims[2], int) else dims[2]), p=pe[2]))

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if zlim is not None:
        ax.set_zlim(zlim)

    if not show_ticks:
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])

    return fig, ax
