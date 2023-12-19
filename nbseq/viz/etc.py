from altair import datum
from nbseq.viz.sample import alt_scale_features
from nbseq.utils import clamp_finite
def figure_to_sparkline(fig):
    import base64
    from io import BytesIO
    import base64    
    
    for ax in fig.axes:
        ax.axis('off')
    
    # squeeze axis to the edges of the figure
    fig.subplots_adjust(left=0)
    fig.subplots_adjust(right=0.99)
    fig.subplots_adjust(bottom=0.1)
    fig.subplots_adjust(top=0.9)
    
    # save the figure to html
    bio = BytesIO()
    plt.savefig(bio)
    plt.close()
    html = """<img src="data:image/png;base64,%s"/>""" % base64.b64encode(bio.getvalue()).decode('utf-8')
    return html

def sparkline(data, plot,
              figsize=(4, 0.25), **kwargs):
    """Create a single HTML image tag containing a base64 encoded
    sparkline style plot
    
    Parameters
    ----------
    data : array-like (list, 1d Numpy array, Pandas Series) sequence of
        data to plot
    figsize : tuple of float, length and height of sparkline plot.  Attribute
        of matplotlib.pyplot.plot.
    """
    from itertools import chain
    from io import BytesIO
    import base64
    
    # fig = plt.figure(figsize=figsize)  # set figure size to be small
    # ax = fig.add_subplot(111)
    fig, ax = plt.subplots(figsize=figsize)
    
    plot(data, ax=ax, fig=fig, **kwargs)
    
    return figure_to_sparkline(fig)


def plot_abundance_sparkline(
    data, ax, fig, 
    point=True, point_color='red', point_marker='.',
    point_fill='red', point_size=6, point_alpha=1.0,
    points=True, fill=True, fill_color='blue', fill_alpha=0.1,**kwargs):
    """
    point : bool, show point marker on last point on right
    point_location : not implemented, always plots rightmost
    point_color : str, matplotlib color code for point, default 'red'
    point_marker : str, matplotlib marker code for point
    point_fill : str, matplotlib marker fill color for point, default 'red'
    point_size : int, matplotlib markersize, default 6
    point_alpha : float, matplotlib alpha transparency for point
    fill : bool, show fill below line
    fill_color : str, matplotlib color code for fill
    fill_alpha : float, matplotlib alpha transparency for fill
    **kwargs : keyword arguments passed to matplotlib.pyplot.plot
    """
    data = list(data)
    
    plot_len = len(data)
    plot_min = min(data)
    point_x = plot_len - 1
    
    plt.plot(data, **kwargs)
    
    # fill between the axes
    plt.fill_between(range(plot_len), data, plot_len*[plot_min],
                     color=fill_color, alpha=fill_alpha)
    
    # plot the right-most point red, probably on makes sense in timeseries
    plt.plot(point_x, data[point_x], color=point_fill,
             marker=point_marker, markeredgecolor=point_color,
             markersize=point_size,
             alpha=point_alpha, clip_on=False)


def plot_abundance_hist_context(data, context, ax, **kwargs):
    ax.hist(np.log10(context), alpha=0.5)
    ymaxes = np.arange(len(data)) / len(data)
    for datum, ymax in zip(np.log10(data), ymaxes):
        ax.axvline(datum, ymin=0, ymax=ymax, color='red')

def plot_log_hist_context(datum, context, ax, **kwargs):
    ax.hist(np.log10(context), alpha=0.5)
    ax.axvline(np.log10(datum), color='red')

def plot_abundance_hists(data, context):
    fig = plot_abundance_hists_context(context)
    fig = mark_abundance_hists(data, fig)
    return g.fig

def plot_abundance_hists_context(context):
    g = sns.displot(data=context, col='r',y='log_abundance', element='step', aspect=0.5, height=1.5)
    for ax in g.axes.flat:
        ax.axis('off')
    g.set_titles('')
    g.fig.subplots_adjust(wspace=0, hspace=0)
    return g.fig

def mark_abundance_hists(data, hists_fig):
    for i, x in enumerate(data):
        hists_fig.axes[i].axhline(y=x, color='red')
    return hists_fig

def _pickle_figure(fig):
    import pickle
    from io import BytesIO
    buf = BytesIO()
    pickle.dump(fig, buf)
    return buf

def plot_abundance_hist_sparklines(df):
    import pickle
    
    context = pd.DataFrame({
        'feature':df['feature'], 'r':df['r'], 
        'log_abundance':np.log10(df['abundance'])}).sort_values(['feature','r'])

    _abundance_hists = plot_abundance_hists_context(context)
    # plt.close()
    with _pickle_figure(_abundance_hists) as _buf:
        group = df.sort_values(['feature','r']).groupby(['sample','feature'])

        def _mark_hist_template(data):
            _buf.seek(0)
            fig = pickle.load(_buf) 
            return figure_to_sparkline(mark_abundance_hists(np.log10(data),fig))

        hists = group['abundance'].apply(_mark_hist_template)
    return hists





df = df_enr.query(f"CDR3ID == '{feature}'").set_index('name')

abundances = df[['R5i','R6i','R7i','R8i']].apply(sparkline, plot=plot_abundance_sparkline, figsize=(2,0.25), marker='.', axis=1).rename('abundance')


df[['enrichment', 'log_enrichment']].join(abundances).style.bar(subset=['enrichment', 'log_enrichment'])












def make_vega_pane(feature, log_scale=False):
    brush = alt.selection_interval(name='brush')  # selection of type "interval"
    
    df = nbseq.ft.fortify(ft[:,feature],obs=True)
    df['log_abundance'] = np.log10(df['abundance'])
    
    if log_scale:
        x = alt.X("log_abundance:Q", bin=True)
    else:
        x = alt.X("abundance:Q", bin=alt.BinParams(anchor=1, maxbins=20))
    
    chart = alt.Chart(df).mark_bar().encode(
        x,
        y=alt.Y('count()', title='# selections'),
        fill=alt.value(nbseq.viz.asv.hash_to_color(feature))
    ).add_selection(
        brush
    ).properties(height=75).facet(row='r')

    # chart = alt.Chart(df).mark_point().encode(
    #     x=alt.X('Beak Length (mm):Q', scale=alt.Scale(zero=False)),
    #     y=alt.Y('Beak Depth (mm):Q', scale=alt.Scale(zero=False)),
    #     color=alt.condition(brush, 'Species:N', alt.value('lightgray'))
    # ).properties(
    #     width=250,
    #     height=250
    # ).add_selection(
    #     brush
    # )

    vega_pane = pn.pane.Vega(chart, debounce=10)

    def filtered_table(selection):
        if selection is None:
            return '## No selection'
        query = ' & '.join(
            f'{crange[0]:.3f} <= `{col}` <= {crange[1]:.3f}'
            for col, crange in selection.items()
        )
        # return query
        sels = set(df.query(query)['selection'])
        df.loc[df.selection.isin(sels), :]
        out = df.groupby(['desc_short','selection'])['abundance'].apply(sparkline, plot=plot_abundance_sparkline, figsize=(2,0.25), marker='.')
        return pn.pane.DataFrame(out, width=400, escape=False)

    return pn.Row(
        vega_pane,
        pn.Column(
            pn.bind(filtered_table, vega_pane.selection.param.brush),
            scroll=True, width=400, height=300
        )
    )

def summarize_clonotype(feature):
    return ex.viz.summarize_clonotypes(feature)

autocomplete = pn.widgets.AutocompleteInput(
    name='CDR3', options=list(ft.var_names.values), value=feature,
    placeholder='Select CDR3')
    
log_scale = pn.widgets.Checkbox(name='log_scale')
pn.Column(
    autocomplete,
    log_scale,
    # pn.bind(summarize_clonotype, autocomplete),
    pn.bind(make_vega_pane, autocomplete, log_scale),
    sizing_mode='stretch_width'
).servable()





from altair_transform import transform_chart
def top_asv_barplot_hist(ex, query, n=100, space='cdr3', select_from_round=8, x='r:O'):
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
    
    
    df2 = nbseq.ft.fortify(ex.fts[space][:,df2_features], feature_col='feature')
    df2['log_abundance'] = np.log10(df2['abundance'])


    # _features = pd.concat([df['feature'], df2['feature']]).unique()
    _features = list(features_by_mean_abundance.index.values)# + ['others']
    feature_scale = alt.Scale(domain=_features, range=[nbseq.viz.utils.hash_to_color(f) for f in _features])
    
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
    

def top_asv_barplot_alt(ex, query, space='cdr3', n=30, select_from_round=8, x='r:O', fancy_sorting=False, phylo=False, **kwargs):

    df = fortify_top_asvs(ex, query, space, n, select_from_round=select_from_round)
    
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
                                  range=[nbseq.viz.utils.hash_to_color(f) for f in _features])
        
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



# ******************


def whittaker_plot(df, n_head=100, n_sample=50, n_tail=100, selector=None):
    dff = pd.concat([
        df.groupby(['ID', 'round']).head(n_head),
        df.groupby(['ID', 'round']).sample(n_sample),
        df.groupby(['ID', 'round']).tail(n_tail)
    ])

    chart = alt.Chart(dff).mark_point().encode(
        x='rank',
        y='log_abundance',
        color=alt.Color('round', scale=alt.Scale(scheme='viridis')),
        tooltip=['CDR3ID', 'round', 'rank', 'abundance']
    ).interactive()

    if selector is not None:
        # chart = chart + alt.Chart().mark_point().encode(
        #     x = alt.Expr()
        #     y
        # )
        pass

    return chart


def feature_traceplot(df, detail='feature', tooltip=['feature', 'name', 'description'], selector=None):
    df_abd = df[list(set(['feature', 'r', 'abundance'] + tooltip))].query("feature != 'others'")
    # feature_scale = alt.Scale(domain=_features, range=[nbseq.viz.utils.hash_to_color(f) for f in _features])

    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.1))
    else:
        _opacity = alt.value(1)

    chart = alt.Chart(df_abd).mark_line(point=True).encode(
        x=alt.X('r:O'),
        y=alt.Y('abundance:Q'),
        color=alt.Color('feature', scale=alt_scale_features(
            df['feature'].unique())),
        opacity=_opacity,
        detail=detail,
        tooltip=['feature', 'name', 'description']
    )

    if selector is not None:
        chart = chart.add_params(selector)
    return chart


def feature_barplot(df, x='r:O', fancy_sorting=False, selector=None, phylo=False, **kwargs):
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
        if selector == True:
            selector = alt.selection_point(fields=['feature'], bind='legend', on='click', clear='dblclick')
        if selector is not None:
            _opacity = alt.condition(selector, alt.value(1), alt.value(0.1))
        else:
            _opacity = alt.value(1)

        feature_scale = alt_scale_features(_features)
        


        chart = (alt.Chart(df)
         .mark_bar()
             .encode(x=x, 
                     y='abundance', 
                     tooltip=[
                         'feature',
                         alt.Tooltip('abundance', format=".2g"),
                         alt.Tooltip('log_abundance', format=".2f")],
                     opacity=_opacity,
                     # order=(alt.Order('mean_log_abundance',sort='ascending') if fancy_sorting else None),
                     color=alt.Color('feature:N',
                                     legend=alt.Legend(columns=n//20,symbolLimit=0,labelExpr="slice(datum.value,0,6)"),
                                     scale=feature_scale,
                                     # sort=(_features if fancy_sorting else None)
                                    ))
         # .facet(column='description')

        if selector is not None:
            chart = chart.add_params(selector)
        return chart
)





# sample_data = df_f.groupby('name').first()


# nbseq.utils.sample_metadata_to_selection_metadata(ex.obs.query("FliC == 1 & io == 'i' & kind == '+'"))




def make_vhh_overview(df, point_selector=None, feature_scale=None):

    selector = alt.selection_interval(name='brush', encodings=['x', 'y'], resolve='intersect', 
                                      on="[mousedown[event.altKey], mouseup] > mousemove",
                                      translate="[mousedown[event.altKey], mouseup] > mousemove",)

    if point_selector is None:
        point_selector = alt.selection_point(
            fields=['feature'], on="click[!event.altKey]", clear="dblclick[!event.altKey]")

    # selector = alt.param(name='brush', select=dict(type='interval', encodings=['x','y'], fields=['feature']))
    # selector = alt.selection_point(fields=['feature'])

    if feature_scale is None:
        feature_scale = alt_scale_features(df['feature'].unique())
    _opacity = alt.condition(selector & point_selector,
                             alt.value(1), alt.value(0.2))
    _size = alt.condition(selector & point_selector,
                          alt.value(100), alt.value(20))

    df['pct_prod_nlogp_neg'] = clamp_finite(df['pct_prod_nlogp_neg'])
    df['pct_prod_nlogp_pos'] = clamp_finite(df['pct_prod_nlogp_pos'])
    df = df.query("binary_pvalue < 1")

    chart1 = alt.Chart(df).mark_point().encode(
        x=alt.X('f_samples_pos_sig_jitter:Q',
                title='VHH enriched in %samples', axis=alt.Axis(format='%')),
        y=alt.Y('binary_nlogp_jitter:Q', title='-log10(binary p-value)'),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'f_samples_pos_sig', 'n_samples_pos_sig',
                 'binary_pvalue', 'pct_prod_nlogp_pos', 'pct_prod_nlogp_neg']
    ).properties(width=200, height=200)

    chart2 = alt.Chart(df).mark_point().encode(
        x=alt.X('pct_prod_nlogp_pos:Q', scale=alt.Scale(clamp=True)),
        y=alt.Y('pct_prod_nlogp_neg:Q', scale=alt.Scale(clamp=True)),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('pct_prod_pos', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_pos', format='.2g'),
                 alt.Tooltip('n_samples_finite_pos', format='.0f'),
                 alt.Tooltip('pct_prod_neg', format='.2g'),
                 alt.Tooltip('pct_prod_pvalue_neg', format='.2g'),
                 alt.Tooltip('n_samples_finite_neg', format='.0f'),
                 ],
    ).properties(width=200, height=200)

    chart3 = alt.Chart(df).mark_point().encode(
        # scale=alt.Scale(clamp=True)),
        x=alt.X('mean_enrichment_pos:Q', scale=alt.Scale(type='log')),
        y=alt.Y('mean_enrichment_neg:Q', scale=alt.Scale(type='log')),
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

    _extents = df[['mean_start', 'mean_end']].min(
    ).min(), df[['mean_start', 'mean_end']].max().max()
    chart4 = alt.Chart(df).mark_point().encode(
        x=alt.X('mean_start:Q', scale=alt.Scale(type='log', domain=_extents,
                padding=+2), axis=alt.Axis(format='.0e')),  # scale=alt.Scale(clamp=True)),
        y=alt.Y('mean_end:Q',  scale=alt.Scale(
            type='log', domain=_extents, padding=+2), axis=alt.Axis(format='.0e')),
        fill=alt.Color('feature:N', scale=feature_scale, legend=None),
        color=alt.value(None),
        opacity=_opacity,
        size=_size,
        tooltip=['feature', 'mn',
                 alt.Tooltip('mean_start', format='.2g'),
                 alt.Tooltip('mean_end', format='.2g'),
                 alt.Tooltip('mean_enrichment_pos', format='.2g')],
    ).properties(width=200, height=200)
    # + alt.Chart(pd.DataFrame(dict(mean_enrichment_pos=[1, _max],mean_enrichment_neg=[1, _max]))).mark_line().encode(
    #     x=alt.X('mean_enrichment_pos:Q',scale=alt.Scale(type='log')),
    #     y=alt.Y('mean_enrichment_neg:Q',scale=alt.Scale(type='log'))
    # )

    return (chart1 | chart2 | chart3 | chart4).add_params(selector, point_selector).transform_calculate(
        f_samples_pos_sig_jitter="datum.f_samples_pos_sig+0.02*random()",
        binary_nlogp_jitter="datum.binary_nlogp+0.1*random()"
    )


def feature_traceplot(
        df, detail='feature', tooltip=['feature', 'name', 'description'],
        features=None,
        selector=None, feature_scale=None):
    df
    # feature_scale = alt.Scale(domain=_features, range=[nbseq.viz.utils.hash_to_color(f) for f in _features])

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
        color=alt.Color('feature:N', scale=feature_scale),
        opacity=_opacity,
        detail=detail,
        tooltip=[
            'feature:N',
            'mn:N',
            'description:N',
            alt.Tooltip('enrichment:Q', format=".2g"),
            alt.Tooltip('log_enrichment:Q', format=".2f")
        ],

        facet=alt.Facet('name:N', title=None, header=alt.Header(
            labels=False, labelPadding=0, labelFontSize=0))
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
            alt.Tooltip('log_abundance:Q', format=".2f")
        ],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.1)),
        color=alt.Color('feature:N',
                        legend=alt.Legend(
                            columns=n//20, symbolLimit=0, labelExpr="slice(datum.value,0,6)"),
                        scale=feature_scale
                        ),
        facet=alt.Facet('name:N', title=None, header=alt.Header(
            labels=False, labelPadding=0, labelFontSize=0))
    ).properties(width=100, height=200)

    return barplot

# chart = alt.vconcat(traceplot, barplot).add_params(selector)

# df = df_samples.rename(columns={feature_col:'feature'})


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
        facet=alt.Facet('name', header=alt.Header(labels=False), title=None)
    ).properties(width=100, height=100).interactive()
    return whittakerplot



def selection_group_dashboard(ex, enr_comparison=('R5i', 'R8i')):

    def setup_selection_group_dashboard(global_query="kind == '+' & io == 'i'", space='cdr3'):

        # setup full dataset
        # --------------------------------------------------------
        # subset feature table according to global_query
        feature_col = nbseq.asvs.get_identifier(space)
        ft = ex.query(global_query, space=space, axis='sample')

        # calculate enrichment df, add p-values and start/end abundance
        df_enr = nbseq.select.enr(ft, method='df', comparison=enr_comparison,
                                dropna=True, add_log=True)
        df_enr = nbseq.select.enrichment_pvalues(
            df_enr, abundance_col=enr_comparison[0], inplace=True)
        df_enr = nbseq.select.enrichment_start_end(df_enr, inplace=True)

        # make several feature tables grouped by selection
        sel_fts = {
            col: nbseq.ft.dataframe_to_anndata(df_enr, obs=ex.selection_metadata, obs_col='name', var_col=feature_col, value_col=col) for col in ['enrichment', 'log_enrichment', 'sig', 'p_value', 'R5i', 'R8i', 'start', 'end']
        }
        ft_enr = sel_fts['enr'] = sel_fts['enrichment']

        def setup_selection_group_dashboard_phenotype(phenotype):
            # setup subset
            pos_query = f"{phenotype} == 1"
            neg_query = f"({phenotype} == 0 | {phenotype} == -1)"

            ft_pos = nbseq.ft.to_relative(
                nbseq.ft.query(ft, pos_query, axis='sample')
            )

            # one row per sample per feature; for Whittaker plots
            df_samples = nbseq.ft.fortify(ft_pos, feature_col=feature_col, obs=True)
            df_selections = nbseq.utils.sample_metadata_to_selection_metadata(df_samples)

            # process features to calculate enrichment statistics
            df_bin = nbseq.select.compare_binary_phenotypes(
                sel_fts['sig'], phenotypes=[phenotype], plot=False, plot_sparkline=True, n_boots=1000)
            df_bin['nlogp'] = -np.log10(df_bin['p_value'])
            df_mean = nbseq.select.calculate_specificities_mean(
                ex.rfts.cdr3, ft_enr, q1=pos_query, q2=neg_query, suffixes=('_pos', '_neg'))
            df_pct = nbseq.select.calculate_specificities_pct_product(
                ft_enr, sel_fts['pct'], q1=pos_query, q2=neg_query, suffixes=('_pos', '_neg'), plot_simulation=False)
            # calculate mean start/end abundance per feature
            df_start_end = (df_enr
                            .loc[df_enr['name'].isin(df_selections.index), [feature_col, 'start', 'end']]
                            .groupby(feature_col)
                            .agg(scipy.stats.mstats.gmean)
                            .rename(columns=lambda x: f'mean_{x}'))

            # join various metrics together; one row per feature
            df_features = (
                df_bin.rename(columns={
                    'p_value': 'binary_pvalue',
                    'nlogp': 'binary_nlogp',
                    'n_samples': 'n_samples_pos_sig',
                    'f_samples': 'f_samples_pos_sig',
                    'plot': 'plot_binary'
                })
                .join(
                    df_pct.rename(columns={
                        'pct_product_pos': 'pct_prod_pos',
                        'p_value_pos': 'pct_prod_pvalue_pos',
                        'nlogp_pos': 'pct_prod_nlogp_pos',
                        'nsamples_pos': 'n_samples_finite_pos',
                        'pct_product_neg': 'pct_prod_neg',
                        'p_value_neg': 'pct_prod_pvalue_neg',
                        'nlogp_neg': 'pct_prod_nlogp_neg',
                        'nsamples_neg': 'n_samples_finite_neg',
                    }),
                    on=feature_col)
                .join(df_mean, on=feature_col)
                .join(df_start_end, on=feature_col)
            )
            df_features['mn'] = df_features[feature_col].apply(hash2mn_short)

            # find relevant features
            features = df_features.query('binary_pvalue < 1')[feature_col]

            # get fortified dataframes
            # one row per sample per feature; relevant features + top n most abundant features
            df_samples_top = (
                nbseq.ft.fortify_top_asvs(ex.fts[space], pos_query, n=100, select_from_round=8,
                                        other_features=features)
                .join(
                    df_enr.set_index(['name', feature_col])[
                        ['enrichment', 'log_enrichment']],
                    on=['name', 'feature']
                )
            )
            df_samples_top['log_abundance'] = np.log10(df_samples_top['abundance'])
            df_samples_top['mn'] = df_samples_top['feature'].apply(hash2mn_short)

            # todo: do we need these p-values and such?
            df_samples['log_abundance'] = np.log10(df_samples['abundance'])
            df_samples['rank'] = (df_samples.groupby(
                'ID')['log_abundance'].rank(ascending=False, pct=False))
            df_samples = (
                df_samples
                .join(
                    df_bin
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
                k: nbseq.ft.query(sel_fts[k], pos_query, axis='obs') for k in sel_fts
            }
            
            def make_samples_table():
                return sample_metadata_to_selection_metadata(df_samples)

            def make_features_table():
                from nbseq.utils import reorder_columns
                from nbseq.viz.utils import plot_mini_histogram, figure_to_sparkline

                df_samples_abundance = df_samples.sort_values([feature_col, 'name', 'round']).pivot(
                    index=[feature_col, 'name'], columns='round', values='abundance')
                _rounds = ['R5i', 'R6i', 'R7i', 'R8i']

                def make_traces():
                    def trace_feature(dfg, feature):
                        fig, ax = plt.subplots(figsize=(2, 0.5))
                        ax.plot(dfg[_rounds].values.T, color='C0', linewidth=0.5)
                        ax.plot(df_samples_abundance.loc[(feature, slice(
                            None)), :].values.T, color='red', marker='.', linewidth=1)
                        return figure_to_sparkline(fig)

                    traces = pd.Series(index=features, dtype='object',
                                    name=f'abundance_traces')
                    for feature, abundances in df_enr_features.groupby(feature_col)[['name'] + _rounds]:
                        traces[feature] = trace_feature(abundances, feature)
                    return traces

                def make_feature_hist(data, feature, column='log_enrichment'):
                    hist_color = nbseq.viz.hash_to_color(feature)
                    fig = plot_mini_histogram(data, color=hist_color)

                    # find abundance of samples in query
                    marks = sel_fts_q[column][:, feature].X.data

                    color = 'red'
                    for m in marks:
                        fig.axes[0].axvline(m, c=color)
                    fig.axes[0].axvline(0, c='black')
                    return figure_to_sparkline(fig, axis_off=False, subplots_adjust=True)

                def make_feature_hists(column):

                    hists = pd.Series(index=features, dtype='object',
                                    name=f'{column}_hist')
                    for feature, abundances in df_enr_features.groupby(feature_col)[column]:
                        hists[feature] = make_feature_hist(
                            abundances, feature, column=column)
                    return hists

                table = (df_features.query('binary_pvalue < 1')
                        .join(make_feature_hists('log_enrichment'), on=feature_col)
                        .join(make_feature_hists('R5i'), on=feature_col)
                        .join(make_feature_hists('R8i'), on=feature_col)
                        .join(make_traces(), on=feature_col)
                        )

                table['color'] = table[feature_col].apply(nbseq.viz.utils.pretty_hex)

                widget = pn.widgets.Tabulator(
                    # table[['log_enrichment',
                    #        'n_samples', 'f_samples', 'p_value']],
                    reorder_columns(table, ['color', 'mn', 'log_enrichment_hist', 'R5i_hist',
                                    'abundance_traces', 'R8i_hist', 'binary_pvalue']).set_index(feature_col),
                    titles={
                    },
                    formatters={
                        'log_enrichment_hist': {'type': 'html'},
                        'R5i_hist': {'type': 'html'},
                        'abundance_traces': {'type': 'html'},
                        'R8i_hist': {'type': 'html'},
                        'plot_binary': {'type': 'html'},
                    },
                    sorters=[
                        {'field': 'binary_pvalue',     'dir': 'asc'},
                        {'field': 'n_samples_pos_sig', 'dir': 'desc'}
                    ],
                    disabled=True,
                    configuration={
                        'columnDefaults': {
                            'titleFormatter': 'html',
                        },
                    }
                )
                # widget.style.applymap(lambda c: f"background-color: #{nbseq.viz.utils.pack_hex(c)}; color:#ffffff", subset=['color'])
                return widget

            def make_overview():
                df = df_samples_top.rename(columns={feature_col: 'feature'})

                selector = alt.selection_point(name='feature_selection', fields=[
                                            'feature'], bind='legend', on='click', clear='dblclick')
                features = df['feature'].unique()
                feature_scale = alt_scale_features(features)

                chart = alt.vconcat(
                    make_vhh_overview(
                        df_features.rename(columns={feature_col: 'feature'}),
                        feature_scale=feature_scale,
                        point_selector=selector
                    ),
                    whittaker_plot(
                        df_samples.rename(columns={feature_col: 'feature'}), selector=selector),
                    alt.vconcat(
                        feature_traceplot(
                            df, selector=selector,
                            features=features, feature_scale=feature_scale).transform_filter(selector),
                        feature_barplot(
                            df, selector=selector,
                            features=features, feature_scale=feature_scale)
                    ),
                    resolve=alt.Resolve(scale={'color': 'independent'})
                ).add_params(selector)  # .configure_header(labels=False, title=None)

                return chart

            phenotype_selector = pn.widgets.AutocompleteInput(
                name='phenotype', options=ex.pheno_names, placeholder='Phenotype to compare'
            )

            overview = pn.pane.Vega(make_overview(), debounce=100)
            table = make_features_table()


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
            return pn.Column(overview, table, sizing_mode='stretch_width')
            
        return pn.Column(
            phenotype_selector, 
            pn.bind(setup_selection_group_dashboard_phenotype, phenotype_selector), 
            sizing_mode='stretch_width'
        )

    global_query = pn.widgets.TextInput(
        name='global_query', value="expt == '027i' & io == 'i' & kind == '+'"
        placeholder='Filter entire dataset')

    space = pn.widgets.Select(
        name='space', options=['cdr3','aa'], value='cdr3'
    )

    # neg_query = pn.widgets.TextInput(
    #     name='neg_query', value="expt == '027i' & io == 'i' & kind == '+'"
    #     placeholder='Select antigen-negative samples')



    pn.Column(
        global_query,
        space,
        # pn.bind(summarize_clonotype, autocomplete),
        # pn.bind(make_overview, autocomplete),
        # pn.bind(make_table, autocomplete),
        pn.bind(setup_selection_group_dashboard, global_query, space),
        sizing_mode='stretch_width'
    ).servable()


def enrichment_binary_volcano_plot(df, selector=None, feature_scale=None, legend=None):
    _opacity = alt.value(1)
    if selector is not None:
        _opacity = alt.condition(selector, alt.value(1), alt.value(0.05))
    if feature_scale is None:
        feature_scale = alt_scale_features(df['feature'].unique())
    return alt.Chart(df.query("binary_pvalue < 1")).mark_point().encode(
        x=alt.X('f_samples_pos_sig_jitter:Q',
                title='VHH enriched in %samples', axis=alt.Axis(format='%')),
        y=alt.Y('binary_nlogp_jitter:Q', title='-log10(binary p-value)'),
        fill=alt.Color('feature:N', scale=feature_scale, legend=legend),
        color=alt.value(None),
        opacity=_opacity,
        tooltip=['feature', 'f_samples_pos_sig', 'n_samples_pos_sig',
                 'binary_pvalue', 'pct_prod_nlogp_pos', 'pct_prod_nlogp_neg']
    ).transform_calculate(
        f_samples_pos_sig_jitter="datum.f_samples_pos_sig+0.02*random()",
        binary_nlogp_jitter="datum.binary_nlogp+0.1*random()"
    )
