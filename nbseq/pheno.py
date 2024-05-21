import pandas as pd

from .ft import (
    sum as ft_sum,
    filter as ft_filter,
    filter_ft_abundance_prevalence,
    get_ids as ft_ids,
    transform as ft_transform,
    query as ft_query,
)
from .utils import sample_metadata_to_selection_metadata


def compare_binary_phenotypes(ft, phenotypes, pos_query=None, neg_query=None, n_boots=100):
    from numpy.random import default_rng
    rng = default_rng()

    p_values = np.zeros(shape=(len(phenotypes), ft.shape[1]))

    for i, phenotype in enumerate(phenotypes):
        if pos_query is not None:
            ph_pos = ft.obs.eval(pos_query, inplace=False).iloc[:,0]
        elif pos_query is None:
            ph_pos = ft.obs[phenotype] == 1

        if neg_query is not None:
            ph_neg = ft.obs.eval(neg_query, inplace=False).iloc[:,0]
        elif neg_query is None:
            ph_neg = ~ph_pos

        # ph_pos = ft.obs[phenotype] == pos
        # ph_neg = ~ph_pos

        # extract feature table containing positive and negative features
        ftf_pos = ft[ph_pos, :]
        ftf_neg = ft[ph_neg, :]

        n_pos = sum(ph_pos)
        n_neg = sum(ph_neg)

        p_features = ft.shape[1]

        # choose `n_pos` random samples from ftf_neg (e.g. the antigen- samples).
        # do this `n_boots` times.
        Xn_sample = rng.integers(low=0, high=n_neg, size=(n_boots, n_pos))

        # store results
        Xn_sums = np.zeros(shape=(n_boots, p_features))

        # for each bootstrap iteration:
        for j in range(n_boots):
            Xn = ftf_neg[Xn_sample[j,:],:]

            # calculate how many of the randomly-chosen Ag- samples meet the criteria
            Xn_sums[j,:] = Xn.X.sum(axis=1)

        # for each feature, calculate how many ag+ samples meet the critera
        Xp_sum = ft_sum(ftf_pos, axis='var')

        # for each feature, calculate the probability that the observed sum exceeds
        # those in the simulation
        p_values[i,:] = (Xp_sum > Xn_sums).mean(axis=1)

    return pd.DataFrame(p_values, index=phenotypes, columns = get_ids(ft, axis='var'))


def compare_phenotypes_resample(ft, phenotypes, n_boots=10, space='cdr3'):
    from numpy.random import default_rng
    rng = default_rng()

    pos = 1

    for phenotype in phenotypes:
        ph_pos = ft.obs[phenotype] == pos
        ph_neg = ~ph_pos

        ftf_pos = ft[ph_pos, :]
        ftf_neg = ft[ph_neg, :]

        for feature in get_ids(ft, axis='var'):
            Xp = ftf_pos[:,feature].X
            Xn = ftf_neg[:,feature].X

            nnz_p = Xp.nnz

            # sample nnz_p random elements from Xn
            Xn_sample = rng.choice(Xn.data, size=(n_boots, nnz_p))




def compare_phenotypes(ft, phenotypes, space='cdr3', min_abundance=10, min_prevalence=2, verbose=True):
    """compare samples between several phenotypes
    """

    results = []

    # transform to log(1+n)
    _ft = ft_transform(ft)

    metadata = ft.obs


    for i, expt, phage_library in metadata[['expt','phage_library']].drop_duplicates().itertuples():

        # subset
        ft = ft_query(_ft, f"expt == '{expt}' & phage_library == '{phage_library}'")

        # filter feature table to consider only features > 1 abundance in this subset
        ft = filter_ft_abundance_prevalence(ft,
            min_abundance=min_abundance,
            min_prevalence=min_prevalence)

        if verbose: print(f"- {repr(ft)}")

        # get metadata for subset
        meta = ft.obs #metadata.loc[metadata['ID'].isin(set(ft.ids('sample'))),:]

        for phenotype in phenotypes:

            ph_plu = meta.loc[meta[phenotype] == +1, :]
            ph_neg = meta.loc[meta[phenotype] == -1, :]
            ph_unk = meta.loc[meta[phenotype] ==  0, :]

            if verbose: print(f"  - {phenotype}: {len(ph_plu)}+, {len(ph_neg)}-, {len(ph_unk)}?")

            ids_plu = ph_plu.ID
            ids_neg = ph_neg.ID
            ids_unk = ph_unk.ID

            ftf_plu = ft_filter(ids_plu, axis='sample', inplace=False)
            ftf_neg = ft_filter(ids_neg, axis='sample', inplace=False)
            ftf_unk = ft_filter(ids_unk, axis='sample', inplace=False)


            def add_results_mannwhitneyu(x, y, comparison):
                test = scipy.stats.mannwhitneyu(x, y)
                # common language effect size
                f = (test.statistic / (len(x) * len(y)))
                rank_biserial_correlation = (f - (1 - f))
                results.append({
                    'expt': expt,
                    'phage_library': phage_library,
                    'phenotype': phenotype,
                    'feature': feature,
                    'pvalue': test.pvalue,
                    'n1': len(x),
                    'n2': len(y),
                    'u1': test.statistic,
                    'u2': len(x)*len(y)-test.statistic,
                    'f': f,
                    'gmean_difference': sparse_gmean_nonzero(x) - sparse_gmean_nonzero(y),
                    'mean_difference': np.average(x) - np.average(y),
                    'median_difference': np.median(x) - np.median(y),
                    'rank_biserial_correlation': rank_biserial_correlation,
                    'comparison': comparison,
                    'test':'mannwhitneyu'
                })

            def add_results_ttest(x, y, comparison):
                test = scipy.stats.ttest_ind(x, y, equal_var=False)
                results.append({
                    'expt': expt,
                    'phage_library': phage_library,
                    'phenotype': phenotype,
                    'feature': feature,
                    'pvalue': test.pvalue,
                    'n1': len(x),
                    'n2': len(y),
                    't': test.statistic,
                    'gmean_difference': sparse_gmean_nonzero(x) - sparse_gmean_nonzero(y),
                    'mean_difference': np.average(x) - np.average(y),
                    'median_difference': np.median(x) - np.median(y),
                    'comparison': comparison,
                    'test':'ttest_ind'
                })

            for feature in ft_ids(ft, axis='observation'):
                try:
                    if ftf_plu.shape[1] > 0:
                        if ftf_neg.shape[1] > 0:
                            add_results_mannwhitneyu(
                                ftf_plu[:,feature].X,
                                ftf_neg[:,feature].X,
                                comparison='+/-')
                        if ftf_unk.shape[1] > 0:
                            add_results_mannwhitneyu(
                                ftf_plu[:,feature].X,
                                ftf_unk[:,feature].X,
                                comparison='+/0')
    #                 print(f"{feature}\t{plus_vs_neg}")
                except(ValueError):
                    continue
        if verbose: print("")

    df = pd.DataFrame(results)
    df['-log(pvalue)'] = -np.log(df['pvalue'])
    df.index.rename('index',inplace=True)

    return df


def summarize_phenotypes(metadata, phenotypes):
    summary = pd.DataFrame(index = phenotypes.name, columns=['selections','rounds','pos','neg'])

    for phenotype in phenotypes.name:
        phenotype_not_na = ~metadata.loc[:,phenotype].isna()
        summary.loc[phenotype,'selections'] = metadata.loc[phenotype_not_na,'selection'].nunique()
        summary.loc[phenotype,'rounds'] = phenotype_not_na.sum()

        # print(f"{phenotype}: {summary.loc[phenotype,'selections']} selections / {summary.loc[phenotype,'rounds']} selection-rounds")

        pos = metadata.query(f"`{phenotype}` == 1")['selection'].nunique()
        neg = metadata.query(f"`{phenotype}` == -1")['selection'].nunique()


        summary.loc[phenotype,'pos'] = list(metadata.query(f"`{phenotype}` == 1").groupby('selection')['genotype_pair'].first())
        summary.loc[phenotype,'neg'] = list(metadata.query(f"`{phenotype}` == -1").groupby('selection')['genotype_pair'].first())

        # print(f"+ {pos:<4} {summary.loc[phenotype,'pos']}")
        # print(f"- {neg:<4} {summary.loc[phenotype,'neg']}")
        # print(f"- {neg:<4}" + str(list(metadata.query(f"`{phenotype}` == -1").groupby('selection')['genotype_pair'].first())))
        # print()
    return summary


def list_phenotypes(metadata, phenotypes):
    from IPython.display import Markdown, display

    summary = pd.DataFrame(index=phenotypes.name, columns=[
                           'selections', 'rounds', 'pos', 'neg'])

    metadata = metadata.copy()

    metadata.loc[metadata['expt'] == '027i.lib','description'] = metadata.loc[metadata['expt'] == '027i.lib', 'sample']
    metadata['description'] = metadata['description'] + \
        metadata['cond_S'].apply(lambda x: f" ({x})" if not pd.isna(x) else "")
    with pd.option_context('display.max_colwidth', None, 'display.max_rows', None, 'display.colheader_justify', 'left'):

        for phenotype in phenotypes.name:

            phenotype_not_na = ~metadata.loc[:, phenotype].isna()
            n_selections = metadata.loc[phenotype_not_na,
                                        'selection'].nunique()
            n_rounds = phenotype_not_na.sum()

            summary.loc[phenotype, 'selections'] = n_selections
            summary.loc[phenotype, 'rounds'] = n_rounds

            # print(f"{phenotype}: {summary.loc[phenotype,'selections']} selections / {summary.loc[phenotype,'rounds']} selection-rounds")

            pos = metadata.query(
                f"(`{phenotype}` == `{phenotype}`) & `{phenotype}` == 1")
            n_pos = pos['selection'].nunique()
            neg = metadata.query(
                f"(`{phenotype}` == `{phenotype}`) & `{phenotype}` != 1")
            n_neg = neg['selection'].nunique()

            display(Markdown(f"## {phenotype}: \n"
                             f"{n_selections} selections /"
                             f"{n_rounds} rounds"
                             ))

            display(
                Markdown(f'### {phenotype}<sup>+</sup> ({n_pos} selections)'))
            display(Markdown(
                pos.groupby('selection')[['description', 'cond_notes']].first().rename(columns={'description': 'Strain pair', 'cond_notes': 'Notes'}).to_markdown())
            )
            # .style.set_properties(**{'text-align': 'left'}))

            display(
                Markdown(f'### {phenotype}<sup>-/0</sup> ({n_neg} selections)'))
            display(Markdown(
                neg.groupby('selection')[['description', 'cond_notes']].first().rename(columns={'description': 'Strain pair', 'cond_notes': 'Notes'}).to_markdown())
            )
            # .style.set_properties(**{'text-align': 'left'}))

            # summary.loc[phenotype,'pos'] = list(metadata.query(f"`{phenotype}` == 1").groupby('selection')['description'].first())
            # summary.loc[phenotype,'neg'] = list(metadata.query(f"`{phenotype}` == -1").groupby('selection')['description'].first())

            # print(f"+ {pos:<4} {summary.loc[phenotype,'pos']}")
            # print(f"- {neg:<4} {summary.loc[phenotype,'neg']}")
            # print(f"- {neg:<4}" + str(list(metadata.query(f"`{phenotype}` == -1").groupby('selection')['genotype_pair'].first())))
            # print()
    return summary


def plot_ag_matrix(matrices, metadata, description_col='description', kws=None, labels=None, ax=None, figsize=None):
    from .viz.utils import trunc_ellipsis
    import matplotlib.pyplot as plt
    
    ag_matrix = matrices[0]
    descriptions = metadata.reindex(ag_matrix.index)[description_col].apply(trunc_ellipsis()).to_frame()
    descriptions['x'] = range(len(descriptions))
    # descriptions


    if ax is None:
        if figsize is None:
            figsize = (0.05*ag_matrix.shape[0], 0.1*ag_matrix.shape[1])
        fig, ax = plt.subplots(figsize=figsize)

    if kws is None:
        kws = [{}] * len(matrices)
    if labels is None:
        labels = [None]*len(matrices)

    ax.grid(True)
    ax.set_axisbelow(True)
    for ag_matrix, kwargs, label in zip(matrices, kws, labels):
        df = ag_matrix.reset_index().melt(id_vars='name', var_name='antigen')

        names = ag_matrix.index.values
        pos = df.query('value == 1')
        neg = df[(~df['value'].isna()) & (df['value'] != 1)]

        # ax.scatter(descriptions.loc[df['name'],'x'], df['antigen'], c=df['value'])
        ax.scatter(descriptions.loc[pos['name'], 'x'], pos['antigen'],
                   label=f'{label}+' if label is not None else '+', color='tab:orange', **kwargs)
        ax.scatter(descriptions.loc[neg['name'], 'x'], neg['antigen'],
                   label=f'{label}-' if label is not None else '-', color='tab:blue', **kwargs)

    ax.set_xticks(descriptions.x, descriptions[description_col], rotation=-90)
    ax.legend()

# fig, ax = plt.subplots(figsize=(30,20))
# plot_ag_matrix([ag_matrix], [dict(s=100, alpha=0.2)], ax=ax)



def plot_selection_phenotype_grid(ft, selection_metadata=None, phenotypes=None, features=None):
    pass

def selection_phenotype_grid(metadata, phenotypes, row_labels=['name','description']):
    from natsort import natsort_keygen
    
    selections = sample_metadata_to_selection_metadata(metadata).reset_index() #metadata.groupby('selection').first()

    selections = selections.sort_values(
        by=row_labels,
        key=natsort_keygen()
    )

    grid = selections.set_index(row_labels).loc[:, phenotypes.name]

    grid.columns = pd.MultiIndex.from_frame(
        phenotypes.set_index('name').loc[grid.columns,['category']].reset_index().rename(columns={'index':'phenotype'})[['category','phenotype']]
    ).map(lambda r: tuple(f"<span>{x}</span>" for x in r))

    grid = grid.sort_index(axis='columns', level=0)


    grid_style = (grid.fillna('').style
                  .format(precision=0)
                  .applymap(lambda x: {1:'background-color:green;',-1:'background-color:red;'}.get(x, None))
                  .set_sticky(axis="index", pixel_size=100)
                  .set_sticky(axis="columns", pixel_size=100)
                 )

    styles = [dict(selector="th.col_heading", props=[('vertical-align', 'bottom'), ('text-align','center')]),
              # dict(selector="th", props=[('width', '40px')]),
              dict(selector="th.col_heading.level0 span",props=[
                ('border-bottom', 'solid 1px silver'),
                ('display','block'),
                ('margin','0px 5px')
              ]),
              dict(selector="th.col_heading.level1 span",props=[
                          # ('text-align', 'left'),
                          ('writing-mode', 'vertical-rl'),
                          ('transform', 'rotate(180deg)'),
                          ('white-space','nowrap'),
                          ('vertical-align', 'top')]),
              dict(selector="th.row_heading", props=[
                          ('white-space','nowrap'),
                          ('text-align','right')
              ]),
              dict(selector='table', props=[
                  ('overflow', 'hidden')
              ]),
              dict(selector='td, th', props=[
                  ('position','relative')
              ]),
              dict(selector='td:hover::after, th:hover::after', props=[
                    ('content', '" "'),
                    ('display', 'inline-block'),
                    ('position', 'absolute'),
                    ('background-color', '#ffa'),
                    ('border', 'solid 1px silver'),
                    ('left', '0'),
                    ('top', '-5000px'),
                    ('height', '10000px'),
                    ('width', '100%'),
                    ('opacity', '40%'),
                    ('z-index', '10000')
              ])


# th:hover::after {
#   content: "";
#   position: absolute;
#   background-color: #ffa;
#   left: 0;
#   top: -5000px;
#   height: 10000px;
#   width: 100%;
#   z-index: -1;
# }
             ]

                   # props=[("writing-mode", "vertical-rl"),
                   #        ('transform', 'rotateZ(180deg)'),
                   #        ('height', '290px'),
                   #        ('vertical-align', 'top')])]
    return grid_style.set_table_styles(styles)

def plot_antigen_label_balance(ag_matrix, orient='h', ax=None):
    import matplotlib.pyplot as plt
    
    control_locs = ag_matrix.index.str.startswith('control_')
    n_controls = control_locs.sum()
    labels = ag_matrix.loc[~control_locs,:]
    
    n_labels = labels.apply(lambda c: c.value_counts())
    n_labels.loc['control',:] = n_controls
    n_labels = n_labels.fillna(0)

    if ax is None:
        fig, ax = plt.subplots(figsize=((16, 8) if orient == 'h' else (8, 16)))

    if orient == 'h':
        ax.bar(x=labels.columns, height= n_labels.loc[1.0,:], label='Antigen increased')
        ax.bar(x=labels.columns, height=-n_labels.loc['control',:], label='Contrived control')
        ax.bar(x=labels.columns, height=-n_labels.loc[0.0,:], bottom=-n_labels.loc['control',:], label='Antigen does not change')
        ax.bar(x=labels.columns, height=-n_labels.loc[-1.0,:], bottom=-n_labels.loc['control',:]-n_labels.loc[0.0,:], label='Antigen decreased')
        ax.axhline(y=0, c='black')
        plt.ylabel('# training examples')
        plt.xlabel('antigen')
        plt.xticks(rotation=-90)
    else:
        ax.barh(y=labels.columns, width= n_labels.loc[1.0,:], label='Antigen increased')
        ax.barh(y=labels.columns, width=-n_labels.loc['control',:], label='Contrived control')
        ax.barh(y=labels.columns, width=-n_labels.loc[0.0,:], left=-n_labels.loc['control',:], label='Antigen does not change')
        ax.barh(y=labels.columns, width=-n_labels.loc[-1.0,:], left=-n_labels.loc['control',:]-n_labels.loc[0.0,:], label='Antigen decreased')
        ax.axvline(x=0, c='black')
        ax.axvline(x=0, c='black')

        plt.xlabel('# training examples')
        plt.ylabel('antigen')

    plt.legend()
    plt.show()

    return n_labels