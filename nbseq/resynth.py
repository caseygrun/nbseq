from . import Experiment
import numpy as np
import pandas as pd
from .viz.utils import pre


def show_or_print(show, text, plain=None):
    if show:
        if hasattr(show, "display"):
            from IPython.display import display
            display(show.display(text))
        else:
            if plain is None:
                plain = text
            print(plain)


def make_consensus(
    fd,
    seq_col="aligned",
    abundance_col="abundance",
    min_fraction=0.15,
    min_count="auto",
    ungap=False,
    plot=True,
    features={},
    clonotype=None,
    show=True,
    library='alpaca',
    figsize=None,
    **kwargs,
):

    from .asvs import marginal_frequencies, consensus

    # TODO: replace with viz.logo.alignment_to_matrix or whatever
    matrix = marginal_frequencies(
        fd[seq_col], counts=fd[abundance_col], characters_to_ignore='X')
    if plot:
        import matplotlib.pyplot as plt

        if figsize is None:
            figsize = (None, 2)
        if figsize[0] is None:
            figsize = (len(fd[seq_col][0]) / 5, figsize[1])

        # if library is not None:
        #     features = ex.get_CDR_positions(library, aa=True)
        #     clonotype = ex.get_clonotype_positions(library, aa=True)
        # else:
        #     features = {}
        #     clonotype = None

        from .viz import logo

        _logo, mat_df = logo.make_logo(
            fd,
            seq_col=seq_col,
            count_col=abundance_col,
            matrix=matrix.drop("-", axis="columns"),
            # matrix=matrix,
            figsize=figsize,
            features=features,
            **kwargs,
        )

        if min_fraction is not None:
            min_count = min_fraction * matrix.sum(axis="columns").max()
            plt.axhline(y=min_count, color="red")

            if clonotype is not None:
                plt.axvline(clonotype[0]-0.5, ymin=0, ymax=1, color='green')
                plt.axvline(clonotype[1]+0.5, ymin=0, ymax=1, color='green')
                # plt.hlines(
                #     y=0,
                #     xmin=clonotype[0],
                #     xmax=clonotype[1],
                #     color="green",
                #     linewidth=10,
                #     alpha=0.5,
                # )
        fig = _logo.ax.figure
        from IPython.display import display
        display(fig)
        plt.close(fig)

    # for positions where abundance < min_fraction, fill with gap characters
    # rather than 'X', to guarantee consensus can be reverse translated
    cons = consensus(matrix, min_fraction=min_fraction, ambiguous="-", gap="-")

    return cons


# def compare_seqs(seq1, seq2, seq1_name, seq2_name, wildcard='X', clonotype=None, features=None, show=None):
#     diff = ''.join(' ' if (seq1[i] == seq2[i] or seq1[i] == wildcard or seq2[i] == wildcard)
#                    else '^' for i in range(len(seq1)))
#     ndiff = sum(seq1[i] != seq2[i] for i in range(len(seq1)))
#     show_or_print(
#         show, f"# {seq1_name} sequence and {seq2_name} sequence differ by {ndiff:.0f} residues:\n"
#         f"{seq1}\n{seq2}\n{diff}\n")
#     return ndiff

def compare_seqs(seq1, seq2, seq1_name, seq2_name, title_format=None, wildcard='X', ignore_initial_gaps=False, clonotype=None, features=None, show=None):
    mask = np.ones(len(seq1))
    if features is not None:
        mask[:] = 0
        for k, pos in features.items():
            mask[pos[0]:pos[1]] = 1

    if ignore_initial_gaps:
        i = 0
        while seq1[i] == '-' and i < len(seq1):
            mask[i] = 0
            i += 1
        i = len(seq1)-1
        while seq1[i] == '-' and i >= 0:
            mask[i] = 0
            i -= 1

    diff = np.array([not ((not mask[i]) or (seq1[i] == seq2[i] or seq1[i]
                    == wildcard or seq2[i] == wildcard)) for i in range(len(seq1))])
    # diff = ''.join(' ' if (not mask[i]) or (seq1[i] == seq2[i] or seq1[i] == wildcard or seq2[i] == wildcard)
    #                else '^' for i in range(len(seq1)))
    # ndiff = sum(seq1[i] != seq2[i] for i in range(len(seq1)))
    ndiff = sum(diff)
    diff_str = ''.join(np.where(diff, '^', ' '))
    if title_format is None:
        title = f"# {seq1_name} sequence and {seq2_name} sequence differ by {ndiff:.0f} residues"
    else:
        title = title_format.format(seq1_name=seq1_name, seq2_name=seq2_name, ndiff=ndiff)
    if features is not None:
        title += f" in {list(features.keys())}"

    show_or_print(show, (
        f"{title}:\n"
        f"{seq1}\n{seq2}\n{diff_str}\n"
    ))
    return ndiff


def compare_to_reference(cons, ref, clonotype=None, features=None, show=None):
    return compare_seqs(cons, ref, 'consensus', 'reference', ignore_initial_gaps=True, features=features, show=show)


def compare_to_mode(cons, mode, clonotype=None, features=None, show=None):
    return compare_seqs(cons, mode, 'consensus', 'mode', show=show)


def decorate_and_trim_consensus(cons, clonotype, features, CDRs, CDR3ID, identify=True, verbose=False, show=None):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from .utils import show_ruler, show_sequence_ranges, md5_ungapped
    from .cdrs import extract_CDRs, identify_clonotype

    if show:
        show_or_print(show,
                      "# consensus AA seq before trimming:\n" +
                      show_sequence_ranges(
                          cons, {**features, 'clonotype': clonotype}) + "\n"
                      # f"{cons}\n"
                      # f"{show_ruler(0,len(cons))}"
                      )

    # if given, trim sequence to clonotype positions and adjust CDRs/features accordingly
    if clonotype is not None:
        features = trim_features_to_clonotype(features, clonotype)
        CDRs = trim_features_to_clonotype(CDRs, clonotype)

        cons = cons[clonotype[0]: clonotype[1]]
    else:
        print("Warning: no clonotype bounds specified; unintended residues may be included")

    cons_features = extract_CDRs(
        [cons], CDRs=features, push_left=False).iloc[0, :].to_dict()

    # extract CDRs and hash to create clonotype ID
    clonotypeID = ""
    if identify:
        if len(CDRs) == 0:
            raise ValueError("identify=True but no CDR positions given.")

        cons_cdrs = extract_CDRs([cons], CDRs=CDRs, push_left=False)
        cons_cdrs = identify_clonotype(
            cons_cdrs, cdrs=CDRs.keys(), id_col="clonotypeID", seq_col=None
        )
        clonotypeID = cons_cdrs.loc[:, "clonotypeID"].iloc[0]
        aaSVID = md5_ungapped(str(cons))

    if show and (verbose > 1):
        show_or_print(show,
                      # f"{cons}\n" +
                      # f"{show_ruler(0,len(cons))}"
                      "# trimmed AA consensus seq:\n" +
                      show_sequence_ranges(cons, features) + "\n"
                      )

    if identify or len(features) > 0:
        from .utils import ungap_annotations
        return (ungap_annotations(cons, features, id=clonotypeID, annotations={
            'CDR3ID': CDR3ID,
            'clonotypeID': clonotypeID,
            'aaSVID': aaSVID,
        }), cons_features)
        # return SeqRecord(Seq(cons), id=_id)
    else:
        return (cons, cons_features)


def back_translate(
    seq, tables=["Standard"], species=["h_sapiens"], overhang5p='', verbose=True, rich=True, **kwargs
):
    """Make, optimize, and compare nucleic acid sequences that code for an amino acid sequence, using multiple translation tables and target species

    Parameters
    ----------
    seq : str or Bio.SeqRecord.SeqRecord
        Amino acid sequence
    tables : list of str
        names of translation tables to use. Will back-translate sequence with translation
        table and warn about any differences if ``verbose=True``
    species : list of str
        names of species to use for sequence optimization.
    verbose : bool
        warn about differences between translation `table`s
    rich : true
    kwargs
        Additional arguments passed to `optimize_sequence`

    Returns
    -------
    Bio.SeqRecord
        Nucleic acid sequence

    """
    import itertools
    import dnachisel
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from .utils import ungapped, ungap_seqrecord, merge_seqrecords, aa_features_to_na

    _id = None
    if isinstance(seq, SeqRecord):
        _rec = ungap_seqrecord(seq)
        seq_aa = str(_rec.seq)
        _id = _rec.id
        _name = _rec.name
        _features = aa_features_to_na(_rec.features, offset=len(overhang5p))
        # seq_aa = ungapped(str(seq.seq))
        # _id = seq.id
        # _name = seq.name
    else:
        seq_aa = ungapped(seq)

    trans = []
    # trans = [dnachisel.reverse_translate(seq, table=table) for table in tables]
    for table, spec in zip(tables, species):
        seq_nt = dnachisel.reverse_translate(seq_aa, table=table)

        seq_nt_o = optimize_sequence(
            seq_nt,
            seq_aa,
            overhang5p=overhang5p,
            codon_optimize_species=spec,
            tables=[table],
            rich=rich,
            verbose=verbose,
            **kwargs,
        )
        trans.append(seq_nt_o)

    if len(tables) > 1:
        tt = zip(tables, trans)
        for ((table1, trans1), (table2, trans2)) in itertools.combinations(tt, 2):
            if isinstance(trans1, SeqRecord):
                trans1 = str(trans1.seq)
            if isinstance(trans2, SeqRecord):
                trans2 = str(trans2.seq)
            if trans1 != trans2:
                mismatches = [c1 != c2 for c1, c2 in zip(trans1, trans2)]
                title = f"{sum(mismatches)} nt differences between {table1} / {table2}"
                body = "\n".join(
                    [
                        f"   {table1:>12} : {trans1}",
                        "                : "
                        + "".join("!" if mm else " " for mm in mismatches),
                        f"   {table2:>12} : {trans2}",
                    ]
                )
                if rich:
                    from .viz.utils import display_accordion
                    from IPython.display import HTML, display

                    display(
                        display_accordion(
                            HTML(
                                f"<pre style='white-space:pre;'>{body}</pre>"),
                            title=title,
                        )
                    )
                else:
                    print(title)
                    print(body)

    out = trans[0]
    if _id is not None:
        if not isinstance(out, SeqRecord):
            out = SeqRecord(Seq(out), id=_id, name=_name)
        else:
            if _rec is not None:
                o1 = merge_seqrecords(_rec, out,
                                      seq=out.seq,
                                      features=_features
                                      #out.features + _features
                                      )
                o1_edits = merge_seqrecords(_rec, out,
                                            seq=out.seq,
                                            # features=_features
                                            features=out.features + _features
                                            )
                return (o1, o1_edits)
                out.features = out.features + _features
            out.id = _id
            out.name = _name

    return (out, out)


def optimize_sequence(
    seq_nt,
    seq_aa,
    codon_optimize_species="e_coli",
    avoid_rare_codon_species=["e_coli", "h_sapiens"],
    min_codon_frequency=0.1,
    tables=["Standard", "Bacterial"],
    constraints=[],
    overhang5p="",
    overhang3p="",
    clamp5p=0,
    clamp3p=0,
    annotations={},
    rich=True,
    verbose=True,
):
    """Optimizes the given nucleic acid sequence subject to several constraints

    Parameters
    ----------
    seq_nt : str or Bio.Seq
        Nucleic acid sequence
    seq_aa : str or Bio.Seq
        Amino acid sequence that should be represented
    codon_optimize_species : str, optional
        Optimize codon usage for this species, by default "e_coli"
    avoid_rare_codon_species : list of str, optional
        For each of these species, avoids codons with frequency < `min_codon_frequency`; by default ["e_coli", "h_sapiens"]
    min_codon_frequency : float
        Minimum codon frequency to use
    tables : str, optional
        Enforce that the coding region of the nucleic acid sequence 
        (i.e. the portion from ``[len(overhang5p):-len(overhang_3p)]``) 
        translates to `seq_aa` using each of these translation tables
    constraints : list of dnachisel.Constraint, optional
        Additional constraints to apply
    overhang5p : str, optional
        Additional sequence to add to beginning of `seq_nt`
    overhang3p : str, optional
        Additional sequence to add to end of `seq_nt`
    clamp5p : int, optional
        Keep the first ``clamp5p`` residues of ``overhang5p + seq_nt + overhang3p`` fixed when optimizing the sequence
    clamp3p : int, optional
        Keep the last ``clamp3p`` residues of ``overhang5p + seq_nt + overhang3p`` fixed when optimizing the sequence
    annotations : dict of tuple : str, optional
        Create `Bio.SeqFeature`s within the coding sequence. Key is a tuble of ``(start, end)``, value is a label for the feature
    rich : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    NotImplementedError
        _description_
    """
    import dnachisel as dc
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    if verbose:
        plan = show_plan(NA=seq_nt, overhang5p=overhang5p,
                         overhang3p=overhang3p, clamp5p=clamp5p, clamp3p=clamp3p)
        plan = "# NA sequence before optimization:\n" + plan + "\n"
        if rich:
            from .viz.syntax import na_highlighter
            from IPython.display import display
            display(na_highlighter.display(plan))
        else:
            print(plan)

    full_seq_nt = overhang5p + seq_nt + overhang3p
    if clamp5p > 0:
        constraints = [dc.AvoidChanges(location=(0, clamp5p))] + constraints
    if clamp3p > 0:
        constraints = [
            dc.AvoidChanges(
                location=(len(full_seq_nt) - clamp3p, len(full_seq_nt)))
        ] + constraints

    translation_frame = (
        len(overhang5p) % 3,
        ((len(full_seq_nt) - (len(overhang5p) % 3)) // 3) * 3,
    )
    constraints = (
        [
            dc.AvoidPattern("NotI_site"),
            dc.AvoidPattern("AscI_site"),
            dc.AvoidPattern("XhoI_site"),
            dc.AvoidPattern("XbaI_site"),
            dc.AvoidPattern("GGGGG"),
            # EnforceGCContent(mini=0.3, maxi=0.7, window=50),
        ]
        + [
            dc.EnforceTranslation(
                translation=seq_aa,
                genetic_table=table,
                location=(len(overhang5p), len(overhang5p) + len(seq_nt)),
            )
            for table in tables
        ]
        + [
            dc.AvoidRareCodons(
                species=species, location=translation_frame, min_frequency=min_codon_frequency
            )
            for species in avoid_rare_codon_species
        ]
        + [
            dc.AvoidStopCodons(location=translation_frame),
        ]
        + constraints
    )

    problem = dc.DnaOptimizationProblem(
        sequence=full_seq_nt,
        constraints=constraints,
        objectives=[
            dc.CodonOptimize(
                species=codon_optimize_species,
                method="match_codon_usage",
                location=(len(overhang5p), len(overhang5p) + len(seq_nt)),
            )
        ],
        logger=None,
    )

    problem.resolve_constraints()
    problem.optimize()

    if rich:
        from nbseq.viz.utils import display_accordion
        from IPython.display import HTML, display

        title = f"Optimizing reverse translation for translation tables {tables}; codon optimizing for '{codon_optimize_species}'; avoiding rare codons (<{min_codon_frequency:.0%}) in species {avoid_rare_codon_species}"
        display(
            display_accordion(
                HTML(
                    f"<pre style='white-space:pre;'>{problem.constraints_text_summary()}\n{problem.objectives_text_summary()}</pre>"
                ),
                title,
            )
        )
    else:
        print(problem.constraints_text_summary())
        print(problem.objectives_text_summary())

    final_sequence = problem.sequence
    final_record = problem.to_record(with_sequence_edits=True)

    final_record.features.append(
        SeqFeature(
            FeatureLocation(0, len(overhang5p)),
            "misc_feature",
            qualifiers={"label": "overhang_5p"},
        )
    )
    final_record.features.append(
        SeqFeature(
            FeatureLocation(len(full_seq_nt) -
                            len(overhang3p), len(full_seq_nt)),
            "misc_feature",
            qualifiers={"label": "overhang_3p"},
        )
    )

    if isinstance(annotations, dict):
        for k, v in annotations.items():
            if (isinstance(v, tuple) or isinstance(v, list)) and isinstance(k, str):
                (v, k) = (k, v)
            final_record.features.append(
                SeqFeature(
                    FeatureLocation(*tuple(i + len(overhang5p) for i in k)),
                    "misc_feature",
                    qualifiers={"label": v},
                )
            )
    elif isinstance(annotations, list):
        raise NotImplementedError()

    return final_record


# cache codon frequency tables
_codon_frequency_tables = {}


def get_codon_frequency_table(species):
    if species in _codon_frequency_tables:
        return _codon_frequency_tables[species]

    import python_codon_tables

    codons = python_codon_tables.get_codons_table(species)
    o = {}
    for d in codons.values():
        for codon in d:
            o[codon] = d[codon]

    _codon_frequency_tables[species] = o
    return o


def trim_features_to_clonotype(features, clonotype_positions):
    length = clonotype_positions[1] - clonotype_positions[0]
    # return {
    #     k: (
    #         v[0] - clonotype_positions[0],
    #         min(v[1] - clonotype_positions[0], clonotype_positions[1] - clonotype_positions[0])
    #     ) for k, v in features.items() if (v[1] > clonotype_positions[0] and v[0] < clonotype_positions[0])
    # }
    return {
        k: (
            max(v[0] - clonotype_positions[0], 0),
            min(v[1] - clonotype_positions[0], length)
        ) for k, v in features.items()
    }


_resynthesize_params = {
    'trim_clonotype'
    'library',
    'global_query',
    'overwrite_consensus',
    'overwrite',
    'overwrite_rules',
    'overhang5p',
    'overhang3p',
    'clamp5p',
    'clamp3p',
}


def resynthesize(
    ex: Experiment,
    CDR3ID=None,
    aaSVID=None,
    library="alpaca",
    global_query=None,
    relative=True,
    trim_clonotype=True,
    overwrite_consensus=None,
    overwrite={},
    overwrite_rules={},
    overhang5p="",
    overhang3p="",
    clamp5p=0,
    clamp3p=0,
    rich=True,
    verbose=True,
):
    """
    resynthesize a given VHH by creating an AA consensus sequence, back-translating, optimizing, and adding cloning overhangs

    Parameters
    ----------
    overwrite : dict of tuple : str
        Overwrite certain parts of the consensus sequence
        each key is a tuple of ``(start, end)`` and each value is a string that should replace that range.
        Raises `ValueError` if ``len(value) != key[1] - key[0]``


    Returns
    -------
    int
        Description of anonymous integer return value.


    """

    from Bio.SeqRecord import SeqRecord
    from .utils import extract_features, show_sequence_ranges, copy_record
    from .ft import fortify_features

    if rich:
        from IPython.display import display, HTML
        from .viz.syntax import aa_highlighter, na_highlighter

    features = None
    if CDR3ID is not None:
        feature = CDR3ID
        # _ft = project(
        #     features=[feature], ft=ft, columns=["aligned", "aaSVID"], library=library
        # )
        if rich:
            display(ex.viz.tag(feature))

        if verbose > 1:
            print("collecting feature sequences:")
        _ft = ex.project([feature], ft=True, relative=relative, query=global_query, from_space='cdr3', to_space='aa', columns=["aligned", "aaSVID"], library=library,
                         show_query=(verbose > 1), verbose=(verbose > 1))
        fd = fortify_features(_ft)

        n_clonotypes = len(fd)
        n_reads = int(fd['abundance'].sum())
        if verbose:
            display(
                f"using {n_clonotypes} features with {n_reads:.0f} total reads")
        cdr3id = CDR3ID  # find_cdr3(CDR3ID, single=True)
        _id = f"CDR3_{cdr3id}"
        name = f"CDR3_{cdr3id[0:6]}"

        features = ex.get_feature_positions(library=library, aa=True)

        mode_sequence = None
        n_modes = 5
        n_modes = min(n_modes, len(fd))

        # fd could be long so avoid sorting the whole thing
        # fds = fd.sort_values('abundance', ascending=False)
        fds = fd.nlargest(n_modes, 'abundance', keep='all')
        seq_len = len(fds['aligned'].iloc[0])

        _mode_summary = ""
        mode_sequence = fds['aligned'][0]  # .iloc[0, 'aligned']

        if verbose:
            for i in range(n_modes):
                mode = fds.iloc[i, :]
                _abundance = f"({mode['abundance']:.2g} abundance)" if relative else f"({mode['abundance']:.0f} reads)"
                _mode_summary += f">{mode.name} {_abundance}\n{mode['aligned']}\n"
            show_or_print(
                aa_highlighter if rich else None,
                (f"# top {n_modes} most abundant sequence(s):\n" +
                 show_sequence_ranges(seq_len, features) +
                 _mode_summary)
            )

        # # adjust indices of features if sequence has been trimmed to clonotype region
        # if trim_clonotype:
        #     cp = ex.get_clonotype_positions(library=library, aa=True)
        #     features = {k: [v[0] - cp[0], v[1] - cp[1]] for k, v in features.items()}
    else:
        raise NotImplementedError("Must provide CDR3ID")

    clonotype = ex.get_clonotype_positions(library=library, aa=True)

    # print(f"{feature} ({library})")
    if overwrite_consensus is not None:
        cons = overwrite_consensus
        show_or_print(aa_highlighter if rich else False,
                      f"# using fixed AA sequence instead of consensus:\n{cons}"
                      )
    else:
        cons = make_consensus(
            fd,
            # library=library,
            clonotype=clonotype,
            features=features,
            show=aa_highlighter if rich else False,
            plot=rich,
        )

    _show = aa_highlighter if rich else None
    n_diffs_mode = None
    if mode_sequence is not None:
        n_diffs_mode = compare_to_mode(cons, mode_sequence, show=_show)

    n_diffs_ref = compare_to_reference(cons, ex.reference('cdr3', library=library, aa=True),
                                       features={
                                           k: v for k, v in features.items() if k in ['FR2', 'FR3']},
                                       show=_show)

    # apply overwrite rules:
    rules_applied = []
    _show = aa_highlighter if rich else None

    from .cdrs import extract_features_dict

    # get all feature sequences for overwrite rules to access
    features_aa = extract_features_dict(cons, features)
    for name, rule in overwrite_rules.items():
        new_cons = rule(cons, features, **features_aa)
        if new_cons != cons:
            rules_applied.append(name)
            if verbose:
                compare_seqs(cons, new_cons, 'cons', 'edited',
                    title_format="# overwrite rule '" + name + "' changed consensus by {ndiff:.0f} residues", 
                    show=_show)
            cons = new_cons
            features_aa = extract_features_dict(new_cons, features)


    cons_rec, features_aa = decorate_and_trim_consensus(cons,
                                                        features=features,
                                                        clonotype=clonotype,
                                                        CDRs=ex.get_CDR_positions(
                                                            library=library, aa=True),
                                                        CDR3ID=CDR3ID,
                                                        identify=True,
                                                        show=aa_highlighter if rich else False,
                                                        verbose=verbose,
                                                        )

    # for name, rule in overwrite_rules.items():
    #     new_cons_rec = copy_record(cons_rec)
    #     new_cons_rec = rule(new_cons_rec, **features_aa)
    #     if new_cons_rec.seq != cons_rec.seq:
    #         rules_applied.append(name)
    #         if verbose:
    #             show_or_print(_show,
    #                 compare_seqs(cons_rec.seq, new_cons_rec, 'cons', 'edited',
    #                 title_format="# overwrite rule '" + name + "' changed consensus by {ndiff:.0f} residues" ))
    #     new_cons_rec = cons_rec
    #     features_aa = extract_features(new_cons_rec)

    # apply overwrites
    cons_rec = overwrite_sequence(
        cons_rec, overwrite, features=features, rich=rich, verbose=rich)

    bt, bt_edits = back_translate(
        cons_rec,
        rich=rich,
        overhang5p=overhang5p,
        overhang3p=overhang3p,
        clamp5p=clamp5p,
        clamp3p=clamp3p,
    )  # , annotations=get_CDR_positions(library, aa=False))

    features_na = extract_features(bt)

    preview_back_translation(feature, cons_rec, bt_edits,
                             overhang5p=overhang5p, overhang3p=overhang3p, show=na_highlighter if rich else None)

    return {
        "clonotypeID": bt.id,
        "aaSVID": cons_rec.annotations['aaSVID'],
        "CDR3ID": CDR3ID,
        "NA": str(bt.seq),
        "AA": str(cons_rec.seq),
        "cons_aa": str(cons),
        "rec_aa": cons_rec,
        "rec_na": bt,
        "mode": mode_sequence,
        "mode_vs_consensus": n_diffs_mode,
        "cons_vs_ref": n_diffs_ref,
        "n_clonotypes": n_clonotypes,
        "n_reads": n_reads,
        "rules_applied": rules_applied,
        **features_aa,
        **{f'{k}_na': v for k, v in features_na.items()}
    }


def show_plan(**rVHH):
    from .utils import complement, rc

    overhang5p = rVHH['overhang5p'].lower() if 'overhang5p' in rVHH else ''
    c_overhang5p = complement(overhang5p)

    c_overhang3p = rc(rVHH['overhang3p'].lower()
                      ) if 'overhang3p' in rVHH else ''
    overhang3p = complement(c_overhang3p)

    if 'NA' in rVHH:
        seq = rVHH['NA']
        c_seq = complement(seq)
    else:
        seq = c_seq = '...VHH...'

    full_seq = f"{overhang5p} {seq} {overhang3p}"
    full_seq_c = f"{c_overhang5p} {c_seq} {c_overhang3p}"

    clamp5p = rVHH['clamp5p'] if 'clamp5p' in rVHH else 0
    clamp3p = rVHH['clamp3p'] if 'clamp3p' in rVHH else 0

    return (
        full_seq + "\n" +
        full_seq_c + "\n" +
        ("^" * clamp5p) + " " + (" " * (len(full_seq) -
                                        clamp5p - clamp3p - 2)) + " " + ("^" * clamp3p)
    )


def preview_back_translation(feature, cons, bt, overhang5p='', overhang3p='', min_codon_frequency=0.1, show=None):
    from Bio.SeqRecord import SeqRecord
    from .utils import ungapped
    from .viz.utils import subplots
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from IPython.display import display

    if isinstance(bt, SeqRecord):
        if isinstance(cons, SeqRecord):
            cons_seq = ungapped(cons.seq)
        else:
            cons_seq = ungapped(cons)

        from dna_features_viewer import (
            GraphicFeature,
            GraphicRecord,
            BiopythonTranslator,
        )

        graphic_record = BiopythonTranslator().translate_record(bt)
        # graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)

        reading_frame = bt.seq[len(overhang5p) % 3:]
        reading_frame = reading_frame[0: ((len(reading_frame) // 3) * 3)]
        translation_check = reading_frame.translate("Standard", to_stop=True)

        species_list = ["e_coli", "h_sapiens"]
        fig, axs = plt.subplots(
            nrows=len(species_list) + 1,
            ncols=1,
            figsize=(len(bt) // 5, (len(species_list) + 5)),
            sharex=True,
            gridspec_kw={"height_ratios": [5] + len(species_list) * [1]},
        )

        xs = (len(overhang5p) % 3) + 1 + \
            np.arange(0, len(reading_frame) // 3) * 3
        for i, species in enumerate(species_list):
            cft = get_codon_frequency_table(species)
            codon_usages = np.zeros(len(reading_frame) // 3)
            for j, c in enumerate(range(0, len(reading_frame), 3)):
                codon = reading_frame[c: c + 3]
                codon_usages[j] = cft[codon]
                y = codon_usages[j]
                if codon_usages[j] <= min_codon_frequency:
                    axs[i+1].annotate("{:.0%}".format(y),
                                      xy=(xs[j], y),
                                      ha='center', va='top', color='red')
            axs[i+1].plot(xs, codon_usages, marker='.')
            axs[i+1].axhline(y=min_codon_frequency, color="red")

            bad_residues = codon_usages < min_codon_frequency
            axs[i+1].scatter(xs[bad_residues],
                             codon_usages[bad_residues], color='red')
            axs[i+1].set_ylabel(species)

        ax = axs[:][0]
        graphic_record.plot(ax=ax, figure_width=len(bt) // 5)
        graphic_record.plot_sequence(ax)
        graphic_record.plot_translation(
            ax,
            translation=cons_seq,
            location=(len(overhang5p), len(overhang5p) + len(bt)),
            fontdict={"weight": "bold"},
        )
        graphic_record.plot_translation(
            ax,
            translation=translation_check,
            location=(
                len(overhang5p) % 3,
                len(overhang5p) % 3 + len(translation_check),
            ),
        )
        display(fig)
        plt.close(fig)

        # check that translation is the right length (i.e. premature stop codon not introduced)
        if len(translation_check)*3 != (len(overhang5p) + (len(cons_seq)*3) + len(overhang3p)):
            raise ValueError(
                "Reverse translated sequence does not produce peptide of correct length; premature stop codon may have been introduced")

        cds = translation_check[(len(overhang5p) // 3)
                                 : -(len(overhang3p) // 3)]
        assert len(cds) == len(cons_seq)
        # cds = cds[0 : len(cons_seq)]

        # check that correct sequence obtained
        assert (
            cds == cons_seq
        ), "Reverse translated sequence does not produce correct amino acid sequence"

        show_or_print(show,
                      f"# final optimized sequence:\n>{feature}\n{bt.seq}\n\n"
                      )
    else:
        print(f"\n{bt}\n")


def overwrite_sequence(cons, overwrite, features={}, rich=True, verbose=True):
    from Bio.SeqRecord import SeqRecord
    import collections.abc
    import copy

    if overwrite is not None and len(overwrite) > 0:
        if isinstance(cons, SeqRecord):
            cons_seq = str(cons.seq)
        else:
            cons_seq = str(cons)
        cons_seq_orig = copy.copy(cons_seq)
        cons_seq = np.array(list(cons_seq))
        changed = np.array([" "] * len(cons_seq_orig))

        for k, v in overwrite.items():
            if isinstance(k, str) and features is not None:
                k = features[k]
            if (
                isinstance(k, collections.abc.Sequence)
                and len(k) == 2
                and isinstance(v, str)
            ):
                l = k[1] - k[0]
                if len(v) != k[1] - k[0]:
                    print(
                        f"Warning: (len('{v}') = {len(v)}) != ({k[1]} - {k[0]} = {l}), this will mess things up. Pad with '-' characters if you want to insert a gap."
                    )
                cons_seq[k[0]: k[1]] = list(v)
                changed[k[0]: k[1]] = "^"
            else:
                raise ValueError(
                    f"Unknown how to process overwrite {k} = '{v}'")
        cons_seq = "".join(cons_seq)
        changed = "".join(changed)
        if cons_seq != cons_seq_orig:
            if verbose:
                if rich:
                    from IPython import display
                    from .viz.syntax import aa_highlighter
                    print("Overwrote to:")
                    display(aa_highlighter.display(f"{cons_seq}\n{changed}"))
                else:
                    print(f"Overwrote to:\n{cons_seq}\n{changed}")

            if isinstance(cons, SeqRecord):
                cons.seq = cons_seq
            else:
                cons = cons_seq
    return cons

def repair_cdr1(seq, features, CDR1, FR1='', **kwargs):
    from .utils import replace_at_index

    if 'GSI' in FR1 and CDR1[0:3] == '---':
        CDR1 = 'GSI' + CDR1[3:]
    elif 'R' in FR1 and CDR1[0:2] == '---':
        CDR1 = 'GRT' + CDR1[3:]
    elif 'R' in FR1 and CDR1[0:2] == '--':
        CDR1 = 'GR' + CDR1[2:]
    elif 'GL' in FR1 and CDR1[0:2] == '--':
        CDR1 = 'GLT' + CDR1[3:]

    if CDR1[0] == '-':
        CDR1 = 'G' + CDR1[1:]

    if CDR1[1] == '-':
        CDR1 = CDR1[0] + 'S' + CDR1[2:]

    if CDR1[2] == '-':
        CDR1 = CDR1[0:2] + 'T' + CDR1[3:]

    return replace_at_index(seq, features['CDR1'], CDR1)

class Cart:
    _queue = None
    _rVHHs = None
    _reports = None

    defaults = dict(
        global_query="kind == '+' & io == 'i'",
        overwrite_rules={
            'repair_cdr1':repair_cdr1
        },
        # overhang5p="AGTGGTGGCGGAGGCTTGGTGCAGCCTGGAGGTTCTCTGAGACTCTCCTGTGCAGCCTCT",
        # overhang5p="TCGGGCGGAGGCTTGGTGCAGCCTGGAGGTTCTCTGAGACTCTCCTGTGCAGCCTCT",  # 2023-02-09
        overhang5p="TCTGGCGGAGGCCTGGTGCAGCCTGGAGGTTCTCTGAGACTCTCCTGTGCAGCCTCT",  # 2023-02-10
        clamp5p=20,
        # overhang3p="CGCGGCCAGGGCACCCAGGTCACCGTCTCCTCCGGCGGTGGTGGCAGC",
        # overhang3p="TGGGGCCAGGGCACCCAGGTCACCGTCTCCTCCGGCGGAGGTGGCTCTGGTGGA",     # 2023-02-09
        overhang3p="TGGGGCCAGGGCACCCAGGTCACCGTCTCCTCCGGCGGAGGT",     # 2023-02-10
        clamp3p=27,
    )

    def save(self, path):
        import pickle
        from .utils import mkdirp_file
        mkdirp_file(path)
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    def __getstate__(self):
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['ex']
        return state

    def load(self, path):
        import pickle
        with open(path, 'rb') as f:
            new_self = pickle.load(f)
        new_self.ex = self.ex
        new_self.ex._cart = new_self
        return new_self

    def reset(self):
        self.ex._cart = None

    def __init__(self, ex=None, **kwargs) -> None:
        self.rVHHs = {}
        self._queue = {}
        self._reports = {}
        self.ex = ex
        self.defaults = {**type(self).defaults, **kwargs}

    def __repr__(self):
        return f"{self.__class__.__name__} with {len(self._queue)} features queued and {len(self.rVHHs)} rVHHs resynthesized"

    def set_defaults(self, **kwargs):
        self.defaults = {**self.defaults, **kwargs}

    def show_plan(self, rVHH=None):
        from .utils import complement, rc
        if rVHH is None:
            rVHH = self.defaults
        return show_plan(**rVHH)

    def add_from_dataframe(
            self, df,
            identifier_col='CDR3ID', 
            antigen_col='antigen', 
            description_col='description', 
            overwrite_consensus_col='overwrite_consensus',
            overwrite_col='overwrite',
            **kwargs):
        for idx, row in df.iterrows():
            if identifier_col is None:
                _id = idx
            else:
                _id = row[identifier_col]

            antigen = row[antigen_col]
            description = row[description_col]
            overwrite_consensus = row[overwrite_consensus_col] if overwrite_consensus_col in row else None
            if overwrite_consensus == '':
                overwrite_consensus = None

            overwrite = row[overwrite_col] if overwrite_col in row else None
            if overwrite is not None:
                import abc
                if isinstance(overwrite, str):
                    import json
                    overwrite = json.loads(overwrite)
                    try:
                        if isinstance(overwrite, dict):
                            overwrite = {tuple(int(l) for l in k.split(',')):v for (k, v) in overwrite}
                        elif isinstance(overwrite, list):
                            overwrite = {tuple(k):v for (k, v) in overwrite}
                    except(ValueError):
                        raise ValueError(f"Cannot interpret 'overwrite' from column {overwrite_col}; "
                        "should be either a list of pairs, e.g. [[[start,end], 'replacement'], [[start,end], 'replacement'], ...] "
                        "or a dict, e.g. {'start,end':'replacement', 'start,end':'replacement', ...} ."
                        f"Instead, recieved {overwrite}."
                        )
                elif not isinstance(overwrite, abc.Mapping):
                    raise ValueError("")

            row = row.drop([identifier_col, antigen_col, description_col, overwrite_consensus_col])

            self.add(
                CDR3ID=_id,
                antigen=antigen,
                description=description,
                overwrite_consensus=overwrite_consensus,
                overwrite=overwrite,
                 **kwargs, **row
            )

    def add(self, CDR3ID, antigen, description="", **kwargs):
        # if self.ex is not None:
        #     CDR3ID = self.ex.find_cdr3(CDR3ID, single=True)

        # _kwargs = {**self.defaults, **kwargs}
        # rVHH = resynthesize(CDR3ID=CDR3ID, ex=self.ex, **_kwargs)
        # if rVHH["CDR3ID"] in self.rVHHs:
        #     print(
        #         f"Warning: overwriting {rVHH['CDR3ID']} as {self.rVHHs[rVHH['CDR3ID']]['clonotypeID']} with {rVHH['clonotypeID']}"
        #     )

        # rVHH["antigen"] = antigen
        # rVHH["description"] = description
        # self.rVHHs[rVHH["CDR3ID"]] = {**_kwargs, **rVHH}
        from .viz.utils import hash_to_mn_short
        var = {}
        if self.ex is not None:
            CDR3ID = self.ex.find_cdr3(CDR3ID, single=True)
            for param in ['CDR3', 'reads', 'nsamples']:
                if param in self.ex.fts.cdr3.var:
                    var[param] = self.ex.fts.cdr3.var.loc[CDR3ID, param]

        rVHH = {
            **self.defaults,
            **kwargs,
            **var,
            'CDR3ID': CDR3ID,
            'CDR3_mn': hash_to_mn_short(CDR3ID),
            'antigen': antigen,
            'description': description
        }
        # rVHH = resynthesize(CDR3ID=CDR3ID, ex=self.ex, **_kwargs)
        if rVHH["CDR3ID"] in self._queue:
            print(
                f"Warning: {rVHH['CDR3ID']} already in queue, overwriting with new params"
            )
        self._queue[CDR3ID] = rVHH

    @property
    def queue(self):
        return pd.DataFrame(self._queue.values())

    def show_rVHHs(
        self,
        columns=['CDR3ID', 'CDR3_mn', 'CDR3', 'NA', 'mode_vs_consensus',
                 'cons_vs_ref', 'n_clonotypes', 'n_reads'],
        **kwargs
    ):
        from .viz.utils import rich_table
        df = self.frame
        df = df[[c for c in columns if c in df.columns]]

        return rich_table(df, **kwargs)

    def show_queue(
        self,
        columns=['CDR3ID', 'CDR3_mn', 'CDR3', 'antigen',
                 'reads', 'samples', 'description'],
        sort=True,
        format={},
        **kwargs
    ):
        from .viz.utils import rich_table

        def sort_cluster(df, col):
            import Bio.Cluster
            strs = df[col].values

            x = np.array(strs, dtype=str)
            y = x.view('U1').reshape((x.size, -1)).view('uint16')
            tree = Bio.Cluster.treecluster(y)

            return df.iloc[tree.sort(), :]

        queue = self.queue
        if columns and columns != True:
            queue = queue[[c for c in columns if c in queue.columns]]
        else:
            drop_cols = ['abundance', 'pick']
            drop_cols = [d for d in drop_cols if d in queue]
            queue = queue.drop(drop_cols, axis='columns')

        if sort and ('CDR3' in queue.columns):
            queue = sort_cluster(
                queue, 'CDR3')

        return rich_table(
            queue,
            format={'reads': '{:.0f}', **format}, **kwargs)

    def check_queue_for_duplicates(self, identifier='CDR3ID'):
        from IPython.display import display
        dff = self.queue
        dups = (dff[identifier].value_counts() > 1)
        dups = dups[dups]
        n_dups = dups.sum()
        if n_dups > 0:
            print(f"Warning: {n_dups:.0f} duplicated sequences!")
            display(dff[dff[identifier].isin(
                dups.index)].sort_values(identifier))
        else:
            print("No duplicated sequences")

    def visualize_queue(self, space='cdr3'):
        from IPython.display import display
        for _id, rVHH in self._queue.items():
            phenotype = rVHH['antigen']
            if phenotype not in self.ex.pheno_names:
                phenotype = None
            display(self.ex.viz.plot_selections_for_feature(
                _id, space=space, phenotype=phenotype))

    def resynthesize(self, items=None, indexes=None, remaining=False, capture=True, verbose=False):
        if self._reports is None:
            self._reports = {}

        if items is not None:
            subqueue = {k: self._queue[k] for k in items}
        elif indexes is not None:
            items = np.array(list(self._queue.keys()))[indexes]
            subqueue = {k: self._queue[k] for k in items}
        elif remaining:
            remaining = set(self._queue.keys()) - set(self.rVHHs.keys())
            _queue = list(self._queue.keys())
            items = []
            items = [k for k in self._queue.keys() if k in remaining]
        else:
            subqueue = self._queue
        n_items = len(subqueue)
        for i, (_id, params) in enumerate(subqueue.items()):
            print(f"Resynthesizing {_id} ({i} / {n_items})")

            # extract relevant arguments to `resynthesize()`
            _kwargs = {k: v for k, v in params.items()
                       if k in _resynthesize_params}
            if verbose:
                print(_kwargs)

            # optionally capture HTML output
            if capture:
                from IPython.utils.capture import capture_output
                from .viz.utils import repr_captured_output
                _context = capture_output()
            else:
                from contextlib import nullcontext
                _context = nullcontext()

            with _context as co:
                rVHH = resynthesize(CDR3ID=_id, ex=self.ex,
                                    verbose=verbose, **_kwargs)

            if capture:
                report = repr_captured_output(co)
                self._reports[_id] = report

            if rVHH["CDR3ID"] in self.rVHHs:
                print(
                    f"Warning: overwriting {rVHH['CDR3ID']} as {self.rVHHs[rVHH['CDR3ID']]['clonotypeID']} with {rVHH['clonotypeID']}"
                )

            self.rVHHs[_id] = {**params, **rVHH}

    def summary(self):
        return pd.DataFrame(self.rVHHs)

    def report(self, id, accordion=False):
        from IPython.display import display, HTML
        from .viz.utils import display_accordion

        row = self.rVHHs[id]
        # title = "{CDR3ID} ({CDR3_mn}): {antigen}".format(**row)
        if accordion:
            title = self.ex.viz.tag(id, short=True, html=False)
            qualifier = ''
            if 'aaSVID' in row:
                qualifier = f" = <code>{row['aaSVID']}</code>"
            return display_accordion(self._reports[id], title=HTML(f"<span>{title}{qualifier}, α-{row['antigen']}</span>"))
        else:
            return self._reports[id]

    def report_all(self, accordion=True):
        from IPython.display import display
        for _id in self._queue:
            if _id in self._reports:
                display(self.report(_id, accordion=accordion))

    @property
    def frame(self):
        rows = []
        index = []
        for k in self._queue:
            if k in self.rVHHs:
                rows.append(self.rVHHs[k])
                index.append(k)
        return pd.DataFrame(data=rows, index=index)

    def __len__(self):
        return len(self.rVHHs)

    # def __repr__(self):

    def preview_plate(self, vendor="twist"):
        pass

    def get_wells(self, vendor="twist"):
        if vendor == "twist":
            pass

    def sequence(self, id, aa=False):
        return self.rVHHs[id]['rec_aa' if aa else 'rec_na']

    @property
    def sequences(self):
        return {id: self.sequence(id) for id, rVHH in self.rVHHs.items()}

    def simulate_cloning(self, method, *args, **kwargs):
        out = {}
        for _id, seq in self.sequences.items():
            out[_id] = method(seq, *args, **kwargs)
        return out

    def save_sequences(self, path, sequences=None, format='gb'):

        from .cloning import write_seqs
        if sequences is None:
            sequences = self.sequences

        write_seqs(sequences, path, format=format)

    def order(self, master=None, order=None, vendor="twist", show=True):
        df = self.frame.drop(columns=['rec_aa', 'rec_na'])
        if len(df) < len(self._queue):
            print("Warning: not all queued sequences have been resynthesized! Order may be incomplete.")

        print(f"Ordering {len(df)} rVHHs...")
        if master is not None:
            print(f" - Saving master file to {master}")
            from .utils import reorder_columns
            reorder_columns(df, ["clonotypeID", "CDR3ID", "NA", "AA", "antigen", "description"]).to_csv(
                master
            )

        if order is not None:
            print(f" - Saving order file for vendor '{vendor}' to {order}")
            if vendor == "twist":
                df[["aaSVID", "NA"]].to_csv(
                    order, index=False, header=False)

        # if show:
        #     display(rich_table(df[['clonotypeID','CDR3ID','NA','AA','antigen','description']], monospace=['NA']))

        return df
