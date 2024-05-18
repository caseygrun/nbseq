from collections.abc import Mapping
from collections import deque
import os
import sys
import pathlib
import itertools
import hashlib
from collections import abc

import numpy as np

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import Bio.SeqIO
    import Bio.Data.CodonTable

    Bio.Data.CodonTable.register_ncbi_table(
            name="Amber suppressing",
            alt_name=None,
            id=999,
            table={
                    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
                    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                    "TAT": "Y", "TAC": "Y",             "TAG": "Q",   # noqa: E241
                    "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
                    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
                    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
            },
            stop_codons=["TAA", "TGA"],
            start_codons=["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
    )


except ImportError:
    print("Warning: BioPython is not installed. Some nbseq functions may not work.")


def md5(s):
    return hashlib.md5(s.encode('utf-8')).hexdigest()


def md5_ungapped(s):
    return md5(ungapped(s))


def complement(seq):
    return str(Bio.Seq.Seq(seq).complement())


def rc(seq):
    return str(Bio.Seq.Seq(seq).reverse_complement())


_rc = rc

def itersplit(data, predicate):
    yes, no = [], []
    for d in data:
        (yes if predicate(d) else no).append(d)
    return [yes, no]


def read_delim_auto(path, **kwargs):
    import pandas as pd

    fn, ext = os.path.splitext(path)
    ext = ext.lower()

    # allow .csv.gz, .tsv.gz to be treated as .csv, .tsv respectively
    # (pandas will detect compression from filename and decompress 
    # automatically)
    if ext == '.gz':
        fn, ext = os.path.splitext(fn)
    if ext == '.csv':
        return pd.read_csv(str(path), **{'sep': ',', **kwargs})
    elif ext == '.tsv' or ext == 'txt':
        return pd.read_csv(str(path), **{'sep': "\t", **kwargs})
    else:
        raise ValueError(
            f"Unable to guess delimiter from file with extension '{ext}'")


def feature_table_to_unique_sequences(feature_table):
    seqs = feature_table.index.unique().values
    print("Unique sequences:")
    print(seqs)
    seq_records = (SeqRecord(Seq(seq), id=md5(seq)) for seq in seqs)
    return seq_records



def sample_metadata_to_expt_metadata(obs, expt_col='expt',selection_col='name')
    """adapt a dataframe of metadata per-sample to one per-experiment (expt), showing the number of selections, number of rounds, and which phage libraries were involved"""

    group = obs.groupby(expt_col)

    df = pd.DataFrame(index=obs[expt_col].unique())
    df['selections'] = group[selection_col].nunique()

    # get number of rounds for selection
    rs = obs.groupby([expt_col,selection_col])['r'].nunique().reset_index()
    df['fewest_rounds_per_selection'] = rs.groupby(expt_col).min()['r']
    df['most_rounds_per_selection'] = rs.groupby(expt_col).max()['r']
    df['rounds'] = group['round'].apply(lambda x: ','.join(sorted(set(x))))
    df['phage_libraries'] = group['phage_library'].apply(lambda x: ','.join(sorted(set(x))))
    return df

def sample_metadata_to_selection_metadata(sample_metadata, selection_col='name'):
    """adapt a dataframe of metadata per-sample to one per-selection. 

    Each selection may have multiple samples from different rounds, replicates, input/output, etc., each of which is represented by a different row in the sample metadata. In the selection metadata, these are grouped together.

    Parameters
    ----------
    sample_metadata : pd.DataFrame
    selection_col : str
            Which column to use to group selections

    Returns
    -------
    pd.DataFrame
    """
    group = sample_metadata.reset_index().groupby(selection_col)
    df = group.first().drop(['round', 'r', 'sample', 'ID'], axis='columns')
    df['samples'] = group['ID'].count()
    df['io'] = group['io'].apply(lambda x: ''.join(sorted(set(x))))
    df['rs'] = group['r'].count()
    df['rounds'] = group['round'].apply(lambda x: ','.join(sorted(set(x))))
    df['replicates'] = group['replicate'].count()

    col_order = ['expt', 'kind', 'io', 'rs', 'rounds',
                 'replicates', 'samples', 'description']
    col_order = col_order + \
        [col for col in list(df.columns) if col not in col_order]
    df = df[col_order]
    return df


def summarize_selection_phenotypes(selection_metadata, phenotypes):
    import pandas as pd

    def summarize_phenotypes(row):
        # print(list(row.items()))
        # print(row)
        return ' '.join(f'{k}+' for k, v in row.items() if v > 0)
    
    good_phenotypes, bad_phenotypes = itersplit(phenotypes, lambda p: p in selection_metadata.columns)

    if len(bad_phenotypes) > 0:
        print(f"Warning: the following phenotypes are not found in the sleection metadata: {bad_phenotypes}.")
    
    return selection_metadata[good_phenotypes].apply(summarize_phenotypes, axis=1)

def get_rounds(obs, round_col='round'):
    """gets a sorted list of rounds (e.g. ['R2i', 'R3i', 'R4i']) from a dataframe of sample metadata"""
    return sorted(obs[round_col].unique())



def seqrecords_to_dataframe(records, transform=None, include_annotations=True):
    """transform iterable of `Bio.SeqRecord.SeqRecord`s into `pd.DataFrame`"""
    import pandas as pd

    rows = []
    for record in records:
        row = dict(
            id=record.id,
            sequence=str(record.seq),
            name=record.name,
            description=record.description
        )

        if include_annotations:
            row = dict(
                **row,
                **record.annotations,
                **record.letter_annotations
            )
        if transform is not None:
            (_id, seq, name, desc) = transform(
                row['id'], row['sequence'], row['name'], row['description'])
            row = {**row, **dict(id=_id, sequence=seq,
                                 name=name, description=desc)}

        rows.append(row)
    return pd.DataFrame(rows)


def dataframe_to_seqrecords(df, seq_col=None, id_col=None, name_col=None, description_col=None, ungap=False, **kwargs):

    if name_col is None:
        names = df.index
    elif callable(name_col):
        names = name_col(df)
    else:
        names = df.loc[:, name_col]

    if id_col is None:
        ids = names
    elif callable(id_col):
        ids = id_col(df)
    else:
        ids = df.loc[:, id_col]

    if seq_col is None:
        seqs = df.iloc[:, 0]
    elif callable(seq_col):
        seqs = seq_col(df)
    else:
        seqs = df.loc[:, seq_col]
    if ungap:
        seqs = seqs.apply(ungapped)

    if description_col is None:
        descs = [''] * len(seqs)
    elif callable(description_col):
        descs = description_col(df)
    else:
        descs = df.loc[:, description_col]

    seq_records = (SeqRecord(Seq(seq), name=str(name), id=str(_id), description=str(desc), **kwargs)
                   for seq, name, _id, desc in zip(seqs, names, ids, descs))
    return seq_records


def dataframe_to_fasta(df, file, **kwargs):
    return Bio.SeqIO.write(dataframe_to_seqrecords(df, **kwargs), file, format='fasta')


def dataframe_to_msa(df, **kwargs):
    from Bio.Align import MultipleSeqAlignment
    return MultipleSeqAlignment(dataframe_to_seqrecords(df, **kwargs))


def dataframe_to_fastas(df, file_col=None, **kwargs):
    if file_col is None:
        files = df.index
    elif callable(file_col):
        files = file_col(df)
    else:
        files = df.loc[:, file_col]
    records = dataframe_to_seqrecords(df, **kwargs)

    for file, record in zip(files, records):
        Bio.SeqIO.write([record], file, format='fasta')


def series_to_seqrecords(sr):
    return dataframe_to_seqrecords(sr.to_frame())


def series_to_fasta(sr):
    return dataframe_to_fasta(sr.to_frame())


def fortify_to_seqrecords(data, **kwargs):
    """accepts a collection and tries to convert it to an iterable of ``SeqRecord``s

    Parameters
    ----------
    data : iterable of SeqRecord, iterable of str, pd.DataFrame, pd.Series
            data in any of the above formats:
            - pd.DataFrame will be passed to ``dataframe_to_seqrecords``
            - pd.Series will be passed to ``series_to_seqrecords``
            - iterable of str, each str is assumed to be a sequence, indexes 0, 1, ...
              are given as the ``name``
    **kwargs
            passed to ``dataframe_to_seqrecords`` or ``series_to_seqrecords``

    Returns
    -------
    iterable of SeqRecord
    """
    import itertools
    import pandas as pd
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if isinstance(data, pd.DataFrame):
        return dataframe_to_seqrecords(data, **kwargs)
    elif isinstance(data, pd.Series):
        return series_to_seqrecords(data, **kwargs)
    elif isinstance(data, abc.Mapping):
        return (SeqRecord(id=name, seq=Seq(seq)) for name, seq in data.items())
    elif isinstance(data, abc.Iterable):
        peek = next(r for r in data)
        # chain = itertools.chain([peek], data)
        if isinstance(peek, SeqRecord):
            return data
        elif isinstance(peek, str):
            return (SeqRecord(id=str(i), seq=Seq(seq)) for i, seq in enumerate(data))


def chromatograms_to_seqrecords(paths, rc=False, names=True, ids=None, descriptions=None, check_quality=True):
    """load chromatograms in scf format to a pd.DataFrame

    Parameters
    ----------
    paths : iterable of str
            paths to chromatograms in .scf format
    rc : bool
            True to reverse-complement chromatogram sequence
    names : bool or iterable or None
            iterable of names for the sequences, or True to take the name as the
            basename sans ext for each of the  ``paths``, or None to add no
            `name` column

    Returns
    -------
    pd.DataFrame
            DataFrame with 2--3 columns:
            - `file` containing the given path,
            - `seq` containing the sequence,
            - `name` containing given or generated according to the ``names`` parameter
    """
    import bioconvert.io.scf
    import pandas as pd
    import numpy as np
    import os

    if names == True:
        names = [os.path.splitext(os.path.basename(x))[0] for x in paths]
    if names is None:
        names = [None for x in paths]

    if ids is None:
        ids = [None] * len(paths)

    if descriptions is None:
        descriptions = [''] * len(paths)

    records = []
    for f, name, _id, description in zip(paths, names, ids, descriptions):
        seq, qualities, comments = bioconvert.io.scf.read_scf(f)
        if check_quality:
            q = np.array(qualities)
            bad_bases = (q < 30).sum()
            if bad_bases > 3:
                print(f"Warning: chromatogram `{f}` has {bad_bases} positions "
                      "with Q < 30. Manually inspect chromatogram to verify "
                      "it is properly trimmed and base calls reflect "
                      "chromatogram traces. Suppress this warning with "
                      "check_quality=False .")

        if rc:
            seq = _rc(seq)
        rec = SeqRecord(seq, name=name, id=_id, description=description, annotations={
                        'comments': comments}, letter_annotations={'phred_quality': qualities})

        records.append(rec)

    return records


def chromatograms_to_dataframe(paths, rc=False, names=True, ids=None, descriptions=None, check_quality=True):
    recs = chromatograms_to_seqrecords(paths)
    return seqrecords_to_dataframe(recs)


def thread_aa_alignment_to_na(aas, nas, gap='-', check_codons=True, verbose=False):
    codon_table = None
    if check_codons == True:
        from Bio.Data.CodonTable import standard_dna_table
        codon_table = standard_dna_table
    elif check_codons != False:
        codon_table = check_codons

    out = []
    for aa, na in zip(aas, nas):
        na_out = []
        na_i = 0
        for r in aa:
            if r == gap:
                na_out.append(gap*3)
            else:
                codon = na[na_i:na_i+3]
                if codon_table is not None:
                    if codon_table.forward_table[codon] != r:
                        raise ValueError("Codon mismatch:\n"
                                         f"Expect {codon} -> {codon_table.forward_table[codon]} != {r}")
                        raise Exception()
                na_out.append(codon)
                na_i += 3
        if na_i < len(na):
            if verbose:
                print("warning: did not consume entire na sequence")
                print("".join(na_out))
                print("..."+na[na_i:])

        out.append("".join(na_out))
    return match_to_series(nas, out)


def make_CDR_features(CDRs):
    """makes a list of `Bio.SeqFeature`s for CDRs/FRs at the given positions

    Parameters
    ----------
    CDRs : dict of str, tuple
            keys are the names of domains (e.g. CDR1, FR2, etc.); values are half-open intervals giving their positions, e.g. (3, 21) = beginning at the fourth base, ending at the 20th base

    Returns
    -------
    list of Bio.SeqFeature.SeqFeature
            features
    """
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    CDR_features = []
    for name, pos in CDRs.items():
        if 'CDR' in name:
            kind = 'variation'
        else:
            kind = 'misc_feature'
        CDR_features.append(SeqFeature(FeatureLocation(
            *pos), type=kind, qualifiers={'note': name}))
    return CDR_features


def print_features_loc(locs, seq=None):
    from nbseq.viz import pre
    sa = np.array([' ', '-'])
    out = ""
    for i, f in enumerate(locs):
        out += f + '\n'
        out += ''.join(sa[locs[:, i]]) + '\n'
    if seq is not None:
        out += seq
    return out


def copy_feature(feature, **kwargs):
    """duplicate a Bio.SeqFeature.SeqFeature
    """
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    return SeqFeature(**{
        **dict(
            location=feature.location,
            type=feature.type,
            id=feature.id,
            qualifiers=dict(feature.qualifiers),
            ref=feature.ref,
            ref_db=feature.ref_db),
        **kwargs}
    )


def merge_seqrecords(rec1, rec2, **kwargs):
    """Merge two seqrecords

    for ``seq``, ``name``, ``id``, ``description``: missing or empty attributes in rec1 will be filled with those from rec2
    for ``features``, ``dbxrefs``, ``annotations``, ``letter_annotations``: elements from both records will be merged, with rec1 preferred over rec2

    Parameters
    ----------
    rec1 : Bio.SeqRecord.SeqRecord
            primary record
    rec2 : Bio.SeqRecord.SeqRecord
            secondary record

    Returns
    -------
    Bio.SeqRecord.SeqRecord
            merged record
    """
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(**{
        **dict(
            seq=rec1.seq,
            name=rec2.name if (
                rec1.name is None or rec1.name == '') else rec1.name,
            description=rec2.description if (
                rec1.description is None or rec1.description == '') else rec1.description,
            id=rec2.id if (rec1.id is None or rec1.id ==
                           '') else rec1.id,
            features=(rec1.features + rec2.features),
            dbxrefs=(rec1.dbxrefs + rec2.dbxrefs),
            annotations={**rec2.annotations, **rec1.annotations},
            letter_annotations={
                **rec2.letter_annotations, **rec1.letter_annotations},
        ),
        **kwargs
    }
    )
    return rec

def index_features_by_name(rec):
    features = {}
    for feature in rec.features:
        if 'note' in feature.qualifiers:
            name = feature.qualifiers['note']
            if name != '':
                features[name] = feature
    return features

def replace_at_index(x, where, repl):
    x = np.array(list(x))
    x[where[0]:where[1]] = list(repl)
    return ''.join(x)

def overwrite_feature(rec, **kwargs):
    from Bio.Seq import Seq
    features = index_features_by_name(rec)
    new_seq = np.array(list(rec.seq))
    for feature_name, seq in kwargs.items():
        if feature_name in features:
            feature = features[feature_name]
            new_seq[int(feature.location.start):int(feature.location.end)] = list(seq)
    rec.seq = Seq(''.join(new_seq))
    return rec


def copy_record(rec):
    import copy
    new_rec = rec[:]
    new_rec.annotations = copy.copy(rec.annotations)
    new_rec.dbxrefs = copy.copy(rec.dbxrefs)
    return new_rec

def extract_features(rec, as_str=True):
    """extract a dict of feature_name:sequence from a SeqRecord

    Parameters
    ----------
    rec : Bio.SeqRecord.SeqRecord
        seq record

    Returns
    -------
    dict
        keys are feature names (e.g. feature.qualifiers['note']); values are `Bio.Seq.Seq` representing that feature
    """
    features = {}
    for feature in rec.features:
        if 'note' in feature.qualifiers:
            name = feature.qualifiers['note']
            if name != '':
                features[name] = feature.extract(rec.seq)
                if as_str:
                    features[name] = str(features[name])
    return features

def ungap_annotations(seq, features=[], name=None, gap='-', **kwargs):
    """takes a gapped seq, applies the list of features to it, then removes gaps and adjusts the feature positions accordingly

    Biopython SeqRecord objects are useful for visualizing or exporting sequences annotated with features.
    This function can assign features to CDR and FR positions, based on their positions in the reference 
    sequence, by matching those positions to those in an aligned sequence. 

    Parameters
    ----------
    seq : Bio.Seq.Seq
            gapped input sequence, e.g. from 
    features : list of Bio.SeqFeature.SeqFeature, optional
            list of input features, by default []
    name : str, optional
            if given, applies a CDS feature to the entire sequence with this `name` as the `note` annotation; will also set the `name` of the returned `SeqRecord`
    gap : str, optional
            which character indicates gaps, by default '-'
    **kwargs : dict
            additional arguments passed to Bio.SeqRecord.SeqRecord

    Return
    ------
    record : Bio.SeqRecord.SeqRecord
            Sequence
    """
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # takes an aligned seq and a dict of CDR positions, create a
    # `Bio.SeqRecord.SeqRecord` with `SequenceFeature`s for each `CDR` position
    # """
    if features is None:
        features = []
    elif isinstance(features, dict):
        features = make_CDR_features(features)
    if len(seq) == 0:
        raise ValueError("Empty sequence")

    # loc[i, j] = 1 if position i is part of feature j
    locs = np.zeros((len(seq), len(features)), dtype=int)

    for j, feature in enumerate(features):
        locs[int(feature.location.start):int(feature.location.end), j] = 1

    # filter `loc` to include non-gap elements
    seq_array = np.array(list(seq))
    non_gaps = seq_array != gap
    floc = locs[non_gaps, :]
    seq_ungapped = ''.join(seq_array[non_gaps])
    if len(seq_ungapped) == 0:
        raise ValueError("Sequence contains only gap characters")

    _features = []
    for j, feature in enumerate(features):

        # get position of nonzero elements of loc, corresponding to indices where feature is present
        # np.nonzero returns a tuple with one array per axis; since this is a 1d array, get 0th element
        # nz[0] = index of first non-zero element
        # nz[-1] = index of last non-zero element
        nz = np.nonzero(floc[:, j])[0]
        if len(nz) > 0:
            _feature = copy_feature(
                feature,
                location=FeatureLocation(
                    int(nz[0]), 
                    # add 1 because locations are half-open intervals
                    int(nz[-1])+1, 
                    features[j].strand)
            )
            _features.append(_feature)

    if name is not None:
        _features.append(SeqFeature(FeatureLocation(
            0, len(seq_ungapped)), type='CDS', qualifiers={'note': name}))
    else:
        name = ''

    return SeqRecord(seq_ungapped, features=_features, name=name, **kwargs)


def ungap_seqrecord(rec, **kwargs):
    return ungap_annotations(rec.seq,
                             features=rec.features,
                             name=rec.name,
                             description=rec.description,
                             id=rec.id,
                             annotations=rec.annotations,
                             dbxrefs=rec.dbxrefs,
                             letter_annotations=rec.letter_annotations,
                             **kwargs)


def aa_features_to_na(features, offset=0):
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    _features = [
        copy_feature(
            feature,
            location=FeatureLocation(
                feature.location.start * 3 + offset,
                feature.location.end * 3 + offset,
                feature.location.strand
            )
        ) for feature in features
    ]
    return _features


def make_cmd(*args, **kwargs):
    # {cmd[0]} {cmd[1]} ... {cmd[-1]} --foo {kwargs.foo} -b {kwargs.b}
    cmd = " ".join([str(x) for x in args] + [f"--{k.replace('_','-')} '{v}'" if len(
        k) > 1 else f"-{k} '{v}'" for k, v in kwargs.items()])
    return cmd


def _cmd(cmd):
    if not isinstance(cmd, str):
        cmd = " ".join([str(x) for x in cmd])
    return cmd


def run_cmd(cmd, verbose=False, conda=None, **kwargs):
    import subprocess
    import tempfile

    if conda is not None:
        cmd = _cmd(cmd)

        # https://github.com/conda/conda/issues/7980
        cmd = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda} && {cmd}'
    else:
        cmd = [str(x) for x in cmd]
    if verbose:
        print(cmd)

    try:
        out = subprocess.run(cmd,
                             check=True,
                             shell=(True if conda else False),
                             capture_output=verbose,
                             stdout=(None if verbose else subprocess.DEVNULL),
                             stderr=(None if verbose else subprocess.DEVNULL),
                             **kwargs)
    except subprocess.CalledProcessError as err:
        if verbose:
            print(err.stdout.decode('utf8'))
            print(err.stderr.decode('utf8'))
        raise err

    if verbose:
        print(out.stdout.decode('utf8'), file=sys.stdout)
        print(out.stderr.decode('utf8'), file=sys.stderr)


def run_script(cmd, verbose=False, conda=None, strict=True, silent=False, capture=False, **kwargs):
    import subprocess
    import tempfile

    if conda is not None:
        cmd = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda}\n' + cmd

    if strict:
        cmd = "set -euo pipefail\n"+cmd

    with tempfile.NamedTemporaryFile('w') as script:
        if verbose:
            print(script.name)
            print(cmd)

        print(cmd, file=script, flush=True)

        subprocess.run(['sh', script.name], **dict(
            check=True,
            capture_output=capture,
            stdout=(subprocess.DEVNULL if silent else sys.stdout),
            stderr=(subprocess.DEVNULL if silent else sys.stderr),
            **kwargs)
        )

        # if capture:
        # 	print(out.stdout.decode('utf8'), file=sys.stdout)
        # 	print(out.stderr.decode('utf8'), file=sys.stderr)


def display_or_print(*args):
    try:
        from IPython.display import display
        display(*args)
    except (ImportError):
        print(*args)

def mkdirp_file(path):
    p = pathlib.Path(path).resolve().parent
    p.mkdir(parents=True, exist_ok=True)


def mkdirp(path):
    p = pathlib.Path(path)
    p.mkdir(parents=True, exist_ok=True)


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)


def ungapped(seq, gap='-'):
    return seq.replace(gap, '')


def len_ungapped(seq, gap='-'):
    return len(seq) - seq.count(gap)


def index_to_column(df, name=None, inplace=False):
    if name is not None:
        if inplace:
            df.index.rename(name, inplace=True)
        else:
            df.index = df.index.rename(name)
    if inplace:
        df.reset_index(inplace=inplace)
        return df
    else:
        return df.reset_index()


def show_ruler(start, end, step_size=10):
    out = ''
    for num in range(start, end, step_size):
        label = str(num)
        out = out + label + ' ' * (step_size - len(label))
    return out


def show_sequence_ranges(seq, positions, number=True):
    out = ''
    if isinstance(seq, str):
        seq_len = len(seq)
    if isinstance(seq, int):
        seq_len = seq
        seq = ''

    if number:
        step_size = 10
        out = out + show_ruler(0, seq_len, step_size) + "\n"
    if seq != '':
        out = out + seq + "\n"
    
    cur_pos = 0
    track = ''
    tracks = []
    for feature_name, pos in positions.items():
        feature_len = pos[1] - pos[0]

        if feature_len > 0:
            if pos[0] >= len(track):
                track = track + (' ' * (pos[0] - len(track)))
            else:
                tracks.append(track)
                track = ' ' * pos[0]

            if feature_len > len(feature_name):
                track = track + (
                    feature_name + 
                    ('-' * (feature_len - len(feature_name)))
                )
                # track = track + ' '*(seq_len - len(track))
                # out = out + track + "\n"
            else:
                track0 = track + feature_name
                #track = track + ' '*(seq_len - len(track))
                track1 = (
                    (' ' * pos[0])
                    + ('-' * feature_len)
                    # + (' ' * (seq_len - pos[0] - feature_len))
                )
                track = track1
                tracks.append(track0)
                # out = out + track + "\n" + track1 + "\n"
    tracks.append(track)
    out += '\n'.join(track for track in tracks if track != '')
    return out + '\n'


def show_sequence_translated(seq, translation=None, suppress_amber=False,
                             ruler=True, positions={}, label=None, show_lengths=True):
    if not isinstance(seq, Seq):
        seq = Seq(seq)
    if translation is None:
        if suppress_amber:
            translation = seq.translate(to_stop=True, table=999)
        else:
            translation = seq.translate(to_stop=True)

    tra = translation
    if label is None:
        label = ''
    else:
        label = label+' '
    _label = ' '*len(label)

    label_na = f"{label }" + (f"({len(seq):>4} nt):" if show_lengths else '')
    label_aa = f"{_label}" + (f"({len(tra):>4} aa):" if show_lengths else '')

    out = (f"{label_na} {str(seq)}\n"
           f"{label_aa} {'  '.join(a for a in str(tra))}")
    l_padding = len(f"({len(seq):>4} nt): ")

    if ruler:
        if not isinstance(ruler, dict):
            ruler = {}
        out = (" " * l_padding) + show_ruler(**
                                             {'start': 0, 'end': len(seq), 'step_size': 10, **ruler}) + "\n" + out

    for feature_name, pos in positions.items():
        feature_name = str(feature_name)
        if isinstance(pos, tuple):
            feature_len = pos[1] - pos[0]
            track = ' ' * (l_padding + pos[0])
            track = track + (feature_name + '-' *
                             (feature_len - len(feature_name)))
            out = out + "\n" + track
        else:
            out = out + "\n" + (" " * (l_padding + pos)) + "| " + feature_name
    return out


def get_reference(reference=None, reference_path=None, reference_frame_start=0):
    if reference is None:
        if reference_path is None:
            raise Exception(
                "Must specify either `reference` or `reference_path`")
        reference = next(Bio.SeqIO.parse(reference_path, "fasta")).seq
    return reference[reference_frame_start:]


def translate_reference(reference, suppress_amber=True):

    # read and translate reference sequence
    if suppress_amber:
        reference_translated = reference.translate(to_stop=True, table=999)
    else:
        reference_translated = reference.translate(to_stop=True)

    return reference_translated


def match_to_series(input, output):
    import pandas as pd

    if isinstance(input, pd.Series):
        return pd.Series(output, index=input.index, name=input.name)
    else:
        return output


def map_generic(mapping, values):
    """
    Uses a function, dict-like, or pd.Series to map one set of ``values`` to another.
    """
    import pandas as pd
    import numpy as np

    if callable(mapping):
        out = np.vectorize(mapping)(values)
    elif isinstance(mapping, dict):
        out = np.fromiter((mapping[x] for x in values), count=len(values))
    elif isinstance(mapping, pd.Series):
        if not mapping.index.is_unique:
            raise Exception("the mapping provided as "
                            "a pd.Series has a non-unique index. This makes it ambiguous how "
                            f"to map old values to new ones. Must provide a mapping where "
                            "mapping.index.is_unique == True. ")
        out = mapping.loc[values]
    else:
        raise Exception("mapping must be callable, dict, or pd.Series")
    return out


def lookup(df, row_labels, col_labels, default=None):
    """
    Label-based "fancy indexing" function for DataFrame.

    For some reason pd.DataFrame.lookup is deprecated, despite replacements being extremely confusing
    https://github.com/pandas-dev/pandas/issues/39171
    Source: https://github.com/pandas-dev/pandas/blob/v1.3.5/pandas/core/frame.py#L4521-L4581


    Given equal-length arrays of row and column labels, return an
    array of the values corresponding to each (row, col) pair.
    .. deprecated:: 1.2.0
            DataFrame.lookup is deprecated,
            use DataFrame.melt and DataFrame.loc instead.
            For further details see
            :ref:`Looking up values by index/column labels <indexing.lookup>`.
    Parameters
    ----------
    row_labels : sequence
            The row labels to use for lookup.
    col_labels : sequence
            The column labels to use for lookup.
    Returns
    -------
    numpy.ndarray
            The found values.
    """
    from pandas.core.dtypes.common import is_object_dtype
    from pandas._libs import lib

    n = len(row_labels)
    if n != len(col_labels):
        raise ValueError("Row labels must have same size as column labels")
    if not (df.index.is_unique and df.columns.is_unique):
        # GH#33041
        raise ValueError("DataFrame.lookup requires unique index and columns")

    thresh = 1000
    if not df._is_mixed_type or n > thresh:
        values = df.values
        ridx = df.index.get_indexer(row_labels)
        cidx = df.columns.get_indexer(col_labels)
        nulls = (ridx == -1) | (cidx == -1)
        if default is None:
            if (ridx == -1).any():
                raise KeyError("One or more row labels was not found")
            if (cidx == -1).any():
                raise KeyError("One or more column labels was not found")

        ridx[ridx == -1] = 0
        cidx[cidx == -1] = 0

        flat_index = ridx * len(df.columns) + cidx
        result = values.flat[flat_index]
        result[nulls] = default
    else:
        result = np.empty(n, dtype="O")
        for i, (r, c) in enumerate(zip(row_labels, col_labels)):
            try:
                result[i] = df._get_value(r, c)
            except Exception as e:
                if default is not None:
                    result[i] = default
                else:
                    raise e

    if is_object_dtype(result):
        result = lib.maybe_convert_objects(result)

    return result


def replace_multiple(values, replacements):
    """perform multiple string replacements

    Parameters
    ----------
    values : pd.Series
            values to edit
    replacements : dict
            keys are strings to search for, values are replacements
    """
    import re
    def replace(x):
        for old, new in replacements.items():
            if isinstance(old, re.Pattern):
                x = old.sub(new, x)
            else:
                x = x.replace(old, new)
        return x

    return values.apply(replace)


def pull_and_rename(df, columns={}):
    return df[list(columns.keys())].rename(columns=columns)


def reorder_columns(df, order):
    cols = list(df.columns)
    _ordered_cols = set(order)
    new_order = order + [c for c in cols if c not in _ordered_cols]
    return df[new_order]


def clamp_finite(series, digits=None):
    from builtins import round as roundn
    _max = series[series != np.inf].max()
    _min = series[series != -np.inf].min()

    if digits is not None:
        _max = round(_max, digits)
        _min = round(_min, digits)
    return series.replace({
        float('inf'): _max,
        float('-inf'): _min
    })


class dotdict(dict):
    """
    dot.notation access to dictionary attributes

    https://stackoverflow.com/a/23689767/4091874
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class lazydict(Mapping):
    def __init__(self, *args, **kw):
        self._recipes = dict(*args, **kw)
        self._cache = dict()

    def __getattr__(self, key):
        return self.__getitem__(key)

    def __getitem__(self, key):
        if key not in self._cache:
            func = self._recipes.__getitem__(key)
            self._cache[key] = func()
        return self._cache[key]

    def __iter__(self):
        return iter(self._recipes)

    def __len__(self):
        return len(self._recipes)


def intersection_ordered(lst, st):
    """ Return items from ``lst`` that are also in ``st``, preserving the order of ``lst``
    """
    st = set(st)
    return [l for l in lst if l in st]


def sparse_coerce(X, Y, **kwargs):
    """ Convert X to the scipy.sparse class of X1 """
    import scipy.sparse

    if scipy.sparse.issparse(X) and scipy.sparse.issparse(Y):
        if scipy.sparse.isspmatrix_csr(Y):
            return X.tocsr(**kwargs)
        elif scipy.sparse.isspmatrix_csc(Y):
            return X.tocsc(**kwargs)
        elif scipy.sparse.isspmatrix_coo(Y):
            return X.tocoo(**kwargs)
    return X


def sparse_drop_na(X1):
    """ Set nan values in a sparse matrix to implicit zeros
    """
    import scipy.sparse

    X = X1.tocoo()
    notna = ~np.isnan(X.data)
    X_notna = scipy.sparse.coo_matrix(
        (X.data[notna], (X.row[notna], X.col[notna])), shape=X.shape)
    return sparse_coerce(X_notna, X1)


def sparse_var(a, axis=None):
    """ Variance of sparse matrix a
    var = mean(a**2) - mean(a)**2
    """
    a_squared = a.copy()
    a_squared.data **= 2
    return a_squared.mean(axis) - np.square(a.mean(axis))


def sparse_std(a, axis=None):
    """ Standard deviation of sparse matrix a
    std = sqrt(var(a))
    """
    return np.sqrt(sparse_var(a, axis))


def sparse_gstd(a, axis=None):
    """ Geometric standard deviation of sparse matrix a

    for A = a_0, ... a_n, where g = geometric mean of A:
    gstd = exp(sqrt( sum( ln(a_i / g)^2 / n ) ))
    """
    # TODO: check this
    # g = sparse_gmean_nonzero(a, axis=!axis)
    g = sparse_gmean_nonzero(a, axis=axis)
    ag = a / g
    ag.data = np.log(ag.data)**2

    # densify to row or column vector
    ag = ag.sum(axis=axis) / a.getnnz(axis=axis)
    ag = np.exp(np.sqrt(ag))

    # if axis=None, ag will be scalar; otherwise it will be np.matrix; convert to an array
    if axis is not None:
        return np.squeeze(ag.A1)
    return ag


def sparse_mean_nonzero(a, axis=None):
    """ Mean of non-zero elements in sparse matrix a

    Parameters
    ----------
    a : scipy.sparse.spmatrix
            sparse matrix
    axis : int, optional
            Mean over which axis: 0 = mean of each row, 1 = mean of each column, None = mean of entire matrix

    Returns
    -------
    numpy.matrix
            matrix of appropriate dimensions, or scalar if axis=None
    """
    return (a.sum(axis=axis).A1 / a.getnnz(axis=axis))


def sparse_gmean_nonzero(a, axis=None):
    import scipy.sparse

    # geometric_mean(x) = e^(arithmetic_mean(ln(x)))
    # c = scipy.sparse.coo_matrix(a)
    c = a.copy()
    c.data = np.log(c.data)

    # this operation densifies the matrix along ~axis
    if axis is not None:
        c = (c.sum(axis=axis).A1 / c.getnnz(axis=axis))
    else:
        c = (c.sum(axis=axis) / c.getnnz(axis=axis))
    return np.exp(c)
    # return sparse_coerce(c, a)


def sparse_rank(X, asc=True, pct=False):
    """ calculates the row-wise rank or percentile of non-zero elements of sparse matrix X
    """
    import scipy.sparse
    from scipy.stats import rankdata

    X = sparse_drop_na(X)
    ml = X.tolil()
    data = ml.data

    # to sort in descending order, negate values before argsort
    # if asc == False:
    #     _dir = -1
    #     # can't get this to work...
    #     raise NotImplementedError()
    # else: _dir = 1

    # np.argsort ranks items in *ascending* order:
    # >>> np.argsort([5,4,3,2,1])
    # array([4, 3, 2, 1, 0])

    def rank_row(x):

        # return list(np.argsort(x))
        # return list(np.argsort(x).argsort())
        return list(rankdata(x))

    ml.data = np.array(list(map(rank_row, data)), dtype='object')

    if pct:
        # divide by row-wise nnz
        # do this weird manipulation to maintain sparsity in X2
        X = ml.tocsr()
        nnz = scipy.sparse.diags((1/X.getnnz(axis=1)).ravel())
        X2 = nnz.dot(X)

        # for descending percentile, pct_descending = 1 - pct_ascending
        if asc == False:
            X2.data = 1 - X2.data
        return X2
    else:
        if asc == False:
            raise NotImplementedError()
        return sparse_coerce(ml, X)


def sparse_product(X, log=False, axis=0):
    """ calculates the product of elements of sparse matrix X along specified axis
    """
    X = X.copy()
    X.data = np.log10(X.data)
    if log:
        return X.sum(axis=axis).A1
    else:
        return 10**(X.sum(axis=axis)).A1


def sparse_min_nz(X, axis=None):
    if axis is None:
        return X.data.min()
    elif axis == 1:
        mins = [X.getrow(i).data.min() for i in range(X.shape[0])]
    elif axis == 0:
        mins = [X.getcol(i).data.min() for i in range(X.shape[1])]
    return np.array(mins)


def ecdf_x_from_p(p, train=None, ecdf=None):
    """ find the lowest value x for which ECDF(x) > p
    """
    from statsmodels.distributions.empirical_distribution import ECDF
    if ecdf is None:
        ecdf = ECDF(train)
    return min(ecdf.x[ecdf.y > p])


# https://stackoverflow.com/questions/5384570/whats-the-shortest-way-to-count-the-number-of-items-in-a-generator-iterator/34404546#34404546

# Avoid constructing a deque each time, reduces fixed overhead enough
# that this beats the sum solution for all but length 0-1 inputs
consumeall = deque(maxlen=0).extend


def ilen(it):
    # Make a stateful counting iterator
    cnt = itertools.count()
    # zip it with the input iterator, then drain until input exhausted at C level
    # cnt must be second zip arg to avoid advancing too far
    consumeall(zip(it, cnt))
    # Since count 0 based, the next value is the count
    return next(cnt)


def display_all(df):
    import pandas as pd
    from IPython.display import display
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        display(df)