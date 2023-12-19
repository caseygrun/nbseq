import pandas as pd
from .utils import *

CDR_positions = {
    'alpaca': {
        'CDR1': (17, 26),
        'CDR2': (46, 53),
        'CDR3': (94, 115)
    },
    'synthetic': {
        'CDR1': (0, 6),
        'CDR2': (23, 29),
        'CDR3': (70, 87),
    }
}


def extract_CDRs(seqs, CDRs=None, library=None, push_left=False):
    """Extract sequences of CDRs at defined positions from a sequence of aligned nucleic acid or peptide sequences

    Parameters
    ----------
    seqs : list, dict, pd.Series
            aligned DNA or translated and aligned AA sequences
    CDRs : dict, optional
            keys are the names of the CDR(s), values are `(start, end)` position tuples
    library : str, optional
            if CDRs is omitted, this can be given instead and values will be taken from `CDR_positions`
    push_left : bool, Default=False
            if True, remove leading dashes ('-') from each CDR sequence, in order to push the aligned sequence to the left

    Returns
    -------
    pd.DataFrame
            one column for each key in `CDRs`, indexed the same as `seqs`

    """
    def get_subsequence_extractor(start, end):
        length = end - start

        def extractor(seq):
            if pd.isna(seq):
                return seq
            else:
                if push_left:
                    seq = seq[start:end].lstrip('-')
                    return seq + ('-' * (length-len(seq)))
                else:
                    return seq[start:end]

        return extractor

    if CDRs is None:
        CDRs = CDR_positions[library]

    if not isinstance(seqs, pd.Series):
        seqs = pd.Series(seqs)

    CDR_df = pd.DataFrame(columns=CDRs.keys(), index=seqs.index)

    assert all(len(cdr_pos) == 2 for cdr, cdr_pos in CDRs.items(
    )), 'Error: CDRs must be specified as a dict mapping CDR names to CDR (start, end) positions'
    for subseq_name, (start, end) in CDRs.items():
        CDR_df[subseq_name] = seqs.apply(get_subsequence_extractor(start, end))

    return CDR_df

def extract_features_dict(seq, features):
    return {k: seq[v[0]:v[1]] for k, v in features.items()}


def identify_clonotype(df, cdrs=['CDR1', 'CDR2', 'CDR3'], id_col='clonotypeID', seq_col='clonotype', inplace=True):
    col_hashes = df[cdrs].apply(lambda c: c.apply(md5_ungapped))
    if inplace:
        out = df
    else:
        out = pd.DataFrame(index=df.index, columns=[id_col])

    out[id_col] = col_hashes.apply(lambda row: md5(
        ''.join(row.values.astype(str))), axis=1)
    if seq_col is not None:
        out[seq_col] = df[cdrs].apply(
            lambda row: ''.join(row.values.astype(str)), axis=1)
    return df


def identify_seq(df, space='CDR3', id_col=None, seq_col=None):
    if space == 'CDR3':
        if id_col is None:
            id_col = 'CDR3ID'
        if seq_col is None:
            seq_col = 'CDR3'
    elif space == 'AA':
        if id_col is None:
            id_col = 'aaSVID'
        if seq_col is None:
            if 'AA_aligned' in df.columns:
                seq_col = 'AA_aligned'
            elif 'AA' in df.columns:
                print("Warning: using column 'AA' to identify sequences. Ensure this is the clonotype-trimmed trimmed amino acid sequence and not the full-length amino acid sequence")
            else:
                raise ValueError(
                    "Could not figure out which column to use to identify sequences. Pass `seq_col` or a different value for `space`")
    df[id_col] = df[seq_col].apply(md5_ungapped)
    return df

def extract_CDRs_dataframe(df, CDRs, seq_col='aligned', library_col='library'):
    pass


