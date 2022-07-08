
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from textwrap import indent
from scipy.spatial.distance import hamming


def remove_barcode_collisions(barcodes, reference,
        query_seq_col='barcode', query_name_col='Name', reference_seq_col = 'Sequence', reference_name_col = 'Sequence Name'):
    """check that candidate query `barcodes` don't overlap with a set of existing `reference` barcodes

    Prints alignment between any query barcodes that overlap with reference barcodes
    """
    problems = []

    for i, row in barcodes.iterrows():
        bc = row[query_seq_col]
        name = row[query_name_col]
        collisions = reference[reference_seq_col].str.contains(bc)

        if sum(collisions) > 0:
            reference_seqs = reference.loc[collisions, :]
            print(f"{name}\t{bc}")
            for j, collision in reference_seqs.iterrows():
                print(f"\t{str(collision[reference_name_col])}")
                for a in pairwise2.align.globalxx(collision[reference_seq_col], bc):
                    print(indent(format_alignment(*a),"\t\t"))

            print()
            problems.append(name)

    return barcodes.loc[~barcodes[query_name_col].isin(problems),:]

def print_sep(char = '*'):
    print('\n'+char*80+'\n')


def enforce_minimum_edit_distance(barcodes, seq_col = 'barcode', min_edit_distance = 2, show_summary=True):
    """check the hamming distance between all pairs of `barcodes`; remove barcodes that differ from other barcodes by < `min_edit_distance` nt

    Parameters
    ----------
    barcodes : pd.DataFrame
    seq_col : str
        Name of column that contains the barcode sequence
    min_edit_distance: int
    show_summary : bool
        True to show a table and histogram of edit distances

    Returns
    -------
    pd.DataFrame
        Barcodes that satisfy the minimum edit distance
    """
    distance_matrix = barcodes[seq_col].apply(lambda bc1: barcodes[seq_col].apply(lambda bc2: len(bc2)*hamming(list(bc1), list(bc2))))

    # omit the main diagonal which will be 0
    tril_ind = np.triu_indices_from(distance_matrix.values, k = 1)
    edit_distances = distance_matrix.values[tril_ind]

    if show_summary:

        display(
            pd.DataFrame(np.array(
                np.unique(edit_distances, return_counts=True),dtype=int), index=['Distance','Count']).astype(int)
        )

        plt.axvline(x=min_edit_distance, color='r', linewidth=1)
        plt.hist(edit_distances)
        plt.show()

    return barcodes.loc[distance_matrix.replace(0, np.nan).min() >= min_edit_distance, :]


def order_primers_IDT(primer_set, wells, seq_col='Order', name_col='Name', rng=None):
    if rng is None:
        rng = slice(0, len(wells))
    df = primer_set.iloc[rng, :].copy()
    df['Well'] = wells
    order = df.loc[:,['Well',name_col,seq_col]].rename(columns={seq_col:'Sequence', name_col: 'Name'})
    return (df, order)
