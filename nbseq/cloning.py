from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
import Bio.SearchIO
import Bio.Blast
import pydna
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.readers import read as read_seq
from pydna.dseqrecord import Dseqrecord

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from pathlib import Path


def load_seq(seq):
    if isinstance(seq, SeqRecord):
        return Dseqrecord(seq)
    elif isinstance(seq, Seq):
        return Dseqrecord(seq)
    elif isinstance(seq, str) or isinstance(seq, Path):
        return read_seq(seq)


def show_seq(dseq, max_len=20):
    if isinstance(dseq, Seq):
        dseq = Dseq(dseq)
    if isinstance(dseq, Dseqrecord):
        dseq = dseq.seq

    wf = dseq.watson
    cf = dseq.crick[::-1]

    overlap_len = len(dseq) - abs(dseq.ovhg)
    if dseq.ovhg > 0:
        wf = ' '*dseq.ovhg + wf
        # al = ' '*dseq.ovhg + '|'*overlap_len
        cf = cf + ' '*dseq.ovhg
    else:
        cf = ' '*(-dseq.ovhg) + cf
        # al = '|'*overlap_len + ' '*(-dseq.ovhg)
        wf = wf + ' '*(-dseq.ovhg)

    length = len(dseq) - abs(dseq.ovhg) - max_len
    if max_len is not None:
        if len(wf) > max_len//2:
            flen = f'({length:d}nt)'
            fdot = '.'*len(flen)
            wf = wf[:max_len//2] + flen + wf[-max_len//2:]
            # al = al[:max_len//2] + fdot + al[-max_len//2:]
            cf = cf[:max_len//2] + fdot + cf[-max_len//2:]

    # return f'{wf}\n{al}\n{cf}'
    return f'{wf}\n{cf}'


def plot_seq(seq, record_class=None, ax=None, multiline=False, **kwargs):
    from dna_features_viewer import BiopythonTranslator
    bpt = BiopythonTranslator()
    if record_class is None:
        record_class = 'linear'
        if isinstance(seq, Dseqrecord):
            if seq.circular:
                record_class = 'circular'
        elif isinstance(seq, SeqRecord):
            if 'toplogy' in seq.annotations:
                if seq.annotations['topology'] == 'circular':
                    record_class = 'circular'
    grec = bpt.translate_record(seq, record_class)

    if multiline:
        if ax is None:
            ax, _ = grec.plot_on_multiple_lines(**kwargs)
        else:
            grec.plot_on_multiple_lines(ax=ax, **kwargs)
    else:
        if ax is None:
            ax, _ = grec.plot(**kwargs)
        else:
            grec.plot(ax=ax, **kwargs)
    return ax, grec


def blast_pairwise(query, subject, scratch_dir=None):
    """Perform local BLAST search of single query sequence against single subject sequence, then read results back into BioPython objects

    Parameters
    ----------
    query : Bio.Seq or Bio.SeqRecord
        Query sequence (e.g. from Sanger sequencing of construct, Plasmidsaurus)
    subject : Bio.Seq or Bio.SeqRecord
        Subject sequence (e.g. )
    scratch_dir : str or pathlib.Path, optional
        directory to save BLAST output files; if omitted, a temporary directory is created

    Returns
    -------
    bl : Bio.SearchIO.QueryResult
        BLAST search results in BioPython format. Contains zero or more `Bio.SearchIO.Hit`s. Each hit contains one or more `Bio.SearchIO.HSP`
    df : pd.DataFrame
        BLAST output in tabular format with the following columns:

        - qseqid means Query Seq-id
        - sseqid means Subject Seq-id
        - pident means Percentage of identical matches
        - length means Alignment length
        - nident means Number of identical matches
        - mismatch means Number of mismatches
        - gapopen means Number of gap openings
        - gaps means Total number of gap
        - qstart means Start of alignment in query
        - qend means End of alignment in query
        - sstart means Start of alignment in subject
        - send means End of alignment in subject
        - evalue means Expect value
        - bitscore means Bit score
        - qcovs means Query Coverage Per Subject (for all HSPs)
        - qcovhsp means Query Coverage Per HSP
        - qcovus is a measure of Query Coverage that counts a position in a subject sequence for this measure only once. The second time the position is aligned to the query is not counted towards this measure.
    """
    import tempfile
    from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline, NcbiblastformatterCommandline
    from Bio import SeqIO, SearchIO
    import pandas as pd

    def assign_id(seq):
        if seq.id == '' or seq.id == 'id' or seq.id == '.':
            seq.id = seq.name
        return seq

    q1 = query
    q2 = subject

    def run(cmd):
        stdout, stderr = cmd()
        if len(stdout) > 0:
            print(stdout)
        if len(stderr) > 0:
            print(stderr)

    with tempfile.TemporaryDirectory() as tmpdir:
        if scratch_dir is None:
            scratch_dir = tmpdir

        q1_fasta = Path(scratch_dir) / 'q1.fasta'
        q2_fasta = Path(scratch_dir) / 'q2.fasta'
        SeqIO.write(assign_id(q1), q1_fasta, 'fasta')
        SeqIO.write(assign_id(q2), q2_fasta, 'fasta')

        # cmd = NcbimakeblastdbCommandline(input_file=q1_fasta, input_type='fasta', dbtype='nucl')
        # stdout, stderr = cmd()
        # print(stdout)
        # print(stderr)

        # asn.1 file contains all data; export to XML and TSV to get some additional datat
        asn1_path = Path(scratch_dir) / 'out.asn'
        cmd = NcbiblastnCommandline(
            subject=q2_fasta, query=q1_fasta, outfmt=11, out=asn1_path)
        run(cmd)

        # xml file contains alignment of each high-scoring pair (HSP) between the query and subject sequence
        xml_path = Path(scratch_dir) / 'out.xml'
        cmd = NcbiblastformatterCommandline(
            out=xml_path, outfmt=5, archive=asn1_path)
        run(cmd)

        # xml format does not give stats on fraction of query sequence that that was covered by HSPs.
        # BioPython cannot tile HSPs to calculate this statistic, so we export to tabular format
        # https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Options_for_the_commandline_a_
        # The supported format specifiers are:
        # qseqid means Query Seq-id
        # qgi means Query GI
        # qacc means Query accesion
        # sseqid means Subject Seq-id
        # sallseqid means All subject Seq-id(s), separated by a ';'
        # sgi means Subject GI
        # sallgi means All subject GIs
        # sacc means Subject accession
        # sallacc means All subject accessions
        # qstart means Start of alignment in query
        # qend means End of alignment in query
        # sstart means Start of alignment in subject
        # send means End of alignment in subject
        # qseq means Aligned part of query sequence
        # sseq means Aligned part of subject sequence
        # evalue means Expect value
        # bitscore means Bit score
        # score means Raw score
        # length means Alignment length
        # pident means Percentage of identical matches
        # nident means Number of identical matches
        # mismatch means Number of mismatches
        # positive means Number of positive-scoring matches
        # gapopen means Number of gap openings
        # gaps means Total number of gap
        # ppos means Percentage of positive-scoring matches
        # frames means Query and subject frames separated by a '/'
        # qframe means Query frame
        # sframe means Subject frame
        # btop means Blast traceback operations (BTOP)
        # staxids means unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
        # sscinames means unique Subject Scientific Name(s), separated by a ';'
        # scomnames means unique Subject Common Name(s), separated by a ';'
        # sblastnames means unique Subject Blast Name(s), separated by a ';' (in alphabetical order)
        # sskingdoms means unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order)
        # stitle means Subject Title
        # salltitles means All Subject Title(s), separated by a '<>'
        # sstrand means Subject Strand
        # qcovs means Query Coverage Per Subject (for all HSPs)
        # qcovhsp means Query Coverage Per HSP
        # qcovus is a measure of Query Coverage that counts a position in a subject sequence for this measure only once. The second time the position is aligned to the query is not counted towards this measure.

        tab_path = Path(scratch_dir) / 'out.tsv'
        formats = 'qseqid sseqid pident length nident mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus'
        cmd = NcbiblastformatterCommandline(
            out=tab_path, outfmt=f'6 {formats}', archive=asn1_path)
        run(cmd)

        bl = list(SearchIO.parse(xml_path, 'blast-xml'))

        # query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

        table = pd.read_csv(
            tab_path, sep="\t", header=None,
            names=formats.split(' ')
            # names=['query', 'subject', 'percent_identity', 'length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bit_score']
        )
    return bl, table


def render_VHH_in_vector(VHH, vector='../../attachments/vectors/pD-VHH-2022b.gb', insert_range=(21, 402)):
    vector = load_seq(vector)
    out = (
        vector[insert_range[1]:insert_range[0]] +
        VHH
    ).looped()
    out.name = f'pD-{VHH.name}'.replace(' ', '-')
    return out


def clone_VHH_pD_to_pJEG3(pD_VHH, insert_name=None, pJEG3='../../attachments/vectors/2021-04-pJEG-3.gb'):
    from Bio.Restriction import NotI, AscI

    pJEG3 = read_seq('../../attachments/vectors/2021-04-pJEG-3.gb')

    # NotI/AscI digest VHH in pD vector; discard backbone
    _, insert = sorted(pD_VHH.cut([NotI, AscI]))

    # NotI/AscI digest pJEG3 vector; discard insert
    pJEG3_cut, _ = sorted(pJEG3.cut([NotI, AscI]))

    # directional clone and ligate
    pJEG3_rVHH = (pJEG3_cut + insert).looped()

    # update to petter name
    if insert_name is None:
        insert_name = insert.name

    pJEG3_rVHH.name = f'pJEG3-{insert_name}'.replace(' ', '-')
    return pJEG3_rVHH


pF4 = SeqRecord(name='pD-pCER-F4',
                seq=Seq('ctgtttagaggcgttcagtCTCAGTTGCAGCTCGTGGAGTCG'))
pR5 = SeqRecord(name='pD-pCER-R5',
                seq=Seq('tccaccagagccacctccgcCaGAGGAGACGGTGACCTGGGTCC'))


def clone_VHH_pD_to_pCER243(
        pD_VHH, insert_name=None,
        pCER243='../../attachments/vectors/2021-11-10_pCER243-empty.gb',
        primer_F=pF4, primer_R=pR5):

    from pydna.amplify import pcr
    from pydna.assembly import Assembly
    from Bio.Restriction import XhoI

    pCER243 = read_seq(pCER243)

    # PCR amplify VHH in pD vector
    amplicon = pcr([primer_F, primer_R, pD_VHH])

    # XhoI digest to linearize pCER243 backbone
    pCER243_XhoI = pCER243.linearize(XhoI)

    # NEBuilder assemble and ligate
    assm = Assembly([pCER243_XhoI, amplicon], limit=20)
    pCER243_VHH = assm.assemble_circular()[0]
    pCER243_VHH.annotations["topology"] = "circular"

    # update to petter name
    if insert_name is None:
        insert_name = amplicon.name

    pCER243_VHH.name = f'pCER243-{insert_name}'.replace(' ', '-')
    return pCER243_VHH


def clone_VHH_gene_to_pCER243(
        rVHH,
        primer_F=SeqRecord(
            name='rV-pCER-F4',
            seq=Seq('CTGTTTAGAGGCGTTCAGTCTCAGGTGCAGCTGGTGGAGTCTGGCGGAGGCCTGGTGC')),
        primer_R=SeqRecord(
            name='rV-pCER-R5',
            seq=Seq('TCCACCAGAGCCACCTCCGCCGGAGGAG')),
        pCER243='../../attachments/vectors/2021-11-10_pCER243-empty.gb', insert_name=None):
    from pydna.amplify import pcr
    from pydna.assembly import Assembly
    from Bio.Restriction import XhoI

    pCER243 = read_seq(pCER243)

    # PCR amplify VHH in pD vector
    amplicon = pcr([primer_F, primer_R, rVHH])

    # XhoI digest to linearize pCER243 backbone
    pCER243_XhoI = pCER243.linearize(XhoI)

    # NEBuilder assemble and ligate
    assm = Assembly([pCER243_XhoI, amplicon], limit=20)
    pCER243_VHH = assm.assemble_circular()[0]

    # update to petter name
    if insert_name is None:
        insert_name = rVHH.name

    pCER243_VHH.name = f'pCER243-{insert_name}'.replace(' ', '-')
    return pCER243_VHH


def write_seqs(sequences, path, format='gb'):
    from Bio import SeqIO
    from pathlib import Path
    from .utils import mkdirp
    from collections.abc import Mapping

    if isinstance(sequences, list):
        sequences = {seq.id:seq for seq in sequences}

    if format=='gb' or format=='gbk':
        mkdirp(path)
        path = Path(path)
        for _id, seq in sequences.items():
            with open(path / f"{_id}.gb", "w") as f:
                SeqIO.write(seq, f, "gb")
    else:    
        with open(path, "w") as f:
            SeqIO.write(sequences.values(), f, format)

def codon_table_to_dataframe(table):
    import pandas as pd
    rows = []
    for aa, codons in table.items():
        for codon, freq in codons.items():
            rows.append([aa, codon, freq])
    return pd.DataFrame(rows, columns=['AA','codon','frequency'])

def compare_codon_usage(species_list = ['h_sapiens', 'e_coli']):
    import pandas as pd
    import python_codon_tables
    dfs = []
    for species in species_list:
        df = codon_table_to_dataframe(python_codon_tables.get_codons_table(species))
        df = df.rename(columns={'frequency':species})
        dfs.append(df)

    while len(dfs) > 1:
        dfs[0] = pd.merge(dfs[0], dfs.pop(), on=['AA','codon'])

    return dfs[0]