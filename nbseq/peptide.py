import os
import Bio.SeqIO
import Bio.Align
from Bio.Align import substitution_matrices

def _setup_alignment(args):
	return _initialize_alignment(**args)

_verbose = False
_min_length = 0

def _initialize_alignment(library,
	reference_path,
	reference_length=101,
	reference_frame_start=0,
	min_length=90,
	verbose = False):
	"""
	perform a peptide alignment between each translated read and the reference peptide sequence, in order to identify CDRs
	"""

	print(f"- Setting up alignment worker {os.getpid()}: {locals()}")

	global aligner
	global reference_translated_seq
	global _verbose
	global _min_length
	_verbose = verbose
	_min_length = min_length

	reference = next(Bio.SeqIO.parse(reference_path,"fasta"))
	reference_translated = reference[reference_frame_start:(reference_frame_start+reference_length*3)].translate(to_stop=True)
	reference_translated_seq = str(reference_translated.seq)

	aligner = Bio.Align.PairwiseAligner()
	aligner.alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ-'

	# ALPACA
	if library == 'alpaca':
		# goal: query sequence aligned to target (reference),
		# + with gaps inserted within CDRs;
		# + aligned query is NOT longer than reference
		# + no gaps are added to left side of target, because want to determine

		aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')

		# don't penalize matches to X, otherwise CDRs align badly
		aligner.substitution_matrix[:,'X'] = 0.5

		# the internal regions of the target (reference) sequence should be
		# generously long enough to accomodate variable length query sequences.
		aligner.target_internal_gap_score = -10000

		# in particular, clamp the alignment to the left side of the reference;
		# otherwise, determining CDRs by position won't work
		aligner.target_left_gap_score = -10000

		# # some reads may be longer than the target sequence; that's OK
		# aligner.target_right_gap_score = 0
		# likewise with the right side of the reference, otherwise it will just
		# push all of CDR3 off the right side
		aligner.target_right_gap_score = -1000

		# still need to penalize gaps somewhat, else
		aligner.query_open_gap_score = -1
		aligner.query_extend_gap_score = -0.5
		aligner.mismatch_score = -1

	# SYNTHETIC
	elif library == 'synthetic':

		# these seem to work better with default settings
		aligner.target_open_gap_score = -10000
		aligner.query_left_gap_score = -1000
		aligner.query_right_gap_score = -10



def _align_reference(seq):
	# leading/trailing '-'s might be introduced in nucleic acid alignment
	# don't waste time aligning anything too short at that point, it's a
	# frameshift/premature stop codon
	seq = seq.strip('-')
	if len(seq) < _min_length: return ''
	alignments = aligner.align(reference_translated_seq,seq)

	# lots of alignments suggests poor alignment quality
	if _verbose: print(len(alignments))
	if len(alignments) < 1000:
		alignments = sorted(alignments)

	alignment = str(alignments[0])
	s = alignment.split("\n")[2]
	if _verbose: print(alignment)
	return s


def pairwise_align_reference_parallel(seqs, processes=None, **kwargs):
	import multiprocessing

	# # table = alignment_table[alignment_table[f'{library}_aligned'] & (alignment_table[f'{library}_length'] > peptide_length_threshold)]
	#
	# return table[f'{library}_translation'].apply(align_reference).rename(f'{library}_translation_align')

	print(f"Aligning {len(seqs)} peptide sequences to reference with {processes} processes...")
	print("Alignment parameters: ")
	print(kwargs)

	with multiprocessing.Pool(processes, initializer = _setup_alignment, initargs = [kwargs]) as pool:
		# TODO: consider using imap_unordered and then sorting at the end?
		aligned = pool.map(_align_reference, seqs)
	return dict(zip(seqs, aligned))

def pairwise_align_reference(seqs, **kwargs):
	_initialize_alignment(**kwargs)
	aligned = map(_align_reference, seqs)
	return dict(zip(seqs, aligned))
