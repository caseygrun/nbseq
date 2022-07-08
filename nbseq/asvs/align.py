import hashlib
import os.path
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO

import Bio.Align

from ..utils import *


def index_reference_bowtie(reference_path, output_path, index_name=None):
	if index_name is None:
		index_name = os.path.splitext(os.path.basename(reference_path))[0]

	# cp {reference_path} {output_path}/{index_name}.fasta
	if Path(reference_path) != (Path(output_path) / f'{index_name}.fasta'):
		shutil.copy(reference_path, Path(output_path) / f'{index_name}.fasta')

	# cd {output_path}
	# bowtie2-build {index_name}.fasta "{index_name}"
	subprocess.run(['bowtie2-build',
		f'{index_name}.fasta', index_name],
		cwd=output_path)

	#
	# """
	# mkdir -p {output.dir}
	# cp {input} {output.fasta}
	# cd {output.dir}
	# bowtie2-build $(basename {output.fasta}) "{wildcards.library}"
	# """

def align_local_bowtie(input_fasta_fwd, input_fasta_rev, input_index, library,
	output_bam_coord_sorted, output_bam_name_sorted,
	retry=True, retry_params="-D 20 -R 3 -N 0 -L 18 -i S,1,0.50",
	n_penalty=64, n_ceil=0, max_fragment_length=350, threads=1, conda=None, cleanup=False, log=None):

	from textwrap import dedent
	import sys

	if conda is not None:
		conda_activation = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda}'
	else:
		conda_activation = ''

	cmd = dedent(f"""\
	set -euox pipefail
	# activate conda env if given
	{conda_activation}

	output_name_dir=$(dirname $(readlink -m {output_bam_name_sorted}))
	mkdir -p $output_name_dir
	output_file_name_sorted=$(basename {output_bam_name_sorted})

	# intermediate .sam and -unsorted.bam files will be produced before the
	# final .bam file. Put them in the same directory as the desired output
	# file.
	output_coord_dir=$(dirname $(readlink -m {output_bam_coord_sorted}))
	mkdir -p $output_coord_dir

	# get absolute path to the input FASTA files before we start changing
	# the directory
	seqs_fwd=$(readlink -f {input_fasta_fwd})
	seqs_rev=$(readlink -f {input_fasta_rev})

	# the index is a set of 6 files in a subdirectory of
	# `resources/references`, named according to the nanobody library.
	# `bowtie2` only looks in the current directory for this index, so we
	# go there first, then tell it with `-x $library` to look for files like
	# `{library}..1.bt2`, etc.
	index={input_index}
	library={library}
	cd "$index"
	>&2 pwd

	retry_path="$output_coord_dir/$library-unaligned-%.fasta.gz"

	# perform alignment, compress to BAM file
	echo "Performing first alignment..."
	bowtie2 -x {library} \\
		-1 $seqs_fwd -2 $seqs_rev -f \\
		--ff -I 0 -X {max_fragment_length} \\
		--local --np {n_penalty} --n-ceil {n_ceil} \\
		--no-mixed --no-discordant \\
		--un-conc-gz "$retry_path" --no-unal \\
		--threads {threads} \\
	| samtools view -b -o "$output_coord_dir/{library}-first-unsorted.bam"

	echo "Performing more sensitive alignment on un-aligned sequences..."
	n_lines=$(zcat "${{retry_path/\%/1}}" | wc -l)
	let n_retry_seqs=n_lines/4
	echo "Aligning $n_retry_seqs sequences..."
	bowtie2 -x {library} \\
		-1 "${{retry_path/\%/1}}" -2 "${{retry_path/\%/2}}" -f \\
		--ff -I 0 -X {max_fragment_length} \\
		--local --np {n_penalty} --n-ceil {n_ceil} \\
		--no-mixed --no-discordant \\
		{retry_params} \\
		--threads {threads} \\
	| samtools view -b -o "$output_coord_dir/{library}-retry-unsorted.bam"

	# make name-sorted files for merging
	samtools sort -o "$output_name_dir/${{output_file_name_sorted}}-first" -n -@ {threads-1} "$output_coord_dir/{library}-first-unsorted.bam"
	samtools sort -o "$output_name_dir/${{output_file_name_sorted}}-retry" -n -@ {threads-1} "$output_coord_dir/{library}-retry-unsorted.bam"

	# merge name-sorted files
	samtools merge -o "$output_name_dir/$output_file_name_sorted" -n -@ {threads-1} "$output_name_dir/${{output_file_name_sorted}}-first" "$output_name_dir/${{output_file_name_sorted}}-retry"

	# sort merged file by coordinates
	output_file_coord_sorted=$(basename {output_bam_coord_sorted})
	samtools sort -o "$output_coord_dir/$output_file_coord_sorted" -@ {threads-1} "$output_name_dir/$output_file_name_sorted"
	samtools index "$output_coord_dir/$output_file_coord_sorted"
	""")

	if cleanup:

		cmd += dedent(f"""\
		# clean up temporary files
		rm "${{retry_path/\%/1}}" "${{retry_path/\%/2}}"
		rm "$output_name_dir/${{output_file_name_sorted}}-first" "$output_name_dir/${{output_file_name_sorted}}-retry"
		rm "$output_coord_dir/{library}-first-unsorted.bam" "$output_coord_dir/{library}-retry-unsorted.bam"
		""")

	with tempfile.NamedTemporaryFile('w') as script:
		print(script.name)
		print(cmd, file=script, flush=True)

		print(cmd)

		# if log is not None:
		# 	with open(log,'w') as logfile:
		# 		print(cmd, file=logfile, flush=True)
		# 		print("",file=logfile, flush=True)
		#
		# 		subprocess.run(['sh', script.name],
		# 			stdout=logfile,
		# 			stderr=logfile)
		# else: subprocess.run(['sh', script.name])
		subprocess.run(['sh', script.name], check=True, stdout=sys.stdout, stderr=sys.stderr)



# def align_local_bowtie(input_fasta_fwd, input_fasta_rev, input_index, library,
# 	output_bam_coord_sorted, output_bam_name_sorted,
# 	log=None, n_penalty=64, n_ceil=0, max_fragment_length=350, threads=1, conda=None):
#
# 	from textwrap import dedent
#
# 	if conda is not None:
# 		conda_activation = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda}'
# 	else:
# 		conda_activation = ''
#
# 	cmd = dedent(f"""\
# 	set -euox pipefail
# 	# activate conda env if given
# 	{conda_activation}
#
# 	# intermediate .sam and -unsorted.bam files will be produced before the
# 	# final .bam file. Put them in the same directory as the desired output
# 	# file.
# 	output_coord_dir=$(dirname $(readlink -m {output_bam_coord_sorted}))
# 	mkdir -p $output_coord_dir
#
# 	# get absolute path to the input FASTA files before we start changing
# 	# the directory
# 	seqs_fwd=$(readlink -f {input_fasta_fwd})
# 	seqs_rev=$(readlink -f {input_fasta_rev})
#
# 	# the index is a set of 6 files in a subdirectory of
# 	# `resources/references`, named according to the nanobody library.
# 	# `bowtie2` only looks in the current directory for this index, so we
# 	# go there first, then tell it with `-x $library` to look for files like
# 	# `{library}..1.bt2`, etc.
# 	index={input_index}
# 	library={library}
# 	cd "$index"
# 	>&2 pwd
#
# 	# perform alignment, compress to BAM file
# 	bowtie2 -x $library \\
# 		-1 $seqs_fwd -2 $seqs_rev -f \\
# 		--ff -I 0 -X {max_fragment_length} \\
# 		--local --np {n_penalty} --n-ceil {n_ceil} \\
# 		--threads {threads} \\
# 	| samtools view -b -o "$output_coord_dir/$library-unsorted.bam"
#
# 	# make a coordinate-sorted and indexed file for viewing
# 	output_file_coord_sorted=$(basename {output_bam_coord_sorted})
# 	samtools sort "$output_coord_dir/$library-unsorted.bam" -o "$output_coord_dir/$output_file_coord_sorted"
# 	samtools index "$output_coord_dir/$output_file_coord_sorted"
# 	""")
#
# 	if output_bam_name_sorted is not None:
#
# 		cmd = dedent(f"""\
# 		output_name_dir=$(dirname $(readlink -m {output_bam_name_sorted}))
# 		mkdir -p $output_name_dir
# 		output_file_name_sorted=$(basename {output_bam_name_sorted})
#
# 		""") + cmd + dedent(f"""\
#
# 		# sort by read name (rather than leftmost coordinates) so mate pairs appear consecutively
# 		samtools sort -n "$output_coord_dir/$library-unsorted.bam" -o "$output_name_dir/$output_file_name_sorted"
# 		""")
#
# 	with tempfile.NamedTemporaryFile('w') as script:
# 		print(script.name)
# 		print(cmd, file=script, flush=True)
#
# 		print(cmd)
#
# 		# if log is not None:
# 		# 	with open(log,'w') as logfile:
# 		# 		print(cmd, file=logfile, flush=True)
# 		# 		print("",file=logfile, flush=True)
# 		#
# 		# 		subprocess.run(['sh', script.name],
# 		# 			stdout=logfile,
# 		# 			stderr=logfile)
# 		# else: subprocess.run(['sh', script.name])
# 		subprocess.run(['sh', script.name], check=True)



def pad_sequence_na(seq, reference_cds_len_bp, first_codon_bp=0):
	"""makes a sequence the right length by padding with gap characters and/or trimming
	"""

	# pad the translated sequence to be same length as the reference
	# add dashes to beginning of sequence if it starts within the reference
	seq_padded = ('-' * first_codon_bp) + str(seq)

	# trim to no longer than the length of the reference sequence
	seq_padded = seq_padded[:reference_cds_len_bp]

	# add dashes to make same length as the reference sequence (if shorter)
	seq_padded = seq_padded + ('-' * (reference_cds_len_bp - len(seq_padded)))
	return seq_padded

def samfile_get_query_names(samfile_path):
	"""gets all read names from a sam/bam file"""
	import pysam
	with pysam.AlignmentFile(samfile_path, "rb") as samfile:
		query_names = [rec.query_name for rec in samfile.fetch(until_eof=True)]
		query_names = set(query_names)
	return query_names


def export_aligned_reads_paired(samfile_path=None,
	reference=None, reference_path=None, reference_frame_start=0, suppress_amber=True,
	min_fwd_end=0, max_rev_start=float("inf"),
	verbose=False, limit_records = float("inf")):
	"""Load paired-end reads from a bam file and assemble full-length aligned VHH gene sequences

	Parameters
	----------
	samfile_path : str
	    Path to a read name-sorted BAM file
	reference : str or Bio.Seq.Seq, optional
	    Reference sequence
	reference_path : str, optional
		Path to the reference sequence in FASTA format
	suppress_amber : bool, default=True
		When translating the reference sequence, True to translate UAG codons as Gln instead of Stop
	reference_frame_start : int
		Nucleotide position of the first in-frame codon in the reference
		sequence. This will define the reference frame used to translate the
		aligned sequences, and assembled sequences will be left-trimmed to this
		position.
	min_fwd_end : int, optional
		Read pairs where the forward read ends before this nucleic acid position
		on the reference sequence will be discarded.
		Use to enforce that the forward read must span all of
		CDR1--2 and end in the conserved region between CDR2 and CDR3, e.g. FR3
	max_rev_start : int, optional
		Read pairs where the reverse read begins after this nucleic acid position
		on the reference sequence will be discarded.
		Use to enforce that the reverse read must span all of
		CDR3 and end in the conserved region between CDR2 and CDR3, e.g. FR3

	Returns
	-------
	pd.DataFrame
	    index: names of the query sequences, e.g. pair_IDs
		columns:
		- seq: full length sequence, padded with gap characters to the length of
		  the translated reference sequence (i.e. until the end of the reference
		  sequence or the occurrence of the first in-frame stop codon, whichever
		  is first)
		- fwd: the part of the forward sequence used to assemble 'seq'; left-
		  trimmed to first in-frame codon
		- int: if fwd and rev overlap, the overlapping part of the sequence that
		  was used; if not, the portion of the
		- rev: the part of the reverse sequence used to assemble 'seq'; right-
		  trimmed to the end of the reference seq
	"""

	import pysam

	# read and translate reference sequence
	# read the full length sequence, because samfile coordinates will be given
	# in terms of reference_full.
	# but also trim the left side of the reference sequence to `reference_frame_start`,
	# so we can figure out the correct reference frame and translated length
	reference_full = get_reference(reference = reference, reference_path = reference_path)
	reference = get_reference(reference = reference_full, reference_frame_start = reference_frame_start)
	reference_translated = translate_reference(reference = reference, suppress_amber=suppress_amber)



	if verbose:
		print("Reference sequence:")
		positions = {}
		if min_fwd_end != 0:
			positions['min_fwd_end'] = min_fwd_end-reference_frame_start
		if max_rev_start < float('inf'):
			positions['max_rev_start'] = max_rev_start-reference_frame_start

		print(show_sequence_translated(reference, reference_translated, ruler={'start':reference_frame_start }, positions = positions))


	# create an empty DataFrame to hold results
	table = pd.DataFrame(index=samfile_get_query_names(samfile_path),
		columns=['seq', 'fwd','int','rev'])
	if verbose:
		print(f"Exporting from samfile with {len(table)} unique records...")
		if limit_records < float("inf"): print(f"(Saving at most {limit_records} records)")
		print("")
	n_records = 0;

	stats = {
		'unpaired': 0,
		'unmapped': 0,
		'fwd_too_short': 0,
		'rev_too_short': 0,
		'non_overlapping': 0,
		'overlapping': 0,
		'overlapping_disagree': 0,
	}

	with pysam.AlignmentFile(samfile_path, "rb") as samfile:

		# `samfile` is read-name sorted, in order
		for i, (fwd, rev) in enumerate(grouper(2, samfile.fetch(until_eof=True))):

			label = f"ID {fwd.query_name} (# {i:>8}, valid # {n_records:>8})"

			# check for valid mate pair
			if not fwd.query_name == rev.query_name:
				raise ValueError(f"  ERROR: adjacent read names do not match: {fwd.query_name} != {rev.query_name}; make sure samfile is sorted by read name `samtools sort -n` and mate pairs have the same read name")
				break
			if not (fwd.is_paired and rev.is_paired):
				if verbose:
					if fwd.is_unmapped and rev.is_unmapped: which = 'both reads'
					else: which = 'fwd read' if fwd.is_unmapped else 'rev read'

					print(f"WARNING: ID {fwd.query_name} / {rev.query_name}; SKIPPING, {which} unpaired")
				stats['unpaired'] += 1
				continue
			if (fwd.is_unmapped or rev.is_unmapped):
				if verbose > 1:
					if fwd.is_unmapped and rev.is_unmapped: which = 'both reads'
					else: which = 'fwd read' if fwd.is_unmapped else 'rev read'
					print(f"WARNING: {label}; SKIPPING, {which} is unmapped")
				stats['unmapped'] += 1
				continue

			# identify first complete codon
			first_codon = (((fwd.reference_start - reference_frame_start) // 3) + 1)
			first_query_base = fwd.query_alignment_start + (first_codon * 3 - (fwd.reference_start - reference_frame_start))
			fwd_seq_from_first_full_codon = fwd.query_sequence[first_query_base:fwd.query_alignment_end]

			# resolve sequence between reads or overlapping sequences

			# case 1: non-overlapping reads
			# fwd: ------------>
			# rev:                  <--------------
			if fwd.reference_end <= rev.reference_start:

				# check for both fwd too short and rev too short. If either (or
				# both) are the case we'll discard the read, but want to check
				# for both independently to get accurate stats.
				bad = False

				# Enforce that the forward read extends far enough into the
				# interior part of the reference (e.g. into FR3)
				if (fwd.reference_end < min_fwd_end):
					if verbose: print(f"WARNING: {label}; SKIPPING, fwd sequence ends at {fwd.reference_end}, before min_fwd_end {min_fwd_end}")
					stats['fwd_too_short'] += 1
					bad = True

				# Enforce that the reverse read extends far enough into the
				# interior part of the sequence (e.g. into FR3)
				if (rev.reference_start > max_rev_start):
					if verbose: print(f"WARNING: {label}; SKIPPING, rev sequence starts at {rev.reference_start}, after max_rev_start {max_rev_start}")
					stats['rev_too_short'] += 1
					bad = True

				if bad: continue

				# take entire forward read (starting in-frame), entire reverse read
				# (excluding soft clipped bases), and internal sequence from
				# reference.
				# CRITICAL: samfile coordinates refer to the full length
				# reference sequence, so ensure that the internal sequence is
				# drawn from there, otherwise it will be off by
				# `reference_frame_start` nt. 
				fwd_seq = fwd_seq_from_first_full_codon
				int_seq = str(reference_full[fwd.reference_end:rev.reference_start])
				rev_seq = rev.query_sequence[rev.query_alignment_start:]

				if verbose > 1:
					print(f"   INFO: {label}: Non-overlapping, filled {len(int_seq):>3} nt")
					# print(f"   INFO: {label}: {len(fwd_seq):>3} + {len(int_seq):>3} + {len(rev_seq):>3}")
				stats['non_overlapping'] += 1

			# case 2:
			# fwd: ------------>
			# rev:         <--------------
			elif fwd.reference_end > rev.reference_start:
				overlap_len = fwd.reference_end - rev.reference_start
				if verbose > 1: print(f"   INFO: {label}: Reads overlap by {overlap_len} nt")

				# remove overlapping sequence from fwd read
				fwd_seq = fwd_seq_from_first_full_codon[:-overlap_len]

				# find and compare overlapping sequences
				int_seq_fwd = fwd_seq_from_first_full_codon[-overlap_len:]
				int_seq_rev = rev.query_sequence[rev.query_alignment_start:][:overlap_len]

				if int_seq_fwd != int_seq_rev:
					stats['overlapping_disagree'] += 1
					if verbose:
						print(f"WARNING: {label}: overlapping sequences do not agree:\n"
						      f"         F> {int_seq_fwd}\n"
						      f"         R> {int_seq_rev}")
				int_seq = int_seq_fwd
				rev_seq = rev.query_sequence[rev.query_alignment_start:][overlap_len:]

				# if verbose: print(f"   ->    {label}: {len(fwd_seq):>3} + {len(int_seq):>3} + {len(rev_seq):>3}")
				stats['overlapping'] += 1

			# assemble full length sequence
			seq = (fwd_seq + int_seq + rev_seq)
			#translated = Seq(seq).translate(to_stop=True, table=999)

			# print debugging info
			if verbose > 2:
				print('         F> '+fwd_seq)
				print('         I> '+int_seq)
				print('         R> '+rev_seq)
				print('         => '+seq)
				print('--------------------------------------------------------------------------------')

			seq_padded = pad_sequence_na(seq, len(reference_translated)*3, first_codon*3)

			if not (len(seq_padded) == len(reference_translated)*3):
				raise ValueError(f"ERROR: {label}: {len(seq_padded)} nt != reference length {len(reference_translated)*3} nt; padded query sequence: {seq_padded}")

			table.loc[fwd.query_name,'fwd'] = fwd_seq
			table.loc[fwd.query_name,'int'] = int_seq
			table.loc[fwd.query_name,'rev'] = rev_seq
			table.loc[fwd.query_name,'seq'] = seq_padded

			n_records = n_records + 1
			if n_records > limit_records: break

	table = table.dropna()

	if verbose:
		print('================================================================================')
		print(f'Read {n_records} matched records from {i} records in samfile')
		print("Statistics: ")
		for k in stats:
			print(f" - {k:<24}: {stats[k]:>10} = {stats[k]/n_records:.1%}")
		print(f'Length of table: {len(table)}')
		print()
		print(table.head())
		print()

	return table


def export_aligned_reads(samfile_path=None,
	reference=None, reference_path=None, reference_frame_start=0, suppress_amber=True,
	verbose=True):
	"""Load alignment records from a samfile/bamfile, find the first codon of each read, and export sequence padded to reference

	Parameters
	----------
	samfile_path : str
	reference : Bio.Seq
	reference_path : str
	reference_frame_start : int, default=0
	suppress_amber : bool, default=True
		Translate TAG/UAG (Amber) stop codons as Gln

	Load data from `samfile_path` for `library` into a DataFrame. Use reference sequence
	at `reference_path` (whose relevant reading frame starts at `reference_frame_start`)
	to find the first complete codon.
	"""
	import pysam

	# read and translate reference sequence
	# read the full length sequence, because samfile coordinates will be given
	# in terms of reference_full.
	# but also trim the left side of the reference sequence to `reference_frame_start`,
	# so we can figure out the correct reference frame and translated length
	reference_full = get_reference(reference = reference, reference_path = reference_path)
	reference = get_reference(reference = reference_full, reference_frame_start = reference_frame_start)
	reference_translated = translate_reference(reference = reference, suppress_amber=suppress_amber)

	# create an empty DataFrame to hold results
	table = pd.DataFrame(index=samfile_get_query_names(samfile_path),
		columns=['seq'])

	if verbose: print(f"Translating {len(table)} sequences from samfile...")

	for row in samfile.fetch():

		# identify the first complete codon
		# the correct reading frame for the reference starts with base `reference_frame_start`

		# codon:       0     | 1     | 2
		# refer: 0 1 | 2 3 4 | 5 6 7 | 8
		#            | G C G | G C C | G C T
		# query: - -   - - G   G C C   G C T
		# reference_start = 4
		# reference_frame_start = 2
		# first_codon = (reference_start - reference_frame_start) // 3 + 1 = (4 - 2) // 3 + 1
		first_codon = (((row.reference_start - reference_frame_start) // 3) + 1)

		# codon:           | 0     | 1     | 2
		# refer: 0 1 2 3 4 | 5 6 7 | 8 9 10
		#                  | G C G | G C C | G C T
		# reference_start:       ^
		#
		# query:       0 1   2 3 4   5 6 7   8 9 10
		#        - - - A A   A A G   G C C   G C T
		# query_alignment_start: ^
		# want first_base =          ^
		#
		# reference_start = 7
		# reference_frame_start = 5
		# query_alignment_start = 4
		# first_codon = (reference_start - reference_frame_start) // 3 + 1 = (7 - 5) // 3 + 1 = 1
		# first_base = query_alignment_start + (first_codon * 3 - (reference_start - reference_frame_start)) = 4 + (1 * 3 - (7-5)) = 4 + 3 - 2 = 5
		first_query_base = row.query_alignment_start + (first_codon * 3 - (row.reference_start - reference_frame_start))

		# extract an in-frame sequence
		seq_from_first_full_codon = row.query_sequence[first_query_base:]

		# pad to length of reference
		seq_padded = pad_sequence_na(seq_from_first_full_codon, len(reference_translated)*3, first_query_base)

		table.loc[row.query_name, 'seq'] = seq_padded

	return table


def translate_aligned(seqs,
	reference=None, reference_path=None, reference_frame_start=0,
	reference_length_aa=None, ungap_internal=False, gap='-',
	suppress_amber=True, threads=1, verbose=False):
	"""translate a series of aligned VHH gene sequences

	Assumes each sequence is left- aligned and optionally right-padded with gap
	characters (-). Will identify the first complete codon in the sequence (that
	is in-frame with the translated reference sequence) and translate starting
	there. If the resulting sequence is shorter than the translated reference,
	the end will be padded with gap characthers; if longer, it will be trimmed.

	Parameters
	----------
	seqs : array_like or pd.Series
	reference : Bio.Seq or strL
	reference_path : str or path_like
	reference_frame_start : int, default=0
		nucleotide position where the reference sequence ORF begins (0-indexed
		position of the first codon)
	ungap_internal : bool, default=False
		True to remove all internal gap characters, in addition to gap
		characters at the beginning and end of the sequence. Set to true if
		``seqs`` are globally aligned (i.e. may contain internal gaps).
	suppress_amber : bool, default=True
		Translate TAG/UAG (Amber) stop codons as Gln
	threads : int
		Set to >1 to use parallel processing or 0 to use all available cores

	Returns
	-------
	list or pd.Series
		Translated sequences. If ``seqs`` is a pd.Series, will return an
		identically-indexed sequence with the same name. Otherwise returns
		translated sequences in a list in the same order
	"""
	#translated = seq.translate(to_stop=True, table=999)

	reference = get_reference(reference = reference, reference_path = reference_path, reference_frame_start = reference_frame_start)
	reference_translated = translate_reference(reference = reference, suppress_amber=suppress_amber)

	if reference_length_aa is not None:
		reference_length_nt = reference_length_aa * 3
		reference = reference[:reference_length_nt]
		reference_translated = reference_translated[:reference_length_aa]


	if verbose:

		print("Reference: ")
		if reference_length_aa is not None:
			print(f"(Trimmed to {reference_length_aa} aa:)")
		print(show_sequence_translated(reference, reference_translated))
		print("================================================================================")

	import re
	import Bio.Data.CodonTable
	from Bio.Data.CodonTable import TranslationError

	if suppress_amber:
		codon_table = Bio.Data.CodonTable.unambiguous_dna_by_id[999]
	else:
		codon_table = 'Standard'

	def translate_seq(seq):

		seq_lungapped = seq.lstrip(gap)
		first_codon = (len(seq) - len(seq_lungapped))//3

		if ungap_internal:
			seq_ungapped = ungapped(seq_lungapped)
		else:
			seq_ungapped = seq_lungapped.rstrip(gap)
			if gap in seq_ungapped:
				raise Exception(f'After stripping left and right gap '
					f'characters ("{gap}"), sequence still contains gaps. Set '
					"ungap_internal=True to remove internal gap characters. \n"
					f"Sequence: {seq}\n"
					f"Ungapped sequence: {seq_ungapped}")

		# trim sequence length to multiple of 3
		seq_ungapped = seq_ungapped[:(len(seq_ungapped)//3)*3]

		# if suppress_amber: translated = Seq(seq_ungapped).translate(to_stop=True, table=999)
		# else: translated = Seq(seq_lungapped).translate(to_stop=True)
		try:
			translated = Seq(seq_ungapped).translate(to_stop=True, table=codon_table)
		except TranslationError as err:
			# if ambiguous nucleotides ended up in the sequence, it will throw a
			# TranslationError. replace these with N's so the resulting codon
			# is X and the sequence can be filtered out later.
			seq_ungapped_fixed = re.sub(r"[^ATCG]", "N", seq_ungapped)

			print(f"ERROR: {seq} :: {err}")
			print(f" - {seq_ungapped} -> {seq_ungapped_fixed}")

			translated = Seq(seq_ungapped_fixed).translate(to_stop=True, table=codon_table)



		# pad the translated sequence to be same length as the reference
		# add dashes to beginning of sequence if it starts within the reference
		translated_padded = (gap * first_codon) + str(translated)

		# trim to no longer than the length of the reference sequence
		translated_padded = translated_padded[:len(reference_translated)]

		# add dashes to make same length as the reference sequence (if shorter)
		translated_padded = translated_padded + (gap * (len(reference_translated) - len(translated_padded)))

		if verbose:
			print(translated)
			print(translated_padded)
			print("================================================================================")

		assert len(translated_padded) == len(reference_translated)

		return(translated_padded)

	from joblib import Parallel, delayed
	if isinstance(seqs, pd.Series):
		_seqs = seqs.values
	else: _seqs = seqs

	_translated = Parallel(n_jobs=threads)(delayed(translate_seq)(seq) for seq in _seqs)

	return match_to_series(seqs, _translated)
	# return seqs.apply(translate_seq)


def get_pairwise_aligner(library, aa=False, mode='local', **kwargs):
	from Bio.Align import substitution_matrices

	aligner = Bio.Align.PairwiseAligner()
	aligner.alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ-'

	if aa:
		aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')

		# don't penalize matches to X, otherwise CDRs align badly
		aligner.substitution_matrix[:,'X'] = 0.5
	else:
		aligner.substitution_matrix = substitution_matrices.load('NUC.4.4')
		aligner.substitution_matrix[:,'N'] = 0.5

	# target = reference

	# ALPACA
	if library == 'alpaca':
		# goal: query sequence aligned to target (reference),
		# + with gaps inserted within CDRs;
		# + aligned query is NOT longer than reference
		# + no gaps are added to left side of target, because want to determine

		# the internal regions of the target (reference) sequence should be
		# generously long enough to accomodate variable length query sequences.
		aligner.target_internal_gap_score = -10000


		if mode == 'global':

			# in particular, clamp the alignment to the left side of the reference;
			# otherwise, determining CDRs by position won't work
			aligner.target_left_gap_score = -10000

			# # some reads may be longer than the target sequence; that's OK
			# aligner.target_right_gap_score = 0
			# likewise with the right side of the reference, otherwise it will just
			# push all of CDR3 off the right side
			aligner.target_right_gap_score = -1000

		elif mode == 'seq':
			aligner.target_left_gap_score  = -1
			aligner.target_right_gap_score = -1

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

	for k, v in kwargs.items():
		setattr(aligner, k, v)

	return aligner

def align_to_reference(seqs, reference_path, library,
	reference=None,
	aa=False,
	suppress_amber=True,
	reference_frame_start = 0,
	reference_length = None,
	mode='global',
	min_length=90,
	verbose=False, threads=None, **kwargs):

	from joblib import Parallel, delayed

	if reference is None:
		reference = get_reference(reference_path = reference_path, reference_frame_start = reference_frame_start)

	aligner = get_pairwise_aligner(library=library, aa=aa, mode=mode, **kwargs)

	if aa:
		reference = translate_reference(reference = reference, suppress_amber=suppress_amber)
		if reference_length is not None:
			reference = reference[:reference_length]
	min_length = min([min_length, len(reference)])

	pool = Parallel(n_jobs=threads)
	out = pool(
		delayed(_align_to_reference)(seq, reference, aligner, min_length=min_length, verbose=verbose)
			for seq in seqs)

	return match_to_series(seqs, out)


def _align_to_reference(seq, reference, aligner, min_length=90, verbose=False, return_alignment=False):
	seq = seq.strip('-')
	if len(seq) < min_length:
		if verbose: print(f'too short, {len(seq)} < {min_length}')
		return ''
	alignments = aligner.align(reference,seq)

	# lots of alignments suggests poor alignment quality
	if len(alignments) < 1000:
		alignments = sorted(alignments)

	if return_alignment:
		return alignments[0]

	alignment = str(alignments[0])

	ref,aln,s,*_ = alignment.split("\n")
	loffset = len(ref) - len(ref.lstrip('-'))
	roffset = len(ref) - len(ref.rstrip('-'))

	# s[i:-0] will produce an empty string, whereas s[i:None] == s[i:]
	st = s[loffset:-roffset or None]
	if verbose:
		print(
f"""{len(alignments)} alignments
{str(alignment)}
trimming {loffset} from left / {roffset} from right > {st}
================================================================================
"""
			)
	return st


def pairwise_align(reference, seq, library,
	aa=False,
	mode='global',
	return_alignment=False,
	**kwargs):

	"""aligns a single query sequence to a reference
	"""

	aligner = get_pairwise_aligner(library=library, aa=aa, mode=mode, **kwargs)
	return _align_to_reference(seq, reference, aligner, verbose=False, min_length=0, return_alignment=return_alignment)


def align_chromatograms(paths, reference_path, CDRs=None, reference_frame_start_nt=0, rc=False):
    """aligns chromatograms to a reference sequence and extracts CDRs
    paths: iterable of str, paths to chromatograms in .scf format
    reference_path: str, path to reference sequence
    reference_frame_start_nt: int, nt where the reading frame of the reference sequence begins
    CDRs: dict of str:tuple, keys are CDR names, values are (start_nt, end_nt) positions for the CDR, 0-indexed, half-open
    rc: bool, True to reverse-complement chromatogram sequence before aligning to reference
    """

    seqs = {}
    for f in paths:
        seq, *_ = bioconvert.io.scf.read_scf(f)
        if rc: seq = nbseq.utils.rc(seq)
        seqs[f] = seq

    out = nbseq.asvs.align.align_to_reference(seqs.values(),reference_path=reference_path,
                 reference_frame_start=reference_frame_start_nt,
                 library='alpaca',
                 mode='local', verbose=False
          )

    df = pd.DataFrame({'file':list(seqs.keys()), 'seq':out})

    if CDRs is not None:
        df = pd.merge(df, nbseq.extract_CDRs(df['seq'], CDRs=CDRs),
                      left_index=True, right_index=True)
    df['name'] = df['file'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    return df
