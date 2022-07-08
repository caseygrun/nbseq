import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from ..utils import lookup

def get_identifier(space):
	identifier = {
		'CDR3':'CDR3ID',
		'CDR':'CDRID',
		'AA':'aaSVID',
		'NA':'ASVID'
	}[space.upper()]
	return identifier

def remove_chimeras(seqs, abundances, ids, chimera_algorithm = "uchime_denovo",
	infile_fasta=None, outfile_chimeras=None, outfile_nonchimeras=None, return_nonchimeras=False):

	with tempfile.TemporaryDirectory() as tmpdir:
		if infile_fasta is None:
			infile_fasta = Path(tmpdir) / 'input.fasta'

		with open(infile_fasta) as input:

			for seq, abundance, id in zip(seqs, abundances, ids):
				# write to fasta file with >size={abundance};{id}
				input.write(f">size={abundance};{id}\n")
				input.write(f"{seq}\n")

		if outfile_nonchimeras is None:
			outfile_nonchimeras = str(Path(tmpdir) / 'nonchimeras.fasta')

		cmd = ["vsearch", "--sortbysize", f"--{chimera_algorithm}", input_fasta]
		if outfile_chimeras is not None:
			cmd += [f"--chimeras {outfile_chimeras}"]
		if outfile_nonchimeras is not None:
			cmd += [f"--nonchimeras {outfile_nonchimeras}"]

		print(cmd)
		subprocess.run(cmd)

		if return_nonchimeras:
			return list(SeqIO.parse(outfile_nonchimeras))


def fill_end_gaps_by_reference(asvs, reference):
	"""replaces 5' and 3' gap characters with those from the reference sequence

	each sequence in `asvs` should already be aligned to the reference
	"""

	def fill_end_gap_by_reference(asv, reference):
		l = len(asv)
		lungapped = asv.lstrip(gap)
		l_gaps = l - len(lungapped)
		ungapped = lungapped.rstrip(gap)
		r_gaps = l - l_gaps - len(ungapped)

		return reference[:l_gaps] + ungapped + reference[-r_gaps:]

	if not isinstance(asvs, pd.Series):
		asvs = pd.Series(asvs)

	return asvs.apply(fill_end_gap_by_reference)

# def collapse_no_mismatch(asvs):
# 	"""groups ASVs into zero-width OTUs"""
#
# 	# asvs: DataFrame(index = names, columns=['aligned','abundance'])
#
# 	# return: DataFrame(index=names, columns=['seq','name'])
#
#
# 	for (seq1, abd1), (seq2, abd2):
# 		# find shorter ungapped sequence
#
#
# 		# crop both sequences to same length
#
# 		# check if sequences are equal
# 			# map both sequences to the higher abundance sequence
# 		# else
# 			# map both sequences to one another
#
#
# 	pass


def find_similar_features(query, db=None, space='CDR3', ft=None, fd=None, verbose=False, identifier=None, **kwargs):
	"""Finds features with similar sequences by searching with mmseqs2, and optionally joins them to feature data

	Parameters
	----------
	query : pd.DataFrame, pd.Series, dict, or iterable
		named feature sequences to query
	space : str
		name of feature space to query; options include: 'AA', 'CDR', 'CDR3'
	db : optional
		str path to a mmseqs2 database; if not given, will be inferred from space
	fd : pd.DataFrame or pd.Series, optional
		if given, this feature data will be joined to the resultint matches;
		should be feature data, where the index uses the same identifer as
		``query``.
	ft : AnnData, optional
		if given and fd is None, feature data from ft.var will be joined to the
		results

	Returns
	-------
	pd.DataFrame
		matching features, joined to fd or ft.var if given. Will include at
		least the following columns:
		- ``query_{identifier}``
		- ``search_{identifier}``
		- ``pident``: percent identity of the match
	"""

	from .db import search_mmseqs2

	if identifier is None:
		identifier = get_identifier(space)


	if isinstance(query, pd.DataFrame) and (identifier in query.columns) and (space in query.columns):
		query = query[[identifier,space]].drop_duplicates().set_index(identifier)

	if db is None:
		db = f'intermediate/{space.lower()}/features_db/features'
	matches = (
		search_mmseqs2(query, db, conda='mmseqs2-vsearch',
				verbose=verbose, **kwargs)
			.rename(columns={
				'qseqid':f'query_{identifier}',
				'sseqid':f'search_{identifier}' })
		)

	if ((fd is None) and not (ft is None)):
		fd = ft.var

	if fd is not None:
		if isinstance(fd, pd.Series):
			fd = fd.to_frame()
		# return (pd.merge(
		# 		left=matches, left_on=f'search_{identifier}',
		# 		right=fd, right_index=True, how='left'
		# 	)
		# )
		# return fd.join(matches, on=f'search_{identifier}',)
		return matches.join(fd, on=f'search_{identifier}')
	else:
		return matches


def project(features, from_space='CDR3', to_space='AA', method='sql', columns='*', mapping=None, fd=None, ft=None, show_query=False, **kwargs):
	"""Projects from one feature space to another, optionally retrieving feature data for the new feature space and joining to an existing feature table

	For large feature datasets (e.g. nucleic acid or amino acid sequences), this
	is best accomplished by querying a pre-populated feature space.

	- method='sql', ft is not None: returns a filtered feature table ``ft``, with feature data joined to ft.var
	- method='sql', ft=None: returns pd.DataFrame containing feature data
	- method='mapping': returns filtered version of the feature table ``ft``
	"""
	identifiers = {
		'NA': 'ASVID',
		'AA': 'aaSVID',
		'CDR3': 'CDR3ID'
	}

	from_space = from_space.upper()
	to_space = to_space.upper()

	if method == 'sql':
		from .db import search_sql, make_sql_where_clause

		src_identifier = identifiers[from_space]
		new_identifier = identifiers[to_space]

		if from_space =='CDR3' and to_space == 'AA':
			if 'aligned' in columns:

				# will require a join to the `aa` table. this is extremely slow
				# unless a proper, ideally covering, index is used.
				# two things that help with this:
				# 1) filter by library='{library}' in the WHERE
				# clause, which will reduce the size of the AA table that needs
				# to be considered
				# 2) filter the CDR3ID using LIKE; this is the default for
				# make_sql_where_clause, so it shouldn't be a problem here.
				# Note that using CDR3ID = '{CDR3ID}' instead of
				# CDR3ID LIKE '{CDR3ID}%' is ~2x slower

				columns = set(columns) - {'aligned'}
				cols = [f"`cdrs`.`{col}`" for col in columns] + ['`aa`.`aligned`']
				cols = ",".join(cols)

				predicates = {}
				if "library" in kwargs:
					# select library from the `aa` table, rather than `cdrs` table.
					# pass value as 1-tuple so that query is built with
					# `aa`.`{library}` = {library}, rather than LIKE '{library}%'.
					# LIKE will not use a covering index, which will be ~5x slower
					predicates["`aa`.`library`"] = (kwargs['library'],)
					del kwargs["library"]
				else:
					print("Warning: project() from CDR3 to AA can be slow if library is not specified as a kwarg")

				predicates[f"`cdrs`.`{src_identifier}`"] = features
				for k, v in kwargs.items():
					predicates[f"`cdrs`.`{k}`"] = v

				predicates = make_sql_where_clause(**predicates)
				query = f"SELECT {cols} FROM `cdrs` INNER JOIN `aa` ON aa.aaSVID == cdrs.aaSVID WHERE ({predicates})"

				fd = search_sql(query=query, show_query=show_query)
			else:
				fd = search_sql(**{'table':'cdrs', 'columns':columns, src_identifier:features, 'show_query':show_query, **kwargs})

			if ft is not None:
				ft = ft[:,fd[new_identifier]]

				fd.set_index(new_identifier, inplace=True)

				ft.var = ft.var.join(fd)
				# ft.var = pd.merge(
				# 	left=ft.var, left_index=True,
				# 	right=fd.set_index(new_identifier),
				# 	right_index=True, how='left')

				return ft
			else:
				return fd

		else:
			raise NotImplementedError(f"No method for projecting from feature space {from_space} to feature space {to_space} with method `sql`.")

	elif method == 'mapping':
		from anndata import AnnData

		# use a dict, Series, or callable mapping to find relevant columns in a target feature table

		if mapping is None:
			raise ValueError("Must provide a mapping for method='mapping'")
		new_features = map_generic(features)

		if fd is not None:
			return fd.loc[new_features]

		elif ft is not None:
			if not isinstance(ft,AnnData):
				raise TypeError("If provided, the feature table must be in AnnData format")
			return ft[:,new_features]
		else:
			return new_features


#
#
# def consensus(seqs, alphabet, table=None):
# 	first = next(iter(seqs)))
# 	if table is None:
# 		if alphabet is None:
# 			alphabet = list(set(first))
# 		if isinstance(alphabet, str):
# 			alphabet = list(alphabet)
#
# 		table = pd.DataFrame(data=0, index=alphabet, columns=range(len(first)))



def check(condition, message):
	if not condition:
		raise Exception(message)

def marginal_frequencies(sequences,
						counts=None,
						abundances=None,
						matrix=None,
						alphabet='ACDEFGHIKLMNPQRSTVWY-',
						characters_to_ignore=''):
	"""
	Generates matrix from a sequence alignment

	From https://github.com/jbkinney/logomaker, Copyright (c) 2019 Ammar Tareen and Justin B. Kinney
	MIT License


	parameters
	----------
	sequences: (list of strings)
		A list of sequences, all of which must be the same length
	counts: (None or list of numbers)
		If not None, must be a list of numbers the same length os sequences,
		containing the (nonnegative) number of times that each sequence was
		observed. If None, defaults to 1.
	to_type: (str)
		The type of matrix to output. Must be 'counts', 'probability',
		'weight', or 'information'
	characters_to_ignore: (str)
		Characters to ignore within sequences. This is often needed when
		creating matrices from gapped alignments.

	returns
	-------
	out_df: pd.DataFrame
		A matrix of the requested type.
	"""

	# validate inputs

	# Make sure sequences is list-like
	check(isinstance(sequences, (list, tuple, np.ndarray, pd.Series)),
		  'sequences must be a list, tuple, np.ndarray, or pd.Series.')
	sequences = list(sequences)

	# Make sure sequences has at least 1 element
	check(len(sequences) > 0, 'sequences must have length > 0.')

	# Make sure all elements are sequences
	check(all(isinstance(seq, str) for seq in sequences),
		  'sequences must all be of type string')

	# validate characters_to_ignore
	check(isinstance(characters_to_ignore, str),
		  'type(seq) = %s must be of type str' % type(characters_to_ignore))

	# Get sequence length
	L = len(sequences[0])

	# Make sure all sequences are the same length
	check(all([len(s) == L for s in sequences]),
		  'all elements of sequences must have the same length.')

	# validate counts as list-like
	if counts is None and abundances is not None:
		counts = abundances
	check(isinstance(counts, (list, tuple, np.ndarray, pd.Series)) or
			(counts is None),
			'counts must be None or a list, tuple, np.ndarray, or pd.Series.')

	# make sure counts has the same length as sequences
	if counts is None:
		counts = np.ones(len(sequences))
	else:
		check(len(counts) == len(sequences),
			  'counts must be the same length as sequences;'
			  'len(counts) = %d; len(sequences) = %d' %
			  (len(counts), len(sequences)))

	# Create a 2D array of characters
	char_array = np.array([np.array(list(seq)) for seq in sequences])

	# Get list of unique characters
	if alphabet is None:
		unique_characters = np.unique(char_array.ravel())
	else:
		unique_characters = np.unique(list(alphabet))
	unique_characters.sort()

	# Remove characters to ignore
	columns = [c for c in unique_characters if not c in characters_to_ignore]
	index = list(range(L))
	counts_df = pd.DataFrame(data=0, columns=columns, index=index)

	# Sum of the number of occurrences of each character at each position
	for c in columns:
		tmp_mat = (char_array == c).astype(float) * counts[:, np.newaxis]
		counts_df.loc[:, c] = tmp_mat.sum(axis=0).T

	return counts_df


def consensus(marginals=None, seqs=None, min_fraction=0.1, ambiguous='X', min_count=None, gap='-'):
	if marginals is None and seqs is not None:
		marginals = marginal_frequencies(seqs)
	elif marginals is None and seqs is None:
		raise ValueError("Must specify marginals or seqs")


	row_sum = marginals.sum(axis=1, skipna=True)


	# # if absolute row_sum < ??? -> convert position to gap
	# if min_count == 'auto':
	# 	min_count = max(row_sum.mean() - (2 * row_sum.std()), 0)
	# if min_count is not None:
	# 	row_sum[row_sum < min_count] = 0

	norm = marginals.div(row_sum, axis=0)

	# 0 / 0 = Inf, apparently. Set to nan so these positions are gapped or
	# dropped
	norm[norm == float('inf')] = float('nan')



	# find most abundant character at position
	cons = norm.idxmax(axis=1, skipna=True).fillna(ambiguous)

	# find abundance of the most abundant character at each position
	frac = lookup(norm, cons.index, cons.values, 0)

	# if the fraction of most abundant character at a position was < min_fraction, mark as ambiguous
	cons[frac < min_fraction] = ambiguous

	if gap is not None:
		cons = cons.fillna(gap)
	else:
		cons = cons.dropna()

	return "".join(cons)
