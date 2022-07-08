from pathlib import Path
import hashlib

import pandas as pd


from .utils import *
from .cdrs import *
from .ft import to_relative, add_feature_data, read as ft_read, query as ft_query, query_ids

class Experiment:
	"""Keep track of various feature tables, databases, etc."""

	def __init__(self, fts={}, trees={}, phenotypes=None, sql_db=None, mmseqs_dbs={}, config={}, name=None):
		self._fts = {k.lower(): v for k, v in fts.items()}
		self._trees = {k.lower(): v for k, v in trees.items()}
		self._phenotypes = phenotypes
		self._sql_db = sql_db
		self._mmseqs_dbs = {k.lower(): v for k, v in mmseqs_dbs.items()}
		self._config = config
		self.name = name

	_rfts = None
	_viz = None
	_phenotype_matrix = None
	_antigen_matrix = None

	def __repr__(self):
		out = [f"Experiment('{self.name}') with layers:"]

		for k, ft in self._fts.items():
			out.append(f"- {k:<8}: {ft.shape[0]} x {ft.shape[1]}, database: {self._mmseqs_dbs.get(k, 'None')}")
			out.append(f"	var: {ft.var.columns.values}")

		out.append(f"	obs: {ft.obs.columns.values}")
		out.append(f"SQL: {self._sql_db}")
		return "\n".join(out)

	@property
	def fts(self):
		return dotdict(self._fts)

	@property
	def obs(self):
		return self._fts[next(iter(self._fts))].obs

	@property
	def tree(self):
		return dotdict(self._trees)

	@property
	def rfts(self):
		if self._rfts is None:
			self._rfts = {k: to_relative(v) for k,v in self._fts.items()}
		return dotdict(self._rfts)

	@property
	def config(self):
		return dict(self._config)

	@property
	def viz(self):
		from .viz import ExperimentVisualizer
		if self._viz is None:
			self._viz = ExperimentVisualizer(self)
		return self._viz

	@property
	def ag_names(self):
		return self.ags.index.values

	@property
	def ags(self):
		return self._phenotypes.query("type == 'antigen'")

	@property
	def ag_matrix(self):
		if self._antigen_matrix is None:
			self._antigen_matrix = self.obs.loc[:,self.ag_names]
		return self._antigen_matrix

	@property
	def pheno_names(self):
		return self._phenotypes.index.values

	@property
	def pheno_matrix(self):
		if self._phenotype_matrix is None:
			self._phenotype_matrix = self.obs.loc[:,self.pheno_names]
		return self._phenotype_matrix

	@staticmethod
	def from_files(directory='.', metadata='config/metadata-phenotypes.csv',
		phenotypes='config/phenotypes.csv',
		fd_cdr3='intermediate/cdr3/features/all/asvs.csv',
		ft_aa='results/tables/aa/feature_table.biom',
		ft_cdr3='results/tables/cdr3/feature_table.biom',
		tree_aa='intermediate/aa/features/top_asvs/alpaca/asvs.nwk',
		tree_cdr3='intermediate/cdr3/features/top_asvs/alpaca/asvs.nwk',
		config='config/config.yaml',
		sql_db='intermediate/aa/asvs.db',
		mmseqs_db_aa='intermediate/aa/features_db/features',
		mmseqs_db_cdr3='intermediate/cdr3/features_db/features',
		verbose=True, **kwargs):

		import yaml
		import time
		from datetime import timedelta
		import humanize
		import skbio.tree

		start_time = time.time()

		realdir = os.path.realpath(directory)
		name = os.path.basename(realdir)

		if verbose:
			print(f"Loading experiment {name} from '{realdir}'...")
			print(f"- Reading metadata from {metadata} ...")


		metadata = read_delim_auto(Path(directory) / metadata).set_index('ID')

		fd_cdr3 = Path(directory) / fd_cdr3
		if not os.path.isfile(fd_cdr3):
			print(f"- Warning: feature data table '{fd_cdr3}' does not exist!")
		else:
			if verbose: print(f"- Reading feature data for table 'cdr3' from {fd_cdr3} ...")
			feature_data_cdr3 = read_delim_auto(fd_cdr3).drop_duplicates('CDR3ID').set_index('CDR3ID')

		if not os.path.isfile(Path(directory) / phenotypes):
			print(f"- Warning: phenotypes table '{phenotypes}' does not exist!")
		else:
			if verbose: print(f"- Reading phenotypes from {phenotypes} ...")
			phenotypes = read_delim_auto(phenotypes).set_index('name', drop=True)


		if verbose:
			print(f"- Reading Config from {config} ...")
		with open(config,'r') as f:
			config = yaml.load(f,yaml.FullLoader)

		sql_db_path = os.path.abspath(Path(directory) / sql_db)
		sql_db = 'sqlite:///' + str(sql_db_path)
		if not os.path.isfile(sql_db_path):
			print(f"- Warning: sqlite database '{sql_db_path}' does not exist")
		else:
			if verbose: print(f"- Using SQL database at '{sql_db}'")

		mmseqs_dbs = {
			'aa': Path(directory) / mmseqs_db_aa,
			'cdr3': Path(directory) / mmseqs_db_cdr3
		}
		for k, v in mmseqs_dbs.items():
			if not os.path.isfile(v):
				print(f"- Warning: mmseqs2 database '{k}' at '{v}' does not exist!")
			elif verbose:
				print(f"- Using mmseqs2 database '{k}' at '{v}'")

		if tree_aa is not None:
			tree_aa = Path(directory) / tree_aa
			if verbose:
				print(f"- Reading AA phylogeny from {tree_aa} ({humanize.naturalsize(tree_aa.stat().st_size)})...")
			tree_aa = skbio.tree.TreeNode.read(str(tree_aa), format="newick")

		if tree_cdr3 is not None:
			tree_cdr3 = Path(directory) / tree_cdr3
			if verbose:
				print(f"- Reading CDR3 phylogeny from {tree_cdr3} ({humanize.naturalsize(tree_cdr3.stat().st_size)})...")
			tree_cdr3 = skbio.tree.TreeNode.read(str(tree_cdr3), format="newick")

		ft_aa = Path(directory) / ft_aa
		if verbose:
			print(f"- Reading AA feature table from {ft_aa} ({humanize.naturalsize(ft_aa.stat().st_size)})...")
		ft_aa   = ft_read(ft_aa,
			to='anndata',
			metadata=metadata
		)

		ft_cdr3 = Path(directory) / ft_cdr3
		if verbose:
			print(f"- Reading CDR3 feature table from {ft_cdr3} ({humanize.naturalsize(ft_cdr3.stat().st_size)})...")
		ft_cdr3 = ft_read(ft_cdr3,
			to='anndata',
			metadata=metadata,
			feature_data=feature_data_cdr3
		)
		ft_cdr3 = add_feature_data(ft_cdr3, 'reads', 'nsamples')

		if verbose:
			print("Finished in " + humanize.precisedelta(timedelta(seconds=start_time-time.time())))

		return Experiment(
			fts = {
				'aa': ft_aa,
				'cdr3': ft_cdr3
			},
			trees = {
				'aa': tree_aa,
				'cdr3': tree_cdr3
			},
			phenotypes = phenotypes,
			sql_db = sql_db,
			mmseqs_dbs = mmseqs_dbs,
			config = config,
			name = name
		)


	def get_CDR_positions(self, library='alpaca', aa=True):
		cdrs = self.config['libraries'][library]['CDRs']

		if aa:
			return cdrs
		else:
			return {k: [v[0] * 3, v[1] * 3] for k,v in cdrs.items()}

	def get_clonotype_positions(self, library='alpaca',aa=True):
		cdrs = self.get_CDR_positions(library, aa)
		start = min(v[0] for k,v in cdrs.items())
		stop  = max(v[1] for k,v in cdrs.items())
		return (start, stop)

	def search_sql(self, **kwargs):
		from .asvs.db import search_sql
		return search_sql(db=self._sql_db, **kwargs)


	def find_similar_features(self, query, space='cdr3', **kwargs):
		"""
		See :py:func:`asvs.find_similar_features`
		"""

		from .asvs import find_similar_features
		space = space.lower()
		return find_similar_features(query,
			ft=self.fts[space],
			db=self._mmseqs_dbs[space],
			space=space,
			**kwargs)

	def project(self, features, from_space='cdr3', to_space='aa', ft=True, **kwargs):
		"""
		See :py:func:`asvs.project`
		"""
		from .asvs import project

		if ft:
			ft = self._fts[to_space.lower()]
		else:
			ft = None
		return project(features, from_space, to_space,
			ft=ft,
			db=self._sql_db)

	def query(self, query, axis='sample', space='cdr3', **kwargs):
		return ft_query(self.fts[space], query, axis=axis, **kwargs)

	def query_ids(self, query, axis='sample', space='cdr3', **kwargs):
		return query_ids(self.fts[space], query, axis=axis, **kwargs)

	def find_feature(self, feature, single=False, space='cdr3'):
		ft = self.fts[space.lower()]
		ids = ft[:,ft.var.index.str.startswith(feature)].var_names.values
		if single:
			if len(ids) > 1:
				raise ValueError(f"More than one {space.upper()} feature found with prefix '{feature}': {list(ids)}")
			else:
				return ids[0]
		return ids

	def find_cdr3(self, CDR3ID, single=False):
		return self.find_feature(feature=CDR3ID, single=single, space='cdr3')

	def find_library(self, feature, space='cdr3'):
		return self.fd(feature, fields='library', space=space)

	def feature_data(self, features, fields, space='cdr3'):
		return self.fts[space.lower()].var.loc[features, fields]

def translate_sam_alignment_paired(samfile_path=None,
	reference=None, reference_path=None, reference_frame_start=0,
	verbose=False, limit_records = float("inf")):

	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	import Bio.SeqIO
	import pysam
	reference = get_reference(reference = reference, reference_path = reference_path, reference_frame_start = reference_frame_start)
	reference_translated = translate_reference(reference = reference, reference_path = reference_path, reference_frame_start = reference_frame_start)

	# find all unique sequence names
	with pysam.AlignmentFile(samfile_path, "rb") as samfile:
		query_names = [rec.query_name for rec in samfile.fetch(until_eof=True)]
		query_names = set(query_names)

	# create an empty DataFrame to hold results
	table = pd.DataFrame(index=query_names, columns=['translated', 'seq', 'fwd','int','rev'])
	if verbose:
		print(f"Translating from samfile with {len(query_names)} unique records...")
		if limit_records < float("inf"): print(f"(Saving at most {limit_records} records)")
		print("")
	n_records = 0;

	with pysam.AlignmentFile(samfile_path, "rb") as samfile:

		# `samfile` is read-name sorted, in order
		for i, (fwd, rev) in enumerate(grouper(2, samfile.fetch(until_eof=True))):

			# check for valid mate pair
			if not fwd.query_name == rev.query_name:
				if verbose: print(f"Error: adjacent read names do not match: {fwd.query_name} != {rev.query_name}; make sure samfile is sorted by read name `samtools sort -n` and mate pairs have the same read name")
				break
			if not (fwd.is_paired and rev.is_paired):
				if verbose: print(f"Warning: skipping {fwd.query_name} / {rev.query_name}; one of the mates is unpaired")
				continue
			if (fwd.is_unmapped or rev.is_unmapped):
				if verbose > 1:
					if fwd.is_unmapped and rev.is_unmapped: which = 'both reads'
					else: which = 'fwd read' if fwd.is_unmapped else 'rev read'
					print(f"Skipping {fwd.query_name} / {rev.query_name}, {which} unpaired")
				continue
			if verbose: print(f"Record {fwd.query_name}, record {i}, valid record {n_records}")
			n_records = n_records + 1

			# identify first complete codon
			first_codon = (((fwd.reference_start - reference_frame_start) // 3) + 1)
			first_query_base = fwd.query_alignment_start + (first_codon * 3 - (fwd.reference_start - reference_frame_start))
			fwd_seq_from_first_full_codon = fwd.query_sequence[first_query_base:fwd.query_alignment_end]

			# resolve sequence between reads or overlapping sequences

			# case 1: non-overlapping reads
			# fwd: ------------>
			# rev:                  <--------------
			if fwd.reference_end <= rev.reference_start:

				# take entire forward read (starting in-frame), entire reverse read
				# (excluding soft clipped bases), and internal sequence from
				# reference
				fwd_seq = fwd_seq_from_first_full_codon
				int_seq = str(reference[fwd.reference_end:rev.reference_start])
				rev_seq = rev.query_sequence[rev.query_alignment_start:]

			# case 2:
			# fwd: ------------>
			# rev:         <--------------
			elif fwd.reference_end > rev.reference_start:
				overlap_len = fwd.reference_end - rev.reference_start
				if verbose: print(f"Reads overlap by {overlap_len} nt")

				# remove overlapping sequence from fwd read
				fwd_seq = fwd_seq_from_first_full_codon[:-overlap_len]

				# find and compare overlapping sequences
				int_seq_fwd = fwd_seq_from_first_full_codon[-overlap_len:]
				int_seq_rev = rev.query_sequence[rev.query_alignment_start:][:overlap_len]

				if int_seq_fwd != int_seq_rev:
					if verbose:
						print("Warning: overlapping sequences do not agree:")
						print("  F> "+int_seq_fwd)
						print("  R> "+int_seq_rev)
				int_seq = int_seq_fwd
				rev_seq = rev.query_sequence[rev.query_alignment_start:][overlap_len:]

			# assemble full sequence, and translate, reading through amber stop codons
			seq = Seq(fwd_seq + int_seq + rev_seq)
			translated = seq.translate(to_stop=True, table=999)

			# print debugging info
			if verbose:
				print('F> '+fwd_seq)
				print('I> '+int_seq)
				print('R> '+rev_seq)


				print('=> '+seq)
				print('T> '+translated)
				print('--------------------------------------------------------------------------------')

			# pad the translated sequence to be same length as the reference
			# add dashes to beginning of sequence if it starts within the reference
			translated_padded = ('-' * first_codon) + str(translated)
			seq_padded = ('-' * (first_codon * 3)) + str(seq)

			# trim to no longer than the length of the reference sequence
			translated_padded = translated_padded[:len(reference_translated)]
			seq_padded = seq_padded[:len(reference_translated) * 3]

			# add dashes to make same length as the reference sequence (if shorter)
			translated_padded = translated_padded + ('-' * (len(reference_translated) - len(translated_padded)))
			seq_padded = seq_padded + ('-' * (3 * (len(reference_translated) - len(translated_padded))))

			assert len(translated_padded) == len(reference_translated)

			table.loc[fwd.query_name,'fwd'] = fwd_seq
			table.loc[fwd.query_name,'int'] = int_seq
			table.loc[fwd.query_name,'rev'] = rev_seq
			table.loc[fwd.query_name,'seq'] = seq_padded
			table.loc[fwd.query_name,'translated'] = translated_padded

			if n_records > limit_records: break

	table = table.dropna()

	if verbose:
		print('================================================================================')
		print(f'Read {n_records} matched records from {i} records in samfile')
		print(f'Length of table: {len(table)}')
		print()
		print(table.head())
		print()

	return table

def translate_sam_alignment(samfile=None, samfile_path=None, reference=None, reference_path=None, reference_frame_start=0):
	"""Load alignment records from a samfile/bamfile, find the first codon of each read, and translate that read.

	Parameters
	----------
	samfile : pysam.AlignmentFile
	samfile_path : str
	reference : Bio.Seq
	reference_path : str
	reference_frame_start : int, default=0

	Load data from `samfile_path` for `library` into a DataFrame. Use reference sequence
	at `reference_path` (whose relevant reading frame starts at `reference_frame_start`)
	to find the first complete codon and translate each aligned sequence.
	"""

	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	import Bio.SeqIO
	import pysam


	if samfile is None:
		if samfile_path is None:
			raise Exception("Must specify either `samfile` or `samfile_path`")
		samfile = pysam.AlignmentFile(samfile_path, "rb")

	# read and translate reference sequence
	# if reference is None:
	# 	if reference_path is None:
	# 		raise Exception("Must specify either `reference` or `reference_path`")
	# 	reference = next(Bio.SeqIO.parse(reference_path,"fasta"))
	# reference_translated = reference[reference_frame_start:].translate(to_stop=True)
	reference_translated = translate_reference(reference = reference, reference_path = reference_path, reference_frame_start = reference_frame_start)

	# allocate empty DataFrame, one row per unique ASV
	alignment_table = pd.Series(index=[row.query_sequence for row in samfile.fetch()], name='translation')

	print(f"Translating {len(alignment_table)} sequences from samfile...")

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
		seq = Seq(seq_from_first_full_codon)

		# translate
		# TODO: handle amber suppression
		translated = seq.translate(to_stop=True)

		# pad the translated sequence to be same length as the reference
		# add dashes to beginning of sequence if it starts within the reference
		translated_padded = ('-' * first_codon) + str(translated)
		# trim to no longer than the length of the reference sequence
		translated_padded = translated_padded[:len(reference_translated)]
		# add dashes to make same length as the reference sequence (if shorter)
		translated_padded = translated_padded + ('-' * (len(reference_translated) - len(translated_padded)))
		assert len(translated_padded) == len(reference_translated)

		# alignment_table.loc[row.seq,f'{library}_aligned'] = (not row.is_unmapped)
		# alignment_table.loc[row.seq,f'{library}_length']  = len(translated)
		#
		# if library == 'alpaca':
		# 	alignment_table.loc[row.seq,f'{library}_inframe'] = len(translated) > 100
		# else:
		# 	alignment_table.loc[row.seq,f'{library}_inframe'] = len(translated) > 70

		alignment_table.loc[row.seq] = translated_padded
		# alignment_table.loc[row.seq,f'{library}_padded'] = translated_padded
	return alignment_table


def get_reference(reference = None, reference_path = None, reference_frame_start = 0):
	if reference is None:
		if reference_path is None:
			raise Exception("Must specify either `reference` or `reference_path`")
		reference = next(Bio.SeqIO.parse(reference_path,"fasta")).seq
	return reference[reference_frame_start:]

def translate_reference(reference = None, reference_path = None, reference_frame_start = 0):
	# read and translate reference sequence
	reference = get_reference(reference, reference_path)
	reference_translated = reference.translate(to_stop=True)
	return reference_translated



def read_feature_table_sparse(directory):
	import scipy.io
	from pathlib import Path
	mtx = scipy.io.mmread(Path(directory) / 'table.mtx')
	samples = pd.read_csv(Path(directory) / 'samples.tsv', header=None)
	features = pd.read_csv(Path(directory) / 'features.tsv', header=None)
	return pd.DataFrame.sparse.from_spmatrix(mtx,
		index=samples.iloc[:,0].values,
		columns=features.iloc[:,0].values)

def align_translated_peptides(seqs, **kwargs):
	from .peptide import pairwise_align_reference_parallel
	return pairwise_align_reference_parallel(seqs, **kwargs)

def trim_align(max_length):
	def do_trim(x):
		if pd.isna(x):
			return x
		else:
			x = x[:max_length]
			if len(x) < max_length: x = x + ('-'*(max_length - len(x)))
			return x
	return do_trim





# alpaca_CDR_positions =
# table_align = extract_CDRs(table_align, 'alpaca', alpaca_CDR_positions)
#
# synthetic_CDR_positions =
#
# table_align = extract_CDRs(table_align, 'synthetic', synthetic_CDR_positions)
