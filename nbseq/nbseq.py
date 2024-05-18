from pathlib import Path
import hashlib

import pandas as pd
# import anndata
# from anndata import AnnData

from .utils import *
from .cdrs import *
from .asvs import get_identifier


from typing import Mapping, TYPE_CHECKING
if TYPE_CHECKING:
	import skbio.tree
	import anndata
	from .viz import ExperimentVisualizer
	from .resynth import Cart

class LibraryConfig(dotdict):
	"""Configuration object for a given phage display library"""

	"""Forward primer sequence"""
	primer_fwd = ""
	
	"""Reverse primer sequence"""
	primer_rev = ""

	"""Path to reference sequence in FASTA format"""
	reference = None

	"""Nucleic acid position (0-based) indicating the first base of the first codon of the reference sequence"""
	reference_frame_start_nt = 0
	
	"""Length of the reference sequence in amino acids; if the reference 
	sequence is longer than (reference_frame_start_nt + (reference_length_aa * 3)), it will be trimmed.
	"""
	reference_length_aa = 0


	"""3' (distal) end of forward read must align to this NA position or later"""
	min_fwd_end = 0

	"""3' (distal) end of reverse read must align to this NA position or earlier"""
	max_rev_start = 0

	"""Features where the aligned amino acid sequence (excluding gap characters)
	are shorter than this length will be dropped"""
	min_aa_length = 69
    
	"""Position of CDR and FR regions within the reference sequence, in amino 
	acid space. In this object, position 0 refers to the amino acid 
	corresponding to  :py:attr:`reference_frame_start_nt`. This object should
	be a dict mapping domain names (e.g. 'CDR1', 'FR2', etc.) to `[start, end]`
	positions, where `start` is inclusive and `end` is exclusive (e.g. 
	half-open) intervals, following the Python convention.

	Example::

		library.CDRs = {
			'FR1':  [0,  21],
			'CDR1': [21, 29],
			'FR2':  [29, 46],
			'CDR2': [46, 54],
			'FR3':  [54, 92],
			'CDR3': [92, 116],
			'FR4':  [116,127],
		}

	"""
	CDRs = {}

	"""Features with CDR or FR regions shorter than this length (in amino acids)
	will be dropped. Dict where keys are domain names (e.g. 'CDR1', 'FR2', etc.
	and should correspond to domains defined in :py:attr:`CDRs`) and values are 
	minimum lengths (in amino acids)"""
	min_CDR_length = {}


class Config(dotdict):
	"""Configuration object for an experiment and the accompanying snakemake workflow.
	"""

	"""Path to folder containing raw input sequences"""
	raw_sequence_dir = ''

	"""Path to scratch directory"""
	scratch = ''

	"""Phage display libraries included in this experiment. Dict where keys are
	library names and values are :py:class:`LibraryConfig` objects."""
	libraries: {}


class Experiment:
	"""Track multiple feature tables in different feature spaces, feature data, databases, etc.
	
	After running the Snakemake pipeline for preprocessing, this class is 
	useful for loading and organizing the results, performing common 
	transformations (e.g. calculating enrichments), and creating visualizations.


	It is recommended to initialize this object using :py:meth:`from_files`.
	"""

	def __init__(self, fts={}, trees={}, phenotypes=None, sql_db=None, mmseqs_dbs={}, enr_models={}, config={}, name=None):
		self._fts = {k.lower(): v for k, v in fts.items()}
		self._trees = {k.lower(): v for k, v in trees.items()}
		self._phenotypes = phenotypes

		if phenotypes is not None:
			_phenotype_names = self._phenotypes.index.values
			for _space in self._fts:		
				for _phenotype in _phenotype_names:
					if _phenotype not in self.fts[_space].obs:
						print(f"Warning: phenotype {_phenotype} not in obs data for feature table {_space}. Adding column of NAs")
						self._fts[_space].obs[_phenotype] = pd.NA

		self._sql_db = sql_db
		self._mmseqs_dbs = {k.lower(): v for k, v in mmseqs_dbs.items()}
		self._config = config

		self._enr_models = enr_models


		self.name = name

	_rfts = None
	_dfs = None
	_rdfs = None
	_enr = None
	_enr_invocations = None
	_enr_models = None
	_references = None

	_viz = None
	_phenotype_matrix = None
	_antigen_matrix = None
	_selection_metadata = None
	_cart = None

	def __repr__(self):
		import textwrap
		out = [f"Experiment('{self.name}') with layers:"]

		for k, ft in self._fts.items():
			out.append(f"- {k:<8}: {ft.shape[0]} samples x {ft.shape[1]} features, database: {self._mmseqs_dbs.get(k, 'None')}")
			# out.append(f"	var: {ft.var.columns.values}")

			out.append(
				textwrap.fill(
					f"var: {ft.var.columns.values}", 
					width=80, initial_indent="  ", subsequent_indent="    ")
			)
			if k in self._trees:
				out.append(f"	tree: {repr(self._trees[k])}")

		# out.append(f"	obs: {ft.obs.columns.values}")
		out.append(
			textwrap.fill(
				f"obs: {ft.obs.columns.values}",
				width=80, initial_indent="  ", subsequent_indent="    ")
		)
		out.append(f"SQL: {self._sql_db}")
		return "\n".join(out)

	@property
	def config(self) -> Config:
		"""configuration object from `config/config.yaml`"""
		return dict(self._config)

	@property
	def viz(self) -> 'ExperimentVisualizer':
		"""Helper to create visualizations from this Experiment"""
		from .viz import ExperimentVisualizer
		if self._viz is None:
			self._viz = ExperimentVisualizer(self)
		return self._viz

	@property
	def cart(self) -> 'Cart':
		"""'Cart' of rVHHs to resynthesize"""
		if self._cart is None:
			from .resynth import Cart
			self._cart = Cart(self)
		return self._cart

	@property
	def fts(self) -> Mapping[str, 'anndata.AnnData']:
		"""Count-based feature tables in each feature space. Keys are names of feature :py:attr:`spaces`."""
		return dotdict(self._fts)

	@property
	def spaces(self):
		"""List of feature space names in order"""
		return list(self._fts.keys())

	@property
	def obs(self):
		"""Sample metadata from the first feature table in :py:attr:`fts`."""
		return self._fts[next(iter(self._fts))].obs

	@property
	def tree(self) -> Mapping[str, 'skbio.tree.TreeNode']:
		"""Phylogenetic trees in each feature space. Keys are feature space names"""
		return dotdict(self._trees)

	@property
	def rfts(self) -> Mapping[str, 'anndata.AnnData']:
		"""Relative abundance feature tables in each feature space. Keys are names of feature :py:attr:`spaces`."""
		if self._rfts is None:
			from .ft import to_relative
			
			self._rfts = {k: to_relative(v) for k,v in self._fts.items()}
		return dotdict(self._rfts)

	@property
	def dfs(self) -> Mapping[str, 'pd.DataFrame']:
		"""Feature tables of counts, converted to :py:func:`.ft.fortify` - fortified :py:`pd.DataFrame`s; keys are space names"""
		from functools import partial
		if self._dfs is None:
			self._dfs = lazydict({
				k: partial(self.fortify, space=k, obs=True)
				for k in self._fts
			})
		return self._dfs

	@property
	def rdfs(self):
		"""Feature tables of relative abundances, converted to :py:func:`.ft.fortify` - fortified :py:`pd.DataFrame`s; keys are space names"""
		from functools import partial
		if self._rdfs is None:
			self._rdfs = lazydict({
				k: partial(self.fortify, space=k, relative=True, obs=True)
				for k in self._fts
			})
		return self._rdfs

	@property
	def enr_models(self):
		"""Models of enrichment probability"""
		return dotdict(self._enr_models)

	@property
	def selection_metadata(self):
		"""Table of metadata with one row per selection, rather than one row per-sample.

		Joins the result of :func:`utils.sample_metadata_to_selection_metadata`  with :func:`utils.summarize_selection_phenotypes` 

		Returns
		-------
		pd.DataFrame
		"""
		if self._selection_metadata is None:
			_md = sample_metadata_to_selection_metadata(self.obs)
			from .utils import summarize_selection_phenotypes
			self._selection_metadata = _md.join(summarize_selection_phenotypes(_md, self.ag_names).rename('antigens'))

		return self._selection_metadata

	def _summarize_columns(self, columns, key_cols):
		import textwrap
		ag_names = set(self.ag_names)
		ph_names = set(self.pheno_names)
		key_cols = set(key_cols)
		
		print(textwrap.fill(
			f"- Antigens: {sorted(ag_names)}",
			subsequent_indent="   ",break_on_hyphens=False))

		print(textwrap.fill(
			f"- Other phenotypes: {sorted(ph_names - ag_names)}", 
			subsequent_indent="   ",break_on_hyphens=False))

		other_cols = [col for col in columns if (col not in ph_names and col not in key_cols)]
		print(textwrap.fill(
			f"- Other columns: {other_cols}", 
			subsequent_indent="   ",break_on_hyphens=False))		

	def summarize():
		pass

	@propery
	def expt_metadata(self):
		"""Print a summary of each sub-experiment within this experiment, showing the number of selections, number of rounds, and which phage libraries were involved"""
		if self._expt_metadata is None:
			self._expt_metadata = sample_metadata_to_expt_metadata(self.obs)
		return self._expt_metadata

	def _summarize_expts(self):
		em = self.expt_metadata

		out = []
		out.append(f"{len(em)} sub-experiment(s):")
		for expt, row in em.iterrows():
			out.append(f"- {expt}: {row['selections']} selections, {row['fewest_rounds_per_selection']} - {row['most_rounds_per_selection']} rounds per selection ({row['rounds']}), using libraries {row['phage_libraries']}")
		return "\n".join(out)

	def summarize_expts(self):
		print(self._summarize_expts())

	def summarize_selections(self):
		"""Print a summary of selection conditions represented by this experiment"""
		key_cols = ['expt', 'phage_library', 'selection',
             'replicate', 'description', 'samples', 'rounds', 'antigens', 'notes']

		display_or_print(self.selection_metadata[key_cols])
		self._summarize_columns(self.selection_metadata.columns, key_cols)

	def summarize_obs(self):
		"""Print a summary of samples represented by this experiment"""
		key_cols = ['expt', 'phage_library', 'selection',
             'replicate', 'description', 'round', 'notes']

		display_or_print(self.obs[key_cols])
		self._summarize_columns(self.obs.columns, key_cols)

	"""Calculate the probability of observing a given ``enrichment`` for some
	starting ``abundance`` in a given feature ``space`` by applying the :func:`enr_model`"""
	def p_enrichment(self, enrichment, abundance, space='cdr3'):
		model = self.enr_models[space.lower()]
		return 1 - model(np.log10(abundance), np.log10(enrichment))


	def join_enrichment(self, space, method, relative=True):
		enr = self.enr(space, method)


	def fortify(self, space, **kwargs):
		from .ft import fortify
		return fortify(self.fts[space], **kwargs)

	def enr(self, space, method=None, update=False, obs=None, *args, **kwargs):
		"""calculate or retrieve memoized enrichment matrix

		Parameters
		----------
		space : str
			one of `self.spaces`
		method : str, optional
			enrichment calculation method accepted by `nbseq.select.enr`; 
			if None, will use any enrichment matrix calculated for this space, otherwise `df`
		*args : list, optional
			Additional arguments will be passed to :func:`.select.enr`
		**kwargs : dict, optional
			Additional arguments will be passed to :func:`.select.enr`

		Returns
		-------
		pd.DataFrame or anndata.Anndata
			enrichment matrix or DataFrame
		"""
		from .select import enr
		from collections import defaultdict

		if self._enr is None:
			self._enr = defaultdict(dict)
			self._enr_invocations = defaultdict(dict)

		if method is None:
			method = 'df'

		if method not in self._enr[space] or update:
			_enr = enr(ft=self.fts[space], method=method, *args, **kwargs)
			# self._enr_invocations[space][method] = tuple([*args, *kwargs.items()])

			# if obs is not None:
			# 	if isinstance(_enr, pd.DatFrame):
			# 		_enr = _enr.join(obs, on='name')
			# 	elif isinstance(_enr, AnnData):
			# 		_enr.obs = _enr.obs.join(obs, on='name')
					
			self._enr[space][method] = _enr

		return self._enr[space][method]

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

	def get_feature_positions(self, library='alpaca', aa=True):
		""" Get start/end position of each configured feature for the given library
		"""
		cdrs = self.config['libraries'][library]['CDRs']

		if aa:
			return cdrs
		else:
			return {k: [v[0] * 3, v[1] * 3] for k,v in cdrs.items()}

	def get_CDR_positions(self, library='alpaca', aa=True):
		""" Get start/end position of each configured CDR for the given library
		"""
		features = self.get_feature_positions(library=library, aa=aa)

		return {k:v for k, v in features.items() if k.startswith('CDR')}
		 

	def get_clonotype_positions(self, library='alpaca', aa=True):
		"""Find the minimum and maximum position of CDR domains"""
		cdrs = self.get_CDR_positions(library, aa)
		start = min(v[0] for k,v in cdrs.items())
		stop  = max(v[1] for k,v in cdrs.items())
		return (start, stop)

	def search_sql(self, **kwargs):
		"""Execute a query against attached SQL database"""
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

	def project(self, features, from_space='cdr3', to_space='aa', ft=True, relative=False, query=None, **kwargs):
		"""
		Project between feature spaces. See :py:func:`asvs.project`
		"""
		from .asvs import project

		if ft:
			if relative:
				ft = self.rfts[to_space.lower()]
			else:
				ft = self.fts[to_space.lower()]
			if query is not None:
				from .ft import query as ft_query
				ft = ft_query(ft, query, axis='obs')
		else:
			ft = None
		return project(features, from_space, to_space,
			ft=ft, db=self._sql_db, 
			**kwargs)

	def query(self, query, axis='sample', space='cdr3', relative=False, **kwargs):
		"""Query a feature table in the given feature space; see :py:func:`ft.query`"""
		from .ft import query as ft_query
		if relative:
			_ft = self.rfts[space]
		else:
			_ft = self.fts[space]
		return ft_query(_ft, query, axis=axis, **kwargs)

	def query_ids(self, query, axis='sample', space='cdr3', **kwargs):
		"""Query a feature table in the given feature space; see :py:func:`ft.query` and return a list of sample or feature IDs"""
		from .ft import query_ids
		return query_ids(self.fts[space], query, axis=axis, **kwargs)

	def find_feature(self, feature=None, mn=None, single=False, space='cdr3'):
		"""locate a feature according to the first few digits of its hash identifier or the first few words of its mnemonicode

		Parameters
		----------
		feature : str, optional
			initial digits of the feature hash; can specify this or `mn`
		mn : str, optional
			initial words of feature hash converted to mnemonicode; can specify this or `feature`
		single : bool, optional
			if True, return only a single feature ID; raise a ValueError if multiple features match; if False, always return a list of matching feature IDs; by default False
		space : str, optional
			which feature space, by default 'cdr3'

		Returns
		-------
		str or List[str]
			matching feature IDs

		Raises
		------
		TypeError
			if neither `feature` nor `mn` are specified
		ValueError
			if `single` = True and more than one feature is found
		"""
		from .viz.utils import mn_to_hash
		ft = self.fts[space.lower()]
		if feature is None and mn is not None:
			feature = mn_to_hash(mn)
		elif feature is None and mn is None:
			raise TypeError("Either `feature` or `mn` must not be None")

		ids = ft[:,ft.var.index.str.startswith(feature)].var_names.values
		if single:
			if len(ids) > 1:
				raise ValueError(f"More than one {space.upper()} feature found with prefix '{feature}': {list(ids)}")
			else:
				return ids[0]
		return ids

	def find_cdr3(self, CDR3ID=None, mn=None, single=False):
		return self.find_feature(feature=CDR3ID, mn=mn, single=single, space='cdr3')

	def find_library(self, feature, space='cdr3'):
		return self.feature_data(feature, fields='library', space=space)

	def feature_data(self, features, fields, space='cdr3'):
		"""get feature data for this sequence"""
		return self.fts[space.lower()].var.loc[features, fields]

	def reference(self, space='cdr3', library='alpaca', aa=False):
		"""use `.config` to locate and read the reference sequence for a given phage display ``library``"""
		if self._references is None:
			self._references = {}
		space = space.lower()	

		if (space, library, aa) in self._references:
			return self._references[(space, library, aa)]

		from .utils import get_reference, translate_reference

		lib_config = self.config['libraries'][library]
		reference_path = lib_config['reference']
		reference_frame_start = lib_config['reference_frame_start_nt']
		reference = get_reference(reference_path = reference_path, reference_frame_start = reference_frame_start)

		if aa:
			if 'suppress_amber' in lib_config:
				suppress_amber = lib_config['lib_config']
			else:
				suppress_amber = True
				
			reference = translate_reference(reference, suppress_amber=suppress_amber)
			if 'reference_length_aa' in lib_config:
				reference_length = lib_config['reference_length_aa']
				reference = reference[:reference_length]

		self._references[(space, library, aa)] = reference
		return self._references[(space, library, aa)]

		


	@staticmethod
	def from_files(directory='.', metadata='config/metadata-phenotypes.csv',
                phenotypes='config/phenotypes.csv',
                config='config/config.yaml',
                sql_db='intermediate/aa/asvs.db',
                verbose=True, strict=False, **kwargs):
		"""Generate an Experiment object from a directory of files produced by the `phage-seq` Snakemake pipelines

		By default, this will read sample metadata as well as feature tables, feature metadata, and phylogenetic trees 
		for the `aa` and `cdr3` feature spaces. You can load data for additional feature spaces by passing additional 
		kwargs, e.g. to load a feature table for space `na`, pass `ft_na='results/tables/na/feature_table.biom'`.
		To skip loading data fro the `aa` or `cdr3` space, pass `ft_aa=None`, `tree_aa=None`, etc. 

		Missing files will note generate an error, but if `verbose=True`, a warning will be printed.

		Parameters
		----------
		directory : str, optional
			root directory, by default '.'
		metadata : str, optional
			path to sample metadata (obs), by default 'config/metadata-phenotypes.csv'
		phenotypes : str, optional
			path to phenotype metadata, by default 'config/phenotypes.csv'
		config : str, optional
			path to configuration file, by default 'config/config.yaml'
		ft_$ : str, optional
			path to feature table in BIOM or AnnData format for space `$`; must be readable by :func:`ft.read`
		fd_$ : str, optional
			path to CSV file feature metadata for space `$`
		tree_$ : str, optional
			path to feature phylogeny for space `$`, in Newick (``.nwk``) format
		mmseqs_db_$ : str, optional
			path to mmseqs2 database for space `$`
		enr_model_$ : str, optional
			path to enrichment model for space `$`, to be read with :func:`.select.load_cond_ecdf`
		sql_db : str, optional
			path to SQLite database, by default 'intermediate/aa/asvs.db'
		verbose : bool, optional
			print paths to files, file size, etc. while reading; by default True
		strict : bool, optional
			True to raise Exception if any indicated files are missing, False to print a warning (if `verbose`) or ignore silently; by default False

		Returns
		-------
		Experiment
			experiment object
		"""

		kwargs = {**dict(
			# fd_cdr3='intermediate/cdr3/features/all/asvs.csv',
			fd_cdr3='results/tables/cdr3/asvs.csv',
			ft_aa='results/tables/aa/feature_table.biom',
			ft_cdr3='results/tables/cdr3/feature_table.biom',
			tree_aa='intermediate/aa/features/top_asvs/alpaca/asvs.nwk',
			tree_cdr3='intermediate/cdr3/features/top_asvs/alpaca/asvs.nwk',
			mmseqs_db_aa='intermediate/aa/features_db/features',
			mmseqs_db_cdr3='intermediate/cdr3/features_db/features',
			enr_model_cdr3='results/tables/cdr3/enrichment/null/ecdf.pickle',
		), **kwargs}


		from .ft import read as ft_read, add_feature_data, query as ft_query, query_ids

		import yaml
		import time
		from datetime import timedelta

		def filesize(size):
			try:
				import humanize
				return humanize.naturalsize(size)
			except ImportError:
				return str(size)

		start_time = time.time()

		realdir = os.path.realpath(directory)
		name = os.path.basename(realdir)

		def warn_or_err(msg):
			if strict:
				raise ValueError(msg)
			else:
				print(f"- Warning: {msg}")

		if verbose:
			print(f"Loading experiment {name} from '{realdir}'...")
			print(f"- Reading metadata from {metadata} ...")

		metadata = read_delim_auto(Path(directory) / metadata).set_index('ID')

		if phenotypes is None:
			print(f"- Warning: no phenotypes table given.")
		else:
			phenotypes_path = Path(directory) / phenotypes
			if not os.path.isfile(phenotypes_path):
				warn_or_err(f"phenotypes table '{phenotypes_path}' does not exist!")
			else:
				if verbose:
					print(f"- Reading phenotypes from {phenotypes_path} ...")
				phenotypes = read_delim_auto(phenotypes_path).set_index('name', drop=True)

		if verbose:
			print(f"- Reading Config from {config} ...")
		with open(config, 'r') as f:
			config = yaml.load(f, yaml.FullLoader)

		sql_db_path = os.path.abspath(Path(directory) / sql_db)
		sql_db = 'sqlite:///' + str(sql_db_path)
		if not os.path.isfile(sql_db_path):
			warn_or_err(f"sqlite database '{sql_db_path}' does not exist")
		else:
			if verbose:
				print(f"- Using SQL database at '{sql_db}'")

		# mmseqs_dbs = {
		# 	k.split('mmseqs_db_',1)[-1]: Path(directory) / v
		# 	for k, v in kwargs.items() if k.startswith('mmseqs_db_')
		# }
		# mmseqs_dbs = {
		# 	'aa': Path(directory) / mmseqs_db_aa,
		# 	'cdr3': Path(directory) / mmseqs_db_cdr3
		# }
		# for k, v in mmseqs_dbs.items():
		# 	if not os.path.isfile(v):
		# 		print(f"- Warning: mmseqs2 database '{k}' at '{v}' does not exist!")
		# 	elif verbose:
		# 		print(f"- Using mmseqs2 database '{k}' at '{v}'")

		# feature_data_cdr3 = None
		# if fd_cdr3 is not None:
		# 	fd_cdr3 = Path(directory) / fd_cdr3
		# 	if not os.path.isfile(fd_cdr3):
		# 		print(f"- Warning: feature data table '{fd_cdr3}' does not exist!")
		# 	else:
		# 		if verbose:
		# 			print(f"- Reading feature data for table 'cdr3' from {fd_cdr3} ...")
		# 		feature_data_cdr3 = read_delim_auto(
		# 			fd_cdr3).drop_duplicates('CDR3ID').set_index('CDR3ID')




		# look for feature data first, since it will be used to construct feature tables
		fd = {}
		for k, v in kwargs.items():
			if k.startswith('fd_'):
				from .asvs import get_identifier
				
				space = k.split('fd_', 1)[-1]

				if v is None:
					if verbose: print(f"- Warning: not loading feature data for space '{space}'; path was None")
					continue

				fd_path = Path(directory) / v
				if not os.path.isfile(fd_path):
					warn_or_err(f"feature data for table '{space}' at '{fd_path}' does not exist!")
				else:
					identifier = get_identifier(space)

					if verbose:
						print(f"- Reading feature data for table '{space}' from {fd_path} ({filesize(fd_path.stat().st_size)})...")
					fd[space] = read_delim_auto(
						fd_path).drop_duplicates(identifier).set_index(identifier)


		trees = {}
		fts = {}
		mmseqs_dbs = {}
		enr_models = {}
		

		for k, v in kwargs.items():
			if v is not None:

				if k.startswith('mmseqs_db_'):
					space = k.split('mmseqs_db_', 1)[-1]
					mmseqs_db_path = Path(directory) / v
					if not os.path.isfile(mmseqs_db_path):
						warn_or_err(f"mmseqs2 database for space '{space}' at '{mmseqs_db_path}' does not exist!")
					else:
						if verbose:
							print(f"- Using mmseqs2 database '{space}' at '{mmseqs_db_path}'")
						mmseqs_dbs[space] = str(mmseqs_db_path)

				elif k.startswith('tree_'):
					space = k.split('tree_', 1)[-1]

					tree_path = Path(directory) / v
					if not os.path.isfile(tree_path):
						warn_or_err(f"phylogeny for space '{space}' at '{tree_path}' does not exist!")
					else:
						import skbio.tree

						if verbose:
							print(
								f"- Reading {space} phylogeny from {tree_path} ({filesize(tree_path.stat().st_size)})...")
						trees[space] = skbio.tree.TreeNode.read(str(tree_path), format="newick")

				# feature tables
				elif k.startswith('ft_'):
					space = k.split('ft_', 1)[-1]

					ft_path = Path(directory) / v
					if not os.path.isfile(ft_path):
						warn_or_err(f"feature table for space '{space}' at '{ft_path}' does not exist!")
					if verbose:
						print(
							f"- Reading {space} feature table from {ft_path} ({filesize(ft_path.stat().st_size)})...")
					fts[space] = ft_read(ft_path,
							to='anndata',
							metadata=metadata,
							feature_data=(fd[space] if space in fd else None)
					)
					fts[space] = add_feature_data(fts[space], 'reads', 'nsamples')

					if fts[space].var.name is None:
						fts[space].var.name = get_identifier(space)

				# enrichment ecdf
				elif k.startswith('enr_model_'):
					space = k.split('enr_model_', 1)[-1]

					ecdf_path = Path(directory) / v
					from .select import load_cond_ecdf

					if not os.path.isfile(ecdf_path):
						warn_or_err(f"enrichment model for space '{space}' at '{ecdf_path}' does not exist!")
					else:
						if verbose:
							print(
								f"- Reading enrichment model (conditional ECDF) for space {space} from {ecdf_path} ({filesize(ecdf_path.stat().st_size)})...")
						enr_models[space] = load_cond_ecdf(ecdf_path)

		# if tree_aa is not None:
		# 	tree_aa = Path(directory) / tree_aa
		# 	if verbose:
		# 		print(
		# 			f"- Reading AA phylogeny from {tree_aa} ({filesize(tree_aa.stat().st_size)})...")
		# 	trees['aa'] = skbio.tree.TreeNode.read(str(tree_aa), format="newick")

		# if tree_cdr3 is not None:
		# 	tree_cdr3 = Path(directory) / tree_cdr3
		# 	if verbose:
		# 		print(
		# 			f"- Reading CDR3 phylogeny from {tree_cdr3} ({filesize(tree_cdr3.stat().st_size)})...")
		# 	trees['cdr3'] = skbio.tree.TreeNode.read(str(tree_cdr3), format="newick")

		# fts = {}
		# if ft_aa is not None:
		# 	ft_aa = Path(directory) / ft_aa
		# 	if verbose:
		# 		print(
		# 			f"- Reading AA feature table from {ft_aa} ({filesize(ft_aa.stat().st_size)})...")
		# 	ft_aa = ft_read(ft_aa,
        #                     to='anndata',
        #                     metadata=metadata
        #            )
		# 	fts['aa'] = ft_aa

		# if ft_cdr3 is not None:
		# 	ft_cdr3 = Path(directory) / ft_cdr3
		# 	if verbose:
		# 		print(
		# 			f"- Reading CDR3 feature table from {ft_cdr3} ({filesize(ft_cdr3.stat().st_size)})...")
		# 	ft_cdr3 = ft_read(ft_cdr3,
        #                     to='anndata',
        #                     metadata=metadata,
        #                     feature_data=feature_data_cdr3
        #              )
		# 	ft_cdr3 = add_feature_data(ft_cdr3, 'reads', 'nsamples')
		# 	fts['cdr3'] = ft_cdr3

		if verbose:
			try:
				import humanize
				print("Finished in " +
	  				humanize.precisedelta(timedelta(seconds=start_time-time.time())))
			except ImportError:
				print(f"Finished at {time.time()}")

		return Experiment(
			fts=fts,
			trees=trees,
			phenotypes=phenotypes,
			sql_db=sql_db,
			mmseqs_dbs=mmseqs_dbs,
			enr_models=enr_models,
			config=config,
			name=name
		)
	




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
