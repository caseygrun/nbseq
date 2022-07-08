"""
Functions to evaluate the process of phage display selection: enrichment, enrichment specificity, etc.
"""
import sys
from pathlib import *
import subprocess
import tempfile

import pandas as pd
import numpy as np

from ..ft import dataframe_to_anndata, fortify, query as ft_query
from ..utils import *

def calculate_enrichment_df(ft, group_by=['name'],
	comparison = ('R2i', 'R5i'), comparison_col='round',
	enrichment_col='enrichment', input_library_query="expt == '027i.lib'", fillna='min'):
	""" Calculate the enrichment of each nanobody in a feature table over rounds of selection.

	Generally, enrichment = relative abundance at last round / relative abundance at first round. Enrichment is calculated for each nanobody, for each selection condition.
	The feature table has one row for each selection condition-round.

	Parameters
	----------
	ft : anndata.AnnData
		feature table
	group_by : list of str, optional
		which column(s) of `ft.obs` should be used to identify sample conditions
	comparison : tuple of str, optional
		which values of ``comparison_col`` should be compared.
	comparison_col : str
		column of `ft.obs` that is used to identify different rounds of selection
	input_library_query : str, optional
		if given, identifies rows of the feature table corresponding to the raw
		(input) library. These rows will be used to fill missing (`NA`) values for
		``comparison_col == 'R1i'``; features with missing abundances for this
		round will be filled with the mean relative abundance in the input library.
	fillna : str or None
		how to fill values for nanobodies where ``comparison_col == comparison[0]`` is NaN
		'min': fills with the minimum relative abundance for that sample
		None : give NaN for the enrichment of such nanobodies
		other: (e.g. 'R3i'): give the value where ``comparison_col == fillna``. For example, if a VHH is not seen in R2i, give the value in R3i

	Returns
	-------
	pd.DataFrame
		Columns are named by values of ``comparison_col``, typically ``R2i``, ``R3i``, etc. ``enrichment_col`` contains the enrichment value.

	"""

	identifier = ft.var_names.name
	df_f1 = None

	# the raw library may have been sequenced as a set of separate samples. If so, it could be considered
	# to contain the "R1 input" abundances for all samples. input_library_query describes which samples
	if input_library_query is not None:
		# if given, identify subset of observations that correspond to the input (un-panned) library.
		# calculate the mean
		df_r1 = (
			fortify(
				ft_query(ft, input_library_query, axis='obs'),
				relative=True)
			.groupby(identifier)['abundance'].mean()
			.rename('R1i').to_frame()
		)

	# fortify to a long DataFrame, one row per sample/nanobody pair
	ft_df = fortify(ft, relative=True)

	# attach minimal metadata about grouping variable and panning round
	# df = pd.merge(ft_df, ft.obs[group_by + ['round']], left_on='ID', right_index=True)
	df = ft_df.join(ft.obs[group_by + [comparison_col]], on='ID')

	# construct table, one row per name/feature, with columns for abundance in R1i, R1o, R2i, R2o, etc.
	dfp = pd.pivot_table(df, index=(group_by + [identifier]), columns=comparison_col,values='abundance')

	# if given, fill (missing) values for R1i with those of the raw library. This allows for some samples
	# to have R1i sequenced separately, if they had a different input library for some reason
	if df_f1 is not None:
		# alternative:
		# dfp.drop(['R1i'], axis='columns').join(df_r1)
		dfp = dfp.combine_first(df_r1)

	# optionally, assume un-observed nanobodies in the comparison round have same relative abundance as
	# the least abundant nanobody that was observed in that sample.
	if fillna == 'min':
		# find the minimum relative abundance for each sample at round == comparison[0]
		# note: must convert .to_frame() to make .combine_first work sensibly below
		c0_min = dfp.reset_index().groupby(group_by)[comparison[0]].min().to_frame()

		# for any nanobodies not observed at round == comparison[0], set their abundances
		# to the minimum relative abundance for that sample
		# important: must keep this a pd.DataFrame (hence dfp[[comparison[0]]] instead of dfp[comparison[0]])
		# for .combine_first to work sensibly, but then grab first column so
		# enrichment calculation by division works
		c0 = dfp[[comparison[0]]].combine_first(c0_min).iloc[:,0]
	elif fillna is None:
		c0 = dfp[comparison[0]]
	else:
		# for any nanobodies not observed at round == comparison[0], set their abundances
		# to round == fillna
		c0 = dfp[[comparison[0]]].combine_first(dfp[[fillna]]).iloc[:,0]

	c1 = dfp[comparison[1]]

	# calculate enrichment
	dfp[enrichment_col] = c1 / c0
	return dfp.reset_index()

def calculate_enrichment_ad(ft, name_col='name', feature_col=None, dropna=True, log=False, **kwargs):
	""" Calculate enrichment matrix as an AnnData table

	See ``calculate_enrichment_df``
	"""
	if feature_col is None:
		feature_col = ft.var_names.name

	obs = ft.obs.groupby('name').head(1).set_index('name')
	var = ft.var
	df = calculate_enrichment_df(ft, group_by=['name'])
	if dropna:
		df = df.dropna(subset=['enrichment'])
	enrichment = dataframe_to_anndata(df, obs_col=name_col, var_col=feature_col, value_col='enrichment', obs=obs, var=var)

	if log:
		enrichment.X.data = np.log10(enrichment.X.data)
	enrichment.X = sparse_drop_na(enrichment.X)

	return enrichment

def summarize_enrichment(enrichment):
	print(f"Calculated enrichment for {enrichment.shape[0]} samples across {enrichment.shape[1]} features; {enrichment.X.nnz} "
		  f"({enrichment.X.nnz/(enrichment.shape[0]*enrichment.shape[1]):.1%}) non-empty enrichment values")

def calculate_specificity(ft, q1, q2):
	from .ft import query as ft_query
	e1 = ft_query(ft, q1, axis='sample')
	e2 = ft_query(ft, q2, axis='sample')

	e1.X = sparse_drop_na(e1.X)
	e2.X = sparse_drop_na(e2.X)

	assert (e1.var_names.values == e2.var_names.values).all()

	df = pd.DataFrame({'enrichment_1': e1.X.mean(axis=0).A1,
					   'enrichment_2': e2.X.mean(axis=0).A1},
					  index=e1.var_names)

	# df = (pd.Series(e1.X.mean(axis=0).A1,
	#                 index=e1.var_names, name='e1').join(
	#       pd.Series(e2.X.mean(axis=0).A1,
	#                 index=e2.var_names, name='e2'))
	#      )
	df['specificity'] = df['e1']/df['e2']
	return df

def calculate_specificities(ft, enrichment, q1, q2, mean='geometric'):
	""" Calculate the enrichment and abundance specificity in two different groups
	"""

	from .ft import query as ft_query

	# if isinstance(partitions, dict):
	#     names = list(partitions.keys())
	#     queries = list(partitions.values())
	# else:
	#     queries = partitions
	#     names = list(range(0, len(queries)))

	a1 = ft_query(ft, q1, axis='sample')
	a2 = ft_query(ft, q2, axis='sample')

	enrichment = enrichment[:,a1.var_names.values]
	e1 = ft_query(enrichment, q1, axis='sample')
	e2 = ft_query(enrichment, q2, axis='sample')

	e1X = sparse_drop_na(e1.X)
	e2X = sparse_drop_na(e2.X)
	a1X = sparse_drop_na(a1.X)
	a2X = sparse_drop_na(a2.X)


	if mean == 'geometric':
		_mean = sparse_gmean_nonzero
		_std = sparse_gstd
	else:
		_mean = sparse_mean_nonzero
		_std = sparse_std

	df = pd.DataFrame({'enrichment_1':           _mean(e1X, axis=0).A1,
					   'enrichment_1_std':        _std(e1X, axis=0).A1,
					   'enrichment_2':           _mean(e2X, axis=0).A1,
					   'enrichment_2_std':        _std(e2X, axis=0).A1,
					   'abundance_1':            _mean(a1X, axis=0).A1,
					   'abundance_1_std':         _std(a1X, axis=0).A1,
					   'abundance_2':            _mean(a2X, axis=0).A1,
					   'abundance_2_std':         _std(a2X, axis=0).A1,
					  },
					  index=a1.var_names)
	df['abundance_specificity']  = df[f'abundance_1']  / df[f'abundance_2']
	df['enrichment_specificity'] = df[f'enrichment_1'] / df[f'enrichment_2']
	return df


def enrichment_rank_products(df):
	pass


def enrichment_lme4(
	df_path=None, ag_col=None, feature_col='CDR3ID', output_path=None,
	conda=None, verbose=True, threads=1, return_output=None):

	with tempfile.TemporaryDirectory() as tmpdir:

		# TODO: accept feature table?
		if df_path is None:
			df_path = Path(tmpdir) / 'df.csv'
			# TODO: save df to path
			raise NotImplementedError()

		if output_path is None:
			tmp_output_path = Path(tmpdir) / 'models.csv'
		else:
			tmp_output_path = output_path
		mkdirp_file(tmp_output_path)

		# call scran
		cmd = ['Rscript', Path(__file__).parent / 'lme4.R',
			str(df_path), str(tmp_output_path), str(feature_col), str(ag_col), str(int(threads))
		]

		if verbose: print("Running R...")
		run_cmd(cmd, conda=conda, verbose=verbose)

		if output_path is None:
			return pd.read_csv(tmp_output_path)
