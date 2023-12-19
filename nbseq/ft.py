from pathlib import Path
import builtins
from functools import singledispatch

import numpy as np
import pandas as pd
import biom
from anndata import AnnData
import anndata

from .utils import map_generic, sparse_var


def read_sparse_dir(directory):
	import scipy.io
	from pathlib import Path
	from scipy.sparse import coo_matrix
	import numpy as np
	import biom

	mtx = scipy.io.mmread(Path(directory) / 'table.mtx')
	samples = pd.read_csv(Path(directory) / 'samples.tsv', header=None).iloc[:,0].values
	observations = pd.read_csv(Path(directory) / 'features.tsv', header=None).iloc[:,0].values
	# return pd.DataFrame.sparse.from_spmatrix(mtx,
	# 	index=samples.iloc[:,0].values,
	# 	columns=features.iloc[:,0].values)

	# # determine a stable order and index positioning for the identifiers
	# obs_order = sorted(set(observations))
	# samp_order = sorted(set(samples))
	# obs_index = {o: i for i, o in enumerate(obs_order)}
	# samp_index = {s: i for i, s in enumerate(samp_order)}
	obs_order = observations
	samp_order = samples

	# # fill the matrix
	# row = np.array([obs_index[obs] for obs in observations], dtype=int)
	# col = np.array([samp_index[samp] for samp in samples], dtype=int)
	# data = np.asarray(data)
	# mat = coo_matrix((data, (row, col)))

	# Parameters
	# ----------
	# data : array_like
	#     An (N,M) sample by observation matrix represented as one of these
	#     types:
	#         An 1-dimensional array of values
	#         An n-dimensional array of values
	#         An empty list
	#         A list of numpy arrays
	#         A list of dict
	#         A list of sparse matrices
	#         A dictionary of values
	#         A list of lists
	#         A sparse matrix of values
	# observation_ids : array_like of str
	#     A (N,) dataset of the observation IDs, where N is the total number
	#     of IDs
	# sample_ids : array_like of str
	#     A (M,) dataset of the sample IDs, where M is the total number of IDs
	return biom.Table(mtx.transpose(), obs_order, samp_order)
read_feature_table_sparse = read_sparse_dir

def read_biom(path):
	import biom
	return biom.load_table(path)
read_feature_table_biom = read_biom

def write_biom(ft, path):
	from biom.util import biom_open
	with biom_open(path, 'w') as f:
		ft.to_hdf5(f, "nbseq")
write_feature_table_biom = write_biom

def write_mtx(ft, path):
	from scipy.io import mmwrite, mmread
	mmwrite(str(path), ft.matrix_data.T)
write_feature_table_mtx = write_mtx


def read(path, to='anndata', format=None, metadata=None, feature_data=None, **kwargs):
	if format is None:
		format = Path(path).suffix
		if len(format) > 0 and format[0] == '.':
			format = format[1:]
	format = format.lower()
	if to is not None:
		to = to.lower()

	if format in ['anndata', 'h5ad', 'csv', 'loom']:
		import anndata

		if format in ('anndata', 'h5ad'):
			ad = anndata.read_h5ad(path, **kwargs)
		elif format == 'csv':
			ad = anndata.read_csv(path, **kwargs)
		elif format == 'loom':
			ad = anndata.read_loom(path, **kwargs)

		if metadata is not None:
			ad.obs = metadata
		if feature_data is not None:
			ad.var = feature_data

		if to is not None and to != 'anndata':
			raise NotImplementedError(f"Don't know how to convert AnnData to {to}")
	elif format == 'biom':
		import biom

		bt = biom.load_table(path)

		if to == 'anndata':
			return to_anndata(bt, metadata=metadata, feature_data=feature_data, **kwargs)
		elif to == 'biom' or to is None:
			return biom
		else:
			raise NotImplementedError(f"Don't know how to convert biom to {to}")
	elif format == '':
		ft = read_sparse_dir(path)

		if to == 'anndata':
			return to_anndata(bt, metadata=metadata, feature_data=feature_data, **kwargs)
		elif to == 'biom' or to is None:
			return biom
		else:
			raise NotImplementedError(f"Don't know how to convert biom to {to}")
	else:
		raise NotImplementedError(f"Don't know how to import feature table of format '{format}'")


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def to_biom(ft):
	import biom

	# in anndata, rows are samples, columns are variables; in biom, it's the transpose
	return biom.Table(data=ft.X.T, observation_ids = ft.var_names.values, sample_ids = ft.obs_names.values)


def to_anndata(ft, metadata=None, feature_data=None, obs=None, var=None, **kwargs):
	"""
	Convert a feature table in biom.Table format to anndata.AnnData format
	"""
	from anndata import AnnData

	reorder = (obs is not None) or (metadata is not None)
	if obs is None:
		if metadata is None:
			try:
				obs = ft.metadata_to_dataframe('sample')
			except KeyError:
				obs = pd.DataFrame(index=ft.ids('sample'))
		else:
			obs = metadata

	# ensure passed obs DataFrame is in the same order as  matrix_data in
	# the biom feature table
	if reorder:
		obs = obs.loc[ft.ids('sample'),:]


	reorder = (var is not None) or (feature_data is not None)
	if var is None:
		if feature_data is None:
			try:
				var = ft.metadata_to_dataframe('observation')
			except KeyError:
				var = pd.DataFrame(index=ft.ids('observation'))
			reorder = False
		else:
			var = feature_data

	# ensure passed var DataFrame is in the same order as  matrix_data in
	# the biom feature table
	if reorder:
		var = var.loc[ft.ids('observation'),:]

	return AnnData(ft.matrix_data.T,
			obs=obs,
			var=var, **kwargs)


def dataframe_to_anndata(df, obs_col, var_col, value_col, obs=None, var=None, layers=[]):
	""" Convert a pd.DataFrame in long format (one row per observation/variable pair) to an AnnData feature table.

	Note: use :py:`anndata.AnnData` directly if df is in wide format (one row per observation, one column per variable)
	"""
	from scipy.sparse import coo_matrix
	from anndata import AnnData
	# TODO: what if a desired obs or var name does not appear in the df? also, should order of obs and/or var be preserved

	i_s, obs_names = pd.factorize(df[obs_col])
	j_s, var_names = pd.factorize(df[var_col])
	values = df[value_col]

	if layers is None:
		layers = [x for x in list(df.columns) if x != value_col]
	layer_matrices = { layer: coo_matrix((df[layer], (i_s, j_s))).tocsr() for layer in layers }

	if obs is not None:
		obs = obs.loc[obs_names.values,:]
	else:
		obs = pd.DataFrame(index=pd.Index(obs_names.values, name=obs_col))
	if var is not None:
		var = var.loc[var_names.values,:]
	else:
		var = pd.DataFrame(index=pd.Index(var_names.values, name=var_col))
	return AnnData(X = coo_matrix((values, (i_s, j_s))).tocsr(),
				   obs = obs, var = var , layers = layer_matrices)

def to_dmatrix(ft, label, **kwargs):
	from xgboost import DMatrix
	return DMatrix(ft.X, label=label, feature_names=ft.var_names.values)


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# TODO: replace with https://docs.python.org/3/library/functools.html#functools.singledispatch
def dispatch(obj, methods, *args, **kwargs):
	"""Simple multiple dispatch. Use type of ``obj`` to call one of several ``methods``
	"""
	for cls, method in methods.items():
		if isinstance(obj,cls):
			return method(obj, *args, **kwargs)
	raise Exception(f"Could not dispatch suitable method for {repr(obj)} (class {type(obj)}) among methods for classes {list(methods.keys())}.")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# def ft_sum_ad(ft_ad, axis='sample'):
# 	axis = get_axis(axis, numeric=True)
#
# 	# ids = ft_ad.obs_names if axis else ft_ad.var_names
# 	X_sum = ft_ad.X.sum(axis=axis).A1
# 	# return pd.Series(X_sum, index = ids)
# 	return X_sum
#
# ft_sum = ft_sum_ad
#
# def ft_sum_biom(ft, axis='sample'):
# 	return ft.sum(axis=axis)
#
# def sum(ft, axis='sample', **kwargs):
# 	return dispatch(ft, {
# 		AnnData: ft_sum_ad,
# 		biom.Table: ft_sum_biom
# 	}, axis=axis, **kwargs)


@singledispatch
def sum(ft, axis='sample', **kwargs):
	raise NotImplementedError()

@sum.register(AnnData)
def _ft_sum_ad(ft_ad, axis='sample'):
	from scipy.sparse import isspmatrix

	axis = get_axis(axis, numeric=True)

	# ids = ft_ad.obs_names if axis else ft_ad.var_names
	if isspmatrix(ft_ad.X):
		X_sum = ft_ad.X.sum(axis=axis).A1
	else:
		X_sum = ft_ad.X.sum(axis=axis)
	
	# return pd.Series(X_sum, index = ids)
	return pd.Series(X_sum, index=get_ids_ad(ft_ad, axis=(1 - axis)))


@sum.register(biom.Table)
def _ft_sum_biom(ft, axis='sample'):
	return ft.sum(axis=axis)

ft_sum = sum

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def get_identifier(ft, axis='var'):
	if axis != 'var':
		raise NotImplementedError()
	identifier = ft.var_names.name
	if identifier is None:
		return 'feature'
	else:
		return identifier

def get_ids_ad(ft, axis):
	return (ft.obs_names if get_axis(axis, True) else ft.var_names).values

def get_ids_biom(ft, axis):
	axis = get_axis_biom(axis)
	return ft.ids(axis=axis)


def get_ids(ft, axis):
	"""
	Returns the identifiers of a feature table across the appropriate axis

	For AnnData, this is var_names.values or obs_names.values; for biom.Table, this is .ids(axis)
	"""
	return dispatch(ft, {
		AnnData: get_ids_ad,
		biom.Table: get_ids_biom
	}, axis=axis)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def group(ft, mapping, axis, **kwargs):
	return dispatch(ft, {
		AnnData: group_ft_ad,
		biom.Table: group_feature_table
	}, mapping=mapping, axis=axis, **kwargs)

def group_ft_ad(ft, mapping, axis='var'):
	from scipy.sparse import coo_matrix

	# approach: convert matrix data back to coo, e.g. triplet format;
	# map triplets in the appropriate axis back according to mapping
	# fix shape, etc. and return new table

	coo = ft.X.tocoo(copy=True)
	orig_ids = get_ids_ad(ft, axis)

	axis = get_axis(axis)

	# identify numerical codes from the appropriate axis
	if axis=='obs': codes = coo.row
	else: codes = coo.col

	# map numerical sparse matrix codes (either rows or columns) to ids
	# across relevant axis (dim)
	dim_ids = orig_ids[codes]

	# map the old IDs to the new IDs
	new_ids = map_generic(mapping, dim_ids)
	assert len(new_ids) == len(dim_ids), f"number of new IDs {len(new_ids)} should be equal to number of old IDs {len(dim_ids)}"

	# convert an array of new id strings (new_ids) to an array of numerical
	# codes (new_codes); new_codes represent indices in uniq_new_ids, such that
	# uniq_new_ids[new_codes] = new_ids
	new_codes, uniq_new_ids = pd.factorize(new_ids)
	assert len(new_codes) == len(dim_ids), "length of the new"


	# overwrite the appropriate axis
	if axis=='obs':
		var_ids = uniq_new_ids
		obs_ids = ft.obs_names.values
		(i, j) = (new_codes, coo.col)
	else:
		var_ids = ft.var_names.values
		obs_ids = uniq_new_ids
		(i, j) = (coo.row, new_codes)

	assert len(coo.data) == len(i) == len(j), f"COO tuple columns should be equal length: len(data) = {len(coo.data)}, len(i) = {len(i)}, len(j) = {len(j)}"

	return AnnData(
		coo_matrix(
			(coo.data, (i, j)), shape=(len(obs_ids),len(var_ids))
		).to_csr(),
		obs=pd.DataFrame(index=obs_ids),
		var=pd.DataFrame(index=var_ids)
	)



def group_ft_biom(ft, mapping, axis='observation'):
	"""groups samples or observations of a feature table according to some mapping.

	much faster than using `Table.collapse`, but only handles the case where
	"""

	from scipy.sparse import coo_matrix

	# approach: convert matrix data back to coo, e.g. triplet format;
	# map triplets in the appropriate axis back according to mapping
	# fix shape, etc. and return new table

	coo = ft.matrix_data.tocoo(copy=True)
	orig_ids = ft.ids(axis)

	# identify numerical codes from the appropriate axis
	if axis=='observation': codes = coo.row
	else: codes = coo.col

	# map numerical sparse matrix codes (either rows or columns) to ids
	# across relevant axis (dim)
	dim_ids = orig_ids[codes]

	# map the old IDs to the new IDs
	new_ids = map_generic(mapping, dim_ids)
	assert len(new_ids) == len(dim_ids), f"number of new IDs {len(new_ids)} should be equal to number of old IDs {len(dim_ids)}"

	# convert an array of new id strings (new_ids) to an array of numerical
	# codes (new_codes); new_codes represent indices in uniq_new_ids, such that
	# uniq_new_ids[new_codes] = new_ids
	new_codes, uniq_new_ids = pd.factorize(new_ids)
	assert len(new_codes) == len(dim_ids), "length of the new"

	# u, c = np.unique(uniq_new_ids, return_counts=True)
	# dup = u[c > 1]
	# print(dup)

	# overwrite the appropriate axis
	if axis=='observation':
		obs_ids = uniq_new_ids
		sam_ids = ft.ids('sample')
		(i, j) = (new_codes, coo.col)

		# other_axis_ids = ft.ids('sample')
		# shape = (len(uniq_new_ids), len(other_axis_ids))
		# new_coo = coo_matrix((coo.data, (new_codes, coo.col)), shape=shape)
		# return Table(new_coo, uniq_new_ids, other_axis_ids)
	else:
		obs_ids = ft.ids('observation')
		sam_ids = uniq_new_ids
		(i, j) = (coo.row, new_codes)

		# other_axis_ids = ft.ids('observation')
		# shape = (len(other_axis_ids), len(uniq_new_ids))
		# new_coo = coo_matrix((coo.data, (coo.row, new_codes)), shape=shape)
		# return Table(new_coo, other_axis_ids, uniq_new_ids)

	assert len(coo.data) == len(i) == len(j), f"COO tuple columns should be equal length: len(data) = {len(coo.data)}, len(i) = {len(i)}, len(j) = {len(j)}"

	return biom.Table(
		coo_matrix((coo.data, (i, j)), shape=(len(obs_ids), len(sam_ids))),
		obs_ids, sam_ids)

group_feature_table = group_ft_biom

def group_feature_table_slow(ft, mapping, aggregate='sum', axis='observation'):
	"""groups a feature table, using collapse method from biom package.
	extremely slow due to lots of copying to handle general case.
	"""

	mapping_axis = {
		0: 'sample',
		1: 'observation',
		'sample':'sample',
		'samples':'sample',
		'rows':'sample',
		'observation':'observation',
		'observations':'observation',
		'features':'observation',
		'columns':'observation',
	}
	axis = mapping_axis[axis]

	if callable(mapping):
		mapper = mapping
	elif isinstance(mapping, dict):
		mapper = lambda x, _md: mapping[x]
	elif isinstance(mapping, pd.Series):
		assert not any(mapping.index.duplicated()), "Series used for mapping must not have duplicate values in index! Each old ID must map to a single new new ID"
		mapper = lambda x, _md: mapping[x]

	return ft.collapse(
		mapper,
		axis=axis,
		include_collapsed_metadata=False,
		one_to_many=False)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def transform(ft, f, **kwargs):
	return dispatch(ft, {
		AnnData: transform_ad,
		biom.Table: transform_biom
	}, f=f, **kwargs)

def transform_biom(ft, f=np.log1p):
	# coo = ft.matrix_data.tocoo(copy=True)
	# coo.data = f(coo.data)

	matrix = ft.matrix_data.copy()
	matrix.data = f(matrix.data)

	return biom.Table(
		matrix,
		ft.ids('observation'), ft.ids('sample'))

def transform_ad(ft, f=np.log1p):
	ft = ft.copy()
	ft.X.data = f(ft.X.data)
	return ft

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def ft_filter(ft, ids_to_keep, axis):
	return dispatch(ft, {
		AnnData: lambda ft: ft[ids_to_keep,:] if get_axis(axis, numeric=True) else ft[:,ids_to_keep],
		biom.Table: lambda ft: ft.filter(ids_to_keep, axis=axis)
	})

filter = ft_filter

def filter_feature_table(ft, ids, axis='observation'):
	return ft.filter(ids,axis='observation')

def filter_ids_min_abundance_prevalence(ft, min_abundance=2, min_prevalence=1, axis='feature'):
	if isinstance(ft,biom.Table):
		axis = get_axis_biom(axis, numeric=True)
		n_items_with_min_abundance = (ft.matrix_data > min_abundance).sum(axis=(not axis))
		indices_to_keep = np.asarray(n_items_with_min_abundance > min_prevalence).squeeze()
		ids_to_keep = get_ids(ft, axis)[indices_to_keep]

	elif isinstance(ft,AnnData):
		axis = get_axis(axis, numeric=True)
		n_items_with_min_abundance = (ft.X > min_abundance).sum(axis=axis)

		# .sum above may return a numpy.matrix; make sure it is an ndarray, then
		# reduce to 1 dimension with squeeze
		indices_to_keep = np.asarray(n_items_with_min_abundance > min_prevalence).squeeze()
		ids_to_keep = get_ids(ft, not axis)[indices_to_keep]

	else:
		raise ValueError(f"ft must be biom.Table or AnnData: do not know how to filter {repr(ft)}")
	
	return ids_to_keep

def filter_abundance_prevalence(ft, min_abundance=2, min_prevalence=1, axis='feature'):
	ids_to_keep = filter_ids_min_abundance_prevalence(ft, min_abundance, min_prevalence, axis=axis)
	return filter(ft, ids_to_keep,axis=axis)

filter_ft_abundance_prevalence = filter_abundance_prevalence


def summarize_filtering(ft_in, ft_out, min_abundance=None, min_prevalence=None):
	print("Retained:")
	print(f"{ft_out.shape[1]:>12.0f} = {ft_out.shape[1] / ft_in.shape[1]:>4.0%} features")
	print(f"{ft_out.X.sum():>12.0e} = {ft_out.X.sum() / ft_in.X.sum():>4.0%} reads")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def summarize_biom(ft, axis='sample', id_name=None, summary_name='reads', method='sum'):
	if id_name is None:
		id_name = {
			'sample': 'sample',
			'observation': 'ASVID'
		}[axis]

	if method == 'sum':
		summary = ft.sum(axis=axis)
	elif method == 'count':
		summary = ft.nonzero_counts(axis=axis)
	else:
		raise Exception(f'Unrecognized feature table summary method "{method}"')

	return pd.DataFrame({
		id_name: ft.ids(axis),
		summary_name: summary
	})

def nonzero_counts_ad(ft, axis):
	axis = get_axis(axis, numeric=True)
	# nonzero = [nonzero_row_index, nonzero_col_index], where len(nonzero[0]) == len(nonzero[1]) == nnz
	nonzero = ft.X.nonzero()
	unique, counts = np.unique(nonzero[axis], return_counts=True)

	summary = pd.Series(np.zeros(ft.X.shape[axis]), index = get_ids_ad(ft, axis))
	summary.iloc[unique] = counts

	return summary

def summarize_ad(ft, axis='obs', id_name=None, summary_name='reads', method='sum'):
	if id_name is None:
		if axis == 'obs':
			if ft.obs_names.name:
				id_name = ft.obs_names.name
			else:
				id_name = 'sample'
		elif axis == 'var':
			if ft.var_names.name:
				id_name = ft.var_names.name
			else:
				id_name = 'ASVID'

	if method == 'sum':
		summary = _ft_sum_ad(ft, axis=axis)
	elif method == 'count':
		summary = nonzero_counts_ad(ft, axis=axis)
	else:
		raise Exception(f'Unrecognized feature table summary method "{method}"')

	# return pd.DataFrame({
	# 	id_name: get_ids(ft, axis),
	# 	summary_name: summary
	# })
	return summary

def summarize(ft, **kwargs):
	return dispatch(ft, {
		AnnData: summarize_ad,
		biom.Table: summarize_biom
	}, **kwargs)

def filter_intersect_biom(ft, ids, axis='sample', remove_empty=True, inplace=False):
	"""filters a feature table to include ids that are both in the feature table and in the provided list of `ids`
	returns (ft, ids_to_keep)
	"""
	ids_to_keep = set(ft.ids(axis)).intersection(ids)
	other_axis = {'sample':'observation','observation':'sample'}[axis]
	if inplace:
		ft.filter(ids_to_keep, axis=axis, inplace=True)
		if remove_empty:
			ft.remove_empty(axis=other_axis,inplace=True)
	else:
		ft = ft.filter(ids_to_keep, axis=axis, inplace=False)
		if remove_empty:
			ft = ft.remove_empty(axis=other_axis, inplace=False)
	return ft, ids_to_keep


def remove_empty_ad(ft, axis, inplace=False):
	if inplace:
		raise NotImplementedError()

	axis = get_axis(axis, True)

	# remove empty samples = sum over features, remove any samples where sum == 0
	# remove empty features = sum over samples, remove any features where sum == 0
	keep = (ft.X.sum(axis=int(not axis)) > 0).A1
	if axis == 0:
		return ft[keep,:]
	else:
		return ft[:,keep]

def remove_empty(ft, axis, **kwargs):
	return dispatch(ft, {
		biom.Table: lambda ft: ft.remove_empty(axis=axis, **kwargs),
		AnnData: remove_empty_ad
	},axis, **kwargs)
ft_remove_empty = remove_empty

def filter_intersect_ad(ft, ids, axis='sample', remove_empty=True, inplace=False):
	"""filters a feature table to include ids that are both in the feature table and in the provided list of `ids`
	returns (ft, ids_to_keep)
	"""
	axis = get_axis(axis)

	# much faster to construct the set from the shorter list of IDs, esp. if
	# the longer list is very long
	# ids_to_keep = list(set(get_ids(ft, axis)).intersection(ids))
	_ids = get_ids(ft, axis)
	if len(_ids) < len(ids):
		ids_to_keep = list(set(_ids).intersection(ids))
	else:
		ids_to_keep = list(set(ids).intersection(_ids))

	other_axis = {'var':'obs','obs':'var'}[axis]

	ft = filter(ft, ids_to_keep, axis=axis)
	if remove_empty:
		ft = ft_remove_empty(ft, axis=other_axis)
	return ft, ids_to_keep


def filter_intersect(ft, ids, axis='sample', **kwargs):
	return dispatch(ft, {
		AnnData: filter_intersect_ad,
		biom.Table: filter_intersect_biom
	}, ids=ids, axis=axis, **kwargs)


def fortify_features_biom(ft, seqs=None, abundance_col='abundance'):
	"""transform a feature table to a DataFrame, one row per ASV, with sequence
	and abundance data

	seqs : pd.Series
		indexed by observation IDs
	"""

	if seqs is None:
		return pd.Series(
			ft.sum(axis='observation'),
			index=ft.ids('observation'),
			name=abundance_col
		)
	else:

		# ids = ft.ids('observation')
		ft, ids = filter_intersect(ft, seqs.index.values, axis='observation', inplace=False)
		counts = ft.sum(axis='observation')

		df = pd.concat([
			# works whether `seqs` is a pd.Series or pd.DataFrame
			seqs.loc[ids],
			pd.Series(counts, index=ids, name=abundance_col)
		],axis=1)
		return df.index.rename(seqs.index.name)

def fortify_features_ad(ft, seqs=None, abundance_col='abundance'):
	if seqs is not None:
		ft, ids = filter_intersect(ft, seqs.index.values, axis='var', inplace=False)
	else:
		ids = ft.var_names.values

	counts = pd.Series(sum(ft, axis='var'), index=ids, name=abundance_col)

	if len(ft.var.columns) == 0 and seqs is None:
		df = counts.to_frame()


	concat_columns = [counts, ft.var]
	if seqs is not None:
		concat_columns.append(seqs)

	df = pd.concat(concat_columns, axis=1)
	name = ft.var_names.name
	if seqs is not None:
		name = seqs.index.name or name

	df.index.rename(name, inplace=True)
	return df

def fortify_features(ft, seqs=None, abundance_col='abundance', **kwargs):
	return dispatch(ft, {
		AnnData: fortify_features_ad,
		biom.Table: fortify_features_biom
	}, seqs=seqs, abundance_col=abundance_col, **kwargs)



@singledispatch
def fortify(ft, sample_col='sample', feature_col='feature', abundance_col='abundance', relative=False, **kwargs):
	"""transform a feature table to a tidy dataframe, one row per abundance observation

	Returned DataFrame has three columns: {sample_col}, {feature_col}, and
	{abundance_col}, containing the sample identifier, the feature identifier,
	and the feature abundance (e.g. # of reads or relative abundance)
	"""
	raise NotImplementedError()

@fortify.register(biom.Table)
def fortify_biom(ft, sample_col='sample', feature_col='feature', abundance_col='abundance', relative=False):
	if relative:
		ft = ft.norm(axis='sample', inplace=False)

	coo = ft.matrix_data.tocoo(copy=True)
	sample_ids = ft.ids('sample')
	obs_ids = ft.ids('observation')

	return pd.DataFrame({
		sample_col: sample_ids[coo.col],
		feature_col: obs_ids[coo.row],
		abundance_col: coo.data
	})

@fortify.register(AnnData)
def fortify_ad(ft, sample_col=None, feature_col=None, abundance_col='abundance', relative=False, obs=False, var=False, index=None):
	if relative:
		ft = to_relative_abundance_ad(ft)

	coo = ft.X.tocoo()
	sample_ids = ft.obs_names.values
	obs_ids = ft.var_names.values

	if sample_col is None:
		sample_col = ft.obs_names.name
		if sample_col is None:
			sample_col = 'sample'
	if feature_col is None:
		feature_col = get_identifier(ft, axis='var')

	df = pd.DataFrame({
		sample_col: sample_ids[coo.row],
		feature_col: obs_ids[coo.col],
		abundance_col: coo.data
	})

	if obs:
		#df = pd.merge(left=df, left_on=sample_col, right=ft.obs, right_index=True, how='left')
		df = df.join(ft.obs, on=sample_col)
	if var:
		# df = pd.merge(left=df, left_on=feature_col, right=ft.var, right_index=True, how='left')
		df = df.join(ft.var, on=feature_col)

	# if index is not None:
	# 	if isinstance(index, list):
	# 		for i in index:
	# 			if i == 'var':
	# 				i = ft.obs_names.name
	return df

fortify_anndata = fortify_ad



def collapse_top_asvs(ft_ad, samples, n=20, top_from_samples=None, relative=False, other_features=None):

	if top_from_samples is None:
		top_from_samples = samples

	ft2 = ft_ad[top_from_samples,:]

	asv_counts = ft_sum(ft2, axis='var').sort_values(ascending=False)
	top_asvs = asv_counts.iloc[0:n].index.values

	if relative:
		ft_ad = to_relative(ft_ad)

	if other_features is not None:
		top_asvs = np.union1d(top_asvs, other_features)

	ft_samples = ft_ad[samples,:]
	ft_top = ft_samples[:,top_asvs]

	# sum top ASVs for each sample
	ft_top_sum = ft_top.X.sum(axis=1)

	# sum all ASVs for each sample
	ft_samples_sum = ft_samples.X.sum(axis=1)

	# calculate sum of all non-top ASVs for each sample
	# interestingly, much faster to get sum of non-top-ASVs by summing everything,
	# then summing top ASVs and subtracting, as opposed to summing non-top-ASVs...
	# the below options are ~10x slower
	#     other_sum = fts[:,other_asvs.index.values].X.sum(axis=1)
	#     other_sum = ft_sum(fts[:,other_asvs.index.values])
	other_asv_sum = ft_samples_sum - ft_top_sum

	ft_other = AnnData(
	   other_asv_sum,
	   obs=pd.DataFrame(index=ft_samples.obs_names),
	   var=pd.DataFrame(index=['others'])
	)
	return anndata.concat([ft_top, ft_other], axis=1, merge="first")

def fortify_top_asvs(ft, query, n=30, select_from_round=4,other_features=None):

	samples = query_ids(ft, f"({query}) & kind =='+' & io == 'i'", axis='obs')
	if select_from_round is not None:
		top_from_samples = query_ids(ft, f"({query}) & r == {select_from_round} & kind =='+' & io == 'i'", axis='obs')
	else:
		top_from_samples = None
	
	ft_top = collapse_top_asvs(ft, samples, top_from_samples= top_from_samples, n=n, other_features=other_features)
	df = fortify(ft_top, obs=True, relative=True)
	return df

@singledispatch
def drop_empty(ft, axis='samples'):
	pass


@drop_empty.register(biom.Table)
def drop_empty_biom(ft, axis='samples', inplace=False):
	ids_to_keep = ft.ids(axis)[ft.sum(axis=axis) > 0]
	if inplace:
		ft.filter(ids_to_keep, axis=axis, inplace=True)
		return ft
	else:
		return ft.filter(ids_to_keep, axis=axis, inplace=False)

@drop_empty.register(AnnData)
def drop_empty_ad(ft, axis='obs'):
	axis = get_axis(axis)
	ids_to_keep = get_ids_ad(ft,axis)[ft_sum(ft,axis=axis) > 0]
	return filter(ft,ids_to_keep, axis=axis)


def drop_zero_var(ft, axis='obs'):
	"""return feature table where zero-variance observations or variables are dropped

	Parameters
	----------
	ft : anndata.AnnData
		feature table
	axis : str, optional
		drop entries from which axis; if None, drop zero-variance variables, then zero-variance observations; by default 'obs'

	Returns
	-------
	anndata.AnnData
		filtered feature table
	"""
	if axis is None:
		obs_variance = sparse_var(ft.X, axis=0)
		var_variance = sparse_var(ft.X, axis=1)
		return ft[obs_variance > 0, var_variance > 0]
	else:
		axis = get_axis(axis, numeric=True)
		axis_variance = sparse_var(ft.X, axis=axis)
		return slice_ft(ft, axis_variance > 0, axis=axis)
	

def slice_ft(ft, _slice, axis):
	axis = get_axis(axis)
	if axis == 'obs':
		return ft[_slice, :]
	elif axis == 'var':
		return ft[:, _slice]


def get_ids_where_nz(ft, axis):
	axis = get_axis(axis)
	return get_ids(ft,axis)[ft_sum(ft,axis=axis) > 0]


def get_samples_where_nonzero(ft, features):
	ft = ft.filter(features, axis='observation', inplace=False)
	return ft.ids('observation')[ft.sum(axis='observation') > 0]



# def fortify_ft_anndata(ft, metadata=None, relative=False):
# 	if relative:
# 		ft = ft_to_relative_abundance(ft)
# 	df = (ft.to_df()
# 		 .reset_index().rename(columns={'index':'ID'})
# 		 .melt(id_vars='ID',var_name='feature',value_name='abundance'))
#
# 	if metadata is not None:
# 		df = pd.merge(left=df, right=metadata, on='ID')
# 	return df


def to_relative(ft):
	return dispatch(ft, {
		AnnData: to_relative_abundance_ad,
		biom.Table: lambda ft: ft.norm(axis='sample', inplace=False)
	})

def to_relative_abundance_ad(ft):
	ft = ft.copy()

	# use .multiply instead of *, as * triggers matrix multiplication
	# for sparse arrays https://github.com/scipy/scipy/issues/5881
	# .multiply weirdly sometimes returns COO which is not understood by
	# anndata, so coerce to csr
	ft.X = ft.X.multiply(1/ft.X.sum(axis=1)).tocsr()

#     Y = (1/ft.X.sum(axis=1))
#     if not isinstance(ft.X, scipy.sparse.csr_matrix):
#         raise ValueError('Matrix must be CSR.')
#     Z = ft.X.copy()
#     # simply repeat each value in Y by the number of nnz elements in each row:
#     Z.data *= Y.repeat(np.diff(Z.indptr))

#     ft.X = Z
	return ft

def get_axis(name, numeric=False):
	name = {
		'obs':'obs',
		'var':'var',
		'sample':'obs',
		'feature':'var',
		0:'obs',
		1:'var',
	}[name]

	if numeric:
		return {
			'obs':1,
			'var':0
		}[name]
	else:
		return name

def get_axis_biom(name, numeric=False):
	name = {
		'observation':'observation',
		'sample':'sample',
		'feature':'observation',
		# biom uses rows = features, columns = samples 
		1: 'sample',
		0: 'observation'
	}[name]

	if numeric:
		return {
			'sample':1,
			'observation':0
		}[name]
	else:
		return name



def add_feature_data(ft, *cols, **kwargs):
	for col in cols:
		if col == 'reads':
			ft.var['reads'] = ft.X.sum(axis=0).A1
		elif col == 'nsamples':
			ft.var['nsamples'] = (ft.X > 0).sum(axis=0).A1
		else:
			raise Exception(f"Do not know how to calculate column {col}. Specify column as a kwarg, column=expr, where expr is a string expression for pd.eval or callable.")

	for col, formula in kwargs:
		if isinstance(formula, str):
			ft.var[col] = ft.var.eval(formula)
		elif callable(formula):
			ft.var[col] = formula(ft)

	return ft

def add_metadata(ft, *cols, **kwargs):
	for col in cols:
		if col == 'reads':
			ft.obs['reads'] = ft.X.sum(axis=1).A1
		elif col == 'nfeatures':
			ft.obs['nfeatures'] = (ft.X > 0).sum(axis=1).A1
		else:
			raise Exception(f"Do not know how to calculate column {col}. Specify column as a kwarg, column=expr, where expr is a string expression for pd.eval or callable.")

	for col, formula in kwargs:
		if isinstance(formula, str):
			ft.obs[col] = ft.obs.eval(formula)
		elif callable(formula):
			ft.obs[col] = formula(ft)

	return ft

def join_feature_data(ft, data, **kwargs):
	return join_metadata(ft, data, axis='var', **kwargs)

def join_metadata(ft, data, axis, index=True, how='left', **kwargs):
	axis = get_axis(axis)
	if axis == 'var':
		ft.var = pd.merge(left=ft.var, right=data, **{'left_index':index, 'right_index':index, 'how':'left', **kwargs})
		if ft.var.index.name is None and data.index.name is not None:
			ft.var.index.name = data.index.name
	elif axis == 'obs':
		ft.obs = pd.merge(left=ft.obs, right=data, **{'left_index':index, 'right_index':index, 'how':'left', **kwargs})
		if ft.obs.index.name is None and data.index.name is not None:
			ft.obs.index.name = data.index.name
	return ft


@singledispatch
def query_ids(ft, query, axis='sample', **kwargs):
	raise NotImplementedError()

@query_ids.register(AnnData)
def _query_ids_ad(ad, query, axis='sample', **kwargs):
	axis = get_axis(axis)
	if axis == 'obs':
		return ad.obs.query(query).index.values
	elif axis == 'var':
		return ad.var.query(query).index.values

def query_ad(ad, query, axis='obs'):
	axis = get_axis(axis)
	if axis == 'obs':
		out = ad[_query_ids_ad(ad, query, axis='obs'),:]
		out.var_names.name = ad.var_names.name
	elif axis == 'var':
		out = ad[:,_query_ids_ad(ad, query, axis='var')]
		out.var_names.name = ad.var_names.name
	return out

def query(ft, query, axis, **kwargs):
	return dispatch(ft, {
		AnnData: query_ad
	}, query, axis, **kwargs)






def fortify_feature_data(ft=None, fd=None, relative=False, features=None, var=False, obs=False, sample_col='ID', **kwargs):
	"""Given either a feature table or feature DataFrame, get a feature DataFrame
	"""
	if fd is not None:
		if features is not None:
			fd = fd.loc[fd['feature'].isin(set(features)),:]
	elif ft is not None:
		if relative:
			ft = nbseq.ft.to_relative(ft)
		if features is not None:
			ft = ft[:,features]
		fd = nbseq.ft.fortify(ft, sample_col=sample_col, obs=obs, var=var, **kwargs)
	else:
		raise ValueError("Must specify either ft or fd")

	return fd


# from anndata import AnnData
# class FeatureTable(AnnData):
# 	def sum(self, axis='obs'):
# 		axis = get_axis(axis, numeric=True)
# 		ids = self.obs_names if axis else self.var_names
# 		X_sum = self.X.sum(axis=axis).A1
# 		return pd.Series(X_sum, index = ids)
#
# 	def nz_obs(self):
# 		X_sum = self.X.sum(axis=1)
# 		return self[X_sum > 0,:]
#
# 	def nz_var(self):
# 		X_sum = self.X.sum(axis=0)
# 		return self[:,X_sum]
#
# 	def query(self, query, axis='obs'):
# 		axis = get_axis(axis)
# 		if axis == 'obs':
# 			return self[self.obs.query(query),:]
# 		elif axis == 'var':
# 			return self[:,self.var.query(query)]
#
# 	@classmethod
# 	def from_anndata(cls, ad):
# 		ad.__class__ = cls
# 		return ad
#
# class MultiFeatureTable():
# 	def __init__(self, **kwargs):
# 		self.spaces = kwargs
#
# 	def __getattr__(self, attr):
# 		return self.spaces[attr]
#
# 	def project_features
