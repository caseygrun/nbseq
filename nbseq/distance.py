
# generate list of unique pairs
# compute edit distance between each pair

from .utils import *
import numpy as np

from skbio import DistanceMatrix

def jaccard(a,b):
	'''
	Implementation of the Jaccard distance for sparse vectors.
	Jaccard distance = 1 - (Jaccard Index)
	Returns a distance in the range [0, 1].
	>>> A = ss.csr_matrix([[0,2,0,1]])
	>>> B = ss.csr_matrix([[1,0,3,0]])
	>>> C = ss.csr_matrix([[0,2,3,5]])
	>>> D = ss.csr_matrix([[0.25, 0.5, 0.25]])
	>>> E = ss.csr_matrix([[0.5, 0.3, 0.2]])
	>>> jaccard(A, A)
	0.0
	>>> jaccard(A, B)
	1.0
	>>> jaccard(B, A)
	1.0
	>>> jaccard(B, C)
	0.75


	https://github.com/zyndagj/sparseDist
	Copyright (c) 2015 Greg Zynda
	MIT license
	'''
	# _validate_svs(a,b)
	intersect = len(np.intersect1d(a.indices, b.indices))
	union = len(np.union1d(a.indices, b.indices))
	# prevent divide by zero
	return 1 - np.float(intersect)/max(np.float(union),1)

def braycurtis(a,b):
	'''
	Computes the Bray-Curtis distance between two 1-D sparse matrices. Implementation follows http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.braycurtis.html#scipy.spatial.distance.braycurtis
	```
	sum(|a-b|)/sum(|a+b|)
	```

	The Bray-Curtis distance lies in the range [0,1].
	>>> a = [0,2,0,1]
	>>> b = [1,0,3,0]
	>>> c = [0,2,3,1]
	>>> braycurtis(ss.csr_matrix([a]), ss.csr_matrix([b]))
	1.0
	>>> braycurtis(ss.csr_matrix([c]), ss.csr_matrix([b]))
	0.40000000000000002

	https://github.com/zyndagj/sparseDist
	Copyright (c) 2015 Greg Zynda
	MIT license
	'''
	return np.add.reduce(np.abs((a-b).data)) / float(np.add.reduce(np.abs((a+b).data)))

def shared_features(a, b):
	return len(np.intersect1d(a.indices, b.indices, assume_unique=True))

def shared_features_frac(a,b):
	return len(np.intersect1d(a.indices, b.indices, assume_unique=True)) / max(a.getnnz(), b.getnnz())

def shared_reads(a, b):
	intersect, comm_a, comm_b = np.intersect1d(a.indices, b.indices,
		assume_unique=True,
		return_indices=True)

	return (a.data[comm_a].sum() + b.data[comm_b].sum())

def shared_reads_frac(a, b):
	intersect, comm_a, comm_b = np.intersect1d(a.indices, b.indices,
		assume_unique=True,
		return_indices=True)

	return (a.data[comm_a].sum() + b.data[comm_b].sum()) / (a.sum()+b.sum())

	# tot = max(a.sum(), b.sum())
	# shared = max(a.data[comm_a].sum(), b.data[comm_b].sum())


sparse_metrics = {
	'braycurtis':braycurtis,
	'jaccard':jaccard,
	'shared_features':shared_features,
	'shared_reads':shared_reads,
	'shared_features_frac':shared_features_frac
}

phylogenetic_metrics = {
	'unweighted_unifrac',
	'weighted_unifrac'
}


def pairwise_comparison(ft, metric='braycurtis', axis='sample', threads='auto', **kwargs):
	import sklearn.metrics.pairwise

	# by default, biom feature table has rows as observations, columns as samples
	if axis == 'sample':
		# rows are samples, columns are observations
		# CSR so data is arranged row-wise
		ftr = ft.matrix_data.T.tocsr()
	elif axis == 'observation':
		ftr = ft.matrix_data.tocsr()
	else: raise Exception(f"Unknown axis '{axis}'")

	if metric in sparse_metrics:
		metric = sparse_metrics[metric]

	if threads == 'auto':
		n_jobs = -1
	else: n_jobs = threads

	dist = sklearn.metrics.pairwise.pairwise_distances(ftr,
		metric=metric, n_jobs=n_jobs,
		**kwargs)

	return dist

def pairwise_distance(ft, metric='braycurtis', axis='sample', threads='auto', **kwargs):
	import skbio

	if metric in sparse_metrics:
		dist = pairwise_comparison(ft, metric=metric, axis=axis, threads=threads, **kwargs)

	elif metric in phylogenetic_metrics:
		if axis != 'sample': raise NotImplementedError()
		if not 'tree' in kwargs:
			raise Exception(f'Must provide phylogenetic tree for phylogenetic distance metric {metric}')
		dist = skbio.diversity.beta_diversity(counts=ft, **{'metric':metric, 'otu_ids':ft.ids('sample'), **kwargs})

	# in case any of the distances ended up being nan, fix them
	dist = repair_distance_matrix(dist)

	return DistanceMatrix(dist, ids=ft.ids('sample'))

def pairwise_distance_phylo(ft_path, tree_path, outfile=None, metric='weighted_unifrac', **kwargs):
	import unifrac

	if not metric in phylogenetic_metrics:
		raise NotImplementedError(f"{metric} not implemented.")

	if metric == 'unweighted_unifrac':
		dist = unifrac.unweighted(table = ft_path,
			phylogeny = tree_path, **kwargs)
		# dist = unifrac.unweighted_to_file(table = ft_path,
		# 	phylogeny = tree_path, out_filename = outfile, **kwargs)
	elif metric == 'weighted_unifrac':
		# dist = unifrac.weighted_normalized_to_file(table = ft_path,
		# 	phylogeny = tree_path, out_filename = outfile, **kwargs)
		dist = unifrac.weighted_normalized(table = ft_path,
			phylogeny = tree_path, **kwargs)

	if outfile is not None:
		dist.write(outfile)

	return dist

def repair_distance_matrix(data, min_distance=0, max_distance=1, make_hollow=True):
	"""replace NaNs in a symmetric distance matrix with valid values
	"""
	for i,j in np.argwhere(np.isnan(data)):
		if i == j:
			data[i,j] = min_distance
		elif not np.isnan(data[j,i]):
			data[i,j] = data[j,i]
		else:
			data[i,j] = max_distance
	if make_hollow:
		np.fill_diagonal(data, 0)
	return data

def save_distance_matrix(dist, path):
	mkdirp_file(path)
	with open(path, 'w') as f:
		dist.write(f)

def load_distance_matrix(path):
	with open(path, 'r') as f:
		return DistanceMatrix.read(f)
