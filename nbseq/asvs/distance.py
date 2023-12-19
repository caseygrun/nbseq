"""calculate distances and distributions of distances between features (VHH sequences)
"""
import abc
import numpy as np
from Bio.Align import _substitution_matrices as _substitution_matrices

def populate_substitution_matrix(_data, invert=False, cls=_substitution_matrices.Array):

	alphabet = tuple(set(s[0] for s in _data.keys()) | {'X','-'})
	mat = _substitution_matrices.Array(alphabet, dims=2)

	for s, v in _data.items():
		mat[s[0],s[1]] = v

	for b in ['X','-']:
		for a in alphabet:
			mat[a,b] = 1
			mat[b,a] = 1
		mat[b,b] = 0

	if invert:
		for a in alphabet:
			for b in alphabet:
				mat[a,b] = 1-mat[a,b]
	return mat

def _get_epstein_data():
	return {
		('F', 'M'): 0.1,
		('F', 'L'): 0.15,
		('F', 'I'): 0.15,
		('F', 'V'): 0.2,
		('F', 'P'): 0.2,
		('F', 'Y'): 0.2,
		('F', 'W'): 0.21,
		('F', 'C'): 0.28,
		('F', 'A'): 0.5,
		('F', 'G'): 0.61,
		('F', 'S'): 0.81,
		('F', 'T'): 0.81,
		('F', 'H'): 0.8,
		('F', 'E'): 1,
		('F', 'Q'): 1,
		('F', 'D'): 1,
		('F', 'N'): 1,
		('F', 'K'): 1,
		('F', 'R'): 1,
		('M', 'F'): 0.05,
		('M', 'L'): 0.05,
		('M', 'I'): 0.05,
		('M', 'V'): 0.1,
		('M', 'P'): 0.1,
		('M', 'Y'): 0.22,
		('M', 'W'): 0.24,
		('M', 'C'): 0.22,
		('M', 'A'): 0.45,
		('M', 'G'): 0.56,
		('M', 'S'): 0.8,
		('M', 'T'): 0.8,
		('M', 'H'): 0.8,
		('M', 'E'): 1,
		('M', 'Q'): 1,
		('M', 'D'): 1,
		('M', 'N'): 1,
		('M', 'K'): 1,
		('M', 'R'): 1,
		('L', 'F'): 0.08,
		('L', 'M'): 0.03,
		('L', 'I'): 0,
		('L', 'V'): 0.05,
		('L', 'P'): 0.05,
		('L', 'Y'): 0.22,
		('L', 'W'): 0.25,
		('L', 'C'): 0.21,
		('L', 'A'): 0.43,
		('L', 'G'): 0.54,
		('L', 'S'): 0.8,
		('L', 'T'): 0.8,
		('L', 'E'): 1,
		('L', 'Q'): 1,
		('L', 'D'): 1,
		('L', 'N'): 1,
		('L', 'K'): 1,
		('L', 'R'): 1,
		('I', 'F'): 0.08,
		('I', 'M'): 0.03,
		('I', 'L'): 0,
		('I', 'V'): 0.05,
		('I', 'P'): 0.05,
		('I', 'Y'): 0.22,
		('I', 'W'): 0.25,
		('I', 'C'): 0.21,
		('I', 'A'): 0.43,
		('I', 'G'): 0.54,
		('I', 'S'): 0.8,
		('I', 'T'): 0.8,
		('I', 'H'): 1,
		('I', 'E'): 1,
		('I', 'Q'): 1,
		('I', 'D'): 1,
		('I', 'N'): 1,
		('I', 'K'): 1,
		('I', 'R'): 1,
		('V', 'F'): 0.1,
		('V', 'M'): 0.1,
		('V', 'L'): 0.03,
		('V', 'I'): 0.03,
		('V', 'P'): 0,
		('V', 'Y'): 0.24,
		('V', 'W'): 0.27,
		('V', 'C'): 0.2,
		('V', 'A'): 0.41,
		('V', 'G'): 0.52,
		('V', 'S'): 0.8,
		('V', 'T'): 0.8,
		('V', 'H'): 0.8,
		('V', 'E'): 1,
		('V', 'Q'): 1,
		('V', 'D'): 1,
		('V', 'N'): 1,
		('V', 'K'): 1,
		('V', 'R'): 1.01,
		('P', 'F'): 0.1,
		('P', 'M'): 0.1,
		('P', 'L'): 0.03,
		('P', 'I'): 0.03,
		('P', 'V'): 0,
		('P', 'Y'): 0.24,
		('P', 'W'): 0.27,
		('P', 'C'): 0.2,
		('P', 'A'): 0.41,
		('P', 'G'): 0.52,
		('P', 'S'): 0.8,
		('P', 'T'): 0.8,
		('P', 'H'): 0.8,
		('P', 'E'): 1,
		('P', 'Q'): 1,
		('P', 'D'): 1,
		('P', 'N'): 1,
		('P', 'K'): 1,
		('P', 'R'): 1.01,
		('Y', 'F'): 0.21,
		('Y', 'M'): 0.25,
		('Y', 'L'): 0.28,
		('Y', 'I'): 0.28,
		('Y', 'V'): 0.32,
		('Y', 'P'): 0.32,
		('Y', 'W'): 0.05,
		('Y', 'C'): 0.25,
		('Y', 'A'): 0.4,
		('Y', 'G'): 0.5,
		('Y', 'S'): 0.62,
		('Y', 'T'): 0.61,
		('Y', 'H'): 0.6,
		('Y', 'E'): 0.8,
		('Y', 'Q'): 0.8,
		('Y', 'D'): 0.81,
		('Y', 'N'): 0.81,
		('Y', 'K'): 0.8,
		('Y', 'R'): 0.8,
		('W', 'F'): 0.25,
		('W', 'M'): 0.32,
		('W', 'L'): 0.36,
		('W', 'I'): 0.36,
		('W', 'V'): 0.4,
		('W', 'P'): 0.4,
		('W', 'Y'): 0.1,
		('W', 'C'): 0.35,
		('W', 'A'): 0.49,
		('W', 'G'): 0.58,
		('W', 'S'): 0.63,
		('W', 'T'): 0.63,
		('W', 'H'): 0.61,
		('W', 'E'): 0.81,
		('W', 'Q'): 0.81,
		('W', 'D'): 0.81,
		('W', 'N'): 0.81,
		('W', 'K'): 0.81,
		('W', 'R'): 0.8,
		('C', 'F'): 0.22,
		('C', 'M'): 0.21,
		('C', 'L'): 0.2,
		('C', 'I'): 0.2,
		('C', 'V'): 0.2,
		('C', 'P'): 0.2,
		('C', 'Y'): 0.13,
		('C', 'W'): 0.18,
		('C', 'A'): 0.22,
		('C', 'G'): 0.34,
		('C', 'S'): 0.6,
		('C', 'T'): 0.6,
		('C', 'H'): 0.61,
		('C', 'E'): 0.8,
		('C', 'Q'): 0.8,
		('C', 'D'): 0.8,
		('C', 'N'): 0.8,
		('C', 'K'): 0.8,
		('C', 'R'): 0.81,
		('A', 'F'): 0.43,
		('A', 'M'): 0.41,
		('A', 'L'): 0.43,
		('A', 'I'): 0.43,
		('A', 'V'): 0.4,
		('A', 'P'): 0.4,
		('A', 'Y'): 0.27,
		('A', 'W'): 0.3,
		('A', 'C'): 0.25,
		('A', 'G'): 0.1,
		('A', 'S'): 0.4,
		('A', 'T'): 0.4,
		('A', 'H'): 0.42,
		('A', 'E'): 0.61,
		('A', 'Q'): 0.61,
		('A', 'D'): 0.61,
		('A', 'N'): 0.61,
		('A', 'K'): 0.61,
		('A', 'R'): 0.62,
		('G', 'F'): 0.53,
		('G', 'M'): 0.42,
		('G', 'L'): 0.51,
		('G', 'I'): 0.51,
		('G', 'V'): 0.5,
		('G', 'P'): 0.5,
		('G', 'Y'): 0.36,
		('G', 'W'): 0.39,
		('G', 'C'): 0.31,
		('G', 'A'): 0.1,
		('G', 'S'): 0.3,
		('G', 'T'): 0.31,
		('G', 'H'): 0.34,
		('G', 'E'): 0.52,
		('G', 'Q'): 0.52,
		('G', 'D'): 0.51,
		('G', 'N'): 0.51,
		('G', 'K'): 0.52,
		('G', 'R'): 0.53,
		('S', 'F'): 0.81,
		('S', 'M'): 0.8,
		('S', 'L'): 0.8,
		('S', 'I'): 0.8,
		('S', 'V'): 0.8,
		('S', 'P'): 0.8,
		('S', 'Y'): 0.62,
		('S', 'W'): 0.63,
		('S', 'C'): 0.6,
		('S', 'A'): 0.4,
		('S', 'G'): 0.32,
		('S', 'T'): 0.03,
		('S', 'H'): 0.1,
		('S', 'E'): 0.22,
		('S', 'Q'): 0.22,
		('S', 'D'): 0.21,
		('S', 'N'): 0.21,
		('S', 'K'): 0.22,
		('S', 'R'): 0.24,
		('T', 'F'): 0.81,
		('T', 'M'): 0.8,
		('T', 'L'): 0.8,
		('T', 'I'): 0.8,
		('T', 'V'): 0.8,
		('T', 'P'): 0.8,
		('T', 'Y'): 0.61,
		('T', 'W'): 0.63,
		('T', 'C'): 0.6,
		('T', 'A'): 0.41,
		('T', 'G'): 0.34,
		('T', 'S'): 0.03,
		('T', 'H'): 0.08,
		('T', 'E'): 0.21,
		('T', 'Q'): 0.21,
		('T', 'D'): 0.2,
		('T', 'N'): 0.2,
		('T', 'K'): 0.21,
		('T', 'R'): 0.22,
		('H', 'F'): 0.8,
		('H', 'M'): 0.8,
		('H', 'L'): 0.81,
		('H', 'I'): 0.81,
		('H', 'V'): 0.81,
		('H', 'P'): 0.81,
		('H', 'Y'): 0.6,
		('H', 'W'): 0.61,
		('H', 'C'): 0.62,
		('H', 'A'): 0.47,
		('H', 'G'): 0.42,
		('H', 'S'): 0.1,
		('H', 'T'): 0.08,
		('H', 'E'): 0.2,
		('H', 'Q'): 0.2,
		('H', 'D'): 0.21,
		('H', 'N'): 0.21,
		('H', 'K'): 0.2,
		('H', 'R'): 0.2,
		('E', 'F'): 1,
		('E', 'M'): 1,
		('E', 'L'): 1,
		('E', 'I'): 1,
		('E', 'V'): 1,
		('E', 'P'): 1,
		('E', 'Y'): 0.8,
		('E', 'W'): 0.81,
		('E', 'C'): 0.81,
		('E', 'A'): 0.63,
		('E', 'G'): 0.56,
		('E', 'S'): 0.21,
		('E', 'T'): 0.21,
		('E', 'H'): 0.2,
		('E', 'Q'): 0,
		('E', 'D'): 0.03,
		('E', 'N'): 0.03,
		('E', 'K'): 0,
		('E', 'R'): 0.05,
		('Q', 'F'): 1,
		('Q', 'M'): 1,
		('Q', 'L'): 1,
		('Q', 'I'): 1,
		('Q', 'V'): 1,
		('Q', 'P'): 1,
		('Q', 'Y'): 0.8,
		('Q', 'W'): 0.81,
		('Q', 'C'): 0.81,
		('Q', 'A'): 0.63,
		('Q', 'G'): 0.56,
		('Q', 'S'): 0.21,
		('Q', 'T'): 0.21,
		('Q', 'H'): 0.2,
		('Q', 'E'): 0,
		('Q', 'D'): 0.03,
		('Q', 'N'): 0.03,
		('Q', 'K'): 0,
		('Q', 'R'): 0.05,
		('D', 'F'): 1,
		('D', 'M'): 1,
		('D', 'L'): 1,
		('D', 'I'): 1,
		('D', 'V'): 1,
		('D', 'P'): 1,
		('D', 'Y'): 0.81,
		('D', 'W'): 0.81,
		('D', 'C'): 0.8,
		('D', 'A'): 0.62,
		('D', 'G'): 0.54,
		('D', 'S'): 0.2,
		('D', 'T'): 0.2,
		('D', 'H'): 0.21,
		('D', 'E'): 0.03,
		('D', 'Q'): 0.03,
		('D', 'N'): 0,
		('D', 'K'): 0.03,
		('D', 'R'): 0.08,
		('N', 'F'): 1,
		('N', 'M'): 1,
		('N', 'L'): 1,
		('N', 'I'): 1,
		('N', 'V'): 1,
		('N', 'P'): 1,
		('N', 'Y'): 0.81,
		('N', 'W'): 0.81,
		('N', 'C'): 0.8,
		('N', 'A'): 0.62,
		('N', 'G'): 0.54,
		('N', 'S'): 0.2,
		('N', 'T'): 0.2,
		('N', 'H'): 0.21,
		('N', 'E'): 0.03,
		('N', 'Q'): 0.03,
		('N', 'D'): 0,
		('N', 'K'): 0.03,
		('N', 'R'): 0.08,
		('K', 'F'): 1,
		('K', 'M'): 1,
		('K', 'L'): 1,
		('K', 'I'): 1,
		('K', 'V'): 1,
		('K', 'P'): 1,
		('K', 'Y'): 0.8,
		('K', 'W'): 0.81,
		('K', 'C'): 0.81,
		('K', 'A'): 0.63,
		('K', 'G'): 0.56,
		('K', 'S'): 0.21,
		('K', 'T'): 0.21,
		('K', 'H'): 0.2,
		('K', 'E'): 0,
		('K', 'Q'): 0,
		('K', 'D'): 0.03,
		('K', 'N'): 0.03,
		('K', 'R'): 0.0,
		('R', 'F'): 1,
		('R', 'M'): 1,
		('R', 'L'): 1.01,
		('R', 'I'): 1.01,
		('R', 'V'): 1.02,
		('R', 'P'): 1.02,
		('R', 'Y'): 0.8,
		('R', 'W'): 0.8,
		('R', 'C'): 0.82,
		('R', 'A'): 0.67,
		('R', 'G'): 0.61,
		('R', 'S'): 0.24,
		('R', 'T'): 0.22,
		('R', 'H'): 0.2,
		('R', 'E'): 0.05,
		('R', 'Q'): 0.05,
		('R', 'D'): 0.08,
		('R', 'N'): 0.08,
		('R', 'K'): 0.05
	}

epstein = populate_substitution_matrix(_get_epstein_data(), invert=False)

substitution_matrices = {
	'epstein': epstein
}

def get_matrix(matrix='BLOSUM62',gap='-'):
	"""loads a substituion matrix from Bio.Align.substitution_matrices or from this package
	matrix: str or Bio.Align.substitution_matrices.Array; if str, will try to load from ``substitution_matrices``, then ``Bio.Align.substitution_matrices``
	gap: for BioPython substitution matrices, change the gap character from '*' to this
	"""
	if (isinstance(matrix, str) and matrix == 'identity') or matrix is None:
		return None
	if not isinstance(matrix, _substitution_matrices.Array):
		if matrix in substitution_matrices:
			matrix = substitution_matrices[matrix]
		else:
			m = _substitution_matrices.load(matrix)
			matrix = _substitution_matrices.Array(alphabet=m.alphabet.replace('*',gap), data=m.data)
	return matrix


def similarity_distribution(seqs, n, matrix='BLOSUM62', gap='-', weights=None, normalize=False, null=False, n_jobs=1):
	"""compute a distribution of a sequence similarity/distance measure, by subsampling from a set of sequences

	Parameters
	----------
	seqs: ndarray of str
		list of ASV sequences
	n : int
		number of sequences to sample, with replacement
	matrix: str or Bio.Align.substitution_matrices.Array or None
		will be passed to ``get_matrix`` to find a substitution Matrix; if None, identity will be used
	weights: ndarray, optional
		abundance or weight of each ASV sequences, will be used to weight the
		random selection towards more abundant ASVs. If not given, ASVs will be
		selected uniformly
	null: bool, default=False
		some substitution matrices (e.g. BLOSUM62) do not give zero for
		identical sequences; to compute the similarity metric on a set of
		identical sequences (e.g. the null distribution), set ``null=True``
	normalize: bool, default=False
		``True`` to divide the similarity measure for each pair by the longer of
		the two sequences
	n_jobs: int, default=1
		how many threads to use
	"""
	from joblib import Parallel
	matrix = get_matrix(matrix, gap)
	from numpy.random import default_rng; rng = default_rng()
	inverse = None

	# weight choice of sequences
	if weights is not None:
		p = weights / weights.sum()

		# if null, make an n x 2 array where pick_idxs[i, 0] == pick_idxs[i, 1]
		if null:
			pick_idxs = rng.choice(np.arange(len(seqs)), size=(n,))
			pick_idxs = np.broadcast_to(pick_idxs[:,np.newaxis],(n,2))
		else:
			pick_idxs = rng.choice(np.arange(len(seqs)), size=(n,2), p=p)

		# filter to unique pairs of sequences to avoid re-computing distances
		pick_idxs_uniq, inverse = np.unique(pick_idxs, return_inverse=True, axis=0)

		# picks.shape = (n, 2)
		picks = np.column_stack( (seqs[pick_idxs_uniq[:,0]], seqs[pick_idxs_uniq[:,1]]) )

	else:
		if null:
			picks = rng.choice(seqs, size=(n,))
			picks = np.broadcast_to(picks[:,np.newaxis],(n,2))
		else:
			picks = rng.choice(seqs, size=(n,2))

	with Parallel(n_jobs=n_jobs) as pool:
		similarities = np.array(calculate_sequence_similarities(picks, matrix, normalize=normalize, pool=pool))

	if inverse is not None:
		return similarities[inverse]
	else:
		return similarities

def calculate_similarity_matrix(seqs, **kwargs):
	from skbio.stats.distance import DissimilarityMatrix
	import numpy as np
	import itertools
	seq_pairs = itertools.combinations(seqs, 2)
	dm = np.zeros((len(seqs),) * 2)
	index = dict(zip(seqs, range(len(seqs))))
	scores = calculate_sequence_similarities(seq_pairs=seqs, **kwargs)
	for (seq1, seq2), score in zip(seq_pairs, scores):
		dm[index[seq1], index[seq2]] = score
		dm[index[seq2], index[seq1]] = score
	return DissimilarityMatrix(dm, seqs)

def calculate_sequence_similarities(seq_pairs, matrix=None, normalize=False, pool=None):
	"""calculates similarity metric on an array of sequence pairs

	Parameters
	----------
	seq_pairs : ndarray
		n x 2 array of sequence strings; will calculate sum of
		matrix[seq_pairs[i,0], seq_pairs[i,1]] for all i in range(n)
	pool : joblib.ParallelPool, optional
		if given, will use an existing worker pool

	Returns
	-------
	ndarray
		vector of similarity scores for all pairs in seq_pairs, len = n
	"""
	from joblib import Parallel, delayed
	from ..utils import len_ungapped

	if pool is None:
		pool = Parallel()
	if matrix is None:
		if normalize:
			_f = lambda s1, s2: sum(c1==c2 for c1, c2 in zip(s1, s2))/len(s1)
		else:
			_f = lambda s1, s2: sum(c1==c2 for c1, c2 in zip(s1, s2))
	else:
		if normalize:
			_f = lambda s1,s2: sum(matrix[c1, c2] for c1, c2 in zip(s1, s2)) / max(len_ungapped(s1),len_ungapped(s2))
		else:
			_f = lambda s1,s2: sum(matrix[c1, c2] for c1, c2 in zip(s1, s2))
	scores = pool(delayed( _f )(s1,s2)
		for s1, s2 in seq_pairs)
	return scores

def null_similarity_distribution(*args, **kwargs):
	return similarity_distribution(*args, null=True, **kwargs)

def identity_distribution(seqs, n, matrix=None, normalize=True, **kwargs):
	return similarity_distribution(seqs, n, matrix=matrix, normalize=normalize, **kwargs)


def plot_similarity_distribution_with_sizes(data_frames,
	ns = [1000,10_000,50_000,100_000],
	labels=['alpaca','synthetic'],
	seq_col='align',
	matrix=epstein,
	normalize=True,n_jobs=4):
	"""plot the ``similarity_distribution`` of one or more sets of sequences, subsampling at different values of ``n``

	Parameters
	----------
	data_frames : dict of str:(pd.Series or list), or list
		 sequences to subsample from
	ns : array_like
		values to use for ``n`` when calling ``similarity_distribution``; will produce one plot per entry
	"""
	import matplotlib.pyplot as plt
	import seaborn as sns

	fig, axs = plt.subplots(nrows=len(ns), figsize=(4,len(ns)*3))

	for i,n in enumerate(ns):
		data = {}
		if isinstance(data_frames, abc.Mapping):
			_items = data_frames.items()
			for label, df in data_frames.items():
				data[label] = similarity_distribution(df[seq_col], n, matrix=matrix, n_jobs=n_jobs, normalize=normalize)
		elif len(labels) == len(data_frames):
			for label, df in zip(labels, data_frames):
				data[label] = similarity_distribution(df[seq_col], n, matrix=matrix, n_jobs=n_jobs, normalize=normalize)

		g = sns.histplot(data=data,
						 common_norm=True,
						 ax=axs[i],
						 element="step");
		if normalize:
			axs[i].set_xlabel('Fraction of differing residues');
		else:
			axs[i].set_xlabel('# of differing residues');
		axs[i].set_title(f'n = {n:.0e} sequence pairs')
		axs[i].set_ylabel('# pairs');
	plt.tight_layout()
	return fig, axs
