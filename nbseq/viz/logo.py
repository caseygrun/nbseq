# import weblogo
import numpy as np
import pandas as pd
import logomaker
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

# chemistry_ambig = weblogo.ColorScheme(
# 		[
# 			weblogo.SymbolColor("BSTYCZ", "green", "polar"),
# 			weblogo.SymbolColor("JNQ", "purple", "neutral"),
# 			weblogo.SymbolColor("KRH", "blue", "basic"),
# 			weblogo.SymbolColor("DE", "red", "acidic"),
# 			weblogo.SymbolColor("PAWFLIMV", "black", "hydrophobic"),
# 			weblogo.SymbolColor("-", "#DDDDDD", "gap")
# 		],
# 		alphabet=weblogo.protein_alphabet
# )
# chemistry_ambig = weblogo.ColorScheme(
# 		[
# 			weblogo.SymbolColor("BSTNQYZ", "green", "polar"),
# 			weblogo.SymbolColor("JNGC", "purple", "neutral"),
# 			weblogo.SymbolColor("KRH", "blue", "basic"),
# 			weblogo.SymbolColor("DE", "red", "acidic"),
# 			weblogo.SymbolColor("AVILMFW", "black", "hydrophobic"),
# 			weblogo.SymbolColor("-", "#DDDDDD", "gap")
# 		],
# 		alphabet=weblogo.protein_alphabet
# )
#
# def seq_abundances_to_count_matrix(seqs, abundances, alphabet):
# 	if not alphabet:
# 		raise ValueError("No alphabet")
#
# 	if len(seqs) != len(abundances):
# 		raise ValueError("Must have exactly one count per sequence")
#
# 	# number of characters in alphabet
# 	N = len(alphabet)
# 	ords = [alphabet.ords(s) for s in seqs]
#
# 	# length of first (= each) sequence
# 	L = len(ords[0])
#
# 	#     counts = [[0, ] * N for l in range(0, L)]
# 	# pre-allocate
# 	counts = np.zeros((L,N))
#
# 	# for each sequence
# 	for o, count in zip(ords, abundances):
#
# 		# check length is correct
# 		if len(o) != L:
# 			raise ValueError("Sequences are of incommensurate lengths. Cannot tally.")
#
# 		# build up counts
# 		for j, n in enumerate(o):
# 			if n < N:
# 				counts[j,n] += count
#
# 	return counts
#
# def seq_abundances_to_motif(seqs, abundances, alphabet):
# 	from weblogo.matrix import Motif
#
# 	counts = seq_abundances_to_count_matrix(seqs, abundances, alphabet)
# 	return Motif(alphabet, counts)
#
# def make_seqlogo(
# 		table, seq_col='sequence', abundance_col='abundance',
# 		logo_start=1, logo_end=None, color_scheme=chemistry_ambig, **kwargs):
#
# 	alphabet = weblogo.protein_alphabet
#
# 	valid_rows = (~table[seq_col].isna() & (table[seq_col] != ''))
# 	sequences = table.loc[valid_rows, seq_col]
#
# 	if logo_end is None:
# 		logo_end = len(sequences[0])
#
# 	if abundance_col is not None:
# 		abundances = table.loc[valid_rows, 'abundance']
# 		counts = seq_abundances_to_motif(sequences, abundances, alphabet)
# 		logodata = weblogo.LogoData.from_counts(alphabet, counts)
# 	else:
# 		sequences = [weblogo.seq.Seq(row, name=ix, alphabet=alphabet)
# 			for (i, (ix, row)) in enumerate(sequences.iteritems())]
# 		logodata = weblogo.LogoData.from_seqs(
# 			weblogo.seq.SeqList(sequences, alphabet=alphabet))
#
# 	logooptions = weblogo.LogoOptions(
# 		color_scheme=color_scheme,
# 		logo_start=logo_start,
# 		logo_end=logo_end,
# 		stacks_per_line=logo_end,
# 		 **{'unit_name':'probability',
# 			'show_boxes':False,
# 			'show_errorbars':False,
# 			**kwargs})
# 	logoformat = weblogo.LogoFormat(logodata, logooptions)
# 	return (logodata, logoformat)
#
# def display_logo(logo, format='png'):
# 	if format == 'png':
# 		from IPython.core.display import display_png
# 		display_png(weblogo.png_formatter(*logo),raw=True)
# 	elif format=='svg':
# 		from IPython.core.display import display_svg
# 		display_svg(weblogo.svg_formatter(*logo),raw=True)
# 	elif format=='pdf':
# 		from IPython.core.display import display_pdf
# 		display_pdf(weblogo.pdf_formatter(*logo),raw=True)

chemistry_ambig = [
	("BSTYCZ"    , "green"  , "polar")       ,
	("JNQ"       , "purple" , "neutral")     ,
	("KRH"       , "blue"   , "basic")       ,
	("DE"        , "red"    , "acidic")      ,
	("GPAWFLIMV" , "black"  , "hydrophobic") ,
	("X"         , "grey"   , "unknown")
]

chemistry_ambig_palette = { chars:color for (chars, color, label) in chemistry_ambig }




def alignment_to_matrix(sequences,
                        counts=None,
                        to_type='counts',
                        background=None,
                        characters_to_ignore='.-',
                        unique_characters=None,
                        center_weights=False,
                        pseudocount=1.0):
	"""
	Generates matrix from a sequence alignment

	adapted from logomaker.src.matrix.alignment_to_matrix, with some patches
	that result in ~4x speedup

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
	background: (array, or df)
		Specification of background probabilities. If array, should be the
		same length as df.columns and correspond to the probability of each
		column's character. If df, should be a probability matrix the same
		shape as df.
	characters_to_ignore: (str)
		Characters to ignore within sequences. This is often needed when
		creating matrices from gapped alignments.
	center_weights: (bool)
		Whether to subtract the mean of each row, but only if to_type=='weight'.
	pseudocount: (number >= 0.0)
		Pseudocount to use when converting from counts to probabilities.
	returns
	-------
	out_df: (dataframe)
		A matrix of the requested type.
	"""
	from ..asvs import marginal_frequencies
	from logomaker.src.matrix import transform_matrix, check, MATRIX_TYPES

	# Define valid types
	valid_types = MATRIX_TYPES.copy()

	# Check that to_type is valid
	check(to_type in valid_types,
            'to_type=%s; must be in %s' % (to_type, valid_types))

	# validate center_weights
	check(isinstance(center_weights, bool),
            'type(center_weights) = %s; must be bool.' % type(center_weights))

	counts_df = marginal_frequencies(
            sequences, counts, background=background, characters_to_ignore=characters_to_ignore,
                                  unique_characters=unique_characters)

	# Convert counts matrix to matrix of requested type
	out_df = transform_matrix(counts_df,
                           from_type='counts',
                           to_type=to_type,
                           pseudocount=pseudocount,
                           background=background)

	# Center values only if center_weights is True and to_type is 'weight'
	if center_weights and to_type == 'weight':
		out_df = transform_matrix(out_df, center_values=True)

	return out_df



# def alignment_to_matrix(sequences,
# 						counts=None,
# 						to_type='counts',
# 						background=None,
# 						characters_to_ignore='.-',
# 						unique_characters=None,
# 						center_weights=False,
# 						pseudocount=1.0):
# 	"""
# 	Generates matrix from a sequence alignment

# 	adapted from logomaker.src.matrix.alignment_to_matrix, with some patches
# 	that result in ~4x speedup

# 	parameters
# 	----------
# 	sequences: (list of strings)
# 		A list of sequences, all of which must be the same length
# 	counts: (None or list of numbers)
# 		If not None, must be a list of numbers the same length os sequences,
# 		containing the (nonnegative) number of times that each sequence was
# 		observed. If None, defaults to 1.
# 	to_type: (str)
# 		The type of matrix to output. Must be 'counts', 'probability',
# 		'weight', or 'information'
# 	background: (array, or df)
# 		Specification of background probabilities. If array, should be the
# 		same length as df.columns and correspond to the probability of each
# 		column's character. If df, should be a probability matrix the same
# 		shape as df.
# 	characters_to_ignore: (str)
# 		Characters to ignore within sequences. This is often needed when
# 		creating matrices from gapped alignments.
# 	center_weights: (bool)
# 		Whether to subtract the mean of each row, but only if to_type=='weight'.
# 	pseudocount: (number >= 0.0)
# 		Pseudocount to use when converting from counts to probabilities.
# 	returns
# 	-------
# 	out_df: (dataframe)
# 		A matrix of the requested type.
# 	"""
# 	from logomaker.src.matrix import transform_matrix, check, MATRIX_TYPES

# 	# validate inputs

# 	# Make sure sequences is list-like
# 	check(isinstance(sequences, (list, tuple, np.ndarray, pd.Series)),
# 		  'sequences must be a list, tuple, np.ndarray, or pd.Series.')
# 	# sequences = list(sequences)

# 	# Make sure sequences has at least 1 element
# 	check(len(sequences) > 0, 'sequences must have length > 0.')

# 	# Make sure all elements are sequences
# 	check(all(isinstance(seq, str) for seq in sequences),
# 		  'sequences must all be of type string')

# 	# validate characters_to_ignore
# 	check(isinstance(characters_to_ignore, str),
# 		  'type(seq) = %s must be of type str' % type(characters_to_ignore))

# 	# validate center_weights
# 	check(isinstance(center_weights, bool),
# 		  'type(center_weights) = %s; must be bool.' % type(center_weights))

# 	# Get sequence length
# 	L = len(sequences[0])

# 	# Make sure all sequences are the same length
# 	check(all([len(s) == L for s in sequences]),
# 		  'all elements of sequences must have the same length.')

# 	# validate counts as list-like
# 	check(isinstance(counts, (list, tuple, np.ndarray, pd.Series)) or
# 		  (counts is None),
# 		  'counts must be None or a list, tuple, np.ndarray, or pd.Series.')

# 	# make sure counts has the same length as sequences
# 	if counts is None:
# 		counts = np.ones(len(sequences))
# 	else:
# 		check(len(counts) == len(sequences),
# 			  'counts must be the same length as sequences;'
# 			  'len(counts) = %d; len(sequences) = %d' %
# 			  (len(counts), len(sequences)))

# 	# validate background
# 	check(isinstance(background, (type([]), np.ndarray, pd.DataFrame)) or
# 		  (background is None),
# 		  'type(background) = %s must be None or array-like or a dataframe.' %
# 		  type(background))

# 	# Define valid types
# 	valid_types = MATRIX_TYPES.copy()

# 	# Check that to_type is valid
# 	check(to_type in valid_types,
# 		  'to_type=%s; must be in %s' % (to_type, valid_types))

# 	# View as a 2D array of characters
# 	# we can do this without copying since we've already checked all strings are same length
# 	# char_array = np.array(sequences, dtype='str').view('U1').reshape((sequences.size, -1))

# 	# view as np.uint32 for faster comparisons (~2x speedup)
# 	char_array = np.array(sequences, dtype='str').view(
# 		np.uint32).reshape((sequences.size, -1))

# 	# Get list of unique characters
# 	if unique_characters is None:
# 		unique_characters = np.unique(char_array.ravel())

# 	unique_characters.sort()

# 	# Remove characters to ignore
# 	columns = np.array([c for c in unique_characters.view(
# 		'U1') if not c in characters_to_ignore])
# 	index = list(range(L))
# 	counts_df = pd.DataFrame(data=0, columns=columns, index=index)

# 	# Sum of the number of occurrences of each character at each position
# 	for c, x in zip(columns, columns.view(np.uint32)):
# 		tmp_mat = (char_array == x).astype(float) * counts[:, np.newaxis]
# 		counts_df.loc[:, c] = tmp_mat.sum(axis=0).T

# 	# Convert counts matrix to matrix of requested type
# 	out_df = transform_matrix(counts_df,
# 							  from_type='counts',
# 							  to_type=to_type,
# 							  pseudocount=pseudocount,
# 							  background=background)

# 	# Center values only if center_weights is True and to_type is 'weight'
# 	if center_weights and to_type == 'weight':
# 		out_df = transform_matrix(out_df, center_values=True)

# 	return out_df



def make_logo(df=None, seq_col='aligned',count_col='abundance', matrix=None,
			  color_scheme='chemistry_ambig', labels=None, legend=True, features={}, **kwargs):
	if matrix is None:
		mat_df = alignment_to_matrix(df[seq_col],df[count_col])
	else:
		mat_df = matrix

	if color_scheme == 'chemistry_ambig':
		color_scheme = chemistry_ambig_palette
		if labels is None:
			labels = { l: c for (chars, c, l) in chemistry_ambig}

	logo = logomaker.Logo(mat_df,
		color_scheme=color_scheme,
		stack_order='small_on_top',
		**kwargs
		)

	if len(features) > 0:
		import matplotlib.transforms as transforms
		ax = logo.ax

		y = -.1
		xs = np.arange(-3, len(mat_df),10)
		ys = y*np.ones(len(xs))

		logo.ax.axhline(y, color='k', linewidth=1)
		
		# blended transform lets us use data coords for x axis and axis coords for y axis
		# so labels are plotted below axis
		# https://stackoverflow.com/questions/63153629/matplotlib-text-use-data-coords-for-x-axis-coords-for-y
		# trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
		trans = ax.get_xaxis_transform() # equivalent to above
		
		for name, (start, stop) in features.items():
			f_start = start - 0.5
			f_stop = stop - 0.5
			ax.plot([f_start, f_stop], [y, y], color='k', linewidth=20, solid_capstyle='butt', 
				transform=trans)
			ax.text(f_start, 2*y, name, verticalalignment='top',
			             horizontalalignment='left', transform=trans)

	if legend and labels is not None:
		legend_elements = [Patch(facecolor=c, label=l) for (l, c) in labels.items()]

		if not isinstance(legend,dict): legend = {}
		plt.legend(handles=legend_elements,
			**{'bbox_to_anchor':(1.05, 1), 'loc':'upper left', **legend }
		)

	return logo, mat_df

# a_mat_df, logo = make_logo(a_df_c, seq_col='clonotype')


def _setup_cdr3_logo_plots(fda, cdrs=['CDR1', 'CDR2', 'CDR3'], nrows=1):
	widths = fda.loc[:, cdrs].iloc[0, :].apply(len)
	fig, axs = plt.subplots(ncols=len(cdrs), gridspec_kw=dict(
		width_ratios=widths), figsize=(widths.sum()/4, 2))
	return fig, axs


def cdr_seq_logo_plot(fda, cdrs=['CDR1','CDR2','CDR3'], title=None, axs=None):
	"""Plot a set of sequence logos, one for each CDR

	Parameters
	----------
	fda : pd.DataFrame
		Columns of sequences named by ``cdrs``
	cdrs : list of str, optional
		Names of columns identifying the CDRs
	title : str
		Title of the sample to show as plt.suptitle

	Returns
	-------
	dict
		Keys are CDRs given in ``cdrs``, values are the sequence logos returned by `make_logo`
	"""
	logos = {}

	if axs is None:
		fig, axs = _setup_cdr3_logo_plots(fda, cdrs)
	else:
		fig = axs[0].figure
		
	for ax, cdr in zip(axs, cdrs):
		cdr_logo, _ = make_logo(fda, seq_col=cdr, ax=ax);
		logos[cdr] = cdr_logo
		ax.get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))
		ax.get_yaxis().set_visible(False)
	if title is not None:
		plt.suptitle(title)
	plt.tight_layout()
	return fig
	# return logos


def partition_logo_plot(ft, sample_partitions, cols=['CDR1','CDR2','CDR3']):
	"""show a sequence logo for each of a set of partitions of a given feature table
	"""
	from ..ft import fortify

	if not isinstance(sample_partitions, dict):
		show_titles = False
		iterator = enumerate(sample_partitions)
	else:
		show_titles = True
		iterator = sample_partitions.items()

	plots = {}
	nrows = len(iterator)
	for title, partition in iterator:
		fda = fortify(ft[partition,:], var=True)
		plots[title] = cdr_seq_logo_plot(fda, cols, title=(title if show_titles else None))

	return plots
