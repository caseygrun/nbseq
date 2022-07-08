# import weblogo
import numpy as np
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

def make_logo(df=None, seq_col='aligned',count_col='abundance', matrix=None,
			  color_scheme='chemistry_ambig', labels=None, legend=True, features={}, **kwargs):
	if matrix is None:
		mat_df = logomaker.alignment_to_matrix(df[seq_col],df[count_col])
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
		y = -.2
		xs = np.arange(-3, len(mat_df),10)
		ys = y*np.ones(len(xs))

		logo.ax.axhline(y, color='k', linewidth=1)
		for name, (start, stop) in features.items():
			f_start = start - 0.5
			f_stop = stop - 0.5
			logo.ax.plot([f_start, f_stop],[y, y], color='k', linewidth=10, solid_capstyle='butt')
			logo.ax.text(f_start, 2.5*y,name, verticalalignment='top', horizontalalignment='left')

	if legend and labels is not None:
		legend_elements = [Patch(facecolor=c, label=l) for (l, c) in labels.items()]

		if not isinstance(legend,dict): legend = {}
		plt.legend(handles=legend_elements,
			**{'bbox_to_anchor':(1.05, 1), 'loc':'upper left', **legend }
		)

	return logo, mat_df

# a_mat_df, logo = make_logo(a_df_c, seq_col='clonotype')


def cdr_seq_logo_plot(fda, cdrs=['CDR1','CDR2','CDR3'], title=None):
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
	widths = fda.loc[:,cdrs].iloc[0,:].apply(len)
	fig, axs = plt.subplots(ncols=len(cdrs), gridspec_kw=dict(width_ratios=widths), figsize=(widths.sum()/4,2))
	for ax, cdr in zip(axs, cdrs):
		cdr_logo, _ = make_logo(fda, seq_col=cdr, ax=ax);
		logos[cdr] = cdr_logo
		ax.get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))
		ax.get_yaxis().set_visible(False)
	if title is not None:
		plt.suptitle(title)
	plt.tight_layout()
	return logos


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
	for title, partition in iterator:
		fda = fortify(ft[partition,:], var=True)
		plots[title] = cdr_seq_logo_plot(fda, cols, title=(title if show_titles else None))

	return plots
