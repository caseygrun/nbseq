import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

VHH_LIBRARIES = ['alpaca','synthetic']

def asv_len_plot(asvs_translated, libraries = None, weighted=True, length_type='Peptide residues (aa)'):
	"""
	asvs_translated: pd.DataFrame with columns:
	- library
	- length
	- abundance (if weighted=True)
	"""
	if libraries is None:
		libraries = pd.unique(asvs_translated['library'])

	for library in libraries:
		asvs = asvs_translated.loc[asvs_translated['library'] == library,:]
		weights = asvs['abundance'] if weighted else None
		plt.hist(asvs['length'],
				 weights=weights, label = library, alpha=0.5);

	if length_type == 'aa':
		plt.xlabel('Peptide residues (aa)');
	elif length_type == 'nt':
		plt.xlabel('(nt)');
	else: plt.xlabel(length_type);

	if weighted: plt.ylabel('Abundance (reads)');
	else: plt.ylabel('Abundance (# of ASVs)');
	plt.legend(loc='upper left')


def segment_len_plot(segment_lens, libraries=None,
	segments=['CDR1','CDR2','CDR3'],
	segment_kws = None, engine='seaborn',
	weighted=True, fig_kw=None, **kwargs):

	if segment_kws is None:
		segment_kws = {}
	for segment in segments:
		if segment not in segment_kws:
			segment_kws[segment] = {}

	if libraries is None:
		if 'library' in segment_lens:
			libraries = pd.unique(segment_lens['library'])
		else:
			libraries = [None]

	widths = segment_lens[segments].max() - segment_lens[segments].min()

	if fig_kw is None:
		fig_kw = dict()

	if 'gridspec_kw' not in fig_kw:
		fig_kw['gridspec_kw'] = {}
	fig_kw['gridspec_kw'] = {'width_ratios': widths, **fig_kw['gridspec_kw']}

	fig_kw = {
		**dict(nrows=1, ncols=len(segments),
	 		figsize=(12, 4), 
			sharey=True), 
		**fig_kw
	}

	fig, axs = plt.subplots(**fig_kw)
	for library in libraries:
		kws = dict(**kwargs) 

		if 'library' in segment_lens:
			segment_lens_lib = segment_lens.query(f"library == '{library}'").dropna()
			kws['label'] = library
		else:
			segment_lens_lib = segment_lens
			
		if weighted: 
			kws['weights'] = segment_lens_lib['abundance']

		
		
		if engine == 'seaborn':
			import seaborn as sns
			for i, segment in enumerate(segments):
				_kws = kws
				if segment in segment_kws:
					_kws = {**kws, **segment_kws[segment]}
				sns.histplot(segment_lens_lib[segment], ax=axs[i], **_kws)
				axs[i].set_title(segment)
		else:
			segment_lens_lib.hist(column=segments, alpha=(1.0/len(libraries)), ax=axs, **kws)

	[ax.set_xlabel('length (nt)') for ax in axs.flat];
	if weighted: axs.flat[0].set_ylabel('abundance (reads)');
	else: axs.flat[0].set_ylabel('# of features');

	if len(libraries) > 1:
		axs.flat[0].legend(loc='upper left')
	plt.tight_layout()
	return fig, axs


from collections import OrderedDict
def plot_asv_stats(*args, subplots_kwargs={}, **kwargs):
	fig, axs = plt.subplots(1, 3, **{'figsize':(12,4), **subplots_kwargs})

	tables = []
	libraries = OrderedDict(list(enumerate(args)) + list(kwargs.items()))

	summary = pd.DataFrame(columns = ['total # of ASVs','total # of reads'], index=libraries.keys())
	summary.index.rename('library', inplace=True)

	for name, counts in libraries.items():
		n_ASVs = len(counts)
		n_reads = sum(counts)

		summary.loc[name, 'total # of ASVs'] = n_ASVs
		summary.loc[name, 'total # of reads'] = n_reads

#         print(f"Total # of ASVs: {n_ASVs}")
#         print(f"Total # of reads: {n_reads}")


		counts.hist(ax=axs[0], label=name)
		axs[0].set_yscale('log')
		axs[0].set_xlabel('ASV abundance');
		axs[0].set_ylabel('# of ASVs');
#         axs[0].ticklabel_format(axis='x',style='sci')
		axs[0].tick_params(axis='x', labelrotation = -90)

		# index is abundance
		# value is the number of times (i.e. number of ASVs) where that abundance occurred
		table = counts.value_counts().sort_index()
		table.rename('n_ASVs', inplace=True)
		table.index.rename('reads_per_ASV', inplace=True)
		table = table.reset_index()
		table['total_reads'] = table['reads_per_ASV'] * table['n_ASVs']
		table = table.sort_values(by='reads_per_ASV', ascending=False)
		table['cum_n_ASVs'] = table['n_ASVs'].cumsum()
		table['frac_n_ASVs'] = table['cum_n_ASVs'] / n_ASVs
		table['cum_total_reads'] = table['total_reads'].cumsum()
		table['frac_total_reads'] = table['cum_total_reads'] / n_reads

		table.plot('reads_per_ASV','frac_total_reads', logx=True, legend=False, ax=axs[1], label=name)
		axs[1].set_xlabel('ASV abundance ≥ $x$ (reads)');
		axs[1].set_ylabel('fraction of reads');

		table.plot('reads_per_ASV','frac_n_ASVs', logx=True, legend=False, ax=axs[2], label=name)
		axs[2].set_xlabel('ASV abundance ≥ $x$ (reads)');
		axs[2].set_ylabel('fraction of ASVs');

		tables.append(table)

	if len(libraries) > 1: axs[0].legend()
	plt.tight_layout()

	if len(libraries) == 1: summary.index = ['library']
	display(summary)

	return table

def abundance_prevalence_plot(ft, count_bins=None, prevalence_bins=None, n_counts=30, n_prevalences=20, log_counts=True, discrete=False,
							  weighted=False, kind='heatmap',interactive=True, ax=None):
	""" Plots what fraction of reads or features have abundance > x in > y samples (prevalence)
	"""

	from matplotlib.ticker import PercentFormatter
	def midpoints(x):
		return (x[1:] + x[:-1]) / 2

	if ax is None: 
		ax = plt.gca()
	
	ft_X = ft.X
	
	if count_bins is None:
		count_bins = np.logspace(0,6,num=n_counts+1)
	if prevalence_bins is None:
		prevalence_bins = np.linspace(0, 0.8*ft_X.shape[0], num=n_prevalences+1)

	counts = midpoints(count_bins)
	prevalences = midpoints(prevalence_bins)

	out = np.zeros((len(prevalences), len(counts)))
	if weighted:
		total = ft_X.sum()
	else:
		total = ft_X.shape[1]
	xx, yy = np.meshgrid(prevalence_bins, count_bins)

	out_rows = []
	for i,prevalence in enumerate(prevalences):
		for j,count in enumerate(counts):
			# find only features whose (abundance is > `count`) in > `prevalence` samples
			features = ((ft_X > count).sum(axis=0) > prevalence)
			if weighted:
				# calculate the percentage of reads preserved under this condition
				out[i,j] = ft_X[:,features.A1].sum() / total
			else:
				out[i,j] = features.sum() / total
			out_rows.append([prevalence, count, out[i,j]])
			
	out_df = pd.DataFrame(out_rows, columns=['prevalence','count','percent'])
				
	if kind == 'heatmap':
		pcm = ax.pcolormesh(count_bins, prevalence_bins, out)
		if log_counts:
			ax.set_xscale('log')
		ax.set_xlabel('count ≥ x')
		ax.set_ylabel('prevalence ≥ y')
		plt.colorbar(pcm, label=('% reads' if weighted else '% features'), format=PercentFormatter(xmax=1), ax=ax)
	elif kind == 'trace':
		if discrete:
			out_df['prevalence'] = out_df['prevalence'].map("{:.0f}".format)
		g = sns.lineplot(data=out_df, x='count',y='percent',hue='prevalence', palette='viridis', ax=ax)
		ax.set_ylabel('% reads' if weighted else '% features')
		ax.set_xlabel('count ≥ x')
		if log_counts:
			ax.set_xscale('log')
		ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))


def abundance_prevalence_plot_grid(ft, read_count_bins=np.linspace(0,1000,num=50)):

	fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
	abundance_prevalence_plot(ft, ax=axs[0,0], weighted=True, kind='heatmap')
	abundance_prevalence_plot(ft, ax=axs[0,1], weighted=True, kind='trace', discrete=False)

	abundance_prevalence_plot(ft, ax=axs[1,0], weighted=False, kind='heatmap', count_bins = read_count_bins, log_counts=False)
	abundance_prevalence_plot(ft, ax=axs[1,1], weighted=False, kind='trace', discrete=False)
	plt.tight_layout()

	return fig, axs