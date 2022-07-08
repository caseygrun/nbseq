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
	weighted=True):

	if libraries is None:
		libraries = pd.unique(segment_lens['library'])

	fig, axs = plt.subplots(nrows=1, ncols=len(segments), figsize=(12,4), sharey=True)
	for library in libraries:
		if 'library' in segment_lens:
			segment_lens_lib = segment_lens.query(f"library == '{library}'").dropna()
			weights = segment_lens_lib['abundance'] if weighted else None
			segment_lens_lib.hist(column=segments,
							  weights = weights,
							  label=library, alpha=0.5,
							  ax=axs);
	[ax.set_xlabel('length (nt)') for ax in axs.flat];
	if weighted: axs.flat[0].set_ylabel('abundance (reads)');
	else: axs.flat[0].set_ylabel('# of ASVs');

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

plot_asv_stats(counts);
