from .utils import *
from ..ft import fortify

import pandas as pd

class ExperimentVisualizer():

	def __init__(self, expt):
		self.ex = expt

	def tag(self, id, space='cdr3'):
		id = self.ex.find_feature(id, space=space)
		return tag(id=id, library=self.ex.find_library(id, space=space))

	def find_and_summarize(self, CDR3ID, **kwargs):
		return self.summarize_cdr3(CDR3ID, **kwargs)

	# def summarize_top_samples(features, ft=None, fd=None, relative=False, table=True, plot=True):
	# 	fd = fortify_feature_data(ft, fd, features=features, relative=relative, obs=True)
	#
	# 	if table:
	# 		display(top_samples_table(features, fd=fd, relative=relative))
	# 	if plot:
	# 		sample_abundance_plot(features, fd=fd, relative=relative).draw(show=True)

	def summarize_top_samples(self, features, space='cdr3', relative=False, table=True, plot=True):
		from .asv import sample_abundance_plot, top_samples_table
		from ..asvs import get_identifier

		space = space.lower()

		if relative:
			ft = self.ex.rfts[space]
		else:
			ft = self.ex.fts[space]

		fd = fortify(ft[:,features], sample_col='ID', feature_col='feature', obs=True)
		# fd['feature'] = fd[get_identifier(space)]

		if table:
			display(top_samples_table(features, fd=fd, relative=relative))
		if plot:
			sample_abundance_plot(features, fd=fd, relative=relative).draw(show=True)


	def summarize_clonotypes(self, cdr3s, sort_by=['CDR3'], ascending=True, other_columns=[], show_query=False, table=True, plot=True):
		from IPython.display import display
		from .logo import cdr_seq_logo_plot

		# find AA sequences contributing to the clonotype
		ft = self.ex.project(cdr3s, from_space='CDR3', to_space='AA', show_query=show_query)
		fda = fortify(ft, var=True)

		if plot:
			# plot sequence logo
			cdr_seq_logo_plot(fda)

		if table:
			fda_gb = fda.groupby('CDR3ID')

			# make summary table of # CDR1, # CDR2, etc.
			summary_df = pd.concat([
				fda_gb[['library','CDR3'] + other_columns].first(),
				fda_gb[['clonotype','CDR1','CDR2']].nunique().rename(columns=lambda c: f"# {c}s"),
				fda_gb['abundance'].sum().astype(int),
			],axis=1,join='inner')
			display(rich_table(summary_df))
			return summary_df

	def summarize_similar_features(self, features, space='CDR3', threshold=0.85):
		from ..asvs import get_identifier
		from IPython.display import display

		fd = self.ex.find_similar_features(
			self.ex.fts[space.lower()].var.loc[features,space],
			space=space.lower()
		)

		identifier = get_identifier(space)

		# filter out hits for the query itself
		fd = fd.loc[~fd[f'search_{identifier}'].isin(features),:]

		# filter out hits with pident < threshold
		fd = fd.loc[fd['pident'] >= threshold,:]

		if len(fd) > 0:
			title = f"{(fd.groupby(f'search_{space}ID')['pident'].first() > threshold).sum()} very similar {space}s (>{threshold:.0%} identity)"
			display(display_accordion(
				rich_table(
					fd[[f'search_{identifier}','library','CDR3','pident','reads','nsamples']]
						.rename(columns={f'search_{identifier}':identifier}),
					sort_by=['reads'],ascending=False),
				title
			))
		else:
			print(f"No very similar {space}s (>{threshold:.0%} identity)")
		return fd

	def summarize_cdr3(self, cdr3id,
					   cdrs=['CDR1','CDR2','CDR3'], other_columns=[],
					   plot_seq_logo=True, list_clonotypes=True, plot_samples=True, list_samples=True, list_similar=True,
					   show_query=False,
					   header=True,
					  **kwargs):

		from IPython.display import Markdown, display
		cdr3id = self.ex.find_cdr3(cdr3id, single=True)

		if header:
			if isinstance(header, bool) or not isinstance(header, int):
				header = '#####'
			else:
				header = ('#' * header)
			display(Markdown(f'{header} `{cdr3id[0:6]}`'))

		cdr3_tag(cdr3id)

		cdr3ids = [cdr3id]

		if plot_seq_logo or list_clonotypes:
			self.summarize_clonotypes(cdr3ids, plot=plot_seq_logo, table=list_clonotypes, show_query=show_query, **kwargs)

		if plot_samples or list_samples:
			self.summarize_top_samples(cdr3ids, relative=True, space='CDR3', plot=plot_samples, table=list_samples)

		if list_similar:
			self.summarize_similar_features(cdr3ids, space='CDR3')

	def summarize_cdr3_and_similar(self, cdr3id,
					   cdrs=['CDR1','CDR2','CDR3'], other_columns=[],  library=None,
					   plot_seq_logo=True, list_clonotypes=True, plot_samples=False, list_samples=False, threshold=0.85,
					   show_query=False,
					   header=True,
					  **kwargs):

		from IPython.display import Markdown, display
		from ..asvs import get_identifier

		# might make this generic eventually
		space = 'CDR3'


		space = space.lower()
		cdr3id = self.ex.find_feature(cdr3id, single=True, space=space)
		identifier = get_identifier(space)

		fd = self.ex.find_similar_features(
			{cdr3id: self.ex.fts[space].var.loc[cdr3id,space]},
			space='cdr3'
		)
		if library is not None:
			fd = fd.loc[fd['library'] == library,:]
		fd = fd.loc[fd['pident'] > threshold,:]
		features = fd[f'search_{identifier}'].unique()

		title = (f'`{cdr3id[0:6]}` and {len(features)} similar {space}s ' + (f'from library {library}' if library is not None else ''))
		if header:
			display(Markdown(f'#### {title}'))
		else:
			print(title)

		cdr3ids = features

		if plot_seq_logo or list_clonotypes:
			self.summarize_clonotypes(cdr3ids, plot=plot_seq_logo, table=list_clonotypes, show_query=show_query, **kwargs)

		if plot_samples or list_samples:
			self.summarize_top_samples(cdr3ids, space=space, relative=True, plot=plot_samples, table=list_samples)


	# def summarize_partition_logo(ft, sample_partitions, cols=['CDR1','CDR2','CDR3']):
	#     """show a sequence logo for each of a set of partitions of a given feature table
	#     """
	#     if not isinstance(sample_partitions, dict):
	#         show_titles = False
	#         iterator = enumerate(sample_partitions)
	#     else:
	#         show_titles = True
	#         iterator = sample_partitions.items()
	#
	#     for title, partition in iterator:
	#         fda = fortify(ft[partition,:], var=True)
	#         cdr_seq_logo_plot(fda, cols, title=(title if show_titles else None))

	def summarize_cdr3_partition(self, cdr3s, sample_partitions, cdrs=['CDR1','CDR2','CDR3'], plot_logo=True, plot_samples=True, **kwargs):
		"""show a sequence logo for CDR1--3 of a particular CDR3 across several partitions of the feature table, e.g. for phenotype +/-/? samples"""
		from .logo import partition_logo_plot
		from .asv import top_samples_table, sample_abundance_plot
		from plotnine import facet_grid

		logos = None
		if plot_logo:
			ft = self.ex.project(cdr3s, from_space='CDR3', to_space='AA')
			logos = partition_logo_plot(ft, sample_partitions, cols=cdrs, **kwargs)

		sample_plots = None
		if plot_samples:
			fd = fortify(ft=self.ex.rfts['cdr3'][:,cdr3s], obs=True, relative=False, sample_col='ID')
			partition_dict = {}
			for (partition, values) in sample_partitions.items():
				for i in values:
					partition_dict[i] = partition

			fd['partition'] = fd['ID'].map(partition_dict)
			sample_plots = sample_abundance_plot(cdr3s,
										  fd=fd.dropna(subset=['partition']),
										  relative='pre') + facet_grid('feature ~ partition')
		return (logos, sample_plots)
