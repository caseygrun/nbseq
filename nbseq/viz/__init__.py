from .utils import *
from ..ft import fortify

import pandas as pd

class ExperimentVisualizer():
	
	global_query = "kind == '+' & io == 'i'"
	_dfs = None
	_df_enr = None
	_enr = None

	def __init__(self, expt):
		self.ex = expt

	def reset(self):
		self.ex._viz = None

	def tag(self, id, space='cdr3', **kwargs):
		ids = self.ex.find_feature(id, space=space)
		if len(ids) > 1:
			print(f"Warning: multiple tags match query `{id}` in space {space}: {ids}!")
		else:
			return tag(id=ids[0], library=self.ex.find_library(ids[0], space=space), **kwargs)

	@property
	def dfs(self):
		from ..utils import lazydict
		from functools import partial
		if self._dfs is None:
			self._dfs = lazydict({
				k: partial(self.fortify, space=k, obs=True)
				for k in self.ex.spaces
			})
		return self._dfs

	def fortify(self, space='cdr3', obs=True, enr=None, **kwargs):
		from ..ft import fortify, query
		from ..asvs import get_identifier
		space = space.lower()
		identifier = get_identifier(space)
		if enr is None:
			enr = (space == 'cdr3')

		ft = self.ex.rfts[space]
		ft = query(ft, self.global_query, axis='obs')
		df = fortify(ft, obs=(obs or enr), **kwargs)

		
		if enr:
			df_enr = self.enr(space)
			cols_to_use = df_enr.columns.difference(df.columns)
			return df.join(
				df_enr.set_index(['name', identifier])[cols_to_use],
				on=['name',identifier]
			)
		return df

	def enr(self, space):
		space = space.lower()
		if self._enr is None:
			self._enr = {}
		
		if space in self._enr:
			return self._enr[space]
		else:
			_enr = self.ex.enr(space,'df', add_start_end=True, add_log=True, add_rank=True, comparison=('start','end'))

			# add p-value from experiment model if available
			if space in self.ex.enr_models:
				from ..select import enrichment_pvalues
				_enr = enrichment_pvalues(_enr, self.ex.enr_models[space], space=space, abundance_col='end')

			# join to selection metadata
			sel_df = self.ex.selection_metadata
			sel_df = shorten_descriptions(sel_df)

			_enr = _enr.join(sel_df, on='name')

			self._enr[space] = _enr

		return self._enr[space]

	def selection_metadata(self):
		if self._sel_metadata is None:
			self._sel_metadata = self.ex.selection_metadata
		return self._sel_metadata

	def find_and_summarize(self, CDR3ID, **kwargs):
		return self.summarize_cdr3(CDR3ID, **kwargs)

	# def summarize_top_samples(features, ft=None, fd=None, relative=False, table=True, plot=True):
	# 	fd = fortify_feature_data(ft, fd, features=features, relative=relative, obs=True)
	#
	# 	if table:
	# 		display(top_samples_table(features, fd=fd, relative=relative))
	# 	if plot:
	# 		sample_abundance_plot(features, fd=fd, relative=relative).draw(show=True)
	
	def abundance_trace_plot(self, feature, space='cdr3', phenotype=None, facet=None, limits=None, title=None):
		from ..ft import fortify
		from ..asvs import get_identifier
		from .utils import hash_to_mn_short, pretty_hex
		import numpy as np
		import plotnine as p9
		
		identifier = get_identifier(space)
		ft = self.ex.rfts[space][:,feature]
		df = fortify(ft, obs=True)

		_aes = dict(x='r',y='abundance')
		_aes_label = dict(label='label')
		if phenotype is not None:
			df[phenotype] = df[phenotype].replace({-1: '-', 0: '~', float('nan'):'?', 1:'+'})
			_aes['color'] = phenotype
			_aes_label['fill'] =phenotype

		
		df = df.join(
			df.groupby('name')['r'].max().rename('last_r'),
			on='name'
		)
		df['label'] = np.where(df['r'] == df['last_r'], df['name'],'')
		df['r'] = df['r'].astype(str)
		df['facet'] = df[identifier] + '\n' + df[identifier].apply(hash_to_mn_short)
		
		gg = (p9.ggplot(df, p9.aes(**_aes))
				+ p9.geom_line(p9.aes(group='name'))
				+ p9.geom_point()
				# + scale_color_md5(df['feature'].unique(), special_values=special_values)
				+ p9.geom_label(p9.aes(**_aes_label), ha='left', color=('white' if phenotype is not None else 'black'), position=p9.position_nudge(x=0.05))
				+ p9.scale_x_discrete(name='Round', expand=(0,0,0,1), drop=True)
				+ p9.scale_y_continuous(name='Abundance')
			)
		if isinstance(feature, str):
			gg += p9.ggtitle(f"{space.upper()} {pretty_hex(feature)}\n({hash_to_mn_short(feature)})")
		if phenotype is not None:
			gg = gg + p9.scale_color_manual(
				breaks=['?', '-', '~', '+'], values=['#bab0ac', '#4c78a8', '#72b7b2', '#f58518'],
				aesthetics=['color','fill']
			)

		if facet is not None:
			gg = gg + p9.facet_wrap(facet)
		if limits is not None:
			gg = gg + p9.coord_cartesian(ylim=limits)
		if title is not None:
			gg = gg + p9.ggtitle(title)
		return gg

	def plot_selections_for_feature(self, feature, space='cdr3', phenotype=None):
		import altair as alt
		from .utils import hash_to_mn_short, pretty_hex

		sample_df = self.dfs[space].rename(columns={'p_value':'enr_p_value'})
		ex = self.ex

		df_abd = sample_df.query(f"CDR3ID == '{feature}' & kind == '+' & io == 'i'")[
                ['CDR3ID', 'r', 'name', 'desc_short', 'abundance', 'enrichment', 'start', 'end', 'enr_p_value'] + list(ex.pheno_names)]


		brush = alt.selection_interval(name='brush',
			on="[mousedown[event.altKey], mouseup] > mousemove",
			translate="[mousedown[event.altKey], mouseup] > mousemove")
		dropdown = alt.binding_select(
			options=['CDR3ID'] + list(ex.pheno_names),
			name='Color by Phenotype'
		)
		phenotype_selector = alt.param(
			name='phenotype_selector', 
			value=(phenotype if phenotype is not None else 'CDR3ID'), 
			bind=dropdown)


		_phenotype_scale = alt.Scale(
			domain=['NaN', -1, 0, 1], range=['#bab0ac', '#4c78a8', '#72b7b2', '#f58518'])
		chart_abd_trace = (
			alt.Chart(df_abd).mark_line(point=True).encode(
				x=alt.X('r:O'),
				y=alt.Y('abundance:Q'),
				color=alt.condition(
					"phenotype_selector != 'CDR3ID'",
					alt.Color('grouping_color:N', title='Phenotype', scale=_phenotype_scale),
					alt.value(hash_to_color(feature)),
					# 'grouping_color:N'
				),
				detail='name',
				tooltip=['name', 'desc_short', 'grouping_color:N']
			)
			.properties(height=200, width=200)
			.transform_filter(brush)
		)#.interactive()
		chart_abd = chart_abd_trace

		
		_opacity = alt.condition(brush, alt.value(1), alt.value(0.05))
		_size = alt.condition(brush, alt.value(100), alt.value(20))

		_tooltips = ['name','desc_short','grouping_color:N',
					alt.Tooltip('start', format=".2g"),
					alt.Tooltip('end', format=".2g"),
					alt.Tooltip('enrichment', format=".2g"),
					]
		if 'enr_p_value' in df_abd.columns:
			_tooltips.append(alt.Tooltip('enr_p_value', format=".2g"))
		if 'log_enrichment' in df_abd.columns:
			_tooltips.append(alt.Tooltip('log_enrichment', format=".2g"))
		if 'enrichment_rank' in df_abd.columns:
			_tooltips.append(alt.Tooltip('enrichment_rank', format=".0f"))

		_scales = alt.selection_interval(bind='scales',
			on="[mousedown[!event.altKey], mouseup] > mousemove",
			translate="[mousedown[!event.altKey], mouseup] > mousemove"
		)

		enrabdplot = alt.Chart(df_abd).mark_point().encode(
			x=alt.X('end', title='end abundance',
					scale=alt.Scale(type='log')),
			# x='end',
			# y='log_enrichment',
			y=alt.Y('enrichment', scale=alt.Scale(type='log')),
			tooltip=_tooltips,
			color=alt.condition(
				"phenotype_selector == 'CDR3ID'",
				alt.value(hash_to_color(feature)),
				alt.Color('grouping_color:N', title='Phenotype',
							scale=_phenotype_scale)
				# 'grouping_color:N'
			),
			opacity=_opacity,
			size=_size,
		).properties(height=200, width=200).add_params(_scales)  # .interactive()

		# title = f"{space} {hash_to_mn_short(feature)} {pretty_hex(feature)}"

		chart_abd_enrs = ((chart_abd_trace | enrabdplot.add_params(brush))
                .transform_calculate(
                    # 'datum[phenotype_selector]'
                    grouping_color='datum[phenotype_selector]'
                )
				.add_params(phenotype_selector)
            ).properties(
				title=alt.TitleParams(
					text=hash_to_mn_short(feature),
					subtitle=f"{space} {feature}",
					color=hash_to_color(feature),
					anchor='start'
				)
			)
		return chart_abd_enrs

	def summarize_top_selections(self, feature, space='cdr3'):
		from .asv import top_selections_table
		space = space.lower()
		return top_selections_table(feature, self.enr(space))

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
			from IPython.display import display
			display(top_samples_table(features, fd=fd, relative=relative))
		if plot:
			sample_abundance_plot(features, fd=fd, relative=relative).draw(show=True)


	def summarize_clonotypes(self, cdr3s, sort_by=['CDR3'], ascending=True, other_columns=[], show_query=False, table=True, plot=True, **kwargs):
		from IPython.display import display
		from .logo import cdr_seq_logo_plot

		# find AA sequences contributing to the clonotype
		ft = self.ex.project(cdr3s, from_space='CDR3', to_space='AA', show_query=show_query)
		fda = fortify(ft, var=True)

		out = []
		if plot:
			# plot sequence logo
			out.append(cdr_seq_logo_plot(fda))

		if table:
			fda_gb = fda.groupby('CDR3ID')

			# make summary table of # CDR1, # CDR2, etc.
			summary_df = pd.concat([
				fda_gb[['library','CDR3'] + other_columns].first(),
				fda_gb[['clonotype','CDR1','CDR2']].nunique().rename(columns=lambda c: f"# {c}s"),
				fda_gb['abundance'].sum().astype(int),
			],axis=1,join='inner')

			out.append(rich_table(summary_df))
		return tuple(out)

	def summarize_similar_features(self, features, space='CDR3', threshold=0.85, include_self=True, **kwargs):
		from ..asvs import get_identifier
		from IPython.display import display

		ft = self.ex.fts[space.lower()]
		seq_col = space.upper()
		if seq_col not in ft.var.columns:
			raise ValueError(
				f"Feature table for space {space} does not have sequence data in column `{seq_col}` of feature metadata (`ft.var`). "
				f"Feature data columns: {ft.var.columns.values}"
			)
		fd = self.ex.find_similar_features(
			self.ex.fts[space.lower()].var.loc[features,seq_col],
			space=space.lower(),
			**kwargs
		)

		identifier = get_identifier(space)

		if not include_self:
			# filter out hits for the query itself
			fd = fd.loc[~fd[f'search_{identifier}'].isin(features),:]

		# filter out hits with pident < threshold
		total_hits = len(fd) - int(include_self)
		fd = fd.loc[fd['pident'] >= threshold,:]

		if len(fd) > 0:
			title = f"{(fd.groupby(f'search_{identifier}')['pident'].first() > threshold).sum()} very similar {space}s (>{threshold:.0%} identity)"
		else:
			title = (f"No very similar {space}s (>{threshold:.0%} identity)")
			if total_hits == 0 and 's' not in kwargs:
				print("Consider increasing sensitivity of kmer prefiltering by setting s=7 to s=8.5")

		fd['mm'] = fd[f'search_{identifier}'].apply(hash_to_mn_short)
		columns = [f'search_{identifier}','mm','library',seq_col,'pident']
		sort_by = ['pident']
		
		if len(features) > 1:
			columns = [f'query_{identifier}'] + columns
			# sort_by = [f'query_{identifier}', 'pident']
			sort_by = ['query', 'pident']

		# handle case where add_feature_data hasn't been called
		if 'nsamples' in fd.columns: 
			columns += ['nsamples'] 
			sort_by = ['nsamples']
		if 'reads' in fd.columns: 
			columns += ['reads'] 
			sort_by = ['reads']
		_df = (fd[columns]
				.rename(columns={f'search_{identifier}':identifier, f'query_{identifier}':'query'}))
		if len(features) > 1:
			_df = _df.set_index(['query',identifier])

		table = rich_table(
			_df,
			sort_by=sort_by, ascending=False, caption=title,
			format=dict(pident='{:.1%}',reads='{:.0f}', nsamples='{:.0f}')
		)
			
		return table

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

	def selection_phenotype_grid(self, query="1"):
		from ..pheno import selection_phenotype_grid

		return selection_phenotype_grid(
			self.ex.obs.query(query), 
			self.ex.ags.reset_index()
		)

	def plot_selection_ag_matrix(self):
		from ..pheno import plot_ag_matrix
		ex = self.ex
		plot_ag_matrix([ex.selection_metadata[ex.ag_names]], ex.selection_metadata)