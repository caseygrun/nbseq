
import matplotlib.colors as mpl_colors
import colorsys

from .utils import *


def cdr3_tag(CDR3ID,ft=None):
	from nbseq.viz.asv import hash_to_color, contrasting_color, pretty_hex
	from IPython.display import HTML
	cdr3id = find_cdr3(CDR3ID, single=True)
	bg = hash_to_color(cdr3id)
	fg = contrasting_color(bg)
	if ft is not None:
		library = ft.var.loc[cdr3id,'library']
	else:
		library = '??'
	display(HTML(f"<span style='background-color:{bg};color:{fg};padding:0.3em;'>{pretty_hex(cdr3id)}</span><small style='font-size:8pt;'><code>{cdr3id}</code> ({library})</small>"))


def top_samples_table(features, ft=None, fd=None, relative=False, index=['feature','method','description','round','replicate','sample']):
	import pandas as pd, numpy as np
	from ..utils import intersection_ordered

	if relative:
		if ft is not None:
			ft = nbseq.ft.to_relative(ft)
		max_abundance = 1
	else:
		if ft is not None:
			max_abundance = ft.X.max()
		else:
			max_abundance = fd['abundance'].max()
	if fd is None and ft is not None:
		fd = nbseq.ft.fortify(ft[:,features], sample_col='ID', obs=True)
	elif fd is not None:
		fd = fd.loc[fd['feature'].isin(set(features)),:]
	else:
		raise ValueError("Must specify either ft or fd")


	index = intersection_ordered(index, list(fd.columns))
	df2 = pd.merge(
		fd.groupby('feature')['abundance'].agg(feature_mean='mean',feature_std='std'),
		fd,
		left_index=True, right_on='feature'
	)
	df2['z'] = (df2['abundance'] - df2['feature_mean']) / df2['feature_std']
	df2['passes'] = np.abs(df2['z']) > 1
	df2 = df2.query('passes')
	title = f"{len(df2)} / {len(fd)} samples have abundance > μ+σ^2:"
	return display_accordion(
		(df2
			 .set_index(index)
			 .sort_index()[['abundance','z']]
			 .style.bar(subset=['abundance'],vmax=max_abundance)
		),
		title
	)

def fortify_feature_data(ft=None, fd=None, relative=False, features=None, var=False, obs=False, sample_col='ID', **kwargs):
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

def sample_abundance_plot(features, ft=None, fd=None, relative=False):
	"""plots absolute or relative abundance of a set of features across all samples
	"""
	import plotnine as p9
	import pandas as pd
	import numpy as np

	fd = fortify_feature_data(ft=ft, fd=fd, relative=(bool(relative) if relative != 'pre' else False), obs=True)

	relative = bool(relative)
	gg = (p9.ggplot(fd, p9.aes(x='ID', y='abundance',color='pd.Categorical(r, sorted(pd.unique(r)))',
			label='np.where(abundance > np.mean(abundance)+np.std(abundance), description, "")')) +
		  p9.geom_point() +
		  p9.geom_label() + #geom_text() +
		  #geom_text(adjust_text=dict(arrowprops=dict(arrowstyle='->'))) +
		  p9.scale_y_continuous(name=('relative abundance' if relative else 'reads')) +
		  p9.scale_x_discrete(name='Sample-Round') +
		  p9.scale_color_cmap_d('viridis') +
		  p9.facet_grid('feature ~ .') +
		  p9.theme(axis_text_x=p9.element_blank())
		 )
	return gg
