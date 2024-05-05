import matplotlib.pyplot as plt
import numpy as np

# Monkey patch the plot method of RocCurveDisplay to prevent it from trying to show the legend when empty (and thus generating lots of unsuppressable warning messages from matplotlib
# Not needed as of sklearn 1.3.0
# def _RocCurveDisplay_plot(self, ax=None, *, name=None, **kwargs):
# 	"""Plot visualization
# 	Extra keyword arguments will be passed to matplotlib's ``plot``.
# 	Parameters
# 	----------
# 	ax : matplotlib axes, default=None
# 		Axes object to plot on. If `None`, a new figure and axes is
# 		created.
# 	name : str, default=None
# 		Name of ROC Curve for labeling. If `None`, use `estimator_name` if
# 		not `None`, otherwise no labeling is shown.
# 	Returns
# 	-------
# 	display : :class:`~sklearn.metrics.plot.RocCurveDisplay`
# 		Object that stores computed values.
# 	"""

# 	name = self.estimator_name if name is None else name

# 	line_kwargs = {}
# 	if self.roc_auc is not None and name is not None:
# 		line_kwargs["label"] = f"{name} (AUC = {self.roc_auc:0.2f})"
# 	elif self.roc_auc is not None:
# 		line_kwargs["label"] = f"AUC = {self.roc_auc:0.2f}"
# 	elif name is not None:
# 		line_kwargs["label"] = name

# 	line_kwargs.update(**kwargs)

# 	import matplotlib.pyplot as plt

# 	if ax is None:
# 		fig, ax = plt.subplots()

# 	(self.line_,) = ax.plot(self.fpr, self.tpr, **line_kwargs)
# 	info_pos_label = (
# 		f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
# 	)

# 	xlabel = "False Positive Rate" + info_pos_label
# 	ylabel = "True Positive Rate" + info_pos_label
# 	ax.set(xlabel=xlabel, ylabel=ylabel)

# 	if "label" in line_kwargs and (len(line_kwargs['label']) > 0) and (line_kwargs['label'][0] != '_'):
# 		ax.legend(loc="lower right")

# 	self.ax_ = ax
# 	self.figure_ = ax.figure
# 	return self

# sklearn.metrics.RocCurveDisplay.plot = _RocCurveDisplay_plot


def plot_roc(cv_out, aggregate_folds=False, ax=None, title=None, verbose=False, legend=True, show_accuracy=False, fig_kw=None):
	"""Plot receiver-operator characteristic curve for cross-validated model using sklearn.metrics.RocCurveDisplay

	Modified from https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py
	© 2007 - 2024, scikit-learn developers (BSD License)
	
	Parameters
	----------
	cv_out : pd.DataFrame
		DataFrame with columns:
		- fold: name or number of the cross-validation fold
		- y: true label
		- p_y: target score, probability predicted for label, etc. 
	aggregate_folds : bool, optional
		True to combine all cross-validation folds into one dataset and plot a single ROC, False to plot each fold as a distinct curve, by default False
	ax : plt.Axes, optional
		Axes on which to plot the figure; if omitted, new Figure and Axes will be created, by default None
	title : str, optional
		Title for the plot, by default None
	verbose : bool, optional
		Print additional debugging information, by default False
	legend : bool or str, optional
		True to show a legend labeling the mean and std ROC curves and print the mean +/- std AUROC across CV folds;
		'text' to only print the mean +/- std AUROC across folds and, if `show_accuracy` = True, the mean accuracy;
		False to show neither; by default True
	show_accuracy : bool, optional
		Print the mean accuracy in the corner of the plot, by default False
	fig_kw : dict, optional
		kwargs to pass to plt.Figure, by default None

	"""
	from sklearn.metrics import RocCurveDisplay, auc
	from scipy.stats import ttest_1samp
	from scipy.stats import mannwhitneyu

	if ax is None:
		if fig_kw is None:
			fig_kw = dict()
		fig, ax = plt.subplots()
	else:
		fig = None

	mean_fpr = np.linspace(0, 1, 100)
	def plot_fold(df, alpha=0.2, label='', **kwargs):

		# with warnings.catch_warnings():
		#     warnings.simplefilter("ignore")
		if label is not None:
			kwargs['label'] = label

		viz = RocCurveDisplay.from_predictions(
			df.y,
			df.p_y, 
			alpha=alpha,
			ax=ax,
			# name="_ROC fold {}".format(i),
			# legend=False,
			**kwargs
		)
		interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		return {'tpr': interp_tpr, 'auc':viz.roc_auc } #, 'p': res.pvalue, 'U': res.statistic}
		# tprs.append(interp_tpr)
		# aucs.append(viz.roc_auc)

	if not aggregate_folds:
		tprs = []
		aucs = []

		curves = cv_out.groupby('fold').apply(plot_fold)
		tprs = [curve['tpr'] for curve in curves.values]
		aucs = [curve['auc'] for curve in curves.values]

		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)
		if verbose:
			print(f"AUC = {mean_auc:0.2f} ± {std_auc:0.2f}")
		ax.plot(
			mean_fpr,
			mean_tpr,
			color="b",
			label="Mean ROC\n(AUC = %0.2f $\\pm$ %0.2f)" % (mean_auc, std_auc),
			lw=2,
			alpha=0.8,
		)

		std_tpr = np.std(tprs, axis=0)
		tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		ax.fill_between(
			mean_fpr,
			tprs_lower,
			tprs_upper,
			color="grey",
			alpha=0.2,
			label=r"$\pm$ 1 std. dev.",
		)

		ax.set(
			xlim=[-0.05, 1.05],
			ylim=[-0.05, 1.05]
		)

	else:
		out = plot_fold(cv_out, label=None, alpha=1)

	pvalue = ttest_1samp(cv_out.score, popmean=0.5, alternative='greater').pvalue

	# calculate P(AUROC > 0.5); equivalent to asking whether the median predicted probability
	# for positive samples (y=1) is greater than the median predicted probability for 
	# negative samples (y=0). This is the question asked by the Mann-Whitney U test, where 
	# AUROC = U1/n1*n2 so we can use `mannwhitneyu` to calculate the p-value
	(pos_idx,) = np.nonzero(cv_out.y.values)
	(neg_idx,) = np.nonzero(~cv_out.y.values)

	res = mannwhitneyu(cv_out.p_y.iloc[pos_idx], cv_out.p_y.iloc[neg_idx], alternative='greater')
	pvalue_auc = res.pvalue

	mean_acc = cv_out['score'].mean()
	std_acc = cv_out['score'].std()
		
	if verbose:
		print(f"  AUC < 0.5, p = {pvalue_auc:0.1g}")
		print(f"Accuracy = {cv_out['score'].mean():0.2g} ± {cv_out['score'].std():0.2g}")
		print(f"  accuracy < 0.5, p = {pvalue:0.1g}")
		

	ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
	ax.set_xlabel("False positive rate")
	ax.set_ylabel("True positive rate")

	label = (f"AUC = {mean_auc:0.2f} ± {std_auc:0.2f}\n"
		 f"p = {pvalue_auc:0.1g}")
	if show_accuracy:
		label += ("\n"
				  f"accuracy = {mean_acc:0.2g} ± {std_acc:0.2f}\n" 
				  f"p = {pvalue:0.1g}")
	
	if title is not None:
		ax.set_title(title)
		
	if legend==True:
		ax.legend(loc="lower right")
	elif legend=='text':
		ax.get_legend().remove()
		ax.text(0.95, 0.05, 
				label, horizontalalignment='right', verticalalignment='bottom', 
				transform = ax.transAxes, bbox=dict(facecolor='white', alpha=0.5, linestyle=''))
	elif legend=='both':
		ax.legend(loc="lower right")
		
		label = f"        AUC < 0.5, p = {pvalue_auc:0.1g}"
		if show_accuracy:
			label += f"\n accuracy < 0.5, p = {pvalue:0.1g}"
		ax.text(0.05, 0.95, 
				label,
				horizontalalignment='left', verticalalignment='top', transform = ax.transAxes)
	
	if fig is not None:
		fig.show()

def plot_feature_importance_summary(out, ag, design, n_features=None):
	import seaborn as sns

	feature_names = design.var_names.values
	
	estimators = out['cv'].groupby('fold')['estimator'].first()
	feature_importances = estimators.apply(lambda model: pd.Series(model.feature_importances_, index=feature_names))
	n_features = (feature_importances > 0).sum(axis=1) #estimators.apply(lambda model: model.feature_importances)
	
	feature_ranks = feature_importances.rank(axis=1, ascending=False)
	
	mfr = feature_ranks.agg(['mean', 'std']).T.sort_values('mean', ascending=True)
	mfi = feature_importances.agg(['mean', 'std']).T.sort_values('mean', ascending=False)
	
	fig, axs = plt.subplots(nrows=3, figsize=(8, 12))
	
	sns.histplot(data=n_features, ax=axs[0])
	axs[0].set_xlabel('# features in model')
	axs[0].set_ylabel('# folds')
	
	if n_features is None:
		from kneed import KneeLocator
		n_features = KneeLocator(x=range(len(mfi)), y=mfi['mean'], curve='convex', direction='decreasing').knee + 10
		
	# xlim = sum(mfi['mean'] > 0)
	
	top_feature_ranks       = mfr.iloc[:n_features,:]
	top_feature_importances = mfi.iloc[:n_features,:]
	
	
	
	axs[1].errorbar(x=range(len(top_feature_ranks['mean'])), y=top_feature_ranks['mean'], yerr=top_feature_ranks['std'])
	
	axs[1].set_xticks(list(range(len(top_feature_ranks['mean']))))
	axs[1].set_xticklabels(top_feature_ranks.index.values, rotation=-45, horizontalalignment='left')
	axs[1].set_ylabel('mean feature rank')
	
	axs[2].errorbar(x=range(len(top_feature_importances)), y=top_feature_importances['mean'], yerr=top_feature_importances['std'])
	axs[2].set_xticks(list(range(len(top_feature_importances))))
	axs[2].set_xticklabels(top_feature_importances.index.values, rotation=-45, horizontalalignment='left')
	axs[2].set_ylabel('mean feature importance')
	
	return mfr, mfi

def plot_cv_indices(cv, X, y, group=None, ax=None, cmap_data = plt.cm.Paired, cmap_cv = plt.cm.coolwarm, lw = 10):
	"""Plot the indices of a cross-validation object.
	https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html

	© 2007 - 2024, scikit-learn developers (BSD License)
	"""
	import numpy as np

	n_splits = cv.get_n_splits(X, y, group)
	n_y = len(y)
	if ax is None:
		fig, ax = plt.subplots(figsize=(len(y)*0.2, n_splits*0.2))

	# Generate the training/testing visualizations for each CV split
	for ii, (tr, tt) in enumerate(cv.split(X=X, y=y, groups=group)):
		# Fill in indices with the training/test groups
		indices = np.array([np.nan] * len(X))
		indices[tt] = 1
		indices[tr] = 0

		# Visualize the results
		ax.scatter(
			range(len(indices)),
			[ii + 0.5] * len(indices),
			c=indices,
			marker="_",
			lw=lw,
			cmap=cmap_cv,
			vmin=-0.2,
			vmax=1.2,
		)
		pos = sum(y[tr])
		neg = len(y[tr])-pos
		
		ax.text(x=n_y, y=ii+0.5, s=f"{pos}+ / {neg}-")

	# Plot the data classes and groups at the end
	ax.scatter(
		range(len(X)), [ii + 1.5] * len(X), c=y, marker="_", lw=lw, cmap=cmap_data
	)

	if group is not None:
		ax.scatter(
			range(len(X)), [ii + 2.5] * len(X), c=group, marker="_", lw=lw, cmap=cmap_data
		)

	# Formatting
	yticklabels = list(range(n_splits)) + ["class"]
	if group is not None:
		yticklabels += ["group"]

	ax.set(
		yticks=np.arange(n_splits + 1 + (group is not None)) + 0.5,
		yticklabels=yticklabels,
		xlabel="Sample index",
		ylabel="CV iteration",
		ylim=[n_splits + (1.1 * (1 + (group is not None))), -0.2],
		xlim=[0, len(y)],
	)
	#ax.set_title("{}".format(type(cv).__name__), fontsize=15)
	ax.set_title(repr(cv))
	return ax
