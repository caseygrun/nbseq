import sklearn.model_selection
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import auc
from sklearn.base import clone
from joblib import Parallel, delayed

import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


from sklearn.model_selection import LeavePOut, KFold, StratifiedKFold, RepeatedStratifiedKFold, GridSearchCV
from itertools import chain, combinations
from sklearn.utils.validation import _num_samples

from .utils import ilen



class LeavePairOut(LeavePOut):
	""" generate cross-validation indices where each test case has p distinct labels. For p=2, this is "leave-pair-out" CV.
	"""
	def _iter_test_indices(self, X, y, groups=None):
		# labels = np.unique(y)
		# self.p = len(labels)

		n_samples = _num_samples(X)
		if n_samples <= self.p:
			raise ValueError(
				"p={} must be strictly less than the number of samples={}".format(
					self.p, n_samples
				)
			)
		for combination in combinations(range(n_samples), self.p):
			combination = np.array(combination)
			if np.unique(y[combination]).size == self.p:
				yield combination

	def get_n_splits(self, X, y=None, groups=None):
		"""Returns the number of splitting iterations in the cross-validator
		Parameters
		----------
		X : array-like of shape (n_samples, n_features)
			Training data, where `n_samples` is the number of samples
			and `n_features` is the number of features.
		y : array-like of shape (n_samples, )
			Label data.
		groups : object
			Always ignored, exists for compatibility.
		"""
		if X is None:
			raise ValueError("The 'X' parameter should not be None.")
		if y is None:
			raise ValueError("The 'y' parameter should not be None.")
		return ilen(self._iter_test_indices(X,y))



import sklearn.metrics

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


def plot_cv_indices(cv, X, y, group=None, ax=None, cmap_data = plt.cm.Paired, cmap_cv = plt.cm.coolwarm, lw = 10):
    """Plot the indices of a cross-validation object.
    https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html
    """

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


def synthesize_input_controls(ft, ag_names=[], input_library_query="expt == '027i.lib'", verbose=True):
	from .ft import query_ids

	ft = ft.copy()
	ids = query_ids(ft, input_library_query, axis='obs')

	n_rounds = (max(ft.obs['r'])-1)
	n_controls = len(ids) // n_rounds
	if verbose:
		print(f"Samples used to create controls: {ids}")
		print(f"Creating {n_controls} synthetic controls with {n_rounds} each...")

	i = 0
	for control in range(n_controls):
		for r in range(2, n_rounds+2):
			_id = ids[i]
			ft.obs.loc[_id,'name'] = f'control_{control}'
			ft.obs.loc[_id,'round'] = f'R{r}i'
			ft.obs.loc[_id,'r'] = r-1

			# set antigen values to 0 for these fake control samples
			ft.obs.loc[_id, ag_names] = 0
			i += 1

	ids_used = ids[:i]
	drop = ids[i:]

	# display(ft.obs.loc[ids,['name','round','r'] + list(ag_names)])
	# print(f"Dropping {sum(ft.obs.index.isin(drop))} samples")
	# return ft[~ft.obs.index.isin(drop),:]
	print(f"Dropping {len(drop)} samples {drop}")

	ftf = ft[~ft.obs.index.isin(drop),:]
	print(f"Final feature table shape: {ftf.shape}")

	return ftf


# from sklearn.utils.metaestimators import _safe_split
#
# def _fit_and_predict_proba(
# 	estimator,
# 	X,
# 	y,
# 	train,
# 	test,
# 	parameters,
# 	candidate_id,
# 	verbose, **kwargs):
#
# 	# todo: clone?
# 	estimator = estimator.set_params(**parameters)
#
# 	X_train, y_train = _safe_split(estimator, X, y, train)
# 	X_test, y_test = _safe_split(estimator, X, y, test, train)
#
#
# 	if y_train is None:
# 		estimator.fit(X_train, **fit_params)
# 	else:
# 		estimator.fit(X_train, y_train, **fit_params)
#
# 	return {
# 		'test_proba': estimator.predict_proba(X_test)[:,1],
# 		'y_test': y_test,
# 		'parameters': parameters,
# 		'candidate_id': candidate_id,
# 		'estimator': estimator
# 	}
#
#
# from sklearn.utils.metaestimators import _safe_split
#
# def parameter_search_aggregate_folds(estimator, X, y, candidate_params, method='predict_proba', metrics={'roc_auc_score': sklearn.metrics.roc_auc_score}, cv = LeaveOneOut(), method='predict_proba'):
# 	""" GridSearchCV and friends apply a `scoring` function to each CV fold, for
# 	each set of `params` in `candidate_params`. This approach is not compatible
# 	with `roc_auc` for `LeaveOneOut` cross-validation, since calculating this
# 	score requires aggregating the results of all folds.
#
# 	This function, for each set of `params` in `candidate_params`, trains a
# 	model on each `cv` fold, then calculates the score using the results of all
# 	folds.
# 	"""
#
# 	base_estimator = clone(self.estimator)
#
# 	parallel = Parallel(n_jobs=self.n_jobs, pre_dispatch=self.pre_dispatch)
#
# 	splits = list(cv.split(X, y, groups))
# 	ys = np.concatenate((y[test] for train, test in splits))
#
# 	results = []
# 	for (cand_idx, parameters) in enumerate(candidate_params):
# 		_estimator = clone(base_estimator).set_params(**parameters)
# 		y_preds = cross_val_predict(
# 			_estimator,
# 			X,
# 			y,
# 			fit_params=candidate_params,
# 			cv=cv
# 		)
#
# 		results.append({
# 			'estimator': _estimator,
# 			'parameters': parameters,
# 			**parameters,
# 			**{ metric_name: metric(y, y_preds) for metric_name, metric in metrics.items()  }
# 		})
#
# 	return pd.DataFrame(results)

# def parameter_search_aggregate_folds(estimator, X, y, candidate_params, cv = LeaveOneOut()):
# 	""" GridSearchCV and friends apply a `scoring` function to each CV fold, for
# 	each set of `params` in `candidate_params`. This approach is not compatible
# 	with `roc_auc` for `LeaveOneOut` cross-validation, since calculating this
# 	score requires aggregating the results of all folds.
#
# 	This function, for each set of `params` in `candidate_params`, trains a
# 	model on each `cv` fold, then calculates the score using the results of all
# 	folds.
# 	"""
#
# 	base_estimator = clone(self.estimator)
#
# 	parallel = Parallel(n_jobs=self.n_jobs, pre_dispatch=self.pre_dispatch)
#
# 	splits = list(cv.split(X, y, groups))
# 	ys = np.concatenate((y[test] for train, test in splits))
#
#
# 	with parallel:
# 		out = parallel(
# 			delayed(_fit_and_predict_proba)(
# 				clone(base_estimator),
# 				X,
# 				y,
# 				train=train,
# 				test=test,
# 				parameters=parameters,
# 				split_progress=(split_idx, n_splits),
# 				candidate_id=cand_idx,#(cand_idx, n_candidates),
# 				**fit_and_score_kwargs,
# 			)
# 			for (cand_idx, parameters), (split_idx, (train, test)) in product(
# 				enumerate(candidate_params), enumerate(cv.split(X, y, groups))
# 			)
# 		)

def list_of_records_to_dict_of_arrays(records, len_indicator, index=None):
	""" transform a list of dicts into a dict of arrays.
	Each dict may represent more than one row in the output array.
	records[i][`len_indicator`] will be used to mark how many rows record[i]
	should generate.
	"""
	columns = list(records[0].keys())
	n_rows = sum(len(rec[len_indicator]) for rec in records)
	
	result = { c: np.zeros(shape=(n_rows,), dtype=records[0][c].dtype if hasattr(records[0][c], 'dtype') else 'object') for c in columns }
	if index is not None:
		if not hasattr(records[0][index], 'index'):
			index=None
		else:
			result['index'] = np.zeros(shape=(n_rows,), dtype='object')
	# result['estimator'] = np.zeros(shape=(n_rows,))
	i = 0
	for record in records:
		record_len = len(record[len_indicator])
		for c in record.keys():
			v = record[c]
			# if hasattr(v, '__len__'):
			# 	if len(v) != record_len:
			# 		result[c][i:i+record_len] = (v,)
			try:
				result[c][i:i+record_len] = v
			except ValueError:
				result[c][i:i+record_len] = (v,)

		if index:
			result['index'][i:i+record_len] = record[index].index
		i += record_len
	return result



def cross_val_multi_predict(estimator, X, y, groups=None, cv=LeaveOneOut(), methods=['predict_proba'], fit_params={}, verbose=False, n_jobs=None, pre_dispatch='2*n_jobs'):
	""" simplified version of `sklearn.model_selection.cross_val_predict` that can give predictions as well as probabilities
	"""
	if isinstance(methods, list):
		methods = {m:m for m in methods}

	# allow methods to be specified as str (e.g. method name) or lambda
	columns = {}
	for name, method in methods.items():
		if isinstance(method, str):
			# method=method captures the value of the local variable `method` at this loop iteration,
			# https://docs.python.org/3/faq/programming.html#why-do-lambdas-defined-in-a-loop-with-different-values-all-return-the-same-result
			columns[name] = lambda estimator, X, y, method=method: getattr(estimator, method)(X)
		else:
			columns[name] = method

	def _fit_and_predict(estimator, X, y, train, test, verbose, fit_params, methods, fold):
		estimator.fit(X[train], y[train], **fit_params)
		return {'fold':fold, 'estimator': estimator, 'y': y[test], **{ name: method(estimator, X[test], y[test]) for name, method in methods.items() } }

	splits = list(cv.split(X, y, groups))
	if verbose:
		print(f"{cv} -> {len(splits)} splits")

	parallel = Parallel(n_jobs=n_jobs, verbose=verbose, pre_dispatch=pre_dispatch)
	predictions = parallel(
		delayed(_fit_and_predict)(
			clone(estimator), X, y, train, test, verbose, fit_params, methods=columns, fold=fold
		)
		for fold, (train, test) in enumerate(splits)
	)

	return list_of_records_to_dict_of_arrays(predictions, len_indicator='y', index='y')



def plot_hyperparam_grid(search, param_grid, limit_x=None, df=None, ax=None):
	if df is None:
		df = pd.DataFrame(search.cv_results_).sort_values('rank_test_score')

	if ax is None:
		fig, ax = plt.subplots()
	if limit_x is not None:
		df = df.iloc[:limit_x]
	ax.errorbar(range(len(df)), y=df['mean_test_score'], yerr=df['std_test_score'])
	ax.get_xaxis().set_visible(False)
	param_names = list(param_grid.keys())
	ax.table(cellText=df[[f"param_{p}" for p in param_names]].applymap('{:.2e}'.format).T.values,
			 rowLabels=param_names)
	ax.set_title(f"{param_names}")
	return ax

def plot_hyperparam_progress(summary, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    x_labels = summary['new_params'].apply(pprint_dict)
    xs = range(len(x_labels))
    ax.errorbar(xs, summary['mean_test_score'], yerr=summary['std_test_score'])
    ax.set_xticks(xs)
    ax.set_xticklabels(x_labels, rotation=-45, horizontalalignment='left');


from collections import namedtuple
SerialHyperparamSearchResults = namedtuple('SerialHyperparamSearchResults', ['summary', 'best_params', 'spaces'])

def serial_hyperparam_search(estimator, X, y, starting_params, param_grids,
	constant_params = dict(n_jobs = 1, seed=1337, use_label_encoder=False),
	cv=RepeatedStratifiedKFold(n_splits=3, n_repeats=8), search_cls=GridSearchCV,
	plot_progress=True, plot=True, verbose=True, ax=None,
	**kwargs):
	""" Perform serial hyperparameter searches; choose best set of hyperparameters at each step and use for training model at next step.
	"""

	from sklearn import clone
	from sklearn.model_selection import cross_val_score
	from statistics import NormalDist
	import pandas as pd
	import datetime

	base_estimator = clone(estimator)
	current_best_params = { **starting_params }
	base_estimator.set_params(**{ **starting_params, **constant_params })

	if verbose:
		print("Beginning stepwise hyperparameter optimization:")
		print(f"- Base estimator {base_estimator}")
		print(f"- Planning {len(param_grids)} optimization steps...")
		print("")

		print(f"Step 0: Baseline:")
		print(f"  starting parameters: {starting_params}")

	step_scores = []
	step_spaces = []

	# evaluate performance oft he estimator with the starting parameters
	start_time = datetime.datetime.now()
	starting_scores = cross_val_score(base_estimator, X, y, scoring='roc_auc', cv=cv, **kwargs)
	mean_starting_score = np.mean(starting_scores)
	std_starting_score =  np.std(starting_scores)
	finish_time = datetime.datetime.now()
	eval_time = (finish_time - start_time).total_seconds()

	# record results
	step_scores.append({
		'params': current_best_params,
		'new_params': {},
		'mean_test_score': mean_starting_score, 'std_test_score': std_starting_score,
		**current_best_params
	})

	step_spaces.append(pd.DataFrame({'step': 'baseline', **step_scores[0]}))

	if verbose:
		print(f"  score = {np.mean(starting_scores)}")
		print(f"  (Cross-validation took {eval_time} sec)")

	for i, param_grid in enumerate(param_grids):
		if verbose:
			print("-" * 80)
			print(f"Step {i+1}: Optimizing: {list(param_grid.keys())}")
			print(f"  Parameter space: {param_grid}")
			print(f"  Best params: so far {current_best_params}")

		params = {**current_best_params}
		for param in param_grid.keys():
			params.pop(param, None)

		estimator = clone(base_estimator)
		estimator.set_params(**{ **params, **constant_params })

		search = search_cls(
			estimator = estimator,
			param_grid = param_grid,
			scoring='roc_auc',
			cv=cv, **kwargs
		)
		search.fit(X,y)

		# update current set of parameters for the next round
		current_best_params = {**params, **search.best_params_}

		# record results
		df = pd.DataFrame(search.cv_results_).sort_values('rank_test_score')
		df['step'] = i
		step_spaces.append(df)
		step_scores.append({
			# cumulative best parameters
			'params':     current_best_params,
			# parameters optimized in this round of searching
			'new_params':  search.best_params_,
			**current_best_params,
			**df.iloc[0, :][['mean_test_score', 'std_test_score']].to_dict()
		})

		if verbose:
			# print results
			print(f"  Best parameters from step: {search.best_params_}\n  score = {search.best_score_}\n")
			if plot_progress:
				plot_hyperparam_grid(search, param_grid, df=df)

			print()

	summary = pd.DataFrame(step_scores, index=['baseline'] + list(range(len(param_grids))))

	if verbose:
		print("Final result:")
		print(f" - Initial score: {mean_starting_score} ± {std_starting_score}")
		print(f" - Final score:   {search.best_score_}")
		print(f" - Change:        {search.best_score_ - mean_starting_score}")
		print(f" - z:             {NormalDist(mu=mean_starting_score, sigma=std_starting_score).zscore(search.best_score_)}")

	if plot:
		plot_hyperparam_progress(summary, ax=ax)

	return SerialHyperparamSearchResults(summary, current_best_params, step_spaces)
