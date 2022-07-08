"""
Functions for working with nbseq data as tidy DataFrames
"""
import pandas as pd
from .utils import *

def segment_lens(table, segments, len_f = len_ungapped):
	"""Calculate the length of one or more segments

	table : pd.DataFrame
	segments : list of str
		list of column names; calculate the length of each of these columns with len_f
	len_f : callable
	"""
	return pd.concat([
		(table[segments]
			 .applymap(len_f, na_action='ignore')),
		table.drop(segments, axis='columns')
	], axis=1)
