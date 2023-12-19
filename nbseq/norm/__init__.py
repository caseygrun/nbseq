import sys
from pathlib import *
import subprocess
import tempfile

import pandas as pd

from ..ft import read_sparse_dir, write_biom, to_anndata
from ..utils import run_cmd, mkdirp


def scran(ft=None, ft_path=None, out_path=None, obs_path=None, var_path=None, conda=None, verbose=True, return_out=None):

	with tempfile.TemporaryDirectory() as tmpdir:

		if ft_path is None:
			ft_path = Path(tmpdir) / 'ft.biom'
			# TODO: save ft to path
			raise NotImplementedError()

		if obs_path is None:
			obs_path = Path(tmpdir) / 'obs.csv'
		if var_path is None:
			var_path = Path(tmpdir) / 'var.csv'

		tmp_out_path = Path(tmpdir) / 'ft_out'
		mkdirp(tmp_out_path)

		# call scran
		cmd = ['Rscript', Path(__file__).parent / 'scran.R',
			str(ft_path), str(tmp_out_path), str(obs_path), str(var_path), int(verbose)
		]

		if verbose: print("Running R...")
		run_cmd(cmd, conda=conda, verbose=verbose)

		obs = pd.read_csv(obs_path).set_index('ID')
		var = pd.read_csv(var_path).set_index('feature')
		ft_out = read_sparse_dir(tmp_out_path)

	if out_path is not None:
		write_biom(ft_out, out_path)

	if return_out is None:
		return_out = out_path is None	

	if return_out:
		if ft is not None:
			obs = ft.obs.join(obs)
			var = ft.var.join(var)
			# out.obs_names.name = ft.obs_names.name
			# out.var_names.name = ft.var_names.name

		out = to_anndata(ft_out, obs=obs, var=var)
		return out
