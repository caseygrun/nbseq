# Copyright (c) 2016-2021, QIIME 2 development team.
# Distributed under the terms of the Modified BSD License.

from .utils import *
import os

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from Bio.Align import MultipleSeqAlignment
import itertools

def run_command(cmd, output_fp, env=None, **kwargs):
	import subprocess

	print(' '.join(cmd))
	with open(output_fp, 'w') as output_f:
		subprocess.run(cmd, stdout=output_f, check=True, env=env, **kwargs)

def mafft(sequences_fp, result_fp, alignment_fp=None, n_threads=0, parttree=False, addfragments=False):
	# MAFFT uses --threads -1 to mean "auto"
	if n_threads == 0:
		n_threads = -1

	cmd = ["mafft", "--preservecase", "--inputorder",
			"--thread", str(n_threads)]

	if parttree == 'dp':
		cmd += ['--dpparttree']
	elif parttree == 'fasta':
		cmd += ['--fastaparttree']
	elif parttree == True:
		cmd += ['--parttree']

	if alignment_fp is not None:
		add_flag = '--addfragments' if addfragments else '--add'
		cmd += [add_flag, sequences_fp, alignment_fp]
	else:
		cmd += [sequences_fp]

	run_command(cmd, result_fp)



def clustalo(sequences_fp, result_fp, log_fp=None, n_threads=0, **kwargs):
	# clustalo --threads {threads} --in {input} --out {output} --log {log}

	cmd = ["clustalo", "--threads", str(n_threads), '--in', str(sequences_fp), '--out', str(result_fp)]

	if log_fp is not None:
		cmd += ['--log', str(log_fp)]

	run_cmd(cmd, **kwargs)

def muscle(sequences_fp, result_fp, log_fp=None, n_threads=0, method="super5", **kwargs):
	# muscle -super5 {input} -output {output} -threads {threads} 2>&1 > {log}

	cmd = ["muscle", f"-{method}", sequences_fp, "-output", str(result_fp), "-threads", str(n_threads)]

	if log_fp is not None:
		cmd += ["2>&1", str(log_fp)]

	run_cmd(cmd, shell=True, **kwargs)


def msa(sequences_fp, result_fp, method='mafft', **kwargs):
	methods = {
		'mafft': mafft,
		'clustalo': clustalo,
		'muscle': muscle,
	}
	method = method.lower()
	if method in methods:
		return methods[method](sequences_fp, result_fp, **kwargs)
	else:
		raise ValueError(f"Unrecognized msa method {method}. Recognized methods are {list(methods.keys())}")


def fasttree(alignment_fp, tree_fp, mode='aa', n_threads = 0):
	env = None
	if n_threads == 1:
		cmd = ['FastTree']
	else:
		env = os.environ.copy()
		# n_threads = 0 if n_threads == 'auto' else n_threads
		env.update({'OMP_NUM_THREADS': str(n_threads)})
		cmd = ['FastTreeMP']

	if mode == 'nt':
		cmd.extend(['-nt'])

	cmd.extend(['-quote', alignment_fp])
	run_command(cmd, tree_fp, env=env)
	return tree_fp


class ParallelDistanceCalculator(DistanceCalculator):
	def __init__(self,n_jobs,*args,**kwargs):
		self.n_jobs = n_jobs
		super().__init__(*args,**kwargs)

	n_jobs = -1
	def _parallel_pairwise(self, seq1, seq2):
		return (seq1.id, seq2.id, self._pairwise(seq1, seq2))

	def get_distance(self, msa):
		from joblib import Parallel, delayed

		if not isinstance(msa, MultipleSeqAlignment):
			raise TypeError("Must provide a MultipleSeqAlignment object.")

		# names = [s.id for s in msa]
		# dm = DistanceMatrix(names)
		parallel = Parallel(n_jobs=self.n_jobs)

		out = parallel(
			delayed(self._parallel_pairwise)(seq1, seq2)
				for seq1, seq2 in itertools.combinations(msa, 2))
		# for seq1, seq2 in itertools.combinations(msa, 2)
			# dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)

		names = [s.id for s in msa]
		dm = DistanceMatrix(names)
		for s1_id, s2_id, dist in out:
			dm[s1_id, s2_id] = dist
		return dm
