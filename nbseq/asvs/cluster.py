import sys
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import nbseq.utils
import numbers
from ..utils import ungapped

def cluster_asvs(input, asvid_col='ASVID', seq_col='seq', read_id_col='pair_id', cluster_path=None,
	cluster_method='linclust', min_seq_id=1, coverage_mode='bidirectional', min_coverage=None,
	threads=1, verbose=True):

	def run_cmd(cmd):
		cmd = [str(x) for x in cmd]
		print(" ".join(cmd))# )
		return subprocess.run(cmd, check=True,
			stdout=(sys.stdout if verbose else subprocess.DEVNULL),
			stderr=(sys.stderr if verbose else subprocess.DEVNULL)) #, stdout=logfile, stderr=logfile)

	# make sure each ASVID maps to a single unique sequence
	asvs_that_are_uniq = input[[asvid_col,seq_col]].groupby(asvid_col)[seq_col].nunique() == 1
	if not all(asvs_that_are_uniq):
		nonuniq_asvs = asvs_that_are_uniq[~asvs_that_are_uniq].index.values

		# it's possible that the same ASVID maps to the same sequence which
		# aligned differently. Detect this and permit it.
		suspects = input.loc[input[asvid_col].isin(nonuniq_asvs),:].copy()
		suspects['seq_ungapped'] = suspects[seq_col].apply(ungapped)
		asvs_that_are_uniq_ungapped = suspects.groupby(asvid_col)['seq_ungapped'].nunique() == 1

		if not all(asvs_that_are_uniq_ungapped):
			nonuniq_asvs = asvs_that_are_uniq_ungapped[~asvs_that_are_uniq_ungapped].index.values
			print(suspects.set_index(asvid_col).loc[nonuniq_asvs,'seq_ungapped'])
			raise Exception("Error: the following ASVIDs are associated with more than one unique sequence: "
			f"{nonuniq_asvs}")

	# contains one row per unique ASVID and therefore per unique sequence
	uniq_seqs = input.drop_duplicates([asvid_col],keep='first')

	if read_id_col is None:
		if len(uniq_seqs) < len(input):
			print("Warning: some ASVIDs are duplicated in the input. This is not an error, since they have the same sequences, but may indicate a problem with an upstream step. ")

	if verbose:
		if read_id_col is None:
			print(f"Clustering {len(uniq_seqs)} unique ASVs (no read pair ID given, all ASVs should be distinct)...")
		else:
			print(f"Clustering {input[read_id_col].nunique()} unique read pairs, {len(uniq_seqs)} unique ASVs...")

	with tempfile.TemporaryDirectory() as tmpdir:

		# dump unique sequences to FASTA file
		input_fasta_path = Path(tmpdir) / 'input.fasta'
		nbseq.utils.dataframe_to_fasta(uniq_seqs,
			file=input_fasta_path,
			seq_col = seq_col, name_col = asvid_col,
			ungap=True)

		tmp_input_db = Path(tmpdir) / 'input_db'
		tmp_cluster_db = Path(tmpdir) / 'cluster_db'
		tmp_scratch = Path(tmpdir) / 'scratch'
		if cluster_path is None: cluster_path = str(Path(tmpdir) / 'clusters.tsv')

		# create database
		run_cmd(['mmseqs','createdb',
			str(input_fasta_path), str(tmp_input_db)])

		# perform clustering
		cmd = ['mmseqs',cluster_method,
			str(tmp_input_db), str(tmp_cluster_db), str(tmp_scratch),
			'--min-seq-id', min_seq_id, # ensure all aligned residues are identical
			'--threads', threads]
		if coverage_mode is not None:
			if not isinstance(coverage_mode, numbers.Number):
				coverage_modes = {
					'bidirectional':0,
					'target':1,
					'query':2,
					'target_length':3
				}
				if coverage_mode not in coverage_modes:
					raise Exception(f"Unknown coverage mode {coverage_mode}, allowed coverage modes: {coverage_modes.keys()}, or integers 0--3")
				coverage_mode = coverage_modes[coverage_mode]
			cmd = cmd + ['--cov-mode', coverage_mode]
		if min_coverage is not None:
			cmd = cmd + ['-c', min_coverage]
		run_cmd(cmd)

		# export cluster identities to TSV
		run_cmd(['mmseqs','createtsv',
			str(tmp_input_db),
			str(tmp_input_db), # not an error, this needs to be entered twice
			str(tmp_cluster_db),
			cluster_path])

		# read back in clusters
		# two columns, no header:
		# cluster-representative 	cluster-member
		clusters = pd.read_csv(cluster_path,
			sep="\t", header=None, names=['cluster_ASVID','member_ASVID'])

		# write clusters with header, and in mapping order
		clusters[['member_ASVID','cluster_ASVID']].to_csv(cluster_path, index=False)

	if verbose:
		print('mmseq2 finished; cleaning up temporary files.')
		n_clusters = clusters['cluster_ASVID'].nunique()
		n_asvs = len(uniq_seqs)
		print(f"Grouped {n_asvs} ASVs into {n_clusters} clusters ({(n_asvs-n_clusters)/n_asvs*100:.2}% of sequences were collapsed)\n")

		print("Re-mapping ASVs to clusters...")
	# cluster_map = clusters.set_index('member')['cluster']
	# print(cluster_map,flush=True)
	#
	#
	# pair_ids = input[read_id_col]
	# # print(input[asvid_col],flush=True)
	# # new_asvids = pd.merge(left=input)
	#
	# asvs_clustered = input[asvid_col].isin(cluster_map)
	# if not all(asvs_clustered):
	# 	raise Exception(f"{sum(asvs_clustered)} ASVID(s) not assigned to a cluster: {input[asvid_col][~asvs_clustered]}")
	#
	# new_asvids = cluster_map.loc[input[asvid_col]].rename(asvid_col) #cluster_map.reindex(input[asvid_col]) #input[asvid_col].apply(lambda x: cluster_map[x])
	#
	# seqs = uniq_seqs.set_index(asvid_col).loc[new_asvids,seq_col].rename(seq_col)
	#
	# print(new_asvids)
	# print(new_asvids.values)
	# print(seqs)
	# print(seqs.values)
	# print(pair_ids)
	# print(pair_ids.values)
	#
	# # output = pd.concat({asvid_col:new_asvids, seq_col:seqs, read_id_col:pair_ids}, axis='columns', copy=False)
	# output = pd.DataFrame({
	# 	asvid_col: new_asvids.values,
	# 	seq_col:seqs.values,
	# 	read_id_col:pair_ids.values
	# })


	if read_id_col is not None:
		# for each read pair (pair_id), find the cluster_ASVID by mapping the
		# original ASVID to the cluster `member_ASVID`
		df = pd.merge(left=input[[asvid_col,read_id_col]], right=clusters,
			how='left', left_on=asvid_col,right_on='member_ASVID', validate='many_to_one')

		assert len(df) == len(input), f"Error: merged dataset should have one entry per read pair, {len(input)}; length was {len(df)}. Check number of clusters {len(clusters)}."

		# find the sequence for the cluster_ASVID, which will become the new ASVID
		# for the read pair
		df2 = pd.merge(left=df[[read_id_col,'cluster_ASVID']],
			right=uniq_seqs[[asvid_col,seq_col]], # search uniq_seqs rather than input to avoid duplicating rows
			how='left', left_on='cluster_ASVID',right_on=asvid_col)

		assert len(df2) == len(input), f"Error: merged dataset should have one entry per read pair, {len(input)}; length was {len(df2)}. Check number of clusters {len(clusters)}."

		output = df2[[read_id_col,'cluster_ASVID',seq_col]].rename(columns={'cluster_ASVID':asvid_col})
	else:
		df2 = pd.merge(left=clusters, right=uniq_seqs[[asvid_col,seq_col]],
			how='left', left_on='cluster_ASVID',right_on=asvid_col)
		assert len(df2) == len(uniq_seqs), f"Error: merged dataset should have one entry per ASV {len(uniq_seqs)}; length was {len(df2)}. Check number of clusters {len(clusters)}."
		output = df2[['cluster_ASVID',seq_col]].rename(columns={'cluster_ASVID':asvid_col})

	if verbose:
		print("Done re-mapping ASVs.")
	return output
