import collections.abc
import sys
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from nbseq.utils import *

def search_mmseqs2(query_sequences, db_path, tmp_path=None, result_path=None, verbose=False, conda=None, **kwargs):

	def run_cmd(cmd, **kwargs):
		# {cmd[0]} {cmd[1]} ... {cmd[-1]} --foo {kwargs.foo} -b {kwargs.b}
		cmd = " ".join([str(x) for x in cmd] + [f"--{k.replace('_','-')} '{v}'" if len(k) > 1 else f"-{k} '{v}'" for k,v in kwargs.items()])
		if verbose: print(cmd)
		if conda is not None:
			# https://github.com/conda/conda/issues/7980
			cmd = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda} && {cmd}'

		out = subprocess.run(cmd,
			check=True,
			shell=(True if conda else False),
			capture_output=verbose,
			stdout = (None if verbose else subprocess.DEVNULL),
			stderr = (None if verbose else subprocess.DEVNULL)) #, stdout=logfile, stderr=logfile)
		if verbose:
			print(out.stdout.decode('utf8'), file=sys.stdout)
			print(out.stderr.decode('utf8'), file=sys.stderr)

	with tempfile.TemporaryDirectory() as tmpdir:

		# dump unique sequences to FASTA file
		tmp_input_fasta_path = Path(tmpdir) / 'input.fasta'
		with open(tmp_input_fasta_path, "w") as infile:
			query_sequences = fortify_to_seqrecords(query_sequences)
			SeqIO.write(query_sequences, infile, "fasta")

		if result_path is None:
			result_path = Path(tmpdir) / 'alnResults.m8'
		if tmp_path is None:
			tmp_path = Path(tmpdir) / 'tmp'
		mkdirp(tmp_path)

		# perform clustering
		cmd = ['mmseqs','easy-search',
			str(tmp_input_fasta_path), str(db_path), str(result_path), str(tmp_path)
		]
		run_cmd(cmd, **kwargs)

		results = pd.read_csv(result_path,
			sep="\t",
			names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

	return results


def db_uri_to_abs(db):
	"""connectionx requires an absolute path to sqlite databases"""
	if db.startswith("sqlite:///"):
		protocol, pth = db.split(":///")
		return f"{protocol}:///{os.path.abspath(pth)}"
	else:
		return db


_sql_engines = {}
def get_sql_engine(db, **kwargs):
	if db in _sql_engines:
		return _sql_engines[db]
	else:
		from sqlalchemy import create_engine
		_sql_engines[db] = create_engine(db, **kwargs)
		return _sql_engines[db]


def make_sql_where_clause(operator='and',**kwargs):
	predicates = []
	for column, pattern in kwargs.items():
		if isinstance(pattern, tuple) and len(pattern) == 1:
			predicates.append(f"({column} = '{pattern[0]}')")
		if isinstance(pattern,collections.abc.Iterable) and not isinstance(pattern, str):
			# lst = ",".join(f"'{p}'" for p in pattern)
			# predicates.append(f"(`{column}` IN ({lst}))")
			predicate = " OR ".join(f"({column} LIKE '{p}%')" for p in pattern)
			predicates.append(f"({predicate})")
		else:
			predicates.append(f"({column} LIKE '{pattern}%')")

	# predicate = "&&".join(f"(`{column}` LIKE '{pattern}%')" for column, pattern in kwargs.items())
	if operator == 'and':
		operator = " AND "
	elif operator == "or":
		operator = " OR "
	else: raise Exception(f"unrecognized operator {operator}")
	predicate = operator.join(predicates)

	return predicate

def search_sql(db='sqlite:///intermediate/aa/asvs.db', single=False, table='cdrs', operator='and', columns='*', query=None, engine='cx', show_query=False, **kwargs):
	"""searches a SQL database for features matching a set of predicates

	Parameters
	----------
	db : str
		SqlAlchemy connection string, e.g. 'sqlite:///intermediate/aa/asvs.db'
	single : bool
		True to collapse results to a single feature
	table : str, default='cdrs'
		name of the SQL table
	operator : str, default='and'
		Return rows that match all of the predicates (`and`) or any predicate (`or`)?
	query : str, optional
		If given, raw SQL query to execute (table, operator, columns, and predicate will be ignored)
	**kwargs : str
		set of predicates, in COLUMN=VALUE format. Control whether to match
		all predicates or any predicate with `operator`.

		- If an argument value is a scalar, rows which begin with ``pattern``
		  will be returned, e.g. the predicate will be ``(`{column}` LIKE `{pattern}%`)``
		- If an argument value is a 1-tuple, rows will be searched for columns
		  with that exact value, e.g. ``(`{column}` = '{pattern}')``
		- If an argument is an iterable, rows will be returned where the column
		  contains any value from the iterable, e.g. ``(`{column}` IN ({pattern}))``

	Returns
	-------
	pd.DataFrame or pd.Series
		DataFrame with one result per row, or Series if ``single=True``

	"""

	if query is None:
		predicate = make_sql_where_clause(operator=operator, **kwargs)

		if isinstance(columns, list):
			columns = ",",join(f"`{col}`" for col in columns)

		query = f"SELECT {columns} FROM `{table}` WHERE ({predicate})"

		failure_msg = f"Feature(s) matching '{kwargs}' not found"
		single_msg = f"Warning: more than one feature matches '{kwargs}':"
	else:
		failure_msg = f"Feature(s) matching '{query}' not found"
		single_msg = f"Warning: more than one feature matches '{query}':"


	if show_query:
		print(query)
	if engine == 'cx':
		import connectorx as cx
		fd = cx.read_sql(db_uri_to_abs(db), query)
	else:
		con = get_sql_engine(db)
		fd = pd.read_sql(query, con=con)

	if (len(fd) > 1 and single):
		print(f"{single_msg}: {fd}")
		fd = fd.iloc[0,:]
	elif (len(fd) == 0):
		print(failure_msg)

	return fd
