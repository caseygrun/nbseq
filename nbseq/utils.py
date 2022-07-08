import os
import sys
import pathlib
import itertools
import hashlib
from collections import abc

import numpy as np

try:
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	import Bio.SeqIO
	import Bio.Data.CodonTable

	Bio.Data.CodonTable.register_ncbi_table(
		name="Amber suppressing",
		alt_name=None,
		id=999,
		table={
			"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
			"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
			"TAT": "Y", "TAC": "Y",             "TAG": "Q",   # noqa: E241
			"TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
			"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
			"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
			"CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
			"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
			"ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
			"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
			"AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
			"AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
			"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
			"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
			"GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
			"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
		},
		stop_codons=["TAA", "TGA"],
		start_codons=["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
	)


	# # BEGIN BIOPYTHON MONKEY PATCH
	# # FROM https://github.com/peterjc/biopython/blob/ffd6cb6d2d38e6243ae573fdd38c22e84505d21b/Bio/SeqRecord.py
	# # MIT LICENSE
	# # MODIFICATIONS BY PETER COCK
	# #
	# # Modifications to Bio.SeqRecord:
	# # Copyright 2000-2002 Andrew Dalke.
	# # Copyright 2002-2004 Brad Chapman.
	# # Copyright 2006-2009 by Peter Cock.
	#
	# # Modifications to Bio.SeqFeature:
	# # Copyright 2000-2003 Jeff Chang.
	# # Copyright 2001-2008 Brad Chapman.
	# # Copyright 2005-2010 by Peter Cock.
	# # Copyright 2006-2009 Michiel de Hoon.
	#
	# import Bio.SeqRecord
	# import Bio.SeqFeature
	#
	# def _map(self, mapping):
	# 	"""Returns a copy of the feature with its location adjusted (PRIVATE).
	#
	# 	The annotation qualifiers are copied."""
	# 	return SeqFeature(location = self.location._map(mapping),
	# 			type = self.type,
	# 			location_operator = self.location_operator,
	# 			strand = self.strand,
	# 			id = self.id,
	# 			qualifiers = dict(self.qualifiers.iteritems()),
	# 			sub_features = [f._map(mapping) for f in self.sub_features],
	# 			ref = self.ref,
	# 			ref_db = self.ref_db)
	# Bio.SeqFeature.SeqFeature._map=_map
	#
	# def _map(self, mapping):
	# 	"""Returns a copy of the location adjusted (PRIVATE)."""
	# 	return FeatureLocation(start = self._start._start_map(mapping),
	# 						   end = self._end._end_map(mapping))
	# Bio.SeqFeature.FeatureLocation._map = _map
	#
	# def _start_map(self, mapping):
	# 	#We want this to maintain the subclass when called from a subclass
	# 	old = self.position
	# 	new = mapping[old]
	# 	while new is None:
	# 		old += 1
	# 		new = mapping[old]
	# 	return self.__class__(new, self.extension)
	#
	# def _end_map(self, mapping):
	# 	#We want this to maintain the subclass when called from a subclass
	# 	old = self.position - 1
	# 	new = mapping[old]
	# 	while new is None:
	# 		old -= 1
	# 		new = mapping[old]
	# 	return self.__class__(new + 1, self.extension)
	#
	# Bio.SeqFeature.AbstractPosition._start_map = _start_map
	# Bio.SeqFeature.AbstractPosition._end_map = _end_map
	#
	# def _start_map(self, mapping, step):
	# 	#We want this to maintain the subclass when called from a subclass
	# 	return self.__class__([position_choice._start_map(mapping, step) \
	# 						   for position_choice in self.position_choices])
	#
	# def _end_map(self, mapping, step):
	# 	#We want this to maintain the subclass when called from a subclass
	# 	return self.__class__([position_choice._end_map(mapping, step) \
	# 						   for position_choice in self.position_choices])
	#
	# Bio.SeqFeature.OneOfPosition._start_map = _start_map
	# Bio.SeqFeature.OneOfPosition._end_map = _end_map
	#
	# def _ungap(self, gap=None):
	# 	"""Return a SeqRecord without the gap character(s) in the sequence.
	# 	The gap character can be specified in two ways - either as an explicit
	# 	argument, or via the sequence's alphabet. For example:
	# 	>>> from Bio.SeqRecord import SeqRecord
	# 	>>> from Bio.Seq import Seq
	# 	>>> from Bio.Alphabet import generic_dna
	# 	>>> my_dna = SeqRecord(Seq("-ATA--TGAAAT-TTGAAAA-", generic_dna), id="X")
	# 	>>> my_dna.seq
	# 	Seq('-ATA--TGAAAT-TTGAAAA-', DNAAlphabet())
	# 	>>> my_dna.ungap("-").seq
	# 	Seq('ATATGAAATTTGAAAA', DNAAlphabet())
	# 	If the gap character is not given as an argument, it will be taken from
	# 	the sequence's alphabet (if defined). For more details, see the Seq
	# 	object's ungap method.
	#
	# 	Any per-letter-annotation is sliced to match how the sequence gets
	# 	sliced to remove the gaps. Other annotations is retained as is.
	# 	SeqFeature locations are adjusted to use the new coordinates:
	# 	>>> from Bio.SeqFeature import SeqFeature, FeatureLocation
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(0,4)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(0,5)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(1,4)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(1,12)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(4,13)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(5,13)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(6,13)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(12,20)))
	# 	>>> my_dna.features.append(SeqFeature(FeatureLocation(12,21)))
	# 	>>> for f in my_dna.features:
	# 	...     print f.location, f.extract(my_dna.seq)
	# 	[0:4] -ATA
	# 	[0:5] -ATA-
	# 	[1:4] ATA
	# 	[1:12] ATA--TGAAAT
	# 	[4:13] --TGAAAT-
	# 	[5:13] -TGAAAT-
	# 	[6:13] TGAAAT-
	# 	[12:20] -TTGAAAA
	# 	[12:21] -TTGAAAA-
	# 	Notice most of these examples deliberately have the features ending on
	# 	a gap character. The start and end positions are adjusted to ensure the
	# 	feature describes the equivalent ungapped sequence:
	#
	# 	>>> ungapped = my_dna.ungap("-")
	# 	>>> for f in ungapped.features:
	# 	...     print f.location, f.extract(ungapped.seq)
	# 	[0:3] ATA
	# 	[0:3] ATA
	# 	[0:3] ATA
	# 	[0:9] ATATGAAAT
	# 	[3:9] TGAAAT
	# 	[3:9] TGAAAT
	# 	[3:9] TGAAAT
	# 	[9:16] TTGAAAA
	# 	[9:16] TTGAAAA
	#
	# 	For example with per-letter-annotation, we'll use Bio.SeqIO to load an
	# 	Ace assembly which includes quality scores but the sequence will be
	# 	padded with any gap characters (for which there is no quality score
	# 	available). You may want to get the ungapped sequence with its quality
	# 	scores (e.g. to output as FASTQ):
	#
	# 	>>> from Bio import SeqIO
	# 	>>> record = SeqIO.read("Ace/consed_sample.ace", "ace")
	# 	>>> print len(record)
	# 	1475
	# 	>>> print len(record) - record.seq.count("-")
	# 	1468
	# 	>>> print record[860:880].format("fastq")
	# 	@Contig1 <unknown description>
	# 	CAGCAGAGAAGGGTTTGAAA
	# 	+
	# 	z{{{{yyyyyy{{{{{{{{{
	# 	<BLANKLINE>
	# 	In the above example we've selected a subsection of the record to show
	# 	in FASTQ format. Now lets remove the gaps:
	#
	# 	>>> ungapped = record.ungap()
	# 	>>> print len(ungapped)
	# 	1468
	#
	# 	Notice below how the coordinates for the region [860:880] have shifted
	# 	by two since there are two gaps before it in the original record:
	# 	>>> record.seq[0:860].count("-")
	# 	2
	# 	>>> record[860:880].format("fastq") == ungapped[858:878].format("fastq")
	# 	True
	# 	So, using this method we can take the gapped consensus records from any
	# 	ACE file and save them as ungapped records in FASTQ, FASTA, QUAL, etc:
	# 	records = (rec.ungap() for rec in SeqIO.parse(in_file, "ace"))
	# 	count = SeqIO.write(records, out_file, "fastq")
	# 	"""
	# 	new_seq = self.seq.ungap(gap)
	# 	if str(new_seq) == str(self.seq):
	# 		#Unchanged, not even the alphabet - don't need to alter annotation
	# 		return self
	# 	if not gap:
	# 		gap = self.seq.alphabet.gap_char
	#
	# 	if not(self.features or self.letter_annotations):
	# 		return SeqRecord(new_seq, id=self.id, name=self.name,
	# 						 description=self.description,
	# 						 dbxrefs=self.dbxrefs[:],
	# 						 annotations=self.annotations.copy())
	#
	# 	slices = []
	# 	new_index = -1
	# 	if self.seq[0]==gap:
	# 		start = None
	# 		in_gap = True
	# 	else:
	# 		start = 0
	# 		in_gap = False
	# 	mapping = []
	# 	for old_index, letter in enumerate(self):
	# 		if letter == gap:
	# 			mapping.append(None)
	# 			if in_gap:
	# 				pass
	# 			else:
	# 				in_gap = True
	# 				if start is not None:
	# 					slices.append((start, old_index))
	# 					start = None
	# 		else:
	# 			new_index += 1
	# 			assert letter == new_seq[new_index]
	# 			mapping.append(new_index)
	# 			if in_gap:
	# 				in_gap = False
	# 				assert start is None
	# 				start = old_index
	# 	if not in_gap:
	# 		slices.append((start, len(self)))
	# 	assert len(mapping) == len(self)
	# 	mapping.append(len(new_seq))
	# 	if str(new_seq) != "".join(str(self.seq[s:e]) for s,e in slices):
	# 		msg = "%s\n%s\n" % (repr(slices), self.seq)
	# 		for s,e in slices:
	# 			msg += " "*s + self.seq[s:e] + "\n"
	# 		assert False, msg
	#
	# 	answer = SeqRecord(new_seq, id=self.id, name=self.name,
	# 					   description=self.description,
	# 					   dbxrefs=self.dbxrefs[:],
	# 					   annotations=self.annotations.copy())
	# 	#Apply the slices to the per letter annotation
	# 	for key, value in self.letter_annotations.iteritems():
	# 		s, e = slices[0]
	# 		new = value[s:e]
	# 		for s, e in slices[1:]:
	# 			new += value[s:e]
	# 		answer.letter_annotations[key] = new
	# 	#Apply the mapping to the feature coordinates
	# 	for f in self.features:
	# 		answer.features.append(f._map(mapping))
	# 	return answer
	#
	# Bio.SeqRecord.SeqRecord.ungap = _ungap
	# # END BIOPYTHON MONKEY PATCH

except ImportError:
	print("Warning: BioPython is not installed. Some nbseq functions may not work.")

def md5(s):
	return hashlib.md5(s.encode('utf-8')).hexdigest()

def md5_ungapped(s):
	return md5(ungapped(s))

def rc(seq):
	return str(Bio.Seq.Seq(seq).reverse_complement())
_rc = rc


def read_delim_auto(path, **kwargs):
	import pandas as pd

	fn, ext = os.path.splitext(path)
	ext = ext.lower()
	if ext == '.csv':
		return pd.read_csv(str(path), **{'sep':',', **kwargs})
	elif ext == '.tsv' or ext == 'txt':
		return pd.read_csv(str(path), **{'sep':"\t", **kwargs})
	else:
		raise ValueError(f"Unable to guess delimiter from file wiht extension '{ext}'")

def feature_table_to_unique_sequences(feature_table):
	seqs = feature_table.index.unique().values
	print("Unique sequences:")
	print(seqs)
	seq_records = (SeqRecord(Seq(seq), id = md5(seq)) for seq in seqs)
	return seq_records

def seqrecords_to_dataframe(records, transform=None):
	"""transform(id, seq, name, desc)"""
	import pandas as pd

	sequences = []
	ids = []
	names = []
	descriptions = []
	for record in records:
		(id, seq, name, desc) = (record.id, str(record.seq), record.name, record.description)
		if transform is not None:
			(id, seq, name, desc) = transform(id, seq, name, desc)
		ids.append(id)
		sequences.append(seq)
		names.append(name)
		descriptions.append(description)
	return pd.DataFrame({
		'sequence': sequences,
		'id': ids,
		'name': names,
		'description':descriptions
	})

def dataframe_to_seqrecords(df, seq_col = None, name_col = None, description_col=None, ungap=False, **kwargs):
	if name_col is None: names = df.index
	elif callable(name_col): names = name_col(df)
	else: names = df.loc[:,name_col]

	if seq_col is None: seqs = df.iloc[:,0]
	elif callable(seq_col): seqs = seq_col(df)
	else: seqs = df.loc[:,seq_col]
	if ungap: seqs = seqs.apply(ungapped)

	if description_col is None: descs = [''] * len(seqs)
	elif callable(description_col): descs = description_col(df)
	else: descs = df.loc[:,description_col]

	seq_records = (SeqRecord(Seq(seq), id = name, description=desc, **kwargs) for seq,name,desc in zip(seqs,names, descs))
	return seq_records

def dataframe_to_fasta(df, file, **kwargs):
	return Bio.SeqIO.write(dataframe_to_seqrecords(df, **kwargs), file, format='fasta' )

def dataframe_to_msa(df, **kwargs):
	from Bio.Align import MultipleSeqAlignment
	return MultipleSeqAlignment(dataframe_to_seqrecords(df, **kwargs))

def dataframe_to_fastas(df, file_col=None, **kwargs):
	if file_col is None: files = df.index
	elif callable(file_col): files = file_col(df)
	else: files = df.loc[:,file_col]
	records = dataframe_to_seqrecords(df, **kwargs)

	for file, record in zip(files, records):
		Bio.SeqIO.write([record], file, format='fasta' )

def series_to_seqrecords(sr):
	return dataframe_to_seqrecords(sr.to_frame())

def series_to_fasta(sr):
	return dataframe_to_fasta(sr.to_frame())

def fortify_to_seqrecords(data, **kwargs):
	"""accepts a collection and tries to convert it to an iterable of ``SeqRecord``s

	Parameters
	----------
	data : iterable of SeqRecord, iterable of str, pd.DataFrame, pd.Series
		data in any of the above formats:
		- pd.DataFrame will be passed to ``dataframe_to_seqrecords``
		- pd.Series will be passed to ``series_to_seqrecords``
		- iterable of str, each str is assumed to be a sequence, indexes 0, 1, ...
		  are given as the ``name``
	**kwargs
		passed to ``dataframe_to_seqrecords`` or ``series_to_seqrecords``

	Returns
	-------
	iterable of SeqRecord
	"""
	import itertools
	import pandas as pd
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord

	if isinstance(data, pd.DataFrame):
		return dataframe_to_seqrecords(data, **kwargs)
	elif isinstance(data, pd.Series):
		return series_to_seqrecords(data, **kwargs)
	elif isinstance(data, abc.Mapping):
		return (SeqRecord(id=name, seq=Seq(seq)) for name, seq in data.items())
	elif isinstance(data, abc.Iterable):
		peek = next(r for r in data)
		# chain = itertools.chain([peek], data)
		if isinstance(peek, SeqRecord):
			return data
		elif isinstance(peek,str):
			return (SeqRecord(id=str(i), seq=Seq(seq)) for i, seq in enumerate(data))

def chromatograms_to_dataframe(paths, rc=False, names=True, check_quality=True):
	"""load chromatograms in scf format to a pd.DataFrame

	Parameters
	----------
	paths : iterable of str
		paths to chromatograms in .scf format
	rc : bool
		True to reverse-complement chromatogram sequence
	names : bool or iterable or None
		iterable of names for the sequences, or True to take the name as the
		basename sans ext for each of the  ``paths``, or None to add no
		`name` column

	Returns
	-------
	pd.DataFrame
		DataFrame with 2--3 columns:
		- `file` containing the given path,
		- `seq` containing the sequence,
		- `name` containing given or generated according to the ``names`` parameter
	"""
	import bioconvert.io.scf
	import pandas as pd, numpy as np

	seqs = []
	for f in paths:
		seq, qualities, *_ = bioconvert.io.scf.read_scf(f)
		if check_quality:
			q = np.array(qualities)
			bad_bases = (q < 30).sum()
			if bad_bases > 3:
				print(f"Warning: chromatogram `{f}` has {bad_bases} positions "
						"with Q < 30. Manually inspect chromatogram to verify "
						"it is properly trimmed and base calls reflect "
						"chromatogram traces. Suppress this warning with "
						"check_quality=False .")

		if rc: seq = _rc(seq)
		seqs.append(seq)
	cols = {
		'file': paths,
		'seq': seqs
	}
	if names == True:
		names = [os.path.splitext(os.path.basename(x))[0] for x in paths]
	if names is not None:
		cols['name'] = names

	return pd.DataFrame(cols)



def make_cmd(*args, **kwargs):
	# {cmd[0]} {cmd[1]} ... {cmd[-1]} --foo {kwargs.foo} -b {kwargs.b}
	cmd = " ".join([str(x) for x in args] + [f"--{k.replace('_','-')} '{v}'" if len(k) > 1 else f"-{k} '{v}'" for k,v in kwargs.items()])
	return cmd

def _cmd(cmd):
	if not isinstance(cmd, str):
		cmd = " ".join([str(x) for x in cmd])
	return cmd

def run_cmd(cmd, verbose=False, conda=None, **kwargs):
	import subprocess
	import tempfile

	if conda is not None:
		cmd = _cmd(cmd)

		# https://github.com/conda/conda/issues/7980
		cmd = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda} && {cmd}'
	else:
		cmd = [str(x) for x in cmd]
	if verbose: print(cmd)

	try:
		out = subprocess.run(cmd,
			check=True,
			shell=(True if conda else False),
			capture_output=verbose,
			stdout = (None if verbose else subprocess.DEVNULL),
			stderr = (None if verbose else subprocess.DEVNULL),
			**kwargs)
	except subprocess.CalledProcessError as err:
		if verbose:
			print(err.stdout.decode('utf8'))
			print(err.stderr.decode('utf8'))
		raise err

	if verbose:
		print(out.stdout.decode('utf8'), file=sys.stdout)
		print(out.stderr.decode('utf8'), file=sys.stderr)


def run_script(cmd, verbose=False, conda=None, strict=True, silent=False, capture=False, **kwargs):
	import subprocess
	import tempfile

	if conda is not None:
		cmd = f'eval "$(conda shell.`basename -- $SHELL` hook)" && conda activate {conda}\n' + cmd

	if strict:
		cmd = "set -euo pipefail\n"+cmd

	with tempfile.NamedTemporaryFile('w') as script:
		if verbose:
			print(script.name)
			print(cmd)

		print(cmd, file=script, flush=True)

		subprocess.run(['sh', script.name], **dict(
			check=True,
			capture_output=capture,
			stdout=(subprocess.DEVNULL if silent else sys.stdout),
			stderr=(subprocess.DEVNULL if silent else sys.stderr),
			**kwargs)
		)

		# if capture:
		# 	print(out.stdout.decode('utf8'), file=sys.stdout)
		# 	print(out.stderr.decode('utf8'), file=sys.stderr)





def mkdirp_file(path):
	p = pathlib.Path(path).resolve().parent
	p.mkdir(parents=True, exist_ok=True)

def mkdirp(path):
	p = pathlib.Path(path)
	p.mkdir(parents=True, exist_ok=True)

def grouper(n, iterable, fillvalue=None):
	"Collect data into fixed-length chunks or blocks"
	# grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
	args = [iter(iterable)] * n
	return itertools.zip_longest(fillvalue=fillvalue, *args)


def ungapped(seq,gap='-'):
	return seq.replace(gap,'')

def len_ungapped(seq, gap='-'):
	return len(seq) - seq.count(gap)

def index_to_column(df, name=None, inplace=False):
	if name is not None:
		if inplace:
			df.index.rename(name, inplace=True)
		else:
			df.index = df.index.rename(name)
	if inplace:
		df.reset_index(inplace=inplace)
		return df
	else:
		return df.reset_index()

def show_ruler(start, end, step_size = 10):
	out = ''
	for num in range(start,end,step_size):
		label = str(num)
		out = out + label + ' ' * (step_size - len(label))
	return out

def show_sequence_ranges(seq, positions, number = True):
	out = ''
	seq_len = len(seq)
	if number:

		step_size = 10
		out = out + show_ruler(0, seq_len, step_size) + "\n"
	out = out + seq + "\n"
	for feature_name, pos in positions.items():
		track = ' ' * pos[0]
		feature_len = pos[1] - pos[0]
		track = track + (feature_name + '-' * (feature_len - len(feature_name)))
		track = track + ' '*(seq_len - len(track))
		out = out + track + "\n"
	return out

def show_sequence_translated(seq, translation=None, suppress_amber=False, ruler=True, positions={}):
	if not isinstance(seq, Seq):
		seq = Seq(seq)
	if translation is None:
		if suppress_amber: translation = seq.translate(to_stop=True, table=999)
		else: translation = seq.translate(to_stop=True)

	tra = translation
	out =  (f"({len(seq):>4} nt): {str(seq)}\n"
			f"({len(tra):>4} aa): {'  '.join(a for a in str(tra))}")
	l_padding = len(f"({len(seq):>4} nt): ")

	if ruler:
		if not isinstance(ruler, dict):
			ruler = {}
		out = (" " * l_padding) + show_ruler(**{ 'start':0, 'end':len(seq), 'step_size':10, **ruler }) + "\n" + out

	for feature_name, pos in positions.items():
		feature_name = str(feature_name)
		if isinstance(pos, tuple):
			feature_len = pos[1] - pos[0]
			track = ' ' * (l_padding + pos[0])
			track = track + (feature_name + '-' * (feature_len - len(feature_name)))
			out = out + "\n" + track
		else:
			out = out + "\n" + (" " * (l_padding + pos)) + "| " + feature_name
	return out

def get_reference(reference = None, reference_path = None, reference_frame_start = 0):
	if reference is None:
		if reference_path is None:
			raise Exception("Must specify either `reference` or `reference_path`")
		reference = next(Bio.SeqIO.parse(reference_path,"fasta")).seq
	return reference[reference_frame_start:]

def translate_reference(reference, suppress_amber=True):

	# read and translate reference sequence
	if suppress_amber: reference_translated = reference.translate(to_stop=True, table=999)
	else: reference_translated = reference.translate(to_stop=True)

	return reference_translated


def match_to_series(input, output):
	import pandas as pd

	if isinstance(input, pd.Series):
		return pd.Series(output, index=input.index, name=input.name)
	else:
		return output

def map_generic(mapping, values):
	"""
	Uses a function, dict-like, or pd.Series to map one set of ``values`` to another.
	"""
	import pandas as pd, numpy as np

	if callable(mapping):
		out = np.vectorize(mapping)(value)
	elif isinstance(mapping, dict):
		out = np.fromiter((mapping[x] for x in value), count=len(value))
	elif isinstance(mapping, pd.Series):
		if not mapping.index.is_unique:
			raise Exception("the mapping provided as "
			"a pd.Series has a non-unique index. This makes it ambiguous how "
			f"to map old values to new ones. Must provide a mapping where "
			"mapping.index.is_unique == True. ")
		out = mapping.loc[values]
	else: raise Exception("mapping must be callable, dict, or pd.Series")
	return out


def lookup(df, row_labels, col_labels, default=None):
	"""
	Label-based "fancy indexing" function for DataFrame.

	For some reason pd.DataFrame.lookup is deprecated, despite replacements being extremely confusing
	https://github.com/pandas-dev/pandas/issues/39171
	Source: https://github.com/pandas-dev/pandas/blob/v1.3.5/pandas/core/frame.py#L4521-L4581


	Given equal-length arrays of row and column labels, return an
	array of the values corresponding to each (row, col) pair.
	.. deprecated:: 1.2.0
		DataFrame.lookup is deprecated,
		use DataFrame.melt and DataFrame.loc instead.
		For further details see
		:ref:`Looking up values by index/column labels <indexing.lookup>`.
	Parameters
	----------
	row_labels : sequence
		The row labels to use for lookup.
	col_labels : sequence
		The column labels to use for lookup.
	Returns
	-------
	numpy.ndarray
		The found values.
	"""
	from pandas.core.dtypes.common import is_object_dtype
	from pandas._libs import lib

	n = len(row_labels)
	if n != len(col_labels):
		raise ValueError("Row labels must have same size as column labels")
	if not (df.index.is_unique and df.columns.is_unique):
		# GH#33041
		raise ValueError("DataFrame.lookup requires unique index and columns")

	thresh = 1000
	if not df._is_mixed_type or n > thresh:
		values = df.values
		ridx = df.index.get_indexer(row_labels)
		cidx = df.columns.get_indexer(col_labels)
		nulls = (ridx == -1) | (cidx == -1)
		if default is None:
			if (ridx == -1).any():
				raise KeyError("One or more row labels was not found")
			if (cidx == -1).any():
				raise KeyError("One or more column labels was not found")

		ridx[ridx == -1] = 0
		cidx[cidx == -1] = 0

		flat_index = ridx * len(df.columns) + cidx
		result = values.flat[flat_index]
		result[nulls] = default
	else:
		result = np.empty(n, dtype="O")
		for i, (r, c) in enumerate(zip(row_labels, col_labels)):
			try:
				result[i] = df._get_value(r, c)
			except Exception as e:
				if default is not None:
					result[i] = default
				else:
					raise e

	if is_object_dtype(result):
		result = lib.maybe_convert_objects(result)

	return result



class dotdict(dict):
	"""
	dot.notation access to dictionary attributes

	https://stackoverflow.com/a/23689767/4091874
	"""
	__getattr__ = dict.get
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__


def intersection_ordered(lst, st):
	""" Return items from ``lst`` that are also in ``st``, preserving the order of ``lst``
	"""
	st = set(st)
	return [l for l in lst if l in st]


def sparse_coerce(X, Y, **kwargs):
	""" Convert X to the scipy.sparse class of X1 """
	import scipy.sparse

	if scipy.sparse.issparse(X) and scipy.sparse.issparse(Y):
		if scipy.sparse.isspmatrix_csr(Y):
			return X.tocsr(**kwargs)
		elif scipy.sparse.isspmatrix_csc(Y):
			return X.tocsc(**kwargs)
		elif scipy.sparse.isspmatrix_coo(Y):
			return X.tocoo(**kwargs)
	return X

def sparse_drop_na(X1):
	""" Set nan values in a sparse matrix to implicit zeros
	"""
	import scipy.sparse

	X = X1.tocoo()
	notna = ~np.isnan(X.data)
	X_notna = scipy.sparse.coo_matrix((X.data[notna], (X.row[notna], X.col[notna])), shape=X.shape)
	return sparse_coerce(X_notna, X1)

def sparse_var(a, axis=None):
	""" Variance of sparse matrix a
	var = mean(a**2) - mean(a)**2
	"""
	a_squared = a.copy()
	a_squared.data **= 2
	return a_squared.mean(axis) - np.square(a.mean(axis))

def sparse_std(a, axis=None):
	""" Standard deviation of sparse matrix a
	std = sqrt(var(a))
	"""
	return np.sqrt(sparse_var(a, axis))

def sparse_gstd(a, axis=None):
	""" Geometric standard deviation of sparse matrix a

	for A = a_0, ... a_n, where g = geometric mean of A:
	gstd = exp(sqrt( sum( ln(a_i / g)^2 / n ) ))
	"""
	# TODO: check this
	#g = sparse_gmean_nonzero(a, axis=!axis)
	g = sparse_gmean_nonzero(a, axis=axis)
	ag = a / g
	ag.data = np.log(ag.data)**2

	# densify to row or column vector
	ag = ag.sum(axis=axis) / a.getnnz(axis=axis)
	ag = np.exp(np.sqrt(ag))
	return ag


def sparse_mean_nonzero(a, axis=None):
	""" Mean of non-zero elements in sparse matrix a

	Parameters
	----------
	a : scipy.sparse.spmatrix
		sparse matrix
	axis : int, optional
		Mean over which axis: 0 = mean of each row, 1 = mean of each column, None = mean of entire matrix

	Returns
	-------
	numpy.matrix
		matrix of appropriate dimensions, or scalar if axis=None
	"""
	return (a.sum(axis=axis).A1 / a.getnnz(axis=axis))

def sparse_gmean_nonzero(a, axis=None):
	import scipy.sparse

	# geometric_mean(x) = e^(arithmetic_mean(ln(x)))
	# c = scipy.sparse.coo_matrix(a)
	c = a.copy()
	c.data = np.log(c.data)

	# this operation densifies the matrix along ~axis
	if axis is not None:
		c = (c.sum(axis=axis).A1 / c.getnnz(axis=axis))
	else:
		c = (c.sum(axis=axis) / c.getnnz(axis=axis))
	return np.exp(c)
	# return sparse_coerce(c, a)



def sparse_rank(X, asc=True, pct=False):
	""" calculates the row-wise rank or percentile of non-zero elements of sparse matrix X
	"""
	import scipy.sparse
	from scipy.stats import rankdata

	X = sparse_drop_na(X)
	ml = X.tolil()
	data = ml.data

	# to sort in descending order, negate values before argsort
	# if asc == False:
	#     _dir = -1
	#     # can't get this to work...
	#     raise NotImplementedError()
	# else: _dir = 1


	# np.argsort ranks items in *ascending* order:
	# >>> np.argsort([5,4,3,2,1])
	# array([4, 3, 2, 1, 0])
	def rank_row(x):

		# return list(np.argsort(x))
		# return list(np.argsort(x).argsort())
		return list(rankdata(x))

	ml.data = np.array(list(map(rank_row, data)), dtype='object')

	if pct:
		# divide by row-wise nnz
		# do this weird manipulation to maintain sparsity in X2
		X = ml.tocsr()
		nnz = scipy.sparse.diags((1/X.getnnz(axis=1)).ravel())
		X2 = nnz.dot(X)

		# for descending percentile, pct_descending = 1 - pct_ascending
		if asc == False:
			X2.data = 1 - X2.data
		return X2
	else:
		if asc == False:
			raise NotImplementedError()
		return sparse_coerce(ml, X)

def sparse_product(X, log=False, axis=0):
	""" calculates the product of elements of sparse matrix X along specified axis
	"""
	X = X.copy()
	X.data = np.log10(X.data)
	if log:
		return X.sum(axis=axis).A1
	else:
		return 10**(X.sum(axis=axis)).A1


def ecdf_x_from_p(p, train=None, ecdf=None):
	""" find the lowest value x for which ECDF(x) > p
	"""
	from statsmodels.distributions.empirical_distribution import ECDF
	if ecdf is None:
		ecdf = ECDF(train)
	return min(ecdf.x[ecdf.y > p])




# https://stackoverflow.com/questions/5384570/whats-the-shortest-way-to-count-the-number-of-items-in-a-generator-iterator/34404546#34404546
from collections import deque

# Avoid constructing a deque each time, reduces fixed overhead enough
# that this beats the sum solution for all but length 0-1 inputs
consumeall = deque(maxlen=0).extend

def ilen(it):
	# Make a stateful counting iterator
	cnt = itertools.count()
	# zip it with the input iterator, then drain until input exhausted at C level
	consumeall(zip(it, cnt)) # cnt must be second zip arg to avoid advancing too far
	# Since count 0 based, the next value is the count
	return next(cnt)
