from pygments.lexer import RegexLexer
from pygments.token import Token
from pygments.style import Style

import pygments
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter

class AminoAcidLexer(RegexLexer):
	name = 'Amino acids clustal colors'

	tokens = {
		'root': [
			# FASTA headers and comment lines
			(r'>.*\n', Token.Name),
			(r'#.*\n', Token.Comment),

			# match amino acids
			# Alanine (ALA)
			(r"[Aa]+", Token.AA_A),
			# Aspartate/Asparagine (ASX)
			(r"[Bb]+", Token.AA_B),
			# Cysteine (CYS)
			(r"[Cc]+", Token.AA_C),
			# Aspartate (ASP)
			(r"[Dd]+", Token.AA_D),
			# Glutamate (GLU)
			(r"[Ee]+", Token.AA_E),
			# Phenylalanine (PHE)
			(r"[Ff]+", Token.AA_F),
			# Glycine (GLY)
			(r"[Gg]+", Token.AA_G),
			# Histidine (HIS)
			(r"[Hh]+", Token.AA_H),
			# Isoleucine (ILE)
			(r"[Ii]+", Token.AA_I),
			# Lysine (LYS)
			(r"[Kk]+", Token.AA_K),
			# Leucine (LEU)
			(r"[Ll]+", Token.AA_L),
			# Methionine (MET)
			(r"[Mm]+", Token.AA_M),
			# Asparagine (ASN)
			(r"[Nn]+", Token.AA_N),
			# Proline (PRO)
			(r"[Pp]+", Token.AA_P),
			# Glutamine (GLN)
			(r"[Qq]+", Token.AA_Q),
			# Arginine (ARG)
			(r"[Rr]+", Token.AA_R),
			# Serine (SER)
			(r"[Ss]+", Token.AA_S),
			# Threonine (THE)
			(r"[Tt]+", Token.AA_T),
			# Valine (VAL)
			(r"[Vv]+", Token.AA_V),
			# Tryptophan (TRP)
			(r"[Ww]+", Token.AA_W),
			# Tyrosine (TYR)
			(r"[Yy]+", Token.AA_Y),
			# Glutamate or Glutamine (GLX)
			(r"[Zz]+", Token.AA_Z),
			(r"[-]+", Token.Gap),

			(r".",     Token.AA_Unk)

		]
	}


class AminoAcidChemistryAmbigStyle(Style):
	default_style = ""
	styles = {
		Token.Name: 'bold #222',
		Token.Comment: '#aaa',
		Token.AA_B: '#0a0',
		Token.AA_S: '#0a0',
		Token.AA_T: '#0a0',
		Token.AA_Y: '#0a0',
		Token.AA_C: '#0a0',
		Token.AA_Z: '#0a0',

		Token.AA_J: '#a0a',
		Token.AA_N: '#a0a',
		Token.AA_Q: '#a0a',

		Token.AA_K: '#00a',
		Token.AA_R: '#00a',
		Token.AA_H: '#00a',

		Token.AA_D: '#f00',
		Token.AA_E: '#f00',

		Token.AA_G: '#000',
		Token.AA_P: '#000',
		Token.AA_A: '#000',
		Token.AA_W: '#000',
		Token.AA_F: '#000',
		Token.AA_L: '#000',
		Token.AA_I: '#000',
		Token.AA_M: '#000',
		Token.AA_V: '#000',

		Token.AA_X: '#333',
		Token.AA_Unk: '#333',
		Token.Gap: '#333',
	}



class NucleicAcidLexer(RegexLexer):
	name = 'IUPAC ambiguous amino acid colors'

	tokens = {
		'root': [
			# FASTA headers and comments
			(r'>.*\n', Token.Name),
			(r'#.*\n', Token.Comment),

			# match nucleic acids
			(r"[Aa]+", Token.NA_A),
			(r"[Cc]+", Token.NA_C),
			(r"[Tt]+", Token.NA_T),
			(r"[Gg]+", Token.NA_G),

			(r"[Uu]+", Token.NA_U),
			(r"[Rr]+", Token.NA_R),
			(r"[Yy]+", Token.NA_Y),
			(r"[Ss]+", Token.NA_S),
			(r"[Ww]+", Token.NA_W),
			(r"[Mm]+", Token.NA_M),
			(r"[Kk]+", Token.NA_K),
			(r"[Dd]+", Token.NA_D),
			(r"[Bb]+", Token.NA_B),
			(r"[Vv]+", Token.NA_V),
			(r"[Hh]+", Token.NA_H),
			(r"[Nn]+", Token.NA_N),
			(r"[Xx]+", Token.NA_X),

			(r"[-]+",  Token.Gap),
			(r".",     Token.NA_Unk)

		]
	}


class NucleicAcidBioSyntaxStyle(Style):
	default_style = ""
	styles = {
		Token.Name: 'bold #222',
		Token.Comment: '#aaa',

		Token.NA_A:   'bg:#47FF19 #000000',
		Token.NA_T:   'bg:#4192FF #000000',
		Token.NA_G:   'bg:#F09000 #000000',
		Token.NA_C:   'bg:#FF4641 #000000',
		Token.NA_U:   'bg:#8A89FF #000000',
		Token.NA_R:   'bg:#FFFE80 #000000',
		Token.NA_Y:   'bg:#E180FF #000000',
		Token.NA_S:   'bg:#FF9B80 #000000',
		Token.NA_W:   'bg:#80FFF2 #000000',
		Token.NA_M:   'bg:#CE8834 #000000',
		Token.NA_K:   'bg:#90B82C #000000',
		Token.NA_D:   'bg:#C7FFB9 #000000',
		Token.NA_B:   'bg:#F8C1C0 #000000',
		Token.NA_V:   'bg:#FFE3B9 #000000',
		Token.NA_H:   'bg:#BFD8F9 #000000',
		Token.NA_N:   'bg:#FFFFFF #000000',
		Token.NA_X:   'bg:#272822 #E6E6E6',
		Token.NA_Gap: 'bg:#272822 #E6E6E6'
	}


class NucleicAcidUnambigStyle(Style):
	default_style = ""
	styles = {
		Token.Name: 'bold #222',
		Token.Comment: '#aaa',

		Token.NA_A:   '#0a0',
		Token.NA_C:   '#00f',
		Token.NA_T:   '#f00',
		Token.NA_G:   '#000',
		Token.Gap:    '#333',
		Token.NA_Unk: '#333'
	}


class Highlighter:
	def __init__(self,lexer, formatter):
		self.lexer = lexer
		self.formatter = formatter

	def highlight(self,text):
		return pygments.highlight(text, self.lexer, self.formatter)

	def display(self, text):
		from IPython.display import HTML, display
		return HTML(self.highlight(text))

	def get_style_defs(self):
		return self.formatter.get_style_defs() + " div.source>pre {white-space:pre;}"

	def get_style_sheet_HTML(self):
		return f"<style>{self.get_style_defs()}</style>"

	def setup_notebook(self):
		from IPython.display import HTML, display
		display(HTML(self.get_style_sheet_HTML()))

def AminoAcid():
	lexer = AminoAcidLexer()
	formatter = HtmlFormatter(cssclass="source", style=AminoAcidChemistryAmbigStyle)
	return Highlighter(lexer, formatter)

def NucleicAcid():
	lexer = NucleicAcidLexer()
	formatter = HtmlFormatter(cssclass="source", style=NucleicAcidUnambigStyle)
	return Highlighter(lexer, formatter)

def NucleicAcidAmbig():
	lexer = NucleicAcidLexer()
	formatter = HtmlFormatter(cssclass="source", style=NucleicAcidBioSyntaxStyle)
	return Highlighter(lexer, formatter)

aa_highlighter = AminoAcid()
na_highlighter = NucleicAcid()