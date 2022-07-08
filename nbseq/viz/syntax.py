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

# class AminoAcidLexer(RegexLexer):
#     name = 'FASTA format with clustal colors'

#     tokens = {
#         'root': [
#             (r'[^/]+', Text),
#             (r'/\*', Comment.Multiline, 'comment'),
#             (r'//.*?$', Comment.Singleline),
#             (r'/', Text)
#         ],
#         'comment': [
#             (r'[^*/]', Comment.Multiline),
#             (r'/\*', Comment.Multiline, '#push'),
#             (r'\*/', Comment.Multiline, '#pop'),
#             (r'[*/]', Comment.Multiline)
#         ]
#     }


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
		return self.formatter.get_style_defs()

	def get_style_sheet_HTML(self):
		return f"<style>{self.get_style_defs()}</style>"

	def setup_notebook(self):
		from IPython.display import HTML, display
		display(HTML(self.get_style_sheet_HTML()))

def AminoAcid():
	lexer = AminoAcidLexer()
	formatter = HtmlFormatter(cssclass="source", style=AminoAcidChemistryAmbigStyle)
	return Highlighter(lexer, formatter)

aa_highlighter = AminoAcid()
