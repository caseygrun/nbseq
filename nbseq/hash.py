# class FormatWords:
# 	"""MILT VARY MAST OK SEES WENT"""
#
# 	def __init__(self, name, dict=DEFAULT_DICT):
# 		self.worddict = dict
# 		self.name = name
#
# 	def encode(self, hash, hashobj):
# 		worddict = self.worddict
# 		if self.worddict is not DEFAULT_DICT:
# 			try:
# 				worddict = self.worddict[hashobj.name]
# 			except KeyError:
# 				raise FormatError, \
# 						"alternative dictionary '%s' doesn't support " \
# 						"hash '%s'" % (self.name, hashobj.name)
# 		n = hash
# 		sum = 0L
# 		for _ in range(32):
# 			sum += n%4
# 			n >>= 2
# 		hash <<= 2
# 		hash |= sum%4
# 		words = []
# 		for _ in range(6):
# 			words.insert(0, worddict[hash%2048])
# 			hash >>= 11
# 		return ' '.join(words)
#
# 	def decode(self, hash, hashobj, detect=0):
# 		validword = 0
# 		words = hash.split()
# 		if len(words) != 6:
# 			return None
# 		worddict = self.worddict
# 		if worddict is not DEFAULT_DICT:
# 			if not detect:
# 				# When not detecting, only "words" matches. Otherwise,
# 				# a failing entry would get through all the 'words'
# 				# alternatives unnecessarily.
# 				return None
# 			try:
# 				worddict = worddict[hashobj.name]
# 			except KeyError:
# 				return None
# 			for word in words:
# 				if word not in worddict:
# 					return None
# 		try:
# 			n = 0L
# 			for word in words:
# 				n <<= 11
# 				n |= DEFAULT_RDICT[word.upper()]
# 				validword = 1
# 		except KeyError:
# 			# Alternative Dictionary Algorithm
# 			#
# 			# The RFC says that *no* words from the standard dictionary
# 			# should be used. We'll be permissive here, and accept
# 			# broken generators.
# 			#
# 			#if validword:
# 			#    return None
# 			n = 0L
# 			for word in words:
# 				# Ditto.
# 				#
# 				#if DEFAULT_RDICT.has_key(word.upper()):
# 				#    return None
# 				n <<= 11
# 				n |= hashobj.hash(word)%2048
# 		hashsum = n%4
# 		n >>= 2
# 		result = n
# 		sum = 0L
# 		for _ in range(32):
# 			sum += n%4
# 			n >>= 2
# 		if hashsum != sum%4:
# 			return None
# 		return result
