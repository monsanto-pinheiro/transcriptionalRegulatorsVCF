'''
Created on Oct 15, 2016

@author: mmp
'''
import socket

class Constants(object):
	'''
	classdocs
	'''

	vect_ambigous = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '*']
	dt_ambigous = { 'R':'[AG]', 'Y':'[CT]', 'K':'[GT]', 'M':'[AC]', 'S':'[CG]', 
			'W':'[AT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]',
			'N':'[ACGT]', '*':'.' }
	dict_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
					'R': 'Y', 
					'Y': 'R', 
					'K': 'M', 
					'M': 'K', 
					'S': 'S', 
					'W': 'W', 
					'B': 'V', 
					'D': 'H',
					'H': 'D',
					'V': 'B',
					'[': ']',
					']': '[',
					'N': 'N'}

	### base to decimal
	WORD_LENGTH = 5		## in bits
	dt_ambigous_bits = { 'A' : 0b0001 , 'C' : 0b0010, 'G' : 0b0100, 'T' : 0b1000,
			'R': 0b0001 | 0b0100, 'Y': 0b0010 | 0b1000, 'K': 0b0100 | 0b1000,
			'M': 0b0001 | 0b0010 , 'S': 0b0010 | 0b0100, 
			'W': 0b0001 | 0b1000 , 'B': 0b0010 | 0b0100 | 0b1000 , 'D': 0b0001 | 0b0100 | 0b1000,
			'H': 0b0001 | 0b0010 | 0b1000, 'V': 0b0001 | 0b0010 | 0b0100,
			'N': 0b0001 | 0b0010 | 0b0100 | 0b1000,
			'*': 0b10000 } ### other than degenerated base
	

	## software paths
	BCFTOOLS = "/home/projects/flu/software/bcftools-1.5/bcftools" if socket.gethostname() == "cs-nb0008" else "bcftools"
	
	## software paths, not USED
	CLUSTALO = "clustalo"
	
	if socket.gethostname() == "cs-nb0008":
		SOFTWARE_PATH = "/usr/local/software/insaflu/mafft-7.453-without-extensions"
		SOFTWARE_SET_ENV_MAFFT = "export MAFFT_BINARIES={}/binaries".format(SOFTWARE_PATH)
		SOFTWARE_MAFFT_name = "{}/scripts/mafft".format(SOFTWARE_PATH)
	else:
		SOFTWARE_MAFFT_name = "mafft"
		
	SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES = "--maxiterate 1000 --localpair --preservecase --leavegappyregion --quiet"
	
	## File name out
	FILE_NAME_OUT_GAIN = "trans_factors_gain.tsv"
	FILE_NAME_OUT_LOST = "trans_factors_lost.tsv"
	FILE_NAME_OUT_EXTENDED_GAIN = "trans_factors_extended_gain.tsv"
	FILE_NAME_OUT_EXTENDED_LOST = "trans_factors_extended_lost.tsv"
	FILE_NAME_RESULTS_STATS = "stats_results.tsv"
	
	## All files
	VECT_FILES = [FILE_NAME_OUT_GAIN, FILE_NAME_OUT_LOST,
		FILE_NAME_OUT_EXTENDED_GAIN, FILE_NAME_OUT_EXTENDED_LOST]
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	
	### complement
	def complement(self, seq):  
		complseq = [self.dict_complement[base] if base in self.dict_complement else base for base in seq]  
		return ''.join(complseq)
	
	#reverse
	def reverse_complement(self, seq):  
		seq = list(seq)  
		seq.reverse()   
		return self.complement(''.join(seq))
	
		
	def __is_integer(self, n_value):
		try:
			int(n_value)
			return True
		except ValueError: 
			return False
	
	
	def _normalize_sequence(self, sequence):
		
		### solve this problem "GGAGGC(G/A)C(T/A)G"
		sequence = sequence.upper().replace('(', '[').replace(')', ']').replace('/', '').replace('-', '')
		
		## need to order the sequences inside the square brackets 
		if sequence.find('[') != -1:
			sz_out, sz_seq = ("", "")
			b_open_bracket = False
			for base in sequence:
				if base == '[': b_open_bracket = True
				elif base == ']':
					if b_open_bracket:
						sz_out += '[' + "".join(sorted(sz_seq)) + base
						b_open_bracket = False
						sz_seq = ""
					else:
						raise Exception("Error: close bracket without open bracket")
				elif b_open_bracket: sz_seq += base 
				else: sz_out += base
			sequence = sz_out
				
		## replace [AC] or other by R
		for key in self.dt_ambigous.keys():
			sequence = sequence.replace(self.dt_ambigous[key], key)
		
		### test if everything went well
		if sequence.find('[') != -1 or sequence.find(']') != -1:
			raise Exception("Error: sequence could not have '[' or ']'")
		
		### solve this problem HVN{3} or AAN7
		sequence_out = ""
		start_bracket = -1
		numbers = ''
		for index, char in enumerate(sequence):
			if index == 0:
				if char == '{': raise("Sequence starts with '{'")
				if self.__is_integer(char): raise("Sequence starts with a number")
				sequence_out += char
			else:
				if char == '{':
					if start_bracket == -1: start_bracket = index
					else: raise("Two open brackets in a row...")
				elif char == '}':
					if start_bracket == -1: raise("Fail, close brackets without open brackets...")
					elif len(numbers) == 0: raise("Fail, no numbers between brackets...")
					sequence_out += sequence_out[-1] * (int(numbers) - 1)
					start_bracket = -1
					numbers = ''
				elif self.__is_integer(char):
					numbers += char
				else:
					if start_bracket != -1:  raise("Fail, base with open brackets...")
					if len(numbers) > 0: sequence_out += sequence_out[-1] * (int(numbers) - 1)
					numbers = ''
					sequence_out += char
		
		if start_bracket != -1:  raise("Fail, base with open brackets...")
		if len(numbers) > 0: sequence_out += sequence_out[-1] * (int(numbers) - 1)
		sequence = sequence_out
		
		return sequence
			
	###
	def ambiguos_to_unambiguous(self, sequence):
		
		## normalization of sequence
		sequence = self._normalize_sequence(sequence)
		
		for ambig in self.vect_ambigous:
			sequence = sequence.replace(ambig, self.dt_ambigous[ambig])
		return sequence

	def ambiguos_to_unambiguous_bit(self, sequence):
		
		## normalization of sequence
		sequence = self._normalize_sequence(sequence)
		
		out_value = 0
		for seq_ in sequence:
			if seq_ in self.dt_ambigous_bits: out_value = (out_value << self.WORD_LENGTH) | self.dt_ambigous_bits[seq_]
			raise Exception("Error: base '{}' not present in dictionary.".format(seq_))
		return out_value
	
	def ambiguos_to_unambiguous_dict_bit(self, sequence):
		
		## normalization of sequence
		sequence = self._normalize_sequence(sequence)
		
		dt_out = {}
		for index in range(len(sequence)):
			try:
				dt_out[index] = self.dt_ambigous_bits[sequence[index]]
			except KeyError as e: 
				raise Exception("Error: base '{}' not present in dictionary.".format(sequence[index]))
		return dt_out


	def get_diff_between_two_seq(self, seq1, seq2):
		if (len(seq1) != len(seq2)): return 0
		n_diff = 0
		for i in range(0, len(seq1)):
			if (seq1[i] != seq2[i]): n_diff += 1
		return n_diff

