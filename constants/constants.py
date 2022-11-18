'''
Created on Oct 15, 2016

@author: mmp
'''
import socket
from Bio.Align.Applications._ClustalOmega import ClustalOmegaCommandline

class Constants(object):
	'''
	classdocs
	'''

	vect_ambigous = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '*']
	dt_ambigous = { 'R':'[AG]', 'Y':'[TC]', 'K':'[GT]', 'M':'[AC]', 'S':'[GC]', 
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

	## software paths
	BCFTOOLS = "/home/projects/flu/software/bcftools-1.5/bcftools" if socket.gethostname() == "cs-nb0008" else "bcftools"
	
	## software paths
	CLUSTALO = "clustalo"
	
	SOFTWARE_PATH = "/usr/local/software/insaflu/mafft-7.453-without-extensions"
	SOFTWARE_SET_ENV_MAFFT = "export MAFFT_BINARIES={}/binaries".format(SOFTWARE_PATH)
	SOFTWARE_MAFFT_name = "{}/scripts/mafft".format(SOFTWARE_PATH)
	
	SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES = "--maxiterate 1000 --localpair --preservecase --leavegappyregion --quiet"
	
	## File name out
	FILE_NAME_OUT_GAIN = "trans_factors_gain.tsv"
	FILE_NAME_OUT_LOST = "trans_factors_lost.tsv"
	
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
		
	###
	def ambiguos_to_unambiguous(self, sequence):
		### solve this problem "GGAGGC(G/A)C(T/A)G"
		sequence = sequence.replace('(', '[').replace(')', ']').replace('/', '')
		
		### solve this problem HVN{3}
		if sequence.find('{') != -1:
			sequence_out = ""
			start_bracket = -1
			numbers = ''
			for index, char in enumerate(sequence):
				if index == 0:
					if char == '{': raise("Sequence starts with '{'")
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
						if start_bracket == -1:  raise("Fail, number without open brackets...")
						numbers += char
					else:
						if start_bracket != -1:  raise("Fail, base with open brackets...")
						sequence_out += char
			
			if start_bracket != -1:  raise("Fail, base with open brackets...")
			sequence = sequence_out

		for ambig in self.vect_ambigous:
			sequence = sequence.replace(ambig, self.dt_ambigous[ambig])
		return sequence

	def get_diff_between_two_seq(self, seq1, seq2):
		if (len(seq1) != len(seq2)): return 0
		n_diff = 0
		for i in range(0, len(seq1)):
			if (seq1[i] != seq2[i]): n_diff += 1
		return n_diff

