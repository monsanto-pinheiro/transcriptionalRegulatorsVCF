'''
Created on 13/11/2018

@author: mmp
'''
from Bio import SeqIO
from constants.constants import Constants
from gff3tool.lib.gff3 import Gff3

try:
	import matplotlib.pyplot as plt
except ModuleNotFoundError as e:
	print("Warning: " + str(e))
	
import getpass, os, random, stat, gzip, re
import sys, time

class Utils(object):
	'''
	classdocs
	'''
	
	TEMP_DIR = os.getenv("TMP", "/tmp")

	def __init__(self, project_name = None, temp_dir = None):
		'''
		param: project_name -> used in temp diretories
		param: temp_dir -> used in temp diretories
		
		'''
		self.project_name = "generic" if project_name is None else project_name 
		self.temp_dir = Utils.TEMP_DIR if temp_dir is None else temp_dir 
	
	def is_integer(self, n_value):
		try:
			int(n_value)
			return True
		except ValueError: 
			return False
	
	def is_float(self, d_value):
		try:
			float(d_value)
			return True
		except ValueError: 
			return False
		
	def is_string(self, sz_value):
		"""
		Test if it is string
		"""
		return (not self.is_integer(sz_value) and not self.is_float(sz_value))
	
	def read_key_from_file(self, file_name, b_convert_value_to_int = False, b_aggregate_same_key = False):
		"""
		Read key file
		# ID
		sdf
		sdfs
		sgsd
		sdgsd
		
		OR
		
		Read key file
		# ID
		sdf xpto1
		sdfs zpt1
		sgsd art1
		sdgsd lrt1
		
		return  { "sdf":1, "sdfs":1, "sgsd":1, "sdgsd":1 }
		or 
		return  { "sdf":xpto1, "sdfs":zpt1, "sgsd":art1, "sdgsd":lrt1 }
		"""
		
		dict_data = {}
		if (not os.path.exists(file_name)): return dict_data
		with open(file_name) as handle_in:
			for line in handle_in:
				sz_temp = line.strip()
				if (len(sz_temp) == 0 or sz_temp[0] == '#'): continue
				lst_data = sz_temp.split()
				if (lst_data[0] in dict_data):
					if (b_aggregate_same_key): dict_data[lst_data[0]].append(0 if len(lst_data) == 1 else (int(lst_data[1]) if b_convert_value_to_int else lst_data[1]))
					else: continue
				else:
					if (b_aggregate_same_key): dict_data[lst_data[0]] = [0 if len(lst_data) == 1 else (int(lst_data[1]) if b_convert_value_to_int else lst_data[1])]
					else: dict_data[lst_data[0]] = 0 if len(lst_data) == 1 else (int(lst_data[1]) if b_convert_value_to_int else lst_data[1])
		return dict_data
	
	def copy_file(self, sz_file_from, sz_file_to):
		""" copy a file, make a directory if does not exist"""
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "cp " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to copy file") 
			
			### set attributes to file 664
			os.chmod(sz_file_to, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
	
	def move_file(self, sz_file_from, sz_file_to):
		""" copy a file, make a directory if does not exist"""
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "mv " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to copy file") 
			
			### set attributes to file 664
			os.chmod(sz_file_to, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)

	def get_temp_file(self, file_name, sz_type):
		"""
		return a temp file name
		"""
		main_path = os.path.join(self.temp_dir, getpass.getuser(), self.project_name)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		else:
			cmd = "touch {}".format(main_path)
			os.system(cmd)
			
		while 1:
			return_file = os.path.join(main_path, file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except FileExistsError:
				pass
	
	def get_temp_file_with_path(self, main_path, file_name, sz_type):
		"""
		return a temp file name
		"""
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			return_file = os.path.join(main_path, file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except FileExistsError:
				pass
			
	def get_temp_dir_passing_root_path(self, root_path):
		"""
		return a temp directory
		"""
		if (root_path != self.get_main_path()): main_path = root_path
		else: main_path = self.get_main_path()

		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			return_path = os.path.join(main_path, "dir_" + str(random.randrange(100000000, 999999999, 10)))
			if (not os.path.exists(return_path)):
				os.makedirs(return_path)
				return return_path

	def get_temp_dir(self):
		"""
		return a temp directory
		"""
		return self.get_temp_dir_passing_root_path(self.get_main_path())

	def get_main_path(self):
		return os.path.join(self.temp_dir, getpass.getuser(), self.project_name)

	def make_path(self, path_name):
		if (len(path_name) > 0 and not os.path.isdir(path_name) and not os.path.isfile(path_name)):
			cmd = "mkdir -p " + path_name
			os.system(cmd)
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to make a path") 


	def test_file_exists(self, file_name):
		if (os.path.exists(file_name)): return file_name
		sys.exit("Error: file does not exist - " + file_name)

	def remove_file(self, file_name):
		if (not file_name is None and os.path.exists(file_name)): os.unlink(file_name)
	
	def remove_dir(self, path_name):
		if (not path_name is None and os.path.isdir(path_name)):
			main_path = self.get_main_path()
			if path_name == main_path or path_name == (main_path + "/"): cmd = "rm -r {}/*".format(path_name)
			else: cmd = "rm -r {}*".format(path_name)
			os.system(cmd)

	def unzip(self, file_from, file_to):
		"""
		unzip files
		"""
		cmd = "gzip -cd {} > {}".format(file_from, file_to)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			raise Exception("Fail to unzip file: {}".format(cmd)) 

	def compress_file(self, file_from, file_to):
		"""
		compress file
		"""
		cmd = "gzip -c {} > {}".format(file_from, file_to)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			raise Exception("Fail to unzip file: {}".format(cmd)) 

	
	def is_fasta(self, sz_file_name):
		"""
		Test Fata file
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		handle = open(sz_file_name)
		b_pass = False
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp[0] == ">"): 
				b_pass = True
				break
			else: 
				handle.close()
				raise IOError(_("Error: the file is not in FASTA format."))
		handle.close()
		if (not b_pass): raise IOError(_("Error: file is not in FASTA format."))

		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
		raise IOError("Error: file is not in FASTA format.")


	def read_text_file(self, file_name):
		"""
		read text file and put the result in an vector
		"""
		if (not os.path.exists(file_name)):
			raise IOError(_("Error: file '" + file_name + "' doens't exist."))
		
		vect_out = []
		with open(file_name) as handle: 
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				vect_out.append(sz_temp)
		return vect_out


	def is_gzip(self, file_name): 
		"""
		test if the file name ends in gzip
		"""
		return True if (file_name.rfind(".gz") == len(file_name) - 3) else False
	
	def get_file_name_without_extension(self, file_name):
		"""
		return file name without extension
		"""
		return os.path.splitext(os.path.basename(file_name))[0]

	
	### return (b_value, b_is_correct, sz_error_message)
	def get_bool_from_string(self, sz_value):
		vect_values_true = ["true", "True", "TRUE", "1", "t"]
		vect_values_false = ["false", "False", "FALSE", "0", "f"]
		for true_temp_ in vect_values_true: 
			if (sz_value == true_temp_): return (True, True, "")
		for false_temp_ in vect_values_false: 
			if (sz_value == false_temp_): return (False, True, "")
		return (False, False, "Error: only this values available for bool type (%s;%s), this is not valid %s" % (",".join(vect_values_true), ",".join(vect_values_false), sz_value))

	
	def str2bool(self, v):
		"""
		str to bool
		"""
		return v.lower() in ("yes", "true", "t", "1", "y")


	def is_supplementary_alignment(self, value):
		return (0x800 & value) > 0
	
	def is_read_reverse_strand(self, value):
		return (0x10 & value) > 0
	
	def test_equal_files(self, file_1, file_2):
		"""
		Test if these two files are equal
		"""
		temp_diff = self.get_temp_file("diff", ".txt")
		cmd = "md5sum {} {} > {}".format(file_1, file_2, temp_diff)
		os.system(cmd)
		vect_result = self.read_text_file(temp_diff)
		return len(vect_result) == 2 and vect_result[0].split()[0] == vect_result[1].split()[0]
	
class NucleotideCodes(object):
	
	def __init__(self):
		self.dt_codes = {'A': ['A'], 
						'C' : ['C'],
						'G' : ['G'],
						'T' : ['T', 'U'],
						'U' : ['T', 'U'],
						'R' : ['A', 'G'],
						'Y' : ['C', 'T'],
						'S' : ['G', 'C'],
						'W' : ['A', 'T'],
						'K' : ['G', 'T'],
						'M' : ['A', 'C'],
						'B' : ['C', 'G', 'T'],
						'D' : ['A', 'G', 'T'],
						'H' : ['A', 'C', 'T'],
						'V' : ['A', 'C', 'G'],
						'N' : ['A', 'C', 'G', 'T', 'U'],
					}
		self.vect_iupac = ['A', 'C', 'G', 'T', 'U', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
		self.vect_iupac_only_two_degenerated_bases = ['R', 'Y', 'S', 'W', 'K', 'M']
		
		self.dt_primary_bases = { 'A': 1, 'C': 1, 'T': 1, 'G': 1, 'U': 1 }
		
		## has two bases for each yupac code
		self.dt_two_codes_iupac = {}
		self.create_two_codes_iupac()

	def has_this_base(self, base_ref, base_test):
		return base_ref.upper() in self.dt_codes and base_test.upper() in self.dt_codes[base_ref.upper()]

	def create_two_codes_iupac(self):
		
		for key in self.dt_codes:
			if key == 'T': self.dt_two_codes_iupac["TT"] = "T"
			elif key == 'U': self.dt_two_codes_iupac["UU"] = "U"
			elif len(self.dt_codes[key]) == 1: self.dt_two_codes_iupac["{}{}".format(key, key)] = key
			elif len(self.dt_codes[key]) == 2: 
				self.dt_two_codes_iupac["{}{}".format(self.dt_codes[key][0], self.dt_codes[key][1])] = key
				self.dt_two_codes_iupac["{}{}".format(self.dt_codes[key][1], self.dt_codes[key][0])] = key
				
				if self.dt_codes[key][1] == "T":
					self.dt_two_codes_iupac["{}U".format(self.dt_codes[key][0])] = key
					self.dt_two_codes_iupac["U{}".format(self.dt_codes[key][0])] = key
				if self.dt_codes[key][0] == "T":
					self.dt_two_codes_iupac["{}U".format(self.dt_codes[key][1])] = key
					self.dt_two_codes_iupac["U{}".format(self.dt_codes[key][1])] = key

	def get_iupac_based_on_bases(self, base1, base2):
		""" try to find in 
		return (Base to pass, True if change to a degenerated base)"""
		if base1 == '-': return (None, False)
		
		base1 = base1.upper()
		if (not base1 in self.dt_primary_bases): return (base1, False)
		base2 = base2.upper()
		if (not base2 in self.dt_primary_bases): return (base1, False)
		if (base2 == 'U'): base2 = 'T'
		return_base = self.dt_two_codes_iupac.get(base1 + base2, base1)
		return (return_base, False if return_base in self.dt_primary_bases else True)

class CigarElement(object):
	"""
	1) Op
	2) BAM
	3) Description 
	4) Consumes query
	5) Consumes reference
	
	1 2 3                                                     4   5
	M 0 alignment match (can be a sequence match or mismatch) yes yes
	I 1 insertion to the reference yes no
	D 2 deletion from the reference no yes
	N 3 skipped region from the reference no yes
	S 4 soft clipping (clipped sequences present in SEQ) yes no
	H 5 hard clipping (clipped sequences NOT present in SEQ) no no
	P 6 padding (silent deletion from padded reference) no no
	= 7 sequence match yes yes
	X 8 sequence mismatch yes yes
	"""
	
	CIGAR_TAG_M = 'M'
	CIGAR_TAG_S = 'S'
	CIGAR_TAG_I = 'I'
	CIGAR_TAG_D = 'D'
	CIGAR_TAG_H = 'H'
	CIGAR_TAG_N = 'N'

	### add all tags that can be parsed
	DICT_CIGAR_TAGS = { CIGAR_TAG_M : 1, CIGAR_TAG_S : 1, CIGAR_TAG_I : 1, CIGAR_TAG_D : 1,\
				CIGAR_TAG_H : 1, CIGAR_TAG_N : 1 }
	
	## soft and hard clip
	DICT_CIGAR_TAGS_CLIP = { CIGAR_TAG_S : 1, CIGAR_TAG_H : 1 }

	def __init__(self, length, tag):
		self.length = length
		self.tag = tag

	def is_M(self): return self.tag == CigarElement.CIGAR_TAG_M
	def is_S(self): return self.tag == CigarElement.CIGAR_TAG_S
	def is_I(self): return self.tag == CigarElement.CIGAR_TAG_I
	def is_D(self): return self.tag == CigarElement.CIGAR_TAG_D
	def is_N(self): return self.tag == CigarElement.CIGAR_TAG_N
	def is_H(self): return self.tag == CigarElement.CIGAR_TAG_H

	def __str__(self):
		return "{} -> {}".format(self.tag, self.length)

class TransFactorSmall(object):
	""" Trans Factor, forward and reverse """
	
	def __init__(self, transfactor_name, sequence_found, position_start):
		self.name = transfactor_name
		self.sequence_found = sequence_found
		self.position_start = position_start
	
#	def __eq__(self, other):
#		return self.name == other.name and self.sequence_found == other.sequence_found \
#				and self.position_start == other.position_start

	def __eq__(self, other):
		return self.name == other.name and \
				self.position_start == other.position_start
				
	def __str__(self):
		return f"Name:{self.name} sequence:{self.sequence_found} position:{self.position_start}"
	
	def get_info(self):
		return f"{self.name}:{self.sequence_found}:{self.position_start}"

class TransFactor(object):
	""" Trans Factor, forward and reverse """
	
	constants = Constants()
	def __init__(self, name, sequence):
		self.name = name
		self.sequence = sequence
		self.forward = self.constants.ambiguos_to_unambiguous(sequence)
		self.reverse = self.constants.reverse_complement(self.constants.ambiguos_to_unambiguous(sequence))
		self.re_forward = re.compile(self.constants.ambiguos_to_unambiguous(sequence))
		self.re_reverse = re.compile(self.constants.reverse_complement(self.constants.ambiguos_to_unambiguous(sequence)))
	
		## create the bit forward and reverse
		self.create_bit_sequence()
		
	def get_pattern(self, b_reverse):
		return self.re_reverse if b_reverse else self.re_forward

	def create_bit_sequence(self):
		"""
		Create a dict with index and bit { 0 : 0b0010,  1 : 0b0001, 2 : 0b01000, 3 : 0b00100, ...}
		"""
		self.dt_bit_sequence = self.constants.ambiguos_to_unambiguous_dict_bit(self.sequence)
	
	def is_match(self, position, bit_base):
		"""
		:out True if match with the base
		"""
		if position >= len(self.dt_bit_sequence): return True
		if int(bit_base & self.dt_bit_sequence[position]) > 0: return True
		return False
		
	def get_match_bit(self, post_start, sequence_fasta):
		"""
		:out vect_data [TransFactorSmall, TransFactorSmall, ...]
		"""
		vect_out = []
		for index in range(len(sequence_fasta) - len(self.dt_bit_sequence) + 1):
			b_found = True
			for index_trans in range(len(self.dt_bit_sequence)):
				try:
					if int(self.dt_bit_sequence[index_trans] & self.constants.dt_ambigous_bits[sequence_fasta[index + index_trans]]) == 0:
						b_found = False
						break
				except KeyError as e:
					b_found = False
					break
			
			if b_found:
				vect_out.append(TransFactorSmall(self.name, sequence_fasta[index: index + index_trans + 1], index + post_start))
		return vect_out
	
class TransFactors(object):
	""" dictonary with Trans Factors """
	def __init__(self, file_name):
		self.file_name = file_name
		self.dt_trans_factors = {}
		self.dt_trans_factors_repeated = {}
		self.read_file()
		
	def read_file(self):
		
		utils = Utils()
		with (gzip.open(self.file_name, mode='rt') if utils.is_gzip(self.file_name) else open(self.file_name, mode='r')) as handle_read:
			try:
				for record in SeqIO.parse(handle_read, "fasta"):
					self.dt_trans_factors[record.id] = TransFactor(record.id, str(record.seq))
			except ValueError:
				pass
		
		### test if there is any equal sequence
		dt_test = {}
		for key in self.dt_trans_factors:
			if self.dt_trans_factors[key].sequence in dt_test:
				dt_test[self.dt_trans_factors[key].sequence].append(key)
			else: dt_test[self.dt_trans_factors[key].sequence] = [key]
		
		### has the repeated
		self.dt_trans_factors_repeated = { key: dt_test[key] for key in dt_test if len(dt_test[key]) > 1}

	def get_length(self):
		return len(self.dt_trans_factors)
	
	def get_names(self):
		return sorted(list(self.dt_trans_factors.keys()))


	def get_match_factors_from_fasta(self, slice_change_file, window, b_reverse = False):
		"""
		Remove * from the sequence
		:param b_reverse, true if reversed
		try to find all trans factors
		:out dt_trans_out: {tras_name : [TransFactor + pos1, TransFactor + pos2], 
						tras_name_2 : [TransFactor + pos1, TransFactor + pos2]
						}
		"""
		
		dt_trans_out = {}
		with open(slice_change_file) as handle_read:
			for record in SeqIO.parse(handle_read, "fasta"):
				pos_start, pos_end = int(record.id.split(':')[1].split('-')[0]), int(record.id.split(':')[1].split('-')[1])
				if b_reverse: sequence_fasta = str(record.seq.reverse_complement()).replace("*", "")
				else: sequence_fasta = str(record.seq).replace("*", "")
				
				for key in self.dt_trans_factors:
					for findRe in re.finditer(self.dt_trans_factors[key].get_pattern(b_reverse), sequence_fasta):
						n_start = findRe.start()
						end = findRe.end()
						trans_factor = TransFactorSmall(self.dt_trans_factors[key].name, str(findRe.group()), 
							pos_start + window - end if b_reverse else pos_start + n_start)
						if key in dt_trans_out: dt_trans_out[key].append(trans_factor)
						else: dt_trans_out[key] = [trans_factor]
		return dt_trans_out
	
	
	def get_match_factors_from_fasta_through_bit(self, slice_change_file):
		"""
		Remove * from the sequence
		:param b_reverse, true if reversed
		try to find all trans factors
		:out dt_trans_out: {tras_name : [TransFactor + pos1, TransFactor + pos2], 
						tras_name_2 : [TransFactor + pos1, TransFactor + pos2]
						}
		"""
		dt_trans_out = {}
		with open(slice_change_file) as handle_read:
			for record in SeqIO.parse(handle_read, "fasta"):
				#pos_start, pos_end = int(record.id.split(':')[1].split('-')[0]), int(record.id.split(':')[1].split('-')[1])
				pos_start = int(record.id.split(':')[1].split('-')[0])
				sequence_fasta = str(record.seq).replace("*", "")
				
				## only works on forward
				for key in self.dt_trans_factors:
					vect_result = self.dt_trans_factors[key].get_match_bit(pos_start, sequence_fasta)
					if len(vect_result) > 0: dt_trans_out[key] = vect_result
		return dt_trans_out


class Gene(object):
	
	def __init__(self, name, chromosome, start, end, strand):
	
		self.name = name
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.strand = strand
		
	def is_forward(self):
		return self.strand == '+'
	
	def __str__(self):
		return "Gene: {}  Position: {}:({}) {}-{}".format(self.name, self.chromosome,
				self.strand, self.start, self.end)

class Genes(object):
	
	## type to process
	vect_type_to_process = ['gene']
	utils = Utils()
	
	def __init__(self, file_name_gff):

		self.dt_chromossomes = {}
		self._process_genes(file_name_gff)
	
	def _process_genes(self, file_name_gff):
		
		with (gzip.open(file_name_gff, mode='rt') if self.utils.is_gzip(file_name_gff) else open(file_name_gff, mode='r')) as handle_read:
			gff = Gff3(handle_read)
		
			for line_gff in gff.lines:
				## {'line_index': 34, 'line_raw': 'chrI\tS01\tTY1/TY2_soloLTR\t36933\t37200\t.\t+\t.\tID=TY1/TY2_soloLTR:chrI:36933-37200:+;Name=TY1/TY2_soloLTR:chrI:36933-37200:+\n', 
				## 'line_status': 'normal', 'parents': [], 'children': [], 'line_type': 'feature', 'directive': '', 'line_errors': [], 'type': 'TY1/TY2_soloLTR', 'seqid': 'chrI', 'source': 'S01', 'start': 36933, 'end': 37200, 'score': '.', 'strand': '+', 'phase': '.', 
				## 'attributes': {'ID': 'TY1/TY2_soloLTR:chrI:36933-37200:+', 'Name': 'TY1/TY2_soloLTR:chrI:36933-37200:+'}}
				if line_gff['line_type'] == 'feature' and (line_gff['type'] in self.vect_type_to_process):

					gene_name = line_gff['attributes']['ID']
					chr_name, position_start, position_end = line_gff['seqid'], line_gff['start'], line_gff['end']
					
					gene = Gene(gene_name, chr_name, position_start, position_end, line_gff['strand'])
					if chr_name in self.dt_chromossomes: self.dt_chromossomes[chr_name].append(gene)
					else: self.dt_chromossomes[chr_name] = [gene]
					
		### sort genes by start position
		for chr_name in self.dt_chromossomes:
			self.dt_chromossomes[chr_name] = sorted(self.dt_chromossomes[chr_name], key= lambda gene : gene.start)
		
	def get_genes(self, chromosome, size_window = -1):
		"""
		:out list of genes
		"""
		vect_return = []
		if chromosome in self.dt_chromossomes:
			## if not has size return all
			if size_window == -1: return self.dt_chromossomes[chromosome]
			## 
			for index, gene in enumerate(self.dt_chromossomes[chromosome]):
				if gene.is_forward():
					if index == 0: continue
					distance = self.dt_chromossomes[chromosome][index].start - self.dt_chromossomes[chromosome][index - 1].end
				else:
					if index == (len(self.dt_chromossomes[chromosome]) - 1): continue
					distance = self.dt_chromossomes[chromosome][index + 1].end - self.dt_chromossomes[chromosome][index].start
				if distance >= size_window: vect_return.append(gene)
			return vect_return
		return []

	def print_all_genes_by_chr(self, chr_name):
		""" print all genes """
		for gene in self.get_genes(chr_name):
			print(gene)
	
	def get_all_chromosomes(self):
		""" return all chromosomes from dictonary """
		return list(self.dt_chromossomes.keys())


	def create_histogram(self, out_file, number_bins = 10, max_distance = 6000, window = 1000):
		"""
		create file with histogram
		"""
		
		vect_distances = []
		vect_distances_till_windown = []
		count_gene, count_genes_overlap = 0, 0
		for key in self.dt_chromossomes:
			for index, gene in enumerate(self.dt_chromossomes[key]):
				count_gene += 1
				if gene.is_forward():
					if index == 0: continue
					distance = self.dt_chromossomes[key][index].start - self.dt_chromossomes[key][index - 1].end
				else:
					if index == (len(self.dt_chromossomes[key]) - 1): continue
					distance = self.dt_chromossomes[key][index + 1].end - self.dt_chromossomes[key][index].start
				if distance > 0: vect_distances.append(distance if distance < max_distance else max_distance)
				if distance < window and distance > 0: vect_distances_till_windown.append(distance)
				if distance < 0:
					count_genes_overlap += 1
					print("Distance:", distance, "##", gene, " -> ",
						self.dt_chromossomes[key][index - 1] if gene.is_forward() \
						else self.dt_chromossomes[key][index + 1])
#		plt.hist(vect_distances, bins=number_bins)
#		plt.title('Distance between genes')
#		plt.ylabel('Number of genes')
#		plt.xlabel('Distance (distance > {} = {})'.format(max_distance, max_distance))
#		plt.savefig(out_file)
		
		plt.hist(vect_distances_till_windown, bins=number_bins)
		plt.title('Distance between genes')
		plt.ylabel('Number of genes')
		plt.xlabel('Distance (distance < window {})'.format(window))
		plt.text(10, 63, 'Occurrences: {}/{}'.format(len(vect_distances_till_windown), count_gene))
		plt.savefig(os.path.splitext(out_file)[0] + "_less_than_window_{}".format(window) + \
				os.path.splitext(out_file)[1])

		## genes overlapped
		print("Genes overlapped: {}/{}".format(count_genes_overlap, count_gene))

	def collect_gene_affected(self, promotor_window):
		"""
		If some gene has less than the promoter window, it could be in 
		"""
		pass
	
class ProcessorThreading(object):
	
	PROCESS_TO_RUN = 8
	SLEEP_TIME_BETWEEN_ALL_END = 20		##
	
	def __init__(self):
		#self.vect_manage_process = []
		self.vect_manage_process = self.PROCESS_TO_RUN * [-1]
		self.number_forks = 0
		self._is_to_fork = True

	def search_end_process(self):
		#test all process in the list
		for i in range(self.PROCESS_TO_RUN):
			if (self.vect_manage_process[i] != -1 and self.vect_manage_process[i] != 111010001):
				# print "find Process %d" % (self.vect_manage_process[i])
				try:
					(p, v) = os.waitpid(self.vect_manage_process[i], os.WNOHANG)
				except ChildProcessError as e:
					if (str(e) == "[Errno 10] No child processes"):
						self.vect_manage_process[i] = -1 # process already finish
						p = 0
					else:
						raise e
				if (p != 0): self.vect_manage_process[i] = -1 # process already finish

	def is_all_end(self):
		"""
		Test if is all end
		"""
		if not self.is_to_fork(): return True
		
		nCount = 0 
		self.search_end_process()
		for nProcess in self.vect_manage_process:
			if (nProcess == -1): nCount += 1
		if (nCount == self.PROCESS_TO_RUN): return True
		return False
	
	def is_all_end_with_time(self):
		""" Test if all end """
		if not self.is_to_fork(): return
		
		while 1:
			time.sleep(self.SLEEP_TIME_BETWEEN_ALL_END)
			if (self.is_all_end()): break

	def get_pos_process(self):
		## not to fork
		if not self.is_to_fork(): return 0
		
		while 1:
			self.search_end_process()
			for i in range(self.PROCESS_TO_RUN):
				if (self.vect_manage_process[i] == -1 and self.vect_manage_process[i] != 111010001):
					self.vect_manage_process[i] = 111010001	#put some data to prevent to others alloc in the same spot
					return i
			time.sleep(10) # in seconds
		return -1
	
	def add_number_forks(self):
		if self.is_to_fork(): self.number_forks += 1
		else: self.number_forks = 1

	def get_number_forks(self):
		if not self.is_to_fork(): return 1
		return self.number_forks
	
	def set_not_forking(self):
		""" set fork off, good for testing """
		self._is_to_fork = False
		
	def is_to_fork(self):
		return self._is_to_fork
