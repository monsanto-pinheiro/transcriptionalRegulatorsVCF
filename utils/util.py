'''
Created on 13/11/2018

@author: mmp
'''
from Bio import SeqIO
from constants.constants import Constants
import getpass, os, random, stat, gzip, re
import sys

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
	
	def __init__(self, transfactor, sequence_found, position_start):
		self.name = transfactor.name
		self.sequence_found = sequence_found
		self.position_start = position_start
	
	def __eq__(self, other):
		return self.name == other.name and self.sequence_found == other.sequence_found \
				and self.position_start == other.position_start

	def __str__(self):
		return f"Name:{self.name} sequence:{self.sequence_found} position:{self.position_start}"

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
	
	def get_pattern(self, b_reverse):
		return self.re_reverse if b_reverse else self.re_forward

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
						trans_factor = TransFactorSmall(self.dt_trans_factors[key], str(findRe.group()), 
							pos_start + window - end if b_reverse else pos_start + n_start)
						if key in dt_trans_out: dt_trans_out[key].append(trans_factor)
						else: dt_trans_out[key] = [trans_factor]
		return dt_trans_out


