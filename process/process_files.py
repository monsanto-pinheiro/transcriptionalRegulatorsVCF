'''
Created on 11/11/2022

@author: mmp
'''

from gff3tool.lib.gff3 import Gff3
import gzip, os, sys
import socket
from utils.util import Utils
from constants.constants import Constants
from Bio import SeqIO
from Bio.Seq import Seq

### merge VCFs
## $ for i in `cat list.txt`; do echo $i; picard_tools.sh MergeVcfs I=${i}_indels.vcf.gz I=${i}_snp.vcf.gz O=merge/${i}.vcf.gz; done;
## Important:
## 1) only first changes are take into account	"Ca22chr1A_C_albicans_SC5314	691	.	A	T,G	4286.02	PASS"
## 2) * pass are deleted, are DEL 				"Ca22chr1A_C_albicans_SC5314	693	.	A	*,G	1157.01	PASS"


class Gene(object):
	
	def __init__(self, gene_name, chr_name):
		self.gene_name = gene_name
		self.chr_name = chr_name
	
	
class ProcessGFF(object):
	'''
	classdocs
	'''

	constants = Constants()
	utils = Utils("trans_factor", "/tmp") 		### make it quick
	vect_type_to_process = ['gene']
	
	def __init__(self, read_factors, sample_list = None):
		'''
		Constructor
		'''
		self.read_factors = read_factors
		
		## has the sample list for the output
		self.vect_sample_list = []
		if not sample_list is None and os.path.exists(sample_list):
			self.vect_sample_list = self.utils.read_text_file(sample_list)
		
	def get_vcf_header(self, file_name):
		"""
		return a VCF file with a slice
		"""
		### set the header
		temp_file = self.utils.get_temp_file("slice", ".vcf")
		cmd = 'zcat {} | grep "^#" > {}'.format(file_name, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file)
			raise Exception("Fail to run zcat")
		
		### get number
		cmd = 'bgzip {}; tabix -p vcf {}.gz'.format(temp_file, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file)
			raise Exception("Fail to run bgzip")
		
		return temp_file + '.gz'
	
	def run_mafft(self, input_file, out_file):
		"""
		run mafft
		out: out_file
		"""
		if socket.gethostname() == "cs-nb0008":
			cmd = "{}; {} {} --thread 2 {} > {}".format(self.constants.SOFTWARE_SET_ENV_MAFFT,\
					self.constants.SOFTWARE_MAFFT_name,\
					self.constants.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES, input_file, out_file)
		else:
			cmd = "{} {} --thread 2 {} > {}".format(self.constants.SOFTWARE_MAFFT_name,\
					self.constants.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES, input_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			raise Exception("Fail to run mafft: " + cmd)
		return out_file
	
	def run_clustalo(self, input_file, out_file, parameters = ""):
		"""
		run clustalo
		out: out_file
		"""
		cmd = "{} --force --infmt=fa --outfmt=fa --seqtype dna --MAC-RAM 8000 {} --threads=2 -i {} -o {}".format(
				self.constants.CLUSTALO,\
				parameters, input_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			raise Exception("Fail to run clustalo: " +  cmd)
		return out_file
		
		
	def correct_postions(self, dt_result_change, fasta_file_ref, fasta_file_change, start_pos,
				out_file_temp = None, input_file_temp = None):
		"""
		Correct positions
		AACATTGTGAATGGG--AACCACCA
		AACATT-TGAATGGGCAAA---CCA
		"""
		REF_NAME = 'ref'
		CHANGE_NAME = 'change'
		
		## first alignment the sequence
		out_file = self.utils.get_temp_file("clustalo_out", ".fasta") if out_file_temp is None else out_file_temp
		input_file = self.utils.get_temp_file("clustalo_in", ".fasta") if input_file_temp is None else input_file_temp
		vect_record = []
		with open(fasta_file_ref) as handle_read:
			for record in SeqIO.parse(handle_read, "fasta"):
				record.id = REF_NAME
				record.seq = Seq(str(record.seq).replace('*', ''))
				vect_record.append(record)
		with open(fasta_file_change) as handle_read:
			for record in SeqIO.parse(handle_read, "fasta"):
				record.id = CHANGE_NAME
				record.seq = Seq(str(record.seq).replace('*', ''))
				vect_record.append(record)
				
		### save files
		with open(input_file, 'w') as handle_out:
			SeqIO.write(vect_record, handle_out, "fasta")
		#self.run_clustalo(input_file, out_file)
		self.run_mafft(input_file, out_file)
		
		## read file with the alignments
		with open(out_file) as handle:
			## get both sequences
			seq_ref = ""
			seq_other = ""
			for record_dict in SeqIO.parse(handle, "fasta"):
				if (record_dict.id == REF_NAME):	## ref seq
					seq_ref = str(record_dict.seq).upper()
				else:
					seq_other = str(record_dict.seq).upper()

		### check the length		
		if len(seq_ref) != len(seq_other) or len(seq_ref) == 0:
			sys.exit("Error: alignment sequences are different sizes or zero." +\
					"\nOutput file: " + out_file +\
					"\nInput file: " + input_file)
		
		## remove file
		if out_file_temp is None: self.utils.remove_file(out_file)
		if input_file_temp is None: self.utils.remove_file(input_file)
		
		### create a dictonary with the translation
		### start on zero
		dt_translation = {}	## { pos_change_1: pos_ref_1, pos_change_2: pos_ref_2, ...}
		pos_ref = 0
		pos_change = 0
		for i in range(len(seq_ref)):
			if seq_ref[i] == '-' and seq_other[i] == '-': continue
			if seq_ref[i] == '-': pos_change += 1
			elif seq_other[i] == '-': pos_ref += 1
			else:
				pos_ref += 1
				pos_change += 1

			if not pos_change in dt_translation:
				dt_translation[pos_change] = pos_ref
			
		for key in dt_result_change:
			for trans in dt_result_change[key]:
				if (trans.position_start - start_pos) in dt_translation:
					trans.position_start = start_pos + dt_translation[trans.position_start - start_pos]
		return dt_result_change
	
	def get_slice_vcf(self, file_name, chr_name, position_start, position_end):
		"""
		return a VCF file with a slice
		"""
		### set the header
		temp_file = self.utils.get_temp_file("slice", ".vcf")
		cmd = 'zcat {} | grep "^#" > {}'.format(file_name, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file)
			raise Exception("Fail to run zcat")
		
		### get slice of VCF
		cmd = "tabix {} {}:{}-{} >> {}".format(file_name,
			chr_name, position_start, position_end, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file)
			raise Exception("Fail to run tabix")

		###
		temp_file_2 = self.utils.get_temp_file("slice", ".vcf") 
		cmd = "awk -v minus={}".format(position_start) +\
				" '{ if ( $0 ~ /^#/ ) { print $0 } else " + '{ printf(\"%s\\t%d\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\" ' +\
				",$1, $2-minus, $3, $4, $5, $6, $7, $8, $9, $10) } }'" +\
				" {} > {}".format(temp_file, temp_file_2)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file)
			self.utils.remove_file(temp_file_2)
			raise Exception("Fail to run awk: " + cmd)

		### remove first temp file
		self.utils.remove_file(temp_file)
		
		### get number
		cmd = 'bgzip {}; tabix -p vcf {}.gz'.format(temp_file_2, temp_file_2)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file_2)
			self.utils.remove_file(temp_file_2 + ".gz")
			self.utils.remove_file(temp_file_2 + ".gz.tbi")
			raise Exception("Fail to run bgzip")
		
		return temp_file_2 + '.gz'

	def make_consensus(self, ref_file, position_to_cut, vcf_file, out_file):
		"""
		:parm position_to_cut <chr>:start:end
		"""
		
		### need to create the index
		if not os.path.exists(ref_file + ".fai"):
			cmd = "samtools faidx {}".format(ref_file)
			os.system(cmd)
			
		cmd = "samtools faidx {} {} | {} consensus {} -o {}".format(ref_file, position_to_cut,
					self.constants.BCFTOOLS, vcf_file, out_file);
		exist_status = os.system(cmd)
		if (exist_status != 0):
			raise Exception("Fail to run samtools: " + cmd)
		
	def get_diff_factors(self, dt_ref_factors, dt_change_factors):
		"""
		:param dt_ref_factors has the factors that appear in reference
		:param dt_change_factors has the factors that appear in reference+variations
		
		dt_ref_factors = {factor: pos, factor1: pos_2, factor2: pos_3}
		"""
		dt_out_lost_ref = {}
		for key in dt_ref_factors:
			## add all changes
			if not key in dt_change_factors: dt_out_lost_ref[key] = dt_ref_factors[key]
			else:	## else
				for factor in dt_ref_factors[key]:
					if not factor in dt_change_factors[key]:
						if key in dt_out_lost_ref: dt_out_lost_ref[key].append(factor)
						else: dt_out_lost_ref[key] = [factor]
		
		### other direction
		dt_out_gain_change = {}
		for key in dt_change_factors:
			## add all changes
			if not key in dt_ref_factors: dt_out_gain_change[key] = dt_change_factors[key]
			else:	## else
				for factor in dt_change_factors[key]:
					if not factor in dt_ref_factors[key]:
						if key in dt_out_gain_change: dt_out_gain_change[key].append(factor)
						else: dt_out_gain_change[key] = [factor]
		return dt_out_lost_ref, dt_out_gain_change

	def save_result_file(self, dt_result, file_out):
		"""
		Save the results
		"""
		b_gain_file = True if file_out.find(Constants.FILE_NAME_OUT_GAIN) != -1 or \
					file_out.find(Constants.FILE_NAME_OUT_EXTENDED_GAIN) != -1 else False
		
		### header Sample, then factors
		self.utils.make_path(os.path.dirname(file_out))
		
		## process file
		with open(file_out, 'w') as handle_out:
			
			### header Sample, then factors
			handle_out.write("\t" + "\t".join([sample_name + "\t" * (self.read_factors.get_length() -1) for sample_name in self.vect_sample_list]) + "\n")
			handle_out.write("Genes" + ("\t" + "\t".join(self.read_factors.get_names())) * len(self.vect_sample_list) + "\n")
		
			## 
			for gene in dt_result:
				has_data = False
				sz_out = "{}".format(gene)
				
				if len(dt_result[gene]) == 0:
					handle_out.write("\n")
					continue
				
				for sample in self.vect_sample_list:
					
					for trans_factor in self.read_factors.get_names():
						if trans_factor in dt_result[gene][sample][1 if b_gain_file else 0]:
							has_data = True
							if file_out.find(Constants.FILE_NAME_OUT_EXTENDED_GAIN) != -1 or \
								file_out.find(Constants.FILE_NAME_OUT_EXTENDED_LOST) != -1:
								sz_out += "\t" + ";".join([trans_factor_temp.get_info() for trans_factor_temp \
										in dt_result[gene][sample][1 if b_gain_file else 0][trans_factor]])
							else:
								sz_out += "\t{}".format(len(dt_result[gene][sample][1 if b_gain_file else 0][trans_factor]))
						else:
							if file_out.find(Constants.FILE_NAME_OUT_EXTENDED_GAIN) != -1 or \
								file_out.find(Constants.FILE_NAME_OUT_EXTENDED_LOST) != -1:
								sz_out += "\t"
							else:	
								sz_out += "\t0"
				sz_out += "\n"
				if has_data: handle_out.write(sz_out)
		print("File saved: " + file_out)
		
		
		
	def process_gff(self, ref_file, gff_file, list_vcf_files, size_window, out_path, dt_gene_names = {}):
		"""
		Read a gff file
		:param list of VCD files
		:param dt_gene_names, list of genes to process, or, if empty, process all
		"""
		
		dt_result = {}		## { SeqID : { sample_name : { lost_ref : [transfactor_1, transfactor_2... ], lost_gain : [transfactor_1, transfactor_2... ]}, {}),
							##   SeqID_2 : { sample_name : [transfactor_1, transfactor_15... ]  
		count_genes = 0
		with (gzip.open(gff_file, mode='rt') if self.utils.is_gzip(gff_file) else open(gff_file, mode='r')) as handle_read:
			gff = Gff3(handle_read)
		
			### fasta files
			slice_change_file = self.utils.get_temp_file("slice_change_file", ".fasta") 
			slice_ref_file = self.utils.get_temp_file("slice_ref_file", ".fasta")
			empty_vcf_file = self.get_vcf_header(list_vcf_files[0])
			
			for line_gff in gff.lines:
				## {'line_index': 34, 'line_raw': 'chrI\tS01\tTY1/TY2_soloLTR\t36933\t37200\t.\t+\t.\tID=TY1/TY2_soloLTR:chrI:36933-37200:+;Name=TY1/TY2_soloLTR:chrI:36933-37200:+\n', 
				## 'line_status': 'normal', 'parents': [], 'children': [], 'line_type': 'feature', 'directive': '', 'line_errors': [], 'type': 'TY1/TY2_soloLTR', 'seqid': 'chrI', 'source': 'S01', 'start': 36933, 'end': 37200, 'score': '.', 'strand': '+', 'phase': '.', 
				## 'attributes': {'ID': 'TY1/TY2_soloLTR:chrI:36933-37200:+', 'Name': 'TY1/TY2_soloLTR:chrI:36933-37200:+'}}
				if line_gff['line_type'] == 'feature' and (line_gff['type'] in self.vect_type_to_process):

					gene_name = line_gff['attributes']['ID']
					### check gene names
					if (len(dt_gene_names) > 0 and not gene_name in dt_gene_names): continue
					count_genes += 1
					
					chr_name, position_start, position_end = line_gff['seqid'], line_gff['start'], line_gff['end']
					if line_gff['strand'] == '+':	## positve
						position_end = position_start - 1
						position_start -= size_window
					else: 
						position_start = position_end + 1
						position_end += size_window 
					
					size_window_temp = size_window
					if position_start < 1:
						position_start = 1
						size_window_temp = position_end - position_start
						
					##  info	
					print("Processing gene ({}): ".format(count_genes) + \
						gene_name + "  {}:{}-{}".format(chr_name, position_start, position_end) + \
						"  Size window: {}".format(size_window_temp))
					
					### for each VCF file, of this gene, check
					dt_result_vcf = {}
					b_add_list_samples = True if len(self.vect_sample_list) == 0 else False

					## fasta file with ref
					self.make_consensus(ref_file, "{}:{}-{}".format(chr_name, position_start, position_end), empty_vcf_file, slice_ref_file)
					### get trans factors for ref file
					dt_ref_factors = self.read_factors.get_match_factors_from_fasta(slice_ref_file, size_window_temp)
						
					for vcf_file in list_vcf_files:
					
						sample_name = os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')
						if b_add_list_samples:
							self.vect_sample_list.append(sample_name)
						elif not sample_name in self.vect_sample_list:
							print("Warning, VCF file not present in sample list names - " +  sample_name)
							continue
							
						## slice_vcf_file_corrected = self.get_slice_vcf(vcf_file, chr_name, position_start, position_end)
						
						## fasta file with changes 
						self.make_consensus(ref_file, "{}:{}-{}".format(chr_name, position_start, position_end), vcf_file, slice_change_file)

						## Important, because the position is based on chromosome there is not necessary to reverse complement
						### get trans factors for change file
						dt_change_factors = self.read_factors.get_match_factors_from_fasta(slice_change_file, size_window_temp)
						## because of insertions and deletions
						dt_change_factors = self.correct_postions(dt_change_factors, slice_ref_file, slice_change_file, position_start)

						### get the difference
						dt_result_vcf[os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')] = \
							self.get_diff_factors(dt_ref_factors, dt_change_factors)
					dt_result[gene_name] = dt_result_vcf
					
					## test limit
					if count_genes > 100: break
					
			### save results...
			self.utils.remove_file(empty_vcf_file)
			self.utils.remove_file(slice_change_file)
			self.utils.remove_file(slice_ref_file)


		### save files
		self.save_result_file(dt_result, os.path.join(out_path, Constants.FILE_NAME_OUT_GAIN))
		self.save_result_file(dt_result, os.path.join(out_path, Constants.FILE_NAME_OUT_LOST))
		self.save_result_file(dt_result, os.path.join(out_path, Constants.FILE_NAME_OUT_EXTENDED_GAIN))
		self.save_result_file(dt_result, os.path.join(out_path, Constants.FILE_NAME_OUT_EXTENDED_LOST))
		

