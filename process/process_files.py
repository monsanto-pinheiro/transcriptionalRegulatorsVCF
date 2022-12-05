'''
Created on 11/11/2022

@author: mmp
'''

import os, sys, time
import socket
from utils.util import Utils, ProcessorThreading, Genes
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
	
	
class ProcessGFF(ProcessorThreading):
	'''
	classdocs
	'''
	
	PREFIX_NAME = "slice"
	constants = Constants()
	utils = Utils("trans_factor", "/tmp") 		### make it quick
	
	def __init__(self, read_factors, sample_list = None):
		'''
		Constructor
		'''
		super(ProcessGFF, self).__init__()
		self.read_factors = read_factors
		
		## has the sample list for the output
		self.vect_sample_list = []
		if not sample_list is None and os.path.exists(sample_list):
			self.vect_sample_list = self.utils.read_text_file(sample_list)
		
		### not to fork
		# self.set_not_forking()
		
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
		
		
	def correct_positions(self, dt_result_change, fasta_file_ref, fasta_file_change, start_pos,
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
		for factor in self.read_factors.get_names():
			file_out_temp = os.path.join(os.path.dirname(file_out), factor, os.path.basename(file_out))
			self.utils.make_path(os.path.dirname(file_out_temp))
			with open(file_out_temp, 'w') as handle_out:
			
				### header Sample, then factors
				handle_out.write("Genes\t" + "\t".join(self.vect_sample_list) + "\n")
				#handle_out.write("Genes" + ("\t" + "\t".join([factor] * len(self.vect_sample_list)) + "\n"))
			
				###
				lines_saved = 0
				for gene in dt_result:
					has_data = False
					sz_out = "{}".format(gene)
					
					if len(dt_result[gene]) == 0:
						handle_out.write("\n")
						continue
					
					for sample in self.vect_sample_list:
						if factor in dt_result[gene][sample][1 if b_gain_file else 0]:
							has_data = True
							lines_saved += 1
							if file_out.find(Constants.FILE_NAME_OUT_EXTENDED_GAIN) != -1 or \
								file_out.find(Constants.FILE_NAME_OUT_EXTENDED_LOST) != -1:
								sz_out += "\t" + ";".join([trans_factor_temp.get_info() for trans_factor_temp \
										in dt_result[gene][sample][1 if b_gain_file else 0][factor]])
							else:
								sz_out += "\t{}".format(len(dt_result[gene][sample][1 if b_gain_file else 0][factor]))
						else:
							if file_out.find(Constants.FILE_NAME_OUT_EXTENDED_GAIN) != -1 or \
								file_out.find(Constants.FILE_NAME_OUT_EXTENDED_LOST) != -1:
								sz_out += "\t"
							else:	
								sz_out += "\t0"
					sz_out += "\n"
					if has_data: handle_out.write(sz_out)
			
			### if none of factor saved it will remove the file
			if lines_saved == 0: self.utils.remove_file(file_out_temp)
				
		#print("File saved: " + file_out_temp)
		
	
	def join_files(self, genes, out_path, file_name, b_remove_slice_files = False):
		"""
		join all files
		"""
		## process file
		dt_out = {}
		for factor in self.read_factors.get_names():
			
			first_file = True
			file_out = os.path.join(out_path, factor, file_name)
			lines_per_factor = 0
			for index_chromosome in range(len(genes.get_all_chromosomes())):
				
				file_out_temp = os.path.join(out_path, factor, "{}_{}_{}".format(self.PREFIX_NAME,
					index_chromosome, file_name))
				if not os.path.exists(file_out_temp): continue
				if first_file: cmd = "cat {} >> {}".format(file_out_temp, file_out)
				else: cmd = "tail -n+2 {} >> {}".format(file_out_temp, file_out)
				os.system(cmd)
				first_file = False
				
				lines_per_factor += len(self.utils.read_text_file(file_out_temp)) - 1

			## set stats			
			dt_out[factor] = lines_per_factor
			
			## remove all slice files, only create a file to no raise exception in remove of all "slice_" files
			if (b_remove_slice_files):
				cmd = "touch {}__".format(os.path.join(out_path, factor, self.PREFIX_NAME))
				os.system(cmd)
				cmd = "rm {}*".format(os.path.join(out_path, factor, self.PREFIX_NAME))
				os.system(cmd)
		return dt_out	## result of factor
		
		
	def process_gff(self, ref_file, gff_file, list_vcf_files, size_window, out_path, dt_gene_names = {},
				only_genes_with_distance = -1):
		"""
		Read a gff file
		:param list of VCD files
		:param dt_gene_names, list of genes to process, or, if empty, process all
		"""
		
		count_genes = 0
		
		### read all genes
		genes = Genes(gff_file)
	
		### empty VCF file, for all threads
		empty_vcf_file = self.get_vcf_header(list_vcf_files[0])
		
		## for each chromosome
		for index_chromosome, chromosome in enumerate(genes.get_all_chromosomes()): 
			
			### set data to forks	
			n_pos_vect = self.get_pos_process()
			if (n_pos_vect == -1):
				raise Exception("Error getting a process ID {}".format(time.ctime()))
		
			### count number of files
			self.add_number_forks()
		
			###
			if self.is_to_fork(): new_ID = os.fork()
			else: new_ID = 0
			
			if new_ID < 0:
				raise Exception("Error forking the main process, time: %s"  % time.ctime())
			elif new_ID == 0: # is the child
				
				### fasta files
				dt_result = {}		## { SeqID : { sample_name : { lost_ref : [transfactor_1, transfactor_2... ], lost_gain : [transfactor_1, transfactor_2... ]}, {}),
							##   SeqID_2 : { sample_name : [transfactor_1, transfactor_15... ]  
				slice_change_file = self.utils.get_temp_file("slice_change_file", ".fasta") 
				slice_ref_file = self.utils.get_temp_file("slice_ref_file", ".fasta")
				out_file = self.utils.get_temp_file("clustalo_out", ".fasta")
				input_file = self.utils.get_temp_file("clustalo_in", ".fasta")
			
				## for all gene
				for gene in genes.get_genes(chromosome, only_genes_with_distance):
					### check gene names
					if (len(dt_gene_names) > 0 and not gene.name in dt_gene_names): continue
					count_genes += 1
					
					position_start, position_end = gene.start, gene.end
					if gene.is_forward():	## positve
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
						gene.name + "  {}:{}-{}".format(chromosome, position_start, position_end) + \
						"  Size window: {}".format(size_window_temp))
					
					### for each VCF file, of this gene, check
					dt_result_vcf = {}
					b_add_list_samples = True if len(self.vect_sample_list) == 0 else False

					## fasta file with ref
					self.make_consensus(ref_file, "{}:{}-{}".format(chromosome, position_start, position_end), empty_vcf_file, slice_ref_file)
					### get trans factors for ref file
					#dt_ref_factors = self.read_factors.get_match_factors_from_fasta(slice_ref_file, size_window_temp)
					dt_ref_factors = self.read_factors.get_match_factors_from_fasta_through_bit(slice_ref_file)
						
					for vcf_file in list_vcf_files:
					
						sample_name = os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')
						if b_add_list_samples:
							self.vect_sample_list.append(sample_name)
						elif not sample_name in self.vect_sample_list:
							print("Warning, VCF file not present in sample list names - " +  sample_name)
							continue
							
						## slice_vcf_file_corrected = self.get_slice_vcf(vcf_file, chr_name, position_start, position_end)
						
						## fasta file with changes 
						self.make_consensus(ref_file, "{}:{}-{}".format(chromosome, position_start, position_end), vcf_file, slice_change_file)

						## Important, because the position is based on chromosome there is not necessary to reverse complement
						### get trans factors for change file
						#dt_change_factors = self.read_factors.get_match_factors_from_fasta(slice_change_file, size_window_temp)
						dt_change_factors = self.read_factors.get_match_factors_from_fasta_through_bit(slice_change_file)
						## because of insertions and deletions
						dt_change_factors = self.correct_positions(dt_change_factors, slice_ref_file, slice_change_file,
									position_start, out_file, input_file)

						### get the difference
						dt_result_vcf[os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')] = \
							self.get_diff_factors(dt_ref_factors, dt_change_factors)
					
					### semaphore
					dt_result[gene.name] = dt_result_vcf
					
					## test limit
					#if count_genes > 5: break
			
				## remove temp files
				self.utils.remove_file(slice_change_file)
				self.utils.remove_file(slice_ref_file)
				self.utils.remove_file(input_file)
				self.utils.remove_file(out_file)
				
				### save files
				self.save_result_file(dt_result, os.path.join(out_path, "{}_{}_".format(self.PREFIX_NAME,
						index_chromosome) + Constants.FILE_NAME_OUT_GAIN))
				self.save_result_file(dt_result, os.path.join(out_path, "{}_{}_".format(self.PREFIX_NAME,
						index_chromosome) + Constants.FILE_NAME_OUT_LOST))
				self.save_result_file(dt_result, os.path.join(out_path, "{}_{}_".format(self.PREFIX_NAME,
						index_chromosome) + Constants.FILE_NAME_OUT_EXTENDED_GAIN))
				self.save_result_file(dt_result, os.path.join(out_path, "{}_{}_".format(self.PREFIX_NAME,
						index_chromosome) + Constants.FILE_NAME_OUT_EXTENDED_LOST))
				
				### kill the thread
				if self.is_to_fork(): sys.exit(0)
				
			### 
			else: self.vect_manage_process[n_pos_vect] = new_ID # is the prent but the ID is from the child
			
		## test number of forks
		if self.get_number_forks() == 0:
			raise Exception("There's no coverage data to process.")

		print("Number of forks: {}".format(self.get_number_forks()))
		### waiting all process
		self.is_all_end_with_time()

		### remove temp file
		self.utils.remove_file(empty_vcf_file)

		### join files
		dt_out = {}
		for index, type_file in enumerate(Constants.VECT_FILES):
			dt_out[type_file] = self.join_files(genes, out_path, type_file,
				index == (len(Constants.VECT_FILES) - 1))
		
		## save stats
		file_name = os.path.join(out_path, Constants.FILE_NAME_RESULTS_STATS)
		with open(file_name, 'w') as handle_out:
			handle_out.write("Number of genes identified with transcription factors in all samples\n" +\
							"Factor\t{}\t{}\t{}\t{}\n".format(
				os.path.splitext(Constants.FILE_NAME_OUT_GAIN)[0],
				os.path.splitext(Constants.FILE_NAME_OUT_LOST)[0],
				os.path.splitext(Constants.FILE_NAME_OUT_EXTENDED_GAIN)[0],
				os.path.splitext(Constants.FILE_NAME_OUT_EXTENDED_LOST)[0]))

			## roll all factors
			for factor in self.read_factors.get_names():
				sz_out = "{}".format(factor)
				for type_file in Constants.VECT_FILES:
					if factor in dt_out[type_file]: sz_out += "\t{}".format(dt_out[type_file][factor])
					else: sz_out += "\t0"
				handle_out.write(sz_out + "\n")

		print("Stats file: " + file_name)

