'''
Created on Oct 15, 2016

@author: mmp
'''
import unittest, os
from constants.constants import Constants
from utils.util import TransFactors, Utils
from process.process_files import ProcessGFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class Test(unittest.TestCase):

	utils = Utils()
	
	def test_transfactors(self):
		
		seq_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/yeastract.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		
		read_factors = TransFactors(seq_file_name)
		self.assertEqual(98, read_factors.get_length())
		self.assertEqual('Ace2p', read_factors.dt_trans_factors['Ace2p'].name)
		self.assertEqual('[AC][AC]CCA[CG]C', read_factors.dt_trans_factors['Ace2p'].forward)
		self.assertEqual('G[CG]TGG[GT][GT]', read_factors.dt_trans_factors['Ace2p'].reverse)
		
		self.assertEqual(5, len(read_factors.dt_trans_factors_repeated))
		self.assertEqual(['Cph2p_2', 'Tye7p_2'], read_factors.dt_trans_factors_repeated['ATCANNTGA'])
		
		sequence = "TAACCAGCCTTTCCCCACCQQTTCC"
		temp_fasta = self.utils.get_temp_file("diff", ".txt")
		### set vect records
		vect_record = [SeqRecord( Seq(sequence), id = 'first:0-20', description="")]
		with open(temp_fasta, 'w') as handle_out:
			SeqIO.write(vect_record, handle_out, "fasta")
			
		dt_result = read_factors.get_match_factors_from_fasta(temp_fasta, 0, False)
		self.assertEqual(1, len(dt_result))
		self.assertEqual(2, len(dt_result['Ace2p']))
		self.assertTrue('Ace2p' in dt_result)
		self.assertEqual('Ace2p', dt_result['Ace2p'][0].name)
		self.assertEqual('AACCAGC', dt_result['Ace2p'][0].sequence_found)
		self.assertEqual(1, dt_result['Ace2p'][0].position_start)
		self.assertTrue('Ace2p' in dt_result)
		self.assertEqual('Ace2p', dt_result['Ace2p'][1].name)
		self.assertEqual('CCCCACC', dt_result['Ace2p'][1].sequence_found)
		self.assertEqual(12, dt_result['Ace2p'][1].position_start)
		self.utils.remove_file(temp_fasta)
		
	def test_get_slice_vcf(self):
		
		seq_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/yeastract.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		vcf_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_1_S31_snps.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name))
		vcf_file_name_result = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/vcf_slice.result.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name_result))
		
		utils = Utils()
		read_factors = TransFactors(seq_file_name)
		process_files = ProcessGFF(read_factors)
		chr_name = "Ca22chr1A_C_albicans_SC5314"
		position_start = 1229	## first value in VCF to appear
		position_end = 2000
		file_out = process_files.get_slice_vcf(vcf_file_name, chr_name, position_start, position_end)
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(self.utils.test_equal_files(file_out, vcf_file_name_result))
		
		utils.remove_file(file_out)
		utils.remove_file(file_out + ".gz")

	def test_get_diff_between_two_seq(self):
		
		seq_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/yeastract.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		ref_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/ref/ref.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		vcf_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_1_S31_snps_1.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name))
		vcf_file_name_1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_1_S31_snps.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name_1))
		gff_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/gff/file.gff3")
		self.assertTrue(os.path.exists(gff_file_name))

		list_vcf_files = [vcf_file_name]
		sample_list = self.utils.get_temp_file("sample_list", ".txt")
		with open(sample_list, 'w') as handle_write:
			handle_write.write(os.path.basename(vcf_file_name).replace('.vcf.gz', ''))
			
		read_factors = TransFactors(seq_file_name)
		process_files = ProcessGFF(read_factors, sample_list)
		process_files.set_not_forking()
		
		out_file = self.utils.get_temp_file("temp_file", ".fasta")
		out_file_ref = self.utils.get_temp_file("temp_file_ref", ".fasta")
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:1-10", vcf_file_name, out_file)
		vect_out = self.utils.read_text_file(out_file)
		self.assertEqual(2, len(vect_out))
		self.assertEqual(">Ca22chr1A_C_albicans_SC5314:1-10", vect_out[0])
		self.assertEqual("GAGTCACGCC", vect_out[1])
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:5-10", vcf_file_name, out_file)
		vect_out = self.utils.read_text_file(out_file)
		self.assertEqual(2, len(vect_out))
		self.assertEqual(">Ca22chr1A_C_albicans_SC5314:5-10", vect_out[0])
		self.assertEqual("CACGCC", vect_out[1])
		
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-700", vcf_file_name, out_file)
		vect_out = self.utils.read_text_file(out_file)
		self.assertEqual(2, len(vect_out))
		self.assertEqual(">Ca22chr1A_C_albicans_SC5314:680-700", vect_out[0])
		self.assertEqual("GCTTTGTCTTTTA*GATTACA", vect_out[1])
		
		vcf_temp_file = process_files.get_vcf_header(vcf_file_name)	## clean VCF file
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-700", vcf_temp_file, out_file)
		vect_out = self.utils.read_text_file(out_file)
		self.assertEqual(2, len(vect_out))
		self.assertEqual(">Ca22chr1A_C_albicans_SC5314:680-700", vect_out[0])
		self.assertEqual("GCTTTGTCCCTAAAGATTACA", vect_out[1])
		
		### get matches
		window = 1000
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-1680", vcf_file_name_1, out_file)
		dt_result_change = read_factors.get_match_factors_from_fasta(out_file, window, False)
		self.assertEqual(19, len(dt_result_change))
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-1680", vcf_temp_file, out_file_ref)
		dt_result_change = process_files.correct_positions(dt_result_change, out_file_ref, out_file, 680)
		self.assertEqual(19, len(dt_result_change))
		dt_result_ref = read_factors.get_match_factors_from_fasta(out_file_ref, window, False)
		self.assertEqual(18, len(dt_result_ref))
		
		dt_out_lost_ref, dt_out_gain_change = process_files.get_diff_factors(dt_result_ref, dt_result_change)
		self.assertEqual(6, len(dt_out_gain_change))
		self.assertEqual("AACCACC", dt_out_gain_change['Ace2p'][0].sequence_found)
		self.assertEqual(813, dt_out_gain_change['Ace2p'][0].position_start)
		self.assertEqual(1, len(dt_out_gain_change['Ace2p']))
		self.assertEqual(1, len(dt_out_gain_change['Skn7p']))
		
		self.assertEqual(2, len(dt_out_lost_ref))
		self.assertEqual(1, len(dt_out_lost_ref['Mrr1p']))
		self.assertEqual(1, len(dt_out_lost_ref['Mrr1p_1']))
		self.assertEqual(1613, dt_out_lost_ref['Mrr1p'][0].position_start)
		self.assertEqual("AAAAAAAAAT", dt_out_lost_ref['Mrr1p'][0].sequence_found)
		
		## with changes
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-1680", vcf_file_name, out_file)
		dt_result_change = read_factors.get_match_factors_from_fasta(out_file, window, False)
		self.assertEqual(19, len(dt_result_change))
		dt_result_change = process_files.correct_positions(dt_result_change, out_file_ref, out_file, 680)
		self.assertEqual(19, len(dt_result_change))
		dt_out_lost_ref, dt_out_gain_change_1 = process_files.get_diff_factors(dt_result_ref, dt_result_change)
		
		self.assertEqual(6, len(dt_out_gain_change_1))
		self.assertEqual(2, len(dt_out_lost_ref))
		
		self.assertEqual('Ace2p', dt_out_gain_change_1['Ace2p'][0].name)
		self.assertEqual('AACCACC', dt_out_gain_change_1['Ace2p'][0].sequence_found)
		self.assertEqual(813, dt_out_gain_change_1['Ace2p'][0].position_start)
		self.assertEqual('Cph2p', dt_out_gain_change_1['Cph2p'][0].name)
		self.assertEqual('ATCA', dt_out_gain_change_1['Cph2p'][0].sequence_found)
		self.assertEqual(1661, dt_out_gain_change_1['Cph2p'][0].position_start)

		## main process
		out_path = self.utils.get_temp_dir()
		size_window = 1000
		process_files.process_gff(ref_file_name, gff_file_name, list_vcf_files, size_window, out_path)
		
		result_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/results/" + Constants.FILE_NAME_RESULTS_STATS)
		self.assertTrue(os.path.exists(result_file_name))
		self.assertTrue(self.utils.test_equal_files(result_file_name, os.path.join(out_path, Constants.FILE_NAME_RESULTS_STATS)))
		
		self.utils.remove_dir(out_path)
		self.utils.remove_file(out_file)
		self.utils.remove_file(sample_list)
		self.utils.remove_file(out_file_ref)
		self.utils.remove_file(vcf_temp_file)
		self.utils.remove_file(vcf_temp_file + ".tbi")

		
	def test_get_diff_between_two_samples(self):
		
		seq_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/yeastract.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		ref_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/ref/ref.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		vcf_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_1_S31.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name))
		vcf_file_name_1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_3_S32.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name_1))
		gff_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/gff/file.gff3")
		self.assertTrue(os.path.exists(gff_file_name))
		sample_list = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/sample_list.txt")
		self.assertTrue(os.path.exists(sample_list))
		
		read_factors = TransFactors(seq_file_name)
		process_files = ProcessGFF(read_factors, sample_list)
		process_files.set_not_forking()
		
		list_vcf_files = [vcf_file_name, vcf_file_name_1]
		## main process
		out_path = self.utils.get_temp_dir()
		size_window = 1000
		process_files.process_gff(ref_file_name, gff_file_name, list_vcf_files, size_window, out_path)
		
		result_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/results_2/" + Constants.FILE_NAME_RESULTS_STATS)
		self.assertTrue(os.path.exists(result_file_name))
		self.assertTrue(self.utils.test_equal_files(result_file_name, os.path.join(out_path, Constants.FILE_NAME_RESULTS_STATS)))
		
		### get number of files inside Mrr1p_1 trans factor
		temp_file = self.utils.get_temp_file_with_path(out_path, "file_name", ".txt")
		cmd = "ls {}/*.tsv | wc -l > {}".format(os.path.join(out_path, "Mrr1p_1"), temp_file)
		os.system(cmd)
		vect_lines = self.utils.read_text_file(temp_file)
		self.assertEqual(1, len(vect_lines))
		self.assertEqual(2, int(vect_lines[0]))
		
		### get number of files inside Mrr1p_1 trans factor
		cmd = "ls {}/*.tsv | wc -l > {}".format(os.path.join(out_path, "Wor1p_1"), temp_file)
		os.system(cmd)
		vect_lines = self.utils.read_text_file(temp_file)
		self.assertEqual(1, len(vect_lines))
		self.assertEqual(0, int(vect_lines[0]))
		
		### get number of files inside Mrr1p_1 trans factor
		cmd = "ls {}/*.tsv | wc -l > {}".format(
			os.path.join(out_path, "Wor1p"),
			temp_file)
		os.system(cmd)
		vect_lines = self.utils.read_text_file(temp_file)
		self.assertEqual(1, len(vect_lines))
		self.assertEqual(2, int(vect_lines[0]))
		
		result_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/results_2/Wor2p/" + Constants.FILE_NAME_OUT_GAIN)
		self.assertTrue(os.path.exists(result_file_name))
		self.assertTrue(self.utils.test_equal_files(result_file_name, os.path.join(out_path, "Wor2p", Constants.FILE_NAME_OUT_GAIN)))
		result_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/results_2/Wor2p/" + Constants.FILE_NAME_OUT_EXTENDED_GAIN)
		self.assertTrue(os.path.exists(result_file_name))
		self.assertTrue(self.utils.test_equal_files(result_file_name, os.path.join(out_path, "Wor2p", Constants.FILE_NAME_OUT_EXTENDED_GAIN)))
		
		result_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/results_2/Mrr1p_1/" + Constants.FILE_NAME_OUT_LOST)
		self.assertTrue(os.path.exists(result_file_name))
		self.assertTrue(self.utils.test_equal_files(result_file_name, os.path.join(out_path, "Mrr1p_1", Constants.FILE_NAME_OUT_LOST)))
		result_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/results_2/Mrr1p_1/" + Constants.FILE_NAME_OUT_EXTENDED_LOST)
		self.assertTrue(os.path.exists(result_file_name))
		self.assertTrue(self.utils.test_equal_files(result_file_name, os.path.join(out_path, "Mrr1p_1", Constants.FILE_NAME_OUT_EXTENDED_LOST)))
		
		self.utils.remove_dir(out_path)
		
	def test_get_diff_between_two_seq_bit(self): ### get matches
		
		seq_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/yeastract.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		ref_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/ref/ref.fasta")
		self.assertTrue(os.path.exists(seq_file_name))
		vcf_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_1_S31_snps_1.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name))
		vcf_file_name_1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/vcf/T0_1_S31_snps.vcf.gz")
		self.assertTrue(os.path.exists(vcf_file_name_1))
		gff_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/gff/file.gff3")
		self.assertTrue(os.path.exists(gff_file_name))

		sample_list = self.utils.get_temp_file("sample_list", ".txt")
		with open(sample_list, 'w') as handle_write:
			handle_write.write(os.path.basename(vcf_file_name).replace('.vcf.gz', ''))
			
		read_factors = TransFactors(seq_file_name)
		process_files = ProcessGFF(read_factors, sample_list)
		process_files.set_not_forking()
		
		out_file = self.utils.get_temp_file("temp_file", ".fasta")
		out_file_ref = self.utils.get_temp_file("temp_file_ref", ".fasta")
		vcf_temp_file = process_files.get_vcf_header(vcf_file_name)	## clean VCF file
		
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-1680", vcf_file_name_1, out_file)
		dt_result_change = read_factors.get_match_factors_from_fasta_through_bit(out_file)
		self.assertEqual(19, len(dt_result_change))
		process_files.make_consensus(ref_file_name, "Ca22chr1A_C_albicans_SC5314:680-1680", vcf_temp_file, out_file_ref)
		dt_result_change = process_files.correct_positions(dt_result_change, out_file_ref, out_file, 680)
		self.assertEqual(19, len(dt_result_change))
		dt_result_ref = read_factors.get_match_factors_from_fasta_through_bit(out_file_ref)
		self.assertEqual(18, len(dt_result_ref))
		
		dt_out_lost_ref, dt_out_gain_change = process_files.get_diff_factors(dt_result_ref, dt_result_change)
		self.assertEqual(6, len(dt_out_gain_change))
		self.assertEqual("AACCACC", dt_out_gain_change['Ace2p'][0].sequence_found)
		self.assertEqual(813, dt_out_gain_change['Ace2p'][0].position_start)
		self.assertEqual(1, len(dt_out_gain_change['Ace2p']))
		self.assertEqual(1, len(dt_out_gain_change['Skn7p']))
		
		self.assertEqual(2, len(dt_out_lost_ref))
		self.assertEqual(1, len(dt_out_lost_ref['Mrr1p']))
		self.assertEqual(1, len(dt_out_lost_ref['Mrr1p_1']))
		self.assertEqual(1613, dt_out_lost_ref['Mrr1p'][0].position_start)
		self.assertEqual("AAAAAAAAAT", dt_out_lost_ref['Mrr1p'][0].sequence_found)
		
		self.utils.remove_file(out_file)
		self.utils.remove_file(sample_list)
		self.utils.remove_file(out_file_ref)
		self.utils.remove_file(vcf_temp_file)
		self.utils.remove_file(vcf_temp_file + ".tbi")

		
if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testConstants']
	unittest.main()

