'''
Created on 25/11/2022

@author: mmp
'''
import unittest
import os
from utils.util import Genes, TransFactor, TransFactorSmall

class Test(unittest.TestCase):


	def setUp(self):
		pass


	def tearDown(self):
		pass


	def test_genes(self):
		
		gff_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files/gff/file.gff3")
		self.assertTrue(os.path.exists(gff_file_name))

		genes = Genes(gff_file_name)

		self.assertEqual(1, len(genes.dt_chromossomes))
		self.assertEqual(19, len(genes.get_genes('Ca22chr1A_C_albicans_SC5314')))
		self.assertEqual(10, len(genes.get_genes('Ca22chr1A_C_albicans_SC5314', 1000)))
		self.assertEqual(8, len(genes.get_genes('Ca22chr1A_C_albicans_SC5314', 2000)))
		#genes.print_all_genes_by_chr('Ca22chr1A_C_albicans_SC5314')
	
	def test_match_bit(self):
		
		trans_factor = TransFactor("xpro1", "TTGG")
		
		post_start, sequence_fasta = 1, "AACCTTGG"
		vect_out = trans_factor.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(1, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro1", "TTGG", 5), vect_out[0])
		
		post_start, sequence_fasta = 1, "AACCBTGV"
		vect_out = trans_factor.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(1, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro1", "BTGV", 5), vect_out[0]) 

		trans_factor_1 = TransFactor("xpro3", "AACC")
		vect_out = trans_factor_1.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(1, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro3", "AACC", 1), vect_out[0])
		
		trans_factor_1 = TransFactor("xpro3", "AACT")
		vect_out = trans_factor_1.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(0, len(vect_out))


		post_start, sequence_fasta = 2, "AACCBTGVAAVV"
		trans_factor_1 = TransFactor("xpro3", "AACC")
		vect_out = trans_factor_1.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(2, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro3", "AACC", 2), vect_out[0])
		self.assertEqual(TransFactorSmall("xpro3", "AAVV", 10), vect_out[1])
		
		
		post_start, sequence_fasta = 2, "AACCBTGVAAVN"
		trans_factor_1 = TransFactor("xpro3", "AACC")
		vect_out = trans_factor_1.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(2, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro3", "AACC", 2), vect_out[0])
		self.assertEqual(TransFactorSmall("xpro3", "AAVN", 10), vect_out[1])
		
		
		post_start, sequence_fasta = 2, "AAC*CBTGVAAVV"
		trans_factor_1 = TransFactor("xpro3", "AACC")
		vect_out = trans_factor_1.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(1, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro3", "AAVV", 11), vect_out[0])
		
		post_start, sequence_fasta = 0, "CCTGCTAATTTAACAAAAACTAATAATTAAA"
		trans_factor_1 = TransFactor("xpro3", "TGCTAATTTAACAAAAACTAATAATTA")
		vect_out = trans_factor_1.get_match_bit(post_start, sequence_fasta)
		self.assertEqual(1, len(vect_out))
		self.assertEqual(TransFactorSmall("xpro3", "TGCTAATTTAACAAAAACTAATAATTA", 2), vect_out[0])

if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.test_genes']
	unittest.main()