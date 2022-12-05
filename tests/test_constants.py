'''
Created on Oct 15, 2016

@author: mmp
'''
import unittest
from constants.constants import Constants 

class Test(unittest.TestCase):


	def testConstants(self):
		
		constants = Constants()
		self.assertEqual(constants.complement("AAATTTCCC"), "TTTAAAGGG", "must be equal")
		self.assertEqual(constants.complement("AAATTGGGTCCC"), "TTTAACCCAGGG", "must be equal")
		self.assertEqual(constants.reverse_complement("AAGGGATTTCCC"), "GGGAAATCCCTT", "must be equal")
		self.assertEqual(constants.reverse_complement("AAATTAGTCCC"), "GGGACTAATTT", "must be equal")
		self.assertEqual(constants.ambiguos_to_unambiguous("YARWATTTCCC"), "[CT]A[AG][AT]ATTTCCC", "must be equal")
		self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN"), "[AG][CT][GT][AC][CG][AT][CGT][AGT][ACT][ACG][ACGT]", "must be equal")
		self.assertEqual(constants.complement("RYKMSWBDHVN"), "YRMKSWVHDBN", "must be equal")

		self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{3}"), "[AG][CT][GT][AC][CG][AT][CGT][AGT][ACT][ACG][ACGT][ACGT][ACGT]", "must be equal")
		self.assertEqual(constants.ambiguos_to_unambiguous("GGAGGC(G/A)C(T/A)G"), "GGAGGC[AG]C[AT]G", "must be equal")
		
		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{3"), "[AG][CT][GT][AC][GC][AT][CGT][AGT][ACT][ACG][ACGT][ACGT][ACGT]", "must be equal")
			self.assertFail("Error, must throw exception")
		except Exception:
			pass
			
		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{3}}"), "")
			self.assertFail("Error, must throw exception")
		except Exception:
			pass
		
		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{}}"), "")
			self.assertFail("Error, must throw exception")
		except Exception:
			pass

		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{A}}"), "")
			self.assertFail("Error, must throw exception")
		except Exception:
			pass
		
		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN}"), "")
			self.assertFail("Error, must throw exception")
		except Exception:
			pass
		
		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("{RYKMSWBDHVN"), "")
			self.assertFail("Error, must throw exception")
		except Exception:
			pass
	
	def test_bit_dict(self):
		
		constants = Constants()
		self.assertEqual({0 : 0b0001}, constants.ambiguos_to_unambiguous_dict_bit("A"))
		self.assertEqual({0 : 0b0011}, constants.ambiguos_to_unambiguous_dict_bit("(A/C)"))
		self.assertEqual({0 : 0b0011, 1 : 0b1111}, constants.ambiguos_to_unambiguous_dict_bit("(A/C)N"))
		self.assertEqual({0 : 0b0011, 1 : 0b1111, 2 : 0b10000 }, constants.ambiguos_to_unambiguous_dict_bit("(A/C)N*"))
		self.assertEqual({0 : 0b0011, 1 : 0b1111, 2 : 0b10000, 3 : 0b1001 }, constants.ambiguos_to_unambiguous_dict_bit("(A/C)N*W"))
		
	def testGet_diff_between_two_seq(self):
		
		constants = Constants()
		self.assertEqual(constants.get_diff_between_two_seq("AAATTTCCC", "TTTAAAGGG"), 9)
		self.assertEqual(constants.get_diff_between_two_seq("AAATTTCCC", "TTTAAAGGGS"), 0)
		self.assertEqual(constants.get_diff_between_two_seq("AAATTTCCC", "AAATTTCCC"), 0)
		self.assertEqual(constants.get_diff_between_two_seq("GAATTTCCC", "AAATTTCCC"), 1)
		self.assertEqual(constants.get_diff_between_two_seq("AAATTTCCC", "AAATTTCCG"), 1)


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testConstants']
	unittest.main()