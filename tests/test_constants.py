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
		self.assertEqual(constants.ambiguos_to_unambiguous("YARWATTTCCC"), "[TC]A[AG][AT]ATTTCCC", "must be equal")
		self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN"), "[AG][TC][GT][AC][GC][AT][CGT][AGT][ACT][ACG][ACGT]", "must be equal")
		self.assertEqual(constants.complement("RYKMSWBDHVN"), "YRMKSWVHDBN", "must be equal")

		self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{3}"), "[AG][TC][GT][AC][GC][AT][CGT][AGT][ACT][ACG][ACGT][ACGT][ACGT]", "must be equal")
		self.assertEqual(constants.ambiguos_to_unambiguous("GGAGGC(G/A)C(T/A)G"), "GGAGGC[GA]C[TA]G", "must be equal")
		
		try:
			self.assertEqual(constants.ambiguos_to_unambiguous("RYKMSWBDHVN{3"), "[AG][TC][GT][AC][GC][AT][CGT][AGT][ACT][ACG][ACGT][ACGT][ACGT]", "must be equal")
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