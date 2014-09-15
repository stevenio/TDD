"""
play with unittest
Author:  Yuhang Wang
Date: Sep-09-2014
"""
#########################################
# Compability with Python-3
#########################################
from __future__ import print_function
#########################################

case = 2

import unittest
import random


if case == 1:
  def f(x):
    return x+1

  class MyTest(unittest.TestCase):
    def test(self):
      self.assertEqual(f(3),4)

elif case == 2:
  class TestSequenceFunctions(unittest.TestCase):
    def setUp(self):
	self.seq = range(10)

    def test_shuffle(self):
	    # make sure the suffled sequence does not lose any elements
	    random.shuffle(self.seq)
	    self.seq.sort()
	    self.assertEqual(self.seq, range(10))

	    # should raise an exception for an immutable sequence
	    self.assertRaises(TypeError, random.shuffle, (1,2,3))

    def test_choice(self):
      element = random.choice(self.seq)
      print("random choice: {0}".format(element))
      self.assertTrue(element in self.seq)


    def test_sample(self):
      with self.assertRaises(ValueError):
        random.sample(self.seq,20)
      for element in random.sample(self.seq, 5):
        self.assertTrue(element in self.seq)

if __name__ == '__main__':
    unittest.main()

