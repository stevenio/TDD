"""
play with unittest
Author:  Yuhang Wang
Date: Sep-09-2014
"""

import unittest

def f(x):
  return x+1

class MyTest(unittest.TestCase):
  def test(self):
    self.assertEqual(f(3),4)



if __name__ == '__main__':
    unittest.main()
