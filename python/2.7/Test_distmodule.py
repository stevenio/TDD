"""
Goal: unit testing for distmodule.py
Author: Yuhang Wang
Date: 2014-09-16
"""

#============================================================
# Compatibility with Python 3
#============================================================
from __future__ import print_function, division
#============================================================

#============================================================
# Dependencies
#============================================================
import nose
import scipy
import distmodule
#============================================================

def test_centmass():
  M = scipy.arange(0,9, dtype=scipy.float64).reshape(3,3)
  b = scipy.array([1,1,1], dtype=scipy.float64)
  result = distmodule.centmass(M,b)
  answer = scipy.array([3., 4., 5.], dtype=scipy.float64)
  print("M=\n",M)
  print("b=\n",b)
  print("center of mass:\n",result)
  assert (result == answer).all()
