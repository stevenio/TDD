"""
Goal: provide a module of small functions for any 
      distance related calculations
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
import scipy
#import scipy.linalg
#============================================================

def centmass(coords_matrix,mass_vector):
  """
  Return the center of mass coordinates for
    a set of coordinates.
  Inputs: coords_matrix << Nx3 matrix(numpy.array object) of x,y,z coordinates
          mass_vector   << Nx1 vector(numpy.array object) of masses 
  """
  totalMass = scipy.sum(mass_vector)
  return 1.0/totalMass*scipy.dot(coords_matrix.T, mass_vector)


if __name__ == "__main__":
  print("main")
  M = scipy.arange(0,9, dtype=scipy.float64).reshape(3,3)
  b = scipy.array([1,1,1], dtype=scipy.float64)
  print("M=\n",M)
  print("b=\n",b)
  print(centmass(M,b))
