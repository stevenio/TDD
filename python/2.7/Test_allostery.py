"""
Goal: test "allostery" module
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
import os
import unittest
import nose
from nose.tools import with_setup
import numpy, scipy
#--------------------------------------------------
# testing target  module
#--------------------------------------------------
import allostery
#============================================================

class TestAllostery(unittest.TestCase):
  def setUp(self):
    print("-----Setting up {0} -------".format(self.__class__.__name__))
    #------------------------------------------------------------
    # Parameters
    #------------------------------------------------------------
    data_dir = "data"
    case = 2
    if case == 1:
      trajectory_name = "sod.pdb"
      psf_name = "sod.psf"
    elif case == 2:
      trajectory_name = "segA3.pdb"
      psf_name = "segA3.psf"
    elif case == 3:
      trajectory_name = "segA.pdb"
      psf_name = "segA.psf"

    data_dir = os.path.realpath(data_dir)
    self.my_trajectory = os.path.join(data_dir, trajectory_name)
    self.my_psf = os.path.join(data_dir, psf_name)
    self.allos = allostery.Allostery(data_dir, self.my_psf, self.my_trajectory)

  def tearDown(self):
    print("-----Tearing down {0} -------".format(self.__class__.__name__))
    del self.allos

  def test_instantiation(self):
    assert self.allos.psf_filename == self.my_psf
    assert self.allos.trajectory_filename == self.my_trajectory

  def test_set_resid_list(self):
    resid_list = range(1,10)
    self.allos.set_resid_list(resid_list)
    self.assertEqual(self.allos.resid_list, resid_list)

  def test_set_pair_dist_matrix(self):
    resid_list = range(1,10)
    N = len(resid_list)
    self.allos.set_resid_list(resid_list)
    self.allos.set_pair_dist_matrix()
    self.assertEqual(numpy.shape(self.allos.pair_dist_matrix), (N,3))


if __name__ == '__main__':
  unittest.main()
