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
      trajectory_name = "segA3.dcd"
      psf_name = "segA3.psf"
    elif case == 3:
      trajectory_name = "segA.pdb"
      psf_name = "segA.psf"

    data_dir = os.path.realpath(data_dir)
    self.my_trajectory = os.path.join(data_dir, trajectory_name)
    self.my_psf = os.path.join(data_dir, psf_name)
    self.allos = allostery.Allostery(data_dir, self.my_psf, self.my_trajectory)
    #-------------------
    # for testing
    #------------------------
    self.resid_list = range(41,43)
    self.segid_list = ['A']
    self.extra_criteria = "not name H*"

  def tearDown(self):
    print("-----Tearing down {0} -------".format(self.__class__.__name__))
    del self.allos

  def test_instantiation(self):
    self.assertEqual(self.allos.psf_filename, self.my_psf)
    self.assertEqual(self.allos.trajectory_filename, self.my_trajectory)

  def test_set_basal_selection(self):
    self.allos.set_basal_selection(self.segid_list, self.extra_criteria)
    self.assertEqual(self.allos.segid_list, self.segid_list)
    selection_str = "(segid {0}) and ({1})".format(self.segid_list[0], self.extra_criteria)
    self.assertEqual(self.allos.basal_selection_str, selection_str)

  def test_build_com_coords(self):
    N = len(self.resid_list)
    self.allos.set_basal_selection(self.segid_list, self.extra_criteria)
    self.allos.build_com_coords(self.resid_list, self.segid_list, self.extra_criteria)
    self.assertEqual(numpy.shape(self.allos.com_coords), (N,3))
  @unittest.skip("")
  def test_build_pair_dist_matrix(self):
    resid_list = range(41,43)
    self.allos.build_com_pair_dist_matrix(resid_list, self.segid_list, self.extra_criteria)
    result = numpy.array([[ 0.        , 4.59339682],
                          [ 4.59339682,  0.        ]], dtype=numpy.float32)
    tol = 1E-5;
    self.assertTrue((numpy.abs(self.allos.com_pair_dist_matrix - result) < tol).all())
  
  def test_get_commute_time(self):
    self.allos.get_commute_time(self.resid_list, self.segid_list, self.extra_criteria)
    print("commute time\n", self.allos.commute_time)
    
    

    
if __name__ == '__main__':
  unittest.main()
