"""
Goal: test "allostery" module
Author: Yuhang Wang
Date: 2014-09-16
"""
#============================================================
# Compatibility with Python 3
#============================================================
from __future__ import print_function, division
from nose.tools import with_setup
import allostery
import nose
import numpy
import scipy
import os
import unittest
import pandas
#============================================================

#============================================================
# Dependencies
#============================================================
#--------------------------------------------------
# testing target  module
#--------------------------------------------------
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
      trajectory_name = "sod.dcd"
      psf_name = "sod.psf"
      self.resid_list = (121, 122)
      self.segid_list = ['O1']
      self.extra_criteria = "not name H*"
      self.answer_key_filename = "sod_comdist.dat"
    elif case == 2:
      trajectory_name = "segA3.pdb"
      psf_name = "segA3.psf"
      self.resid_list = (41, 42, 43)
      self.segid_list = ['A']
      self.extra_criteria = "not name H*"
      self.answer_key_filename = "segA3_comdist.dat"
    elif case == 3:
      trajectory_name = "segA.pdb"
      psf_name = "segA.psf"

    self.data_dir = os.path.realpath(data_dir)
    self.my_trajectory = os.path.join(self.data_dir, trajectory_name)
    self.my_psf = os.path.join(self.data_dir, psf_name)
    self.allos = allostery.Allostery(self.data_dir, self.my_psf, self.my_trajectory)
    self.tol = 1E-5;
    #-------------------
    # for testing
    #------------------------


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
      
  def test_build_pair_dist_matrix(self):
    self.allos.build_com_pair_dist_matrix(self.resid_list, self.segid_list, self.extra_criteria)
    filename = os.path.join(self.data_dir, self.answer_key_filename)
    key = numpy.loadtxt(filename)
    self.assertTrue((numpy.abs(self.allos.com_pair_dist_matrix - key) < self.tol).all())
    
  def test_get_commute_time(self):
    self.allos.get_commute_time(self.resid_list, self.segid_list, self.extra_criteria)
    print("commute time\n", self.allos.commute_time)
    
  def test_save_commute_time(self):
    self.test_get_commute_time()
    filename = "comm_time"
    cwd = os.getcwd()
    myfile = os.path.join(cwd, filename+".h5")

    if os.path.isfile(myfile): # remove previous version
      os.remove(myfile)

    self.allos.save_commute_time(filename)
    self.assertTrue(os.path.isfile(myfile))
    data = pandas.read_hdf(myfile, 'data')
    print("file content: {0}\n".format(myfile),data.values)
    

    
if __name__ == '__main__':
  unittest.main()
