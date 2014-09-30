"""
Goal: test "allostery" module
Author: Yuhang Wang
Date: 2014-09-16
"""
#============================================================
# Compatibility with Python 3D
#============================================================
from __future__ import print_function, division
import allostery
import numpy
# import scipy
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
      self.comdist_data_filename = "sod_comdist.dat"


    elif case == 2:
      trajectory_name = "segA3.pdb"
      psf_name = "segA3.psf"
      self.resid_list = (41, 42, 43)
      self.segid_list = ['A']
      self.extra_criteria = "not name H*"
      self.expected_selection_str = \
        "({0}) and (resid {1} or resid {2} or resid {3}) and (segid {4})".format(self.extra_criteria,
                                                                                 self.resid_list[0],
                                                                                 self.resid_list[1],
                                                                                 self.resid_list[2],
                                                                                 self.segid_list[0])
      self.comdist_data_filename = "segA3_comdist.dat"
      self.com_data_filename = "segA3_com.dat"
      self.expected_numAtomsSelected = 18
      self.coord_filename = "segA3_coord.dat"
      self.mass_filename = "segA3_mass.dat"

    elif case == 3:
      trajectory_name = "segA.pdb"
      psf_name = "segA.psf"

    self.data_dir = os.path.realpath(data_dir)
    self.my_trajectory = os.path.join(self.data_dir, trajectory_name)
    self.my_psf = os.path.join(self.data_dir, psf_name)
    self.allos = allostery.Allostery(self.data_dir, self.my_psf, self.my_trajectory)
    self.tol = 1E-2;
    #-------------------
    # for testing
    #------------------------


  def tearDown(self):
    print("-----Tearing down {0} -------".format(self.__class__.__name__))
    del self.allos

  def test_instantiation(self):
    print("\t[[[ Testing Allostery() object instantiation ]]]")
    self.assertEqual(self.allos.psf_filename, self.my_psf)
    self.assertEqual(self.allos.trajectory_filename, self.my_trajectory)


  def test_atom_selection_string(self):
    print("\t[[[ Testing atom_selection_string() ]]]")
    self.allos.select(self.resid_list, self.segid_list, self.extra_criteria)
    self.assertEqual(self.allos.selection_string, self.expected_selection_str)
    print("Atom selection: \t",self.allos.selection_string)
    print("expected string:\t",self.expected_selection_str)

  def test_number_selected_atoms(self):
    self.test_atom_selection_string()
    print("\t[[[ Testing get_number_selected_atoms() ]]]")
    self.assertEqual(self.allos.get_number_selected_atoms(), self.expected_numAtomsSelected)
    print("number of atoms selected: ", self.allos.get_number_selected_atoms())


  def test_raise_zero_atom_selection_error(self):
    print("\t[[[ Testing raise_zero_atom_selection_error() ]]]")
    nonexistent_resid = [-1]
    nonexistent_segid = [-1]
    self.assertRaises(UserWarning, self.allos.select, nonexistent_resid, nonexistent_segid)

  def test_reading_atom_coords(self):
    self.test_atom_selection_string()
    print("\t[[[ Testing reading atom coordinates ]]]")
    coords = self.allos.get_selected_coords()
    filename = os.path.join(self.data_dir, self.coord_filename)
    expected = numpy.loadtxt(filename)
    self.assertTrue((numpy.abs(coords-expected)<self.tol).all())

  def test_get_selected_mass(self):
    self.test_atom_selection_string()
    print("\t[[[ Testing get mass of selected atoms ]]]")
    masses = self.allos.get_selected_mass()
    filename = os.path.join(self.data_dir, self.mass_filename)
    expected = numpy.loadtxt(filename)
    self.assertTrue(((masses - expected)<self.tol).all())


  def test_build_com_matrix(self):
    print("\t[[[ Testing build_com_matrix() ]]]")
    N = len(self.resid_list)
    self.test_atom_selection_string()
    self.allos._build_com_matrix()
    filename = os.path.join(self.data_dir, self.com_data_filename)
    expected = numpy.loadtxt(filename)
    print("comMatrix\n",self.allos.comMatrix)
    print("expected comMatrix\n",expected)
    self.assertEqual(numpy.shape(self.allos.comMatrix), (N,3))
    self.assertTrue(numpy.abs((self.allos.comMatrix - expected) < self.tol).all())

  @unittest.skip("paircom")
  def test_build_pair_dist_com_matrix(self):
    self.test_atom_selection_string()
    print("\t[[[ Testing build_pairwise_dist_com_matrix() ]]]")
    self.allos._build_pairwise_distance_com_matrix()
    filename = os.path.join(self.data_dir, self.comdist_data_filename)
    key = numpy.loadtxt(filename)
    print("Expected pairComMatrix\n", key)
    self.assertTrue((numpy.abs(self.allos.pairComMatrix - key) < self.tol).all())

  @unittest.skip("tmp")
  def test_get_commute_time_matrix(self):
    self.test_atom_selection_string()
    print("\t[[[ Testing get_commute_time_matrix() ]]]")
    self.allos.get_commute_time_matrix()
    print("commute time\n", self.allos.commute_time_matrix)

  @unittest.skip("tmp")
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
