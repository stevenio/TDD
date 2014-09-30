"""
Goal: compute residue pairwise center of mass distances
Author: Yuhang Wang
Date: 09-17-2014 
"""
#============================================================
# Compatibility with python 3
#============================================================
from __future__ import print_function, division
#============================================================

#============================================================
# External Module Dependencies
#============================================================
import prody
import pandas
import MDAnalysis as mda
import MDAnalysis.core.parallel.distances as mda_dist
import numpy, scipy
import os
#--------------------
# my modules
#--------------------
import distmodule
#============================================================


#============================================================================
#               [[[ Interface: read/manipulate user data]]]
#============================================================================
class UserData(object):
  """
  Role: an universal interface between user input data and analysis tools
  Friend: :class: 'CreateUserDataObj'
  Job:  1). create an instance of user data
        2). make atom selections
  """
  def __init__(self, psf_filename, trajectory_filename):
    """
    Role: use delegation to create an user data object
    """
    self.userDataObj = Proxy_MDAnalysis(psf_filename,trajectory_filename)

  def numResiduesSelected(self):
    """
    Role count the number of residues selected
    :return: scalar (number of residues in self.selectedAtoms
    """
    return self.userDataObj.numberResiduesSelected()

  def select(self, selection_string):
    """
    Role: use delegation to select part of the MD system
    :return None (left the delegated class to store the selected atoms
            in whatever data type it likes)
    """
    self.userDataObj.select(selection_string)


  def build_com_matrix(self):
    """
    Role: build a matrix (Nx3) of center of mass coordinates
    :return: numpy matrix (shape: Nx3)
    """
    self.comMatrix = self.userDataObj.selected_centerOfMasses()
    return self.comMatrix


#---------------------------------------------------------------------------
#                 [[[ Proxy for MDAnalysis ]]]
#---------------------------------------------------------------------------
class Proxy_MDAnalysis(object):
  """
    Title: A particular type of Python module for creating user data object
    Purpose: act as a proxy for using "MDAnalysis" module
    Friend: MDAnalysis module (so Proxy_MDAnslysis is expected to
            know everthing about MDAnalysis module)

    API: select(self, selection_string): select a subset of atoms
    Rule: all proxy must have the following methods:
         1). select(self, selection_string)

    note: Proxy_MDAnalysis should be the only class that interacts with MDAnalysis module
    """
  def __init__(self, psf_filename, trajectory_filename):
    """
    Goal: create an universe object based on user's files
    :rtype : MDAnalysis universe object
    """
    self.userData = mda.Universe(psf_filename, trajectory_filename)

  def select(self, selection_string):
    """
    Role: select part of the MD system
    :return a data object that should know the coordinates of selected atoms
    """
    self.selectedAtoms = self.userData.selectAtoms(selection_string)
    if self.selectedAtoms.numberOfAtoms() == 0:
      msg = "atom selection [{0}] yields zero atoms".format(selection_string)
      raise UserWarning(msg)
    else:
      return self.selectedAtoms

  def numberResiduesSelected(self):
    """
    Role: count the number of residues in the current atom selection
    :return scalar (number of residues)
    """
    return self.selectedAtoms.numberOfResidues()

  def selected_centerOfMasses(self,pbc=True):
    """
    Role: calculate the center of mass of residueObj
    :return numpy array [shape: (3,)]
    """
    self.centerOfMasses = numpy.zeros((self.numberResiduesSelected(), 3))
    for _ccc in range(0,self.numberResiduesSelected()):
      _residue = self.selectedAtoms.residues[_ccc]
      self.centerOfMasses[_ccc,:] = _residue.centerOfMass()
    return self.centerOfMasses





#============================================================================
#               [[[Interface: make atom selections]]]
#============================================================================
class Selection(object):
  """
  Role: an interface for make atom selections, which delegate other
        classes to do the job, in order to match the style
        of different user data object creation modules (like :class:'MDAnalysis')
  """
  @classmethod
  def create(cls, generic_keyword, selection_member_list, extra_criteria=None):
    """
    Role: create an atom selection string
    """
    selection_keyword = MDAnalysis_Style_Selection.keyword_style_convert(generic_keyword)
    return MDAnalysis_Style_Selection.create(selection_keyword,
                                             selection_member_list,
                                             extra_criteria)


class GenericKeywords(object):
  """
  Role: provide a list of generic symbol for atom selection keywords
  """
  #
  ResId = "_resid"
  SegId = "_segid"


class MDAnalysis_Style_Selection(object):
  """
  Title: A particular type of atom selection style
  Role: create a selection string that match the atom selection style of [[MDAnalysis]] python module
  """
  keyword_conversion_dict = dict({GenericKeywords.ResId:"resid",
                                  GenericKeywords.SegId:"segid"})
  @classmethod
  def keyword_style_convert(cls, keyword):
    """
    Role: covert keyword into a style that conforms to the
          style of [[MDAnalysis]] module
    :param keyword:
    :return:
    """
    return MDAnalysis_Style_Selection.keyword_conversion_dict[keyword]

  @classmethod
  def _make_selection_str(cls, keyword, member_list):
    """
    return of string like this "(segid 1 or segid 2)",
    where keyword is "segid" and member_list is [1,2]
    """
    S = []
    for member in member_list:
      S.append("{0} {1}".format(keyword, member))
    return "({0})".format(" or ".join(S))


  @classmethod
  def create(cls, keyword, member_list, extra_criteria=None):
    """
    Role: make an atom selection string according to the [[MDAnalysis]] module style
    """
    if extra_criteria == None:
      S = [] # list to store intermediate atom selections
    else:
      S = ["({0})".format(extra_criteria)] # list to store intermediate atom selections

    S.append(MDAnalysis_Style_Selection._make_selection_str(keyword, member_list))
    return "{0}".format(" and ".join(S))



#============================================================================
#               [[[ Module: tools for data analysis ]]]
#============================================================================
class AnalysisTools(object):
  """
  Role: provide an interface to all analysis tools
  """
  def set_selection_str(self, keyword, member_list):
    """
    return of string like this "(segid 1 or segid 2)",
    where keyword is "segid" and member_list is [1,2]
    """
    S = []
    for member in member_list:
      S.append("{0} {1}".format(keyword, member))
    return "({0})".format(" or ".join(S))

  def set_basal_selection(self, segid_list, extra_criteria="all"):
    """
    Set up base selection string using segid's, [and extra_criteria],
    to be used for residue pairwise distance
    calculation.
    """

    self.segid_list = segid_list
    self.extra_criteria = extra_criteria
    S  = []
    S.append(self.set_selection_str("segid", segid_list))
    self.basal_selection_str = " and ".join(S)
    self.basal_selection_str += " and ({0})".format(extra_criteria)



  def build_com_coords(self, resid_list, segid_list, extra_criteria="all"):
    """
    build an Nx3 matrix of center of mass coordinates for selected residues in resid_list
    """
    self.set_basal_selection(segid_list, extra_criteria)
    self.resid_list = resid_list
    N = len(self.resid_list)
    self.com_coords = numpy.zeros((N,3), dtype=numpy.float32)
    ccc = 0 # counter
    for resid in self.resid_list:
      selection_str = self.basal_selection_str + " and resnum {0}".format(resid)
      tmp_selection = self.universe.selectAtoms(selection_str)
      #------------------------------------------------
      # in case zero atoms got selected
      #------------------------------------------------
      if tmp_selection.numberOfAtoms() == 0:
        msg = "atom selection [{0}] yields zero atoms".format(selection_str)
        raise UserWarning(msg)
      #--------------------------------------------------
      X = tmp_selection.get_positions()
      M = tmp_selection.masses()
      self.com_coords[ccc,:] = tmp_selection.centerOfMass(pbc=True)
      ccc += 1

  def build_com_pair_dist_matrix(self, resid_list, segid_list, extra_criteria="all"):
    """
    build pair-wise distance matrix for center of mass coordinates
    """
    self.build_com_coords(resid_list, segid_list, extra_criteria)
    #-----------------
    # CALL MDAnalysis.core.distances.distance_array
    # instead of the MDAnalysis.core.parallel.distances (which requires inputs to be of Cython DTYPE_t type)
    self.com_pair_dist_matrix = mda.core.distances.distance_array(self.com_coords, self.com_coords)

  def get_commute_time(self, resid_list, segid_list, extra_criteria="all"):
    """
    build commute time matrix based on pair-wise com distance matrix
    """

    ccc = 0 # frame counter
    for fr in self.universe.trajectory:
      frameId = fr.frame
      ccc += 1
      print("frame: {0}".format(fr.frame))
      #---------------------------------------------------------
      # build self.com_pair_dist_matrix for the current matrix
      #----------------------------------------------
      self.build_com_pair_dist_matrix(resid_list, segid_list, extra_criteria)
      #------------------
      #  initialize
      #-----------------------
      if frameId == 1:
        self.commute_time = numpy.zeros(numpy.shape(self.com_pair_dist_matrix), dtype=numpy.float32)
      #-----------------------------------
      self.commute_time += self.com_pair_dist_matrix
    #----------------------------
    # compute the average
    #-----------------------------------
    self.commute_time /= ccc

  def save_commute_time(self, filename, output_dir="./"):
    """
    Save commute time into [filename].h5
    """
    output_data = pandas.DataFrame(self.commute_time, columns=self.resid_list)
    filename += ".h5"
    current_dir = os.path.realpath(output_dir) # get current path
    output_filename = os.path.join(current_dir, filename)
    output_data.to_hdf(output_filename,'data')







class Allostery(object):
  """
  Goal: get the allosteric relations between residues
  Role: This "Allostery" class is a front-end interface
        for reading data and doing analyses
  """

  def __init__(self, data_src_dir, psf_filename, trajectory_filename):
    """
    Role: create an universe object to store user data
    :param data_src_dir: directory where user data live
    :param psf_filename: file name for the *.psf
    :param trajectory_filename: file name for the trajectory files, e.g. *.dcd, *.xtc
    :return: None
    """
    self.psf_filename = os.path.join(data_src_dir, psf_filename)
    self.trajectory_filename = os.path.join(data_src_dir, trajectory_filename)
    # use delegation to read user data
    self.userData = UserData(self.psf_filename, self.trajectory_filename)


  def select(self, resid_list=None, segid_list=None, extra_criteria=None):
    """
    Role: create an user defined atom selection object
    :param resid_list: a list of residue IDs, e.g. [1, 2, 3]
    :param segid_list: a list of segment IDs, e.g. ['A', 'B', 'C']
    :param extra_criteria: some extra restraining criteria, e.g. "not name H*"
    :return: None (let the delegated class to store the selected atoms)
    """
    # add extra criteria
    # note: S is a list which stores intermediate atom selection strings
    if extra_criteria == None:
      S = []
    else:
      S = ["({0})".format(extra_criteria)]

    # select residues
    if resid_list != None:
      S.append(Selection.create(GenericKeywords.ResId, resid_list))
      self.ResIdList = resid_list
    else:
      self.ResIdList = None

    # select segments
    if segid_list != None:
      S.append(Selection.create(GenericKeywords.SegId, segid_list))
      self.SegIdList = segid_list
    else:
      self.SegIdList = None

    self.selection_string = " and ".join(S)
    self.userData.select(self.selection_string)

  def build_com_matrix(self):
    """
    Role: build the center of mass matrix for selected atoms
    :return: numpy matrix (shape: Nx3)
    """
    self.comMatrix = self.userData.build_com_matrix()
    print("comMatrx\n",self.comMatrix)




    



