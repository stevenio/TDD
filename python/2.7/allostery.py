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
# import prody
import pandas
import MDAnalysis
import MDAnalysis.core.parallel.distances as distance_module
import numpy
import os
#--------------------
# my modules
#--------------------
# import distmodule
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



  def get_number_of_selected_residues(self):
    """
    Role: count the number of residues selected
    :return: scalar (number of residues in self.selectedAtoms
    """
    return self.userDataObj.get_number_of_selected_residues()

  def get_numAtomsSelected(self):
    """
    Role: count the number of atoms selected
    :return: Int scalar
    """
    return self.userDataObj.get_numberAtomsSelected()

  def get_selected_mass(self):
    """
    Role: get the mass of selected atoms
    :return: numpy array
    """
    return  self.userDataObj.get_selected_mass()

  def select(self, selection_string):
    """
    Role: use delegation to select part of the MD system
    :return None (left the delegated class to store the selected atoms
            in whatever data type it likes)
    """
    self.userDataObj.select(selection_string)


  def get_selected_com_coords(self):
    """
    Role: build a matrix (Nx3) of center of mass coordinates
    :return: numpy matrix (shape: Nx3)
    """
    self.comMatrix = self.userDataObj.selected_centerOfMasses()
    return self.comMatrix

  def next_frame(self):
    """
    Role: move on to the next traectory frame
    :return: None
    """
    self.userDataObj.next_frame()

  def get_current_frameId(self):
    return self.userDataObj.get_current_frameId()

  def get_total_numFrames(self):
    return self.userDataObj.get_total_numFrames()

  def get_selected_coords(self):
    """
    Role: get coordinates of selected atoms
    :return: numpy array
    """
    return self.userDataObj.get_selected_coords()

  def iter_get_selected_coords(self, frame_start, frame_end, frame_skip=1):
    """
    Role: a generator function that returns the coordinates of
          user-selected atoms iteratively for all trajectory frames
    :return: numpy array (shape: Nx3)
    """
    return self.userDataObj.iter_get_selected_coords(frame_start, frame_end, frame_skip)

  def iter_get_selected_com_coords(self, frame_start=0, frame_end=None, frame_skip=1):
    """
    Role: a generator function that returns the center of mass coordinates
          of selected atoms
    :return: numpy array (shape: Nx3, N is the number of selected residues)
    """
    if frame_end == None:
      frame_end = self.get_total_numFrames()
    for frameId in range(frame_start, frame_end, frame_skip):
      self.userDataObj.go_to_frame(frameId) # update frame
      yield self.get_selected_com_coords()




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
    self.userData = MDAnalysis.Universe(psf_filename, trajectory_filename)
    self.current_frameId = 0
    self.total_numFrames = self.userData.trajectory.numframes

  def select(self, selection_string):
    """
    Role: select part of the MD system
    :return a data object that should know the coordinates of selected atoms
    """
    self.selection_string = selection_string
    self.selectedAtoms = self.userData.selectAtoms(selection_string)
    if self.selectedAtoms.numberOfAtoms() == 0:
      msg = "atom selection [{0}] yields zero atoms".format(selection_string)
      raise UserWarning(msg)
    else:
      return self.selectedAtoms

  def get_selected_mass(self):
    """
    Role: get the mass of selected atoms
    :return: numpy array
    """
    return self.selectedAtoms.masses()

  def get_number_of_selected_residues(self):
    """
    Role: count the number of residues in the current atom selection
    :return scalar (number of residues)
    """
    return self.selectedAtoms.numberOfResidues()

  def get_numberAtomsSelected(self):
    """
    Role: count the number of atoms selected
    :return: Int scalar
    """
    return self.selectedAtoms.numberOfAtoms()

  def get_selected_coords(self):
    """
    Role: get the coordinates of selected atoms
    :return: numpy array
    """
    return self.selectedAtoms.coordinates()

  def selected_centerOfMasses(self,numericType=numpy.float32):
    """
    Role: calculate the center of mass of residueObj
    :return numpy array [shape: (3,)]
    """

    self.centerOfMasses = numpy.zeros((self.get_number_of_selected_residues(), 3), dtype=numericType)
    for _ccc in range(0,self.get_number_of_selected_residues()):
      # note: due to a bug in the implementation of MDAnalysis,
      #       self.selectedAtoms.residues[0] will return all the atoms in residue0
      #       of the original PDB, not from the selected atom set self.selectedAtoms
      _residue = self.selectedAtoms.residues[_ccc].selectAtoms(self.selection_string)
      self.centerOfMasses[_ccc,:] = _residue.centerOfMass()
    return self.centerOfMasses

  def next_frame(self):
    """
    Role: move on to the next trajectory frame
    :return None
    """
    if self.current_frameId + 1 < self.total_numFrames:
      self.current_frameId += 1
      self.userData.trajectory[self.current_frameId]
    else:
      msg = "Already reached the last frame!"
      msg += "(current frameId = {0}/{1})".format(self.current_frameId, self.total_numFrames)
      raise UserWarning(msg)

  def go_to_frame(self, frameId):
    """
    Role: go to a specific trajectory frame
    :param frameId: Int scalar
    :return: None
    """
    if frameId < 0 or frameId >= self.total_numFrames:
      msg = "frameId must be within [0,{0})".format(self.total_numFrames)
      msg += ", but you specified frameId = {0}".format(frameId)
      raise UserWarning(msg)
    self.userData.trajectory[frameId]

  def iter_get_selected_coords(self, frame_start, frame_end, frame_skip=1):
    """
    Role: a generator function that returns the coordinates
          of selected atoms for all trajectory frames
    :param frame_start: starting frame (0-based index)
           frame_end: last frame
           frame_skip: interval between frames
    :return: numpy array
    """
    for frameId in xrange(frame_start, frame_end, frame_skip):
      self.userData.trajectory[frameId]
      yield self.selectedAtoms.coordinates()



  def get_total_numFrames(self):
    """
    Role: get total number of frames in the trajectory
    :return: Int scalar
    """
    return self.total_numFrames

  def get_current_frameId(self):
    """
    Role: get the current frame ID
    :return Int scalar
    """
    return self.current_frameId






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
#               [[[ Interface: tools for data analysis ]]]
#============================================================================
class AnalysisTools(object):
  """
  Role: provide an interface to all analysis tools
  """
  @classmethod
  def pairwise_distances(cls, com_matrix):
    """
    build pair-wise distance matrix for center of mass coordinates
    """
    #-----------------
    # CALL MDAnalysis.core.distances.distance_array
    # instead of the MDAnalysis.core.parallel.distances (which requires inputs to be of Cython DTYPE_t type)
    return distance_module.distance_array(com_matrix, com_matrix)

  @classmethod
  def commute_time(cls,com_coords_iterator, number_of_rows):
    """
    Role: compute the commute time matrix
    :return: numpy matrix
    """
    commute_time_matrix = numpy.zeros((number_of_rows, number_of_rows), dtype=numpy.float32)
    numberFrames = 0
    for com_coords in com_coords_iterator:
      commute_time_matrix += com_coords
      numberFrames += 1
    commute_time_matrix /= numberFrames
    return commute_time_matrix


#============================================================================
#               [[[ Interface: Analysis Data IO]]]
#============================================================================
class DataIO(object):
  """
  Role: Write/Read data generated by the analyses
  """
  @classmethod
  def write(cls,user_data,column_header_array,output_filename, output_data_tag):
    """
    Role: write data into an output file
    """
    Proxy_Pandas.write(user_data,column_header_array,output_filename, output_data_tag)


class Proxy_Pandas(object):
  """
  Role: a proxy interface to the Pandas module
  :return None
  """
  @classmethod
  def write(cls,user_data,column_header_array,output_filename, output_data_tag="data"):
    """
    Role: write data into an output file
    """
    output_filename += ".h5"
    output_data = pandas.DataFrame(user_data, columns=column_header_array)
    output_data.to_hdf(output_filename, output_data_tag)







#============================================================================
#               [[[ Client(View): User Interface Frontend ]]]
#============================================================================
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
    self.total_numFrames = self.userData.get_total_numFrames()


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

  def get_selected_coords(self):
    """
    Role: get the coordinates of selected atoms
    :return: numpy array
    """
    return self.userData.get_selected_coords()

  def get_selected_mass(self):
    """
    Role: get the masses of selected atoms
    :return: numpy array
    """
    return self.userData.get_selected_mass()

  def _build_com_matrix(self):
    """
    Role: build the center of mass matrix for selected atoms
    :return: numpy matrix (shape: Nx3)
    """
    self.comMatrix = self.userData.get_selected_com_coords()
    # print("comMatrx\n",self.comMatrix)

  def _build_pairwise_distance_com_matrix(self):
    """
    Role: build the pair-wise distance matrix using center of mass
          matrix
    :return: None
    """
    self._build_com_matrix()
    self.pairComMatrix = AnalysisTools.pairwise_distances(self.comMatrix)

  def next_frame(self):
    """
    Role: move on to next trajectory frame
    :return: None
    """
    self.userData.next_frame()

  def iter_get_selected_coords(self, frame_start=0, frame_end=None, frame_skip=1):
    """
    Role: a generator function that returns the coordinates of
          user-selected atoms iteratively for all trajectory frames
    :return: numpy array (shape: Nx3)
    """
    if frame_end == None:
      frame_end = self.total_numFrames
    return self.userData.iter_get_selected_coords(frame_start, frame_end, frame_skip)

  def iter_get_selected_com_coords(self, frame_start=0, frame_end=None, frame_skip=1):
    """
    Role: a generator function that returns the coordinates of
          user-selected atoms iteratively for all trajectory frames
    :return: numpy array (shape: Nx3)
    """
    if frame_end == None:
      frame_end = self.total_numFrames
    return self.userData.iter_get_selected_com_coords(frame_start, frame_end, frame_skip)

  def get_number_selected_residues(self):
    return self.userData.get_number_of_selected_residues()

  def get_number_selected_atoms(self):
    return self.userData.get_numAtomsSelected()

  def get_commute_time_matrix(self):
    """
    Role: compute the commute time matrix
    :return: None
    """
    self.commute_time_matrix = AnalysisTools.commute_time(self.iter_get_selected_com_coords(),
                                                          self.get_number_selected_residues())

    
  def save_commute_time(self, filename, output_dir="./"):
    """
    Save commute time into [filename]
    """
    data_column_header = self.ResIdList
    current_dir = os.path.realpath(output_dir) # get current path
    output_filename = os.path.join(current_dir, filename)
    DataIO.write(self.commute_time_matrix, data_column_header, output_filename, 'data')


