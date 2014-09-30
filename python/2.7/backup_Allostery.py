__author__ = 'stevenw'



class Allostery(object):
  """
  Goal: get the allosteric relations between residues
  """

  def __init__(self, data_src_dir, psf_filename, trajectory_filename):
    self.psf_filename = os.path.join(data_src_dir, psf_filename)
    self.trajectory_filename = os.path.join(data_src_dir, trajectory_filename)
    self.universe = mda.Universe(self.psf_filename, self.trajectory_filename)
    self.resid_list   = []
    self.segid_list = False
    self.com_coords = False
    self.commute_time = numpy.array([])


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
