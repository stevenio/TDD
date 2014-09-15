"""
Goal: compute residue pairwise center of mass distances
Author: Yuhang Wang
Date: 09-15-2014
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
import MDAnalysis as mda
import MDAnalysis.core.parallel.distances as mda_dist
import numpy as np
import os
#============================================================


#------------------------------------------------------------
# Parameters
#------------------------------------------------------------
data_dir = "data"
case = 2
if case == 1:
  pdb_name = "sod.pdb"
  psf_name = "sod.psf"
elif case == 2:
  pdb_name = "segA.pdb"
  psf_name = "segA.psf"

data_dir = os.path.realpath(data_dir)
my_pdb = os.path.join(data_dir, pdb_name)
my_psf = os.path.join(data_dir, psf_name)
print(my_psf)


#------------------------------------------------------------
# Read Input 
#------------------------------------------------------------

univ = mda.Universe(my_psf,my_pdb)
SL = univ.selectAtoms("segid A and resid 42:43 and not (name H*)")
print(univ.atoms)
print(SL)

for fr in univ.trajectory:
  print("\nframe: {0}\n".format(fr.frame))
  print(np.shape(SL.get_positions()))
  print((SL.get_positions()))
  pairDist = mda_dist.distance_array(SL.get_positions(),SL.get_positions())
  print(pairDist.shape)
  print(SL.masses())
  print("indices:\n",SL.indices())

