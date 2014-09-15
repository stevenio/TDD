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
import MDAnalysis
import os
#============================================================


#------------------------------------------------------------
# Parameters
#------------------------------------------------------------
data_dir = "data"
pdb_name = "sod.pdb"
psf_name = "sod.psf"

data_dir = os.path.realpath(data_dir)
my_pdb = os.path.join(data_dir, pdb_name)
my_psf = os.path.join(data_dir, psf_name)
print(my_psf)


#------------------------------------------------------------
# Read Input 
#------------------------------------------------------------


print("coordinates")
print(pdbObj.getCoords())
print("masses")
print(pdbObj.getMasses())
