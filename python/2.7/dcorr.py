"""
Goal: calculate the cross correlation of inter-residue distances
Author: Steven(Yuhang) Wang
Date: Sep-09-2014

Blueprint:
1. extract all-atom coordinate set using [[ProDy]] module
2. calculate the inter-residue distance matrix, and save as data frames using [[Pandas]] module

External Dependencies:
Name     |     Link
-------- | -------------
Prody    | http://prody.csb.pitt.edu/downloads
Pandas   | http://pandas.pydata.org
-------- |  
"""

#--------------------#
# Python Builtins    #
#--------------------#
import unittest 

