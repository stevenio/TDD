# calculate center of mass distance between residues
# in order to generate  comparison data set for debugging
# the python code

set base_filename "segA3"
set psf_filename ${base_filename}.psf
set trajectory_filename ${base_filename}.pdb
set output_filename "${base_filename}_comdist.dat"


set OUT [open $output_filename w]


set molId [mol new $psf_filename ]
mol addfile $trajectory_filename watifor all

set S [atomselect top "name CA"]
set resids [$S get resid]

set extra "noh"
foreach resid_1 $resids {
  foreach resid_2 $resids {
      set S1 [atomselect top "resid $resid_1 and $extra"]
      set S2 [atomselect top "resid $resid_2 and $extra"]
      set com1 [measure center $S1 weight mass]
      set com2 [measure center $S2 weight mass]
      set comdist [veclength [vecsub $com1 $com2]]
      puts -nonewline $OUT "$comdist "
  }
  puts $OUT ""
}

close $OUT
exit
