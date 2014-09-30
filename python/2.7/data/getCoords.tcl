# extract coordinates of selected atoms

set base_filename "segA3"
set psf_filename ${base_filename}.psf
set trajectory_filename ${base_filename}.pdb
set output_filename "${base_filename}_coord.dat"


set OUT [open $output_filename w]


set molId [mol new $psf_filename ]
mol addfile $trajectory_filename watifor all

set S [atomselect top "name CA"]
set resids [$S get resid]

set extra "noh"
foreach resid_1 $resids {
      set S1 [atomselect top "resid $resid_1 and $extra"]
      set coords [$S1 get {x y z}]
      foreach xyz $coords {
	puts $OUT "$xyz"
      }
}

close $OUT
exit
