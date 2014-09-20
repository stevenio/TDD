
set base_filename "sod"
set psf_filename ${base_filename}.psf
set pdb_filename ${base_filename}.pdb
set output  "${base_filename}.dcd"

mol new $psf_filename
mol addfile $pdb_filename

set numFrames [molinfo top get numframes]
for {set fr 0 } {$fr < $numFrames} {incr fr} {
  molinfo top set a 100. frame $fr
  molinfo top set b 100. frame $fr
  molinfo top set c 100. frame $fr
}

set S [atomselect top all]

animate write dcd $output $S 
exit
