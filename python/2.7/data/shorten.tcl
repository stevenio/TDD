# select only a few residues
set psf_filename "segA.psf"
set pdb_filename "segA.pdb"
set dcd_filename "segA.dcd"
set output_filename "segA3"

set molId [mol new $psf_filename]

mol addfile $pdb_filename
mol addfile $dcd_filename waitfor all

set S [atomselect top "resid 41 to 43"]
$S frame 0
$S writepdb ${output_filename}.pdb
$S writepsf ${output_filename}.psf
animate write dcd ${output_filename}.dcd beg 1  waitfor all sel $S

exit
