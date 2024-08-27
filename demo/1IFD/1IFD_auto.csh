#!/bin/csh -f
#helical virus

#to process:
runpt.csh 1IFD.cif  1IFD_biomt 

#then make pictures :
env CHIMERAPDB=`pwd`/build_auth.pdb chimera $PTSUITE/chimera/helix_biomt.py
env CHIMERAPDB=`pwd`/build_pointsuite.pdb chimera $PTSUITE/chimera/helix_biomt.py
display build_auth.pdb.jpg &
display build_pointsuite.pdb.jpg &



