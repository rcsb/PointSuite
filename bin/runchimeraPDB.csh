#!/bin/csh -f
#  script to validate matrix representation processing
#  C Lawson last updated Jan 2013
##########################################################################
if (! -e build_auth.pdb) then
    echo "missing build_auth.pdb file"
    exit
endif
if (! -e build_pointsuite.pdb) then
    echo "missing build_pointsuite.pdb file"
    exit
endif
# run chimera to check biomt in, biomt out
cp $PTSUITE/chimera/icos_biomt.py .
echo " "
echo "making picture of author build"
env CHIMERAPDB=`pwd`/build_auth.pdb chimera `pwd`/icos_biomt.py  #before
echo "making picture of PointSuite build"
env CHIMERAPDB=`pwd`/build_pointsuite.pdb chimera `pwd`/icos_biomt.py  #after
echo " "
echo "displaying pictures"
echo "everything is okay if structures look the same (colors can be different)"
echo " "
display build_auth.pdb.jpg  &
display build_pointsuite.pdb.jpg  &
rm -f icos_biomt.py
rm -f icos_biomt.pyc
##########################################################################


