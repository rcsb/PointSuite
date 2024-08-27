#!/bin/csh -f
#  script to validate matrix representation processing
#  C Lawson last updated Jan 2013
##########################################################################
if (! -e build_auth.cif) then
    echo "missing build_auth.cif file"
    exit
endif
if (! -e build_pointsuite.cif) then
    echo "missing build_pointsuite.cif file"
    exit
endif
# run chimera to check biomt in, biomt out
cp $PTSUITE/chimera/icos_cif.py .
echo " "
echo "making picture of author build"
env CHIMERAPDB=`pwd`/build_auth.cif chimera `pwd`/icos_cif.py  #before
echo "making picture of PointSuite build"
env CHIMERAPDB=`pwd`/build_pointsuite.cif chimera `pwd`/icos_cif.py  #after
echo " "
echo "displaying pictures"
echo "everything is okay if structures look the same (colors can be different)"
echo " "
display build_auth.cif.jpg  &
display build_pointsuite.cif.jpg  &
rm -f icos_cif.py
rm -f icos_cif.pyc
##########################################################################


