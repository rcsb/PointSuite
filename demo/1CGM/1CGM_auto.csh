#!/bin/csh -f
# this demo describes how to make an image of a helical structure with an arbitrary number of matrices
# based on  defined _pdbx_helical_symmetry parameters
#
# cif file with pdbx_helical symmetry parameters
set ciffile = 1CGM.cif
# 
# change value of _pdbx_helical_symmetry.number_of_operations from 49 to 900
sed -e "/number_of_operations/s/49/900/1" $ciffile >> 900mats.cif
#
# generate the matrices:
pointmats 900mats.cif
#
cat pointmats.biomt > 1cgm_900mats.pdb
cif2pdb $ciffile 1cgm_900mats.pdb_TMP >/dev/null
grep -v "REMARK 350" 1cgm_900mats.pdb_TMP >>1cgm_900mats.pdb
rm -f 1cgm_900mats.pdb_TMP
#
#
# use chimera multiscale model to make a pretty picture:
env CHIMERAPDB=`pwd`/1cgm_900mats.pdb chimera $PTSUITE/chimera/helix_biomt.py
display 1cgm_900mats.pdb.jpg &


