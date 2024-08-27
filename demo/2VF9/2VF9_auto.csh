#!/bin/csh -f
# generate assembly cif for 2 particles in the crystal asymmetric unit
# "ident" flag : 1st particle position in crystal frame
# 2VF9.x1.mat : transformation required to place coordinates in the 2nd particle position
runpt.csh 2VF9.cif  2VF9.biomt1 ident 2VF9.x1.mat

#make chimera pictures:
runchimera.csh

