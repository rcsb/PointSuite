#!/bin/csh -f
#
# 3N7X X-ray, 1/2 particle/crystal a.u., coordinates deposited in non-crystal frame
#
#to generate assembly cif :
runpt.csh 3N7X.cif 3N7X.biomt 3N7X.x0.mat

#make chimera pictures:
runchimera.csh

