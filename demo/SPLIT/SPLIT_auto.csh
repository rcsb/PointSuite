#!/bin/csh -f
# example  -- EM virus with split coordinates -- 3IYU/3N09

#create assembly cif for 1st entry
runpt.csh SPLIT.3iyu.cif SPLIT.biomt

#create simple pdb 1 to test image script below
mv build_pointsuite.cif  3iyu.cif

#create assembly cif for 2nd entry
runpt.csh SPLIT.3n09.cif SPLIT.biomt

#create simple pdb 2 to test image script below
mv build_pointsuite.cif  3n09.cif

#would normally run this on maxit-generated pdb from cif
RCSBvirusimages.csh 3iyu.cif 3n09.cif

