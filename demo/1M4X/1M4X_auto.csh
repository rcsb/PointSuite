#!/bin/csh -f
# generate pdbs with BIOMT records for different 1M4X assemblies
cif2pdb 1M4X.cif 1M4X-tmp.pdb >/dev/null
#
# assembly #1 whole virus
multiplymats 1M4X.cif "(1-60)(61-88)" > mult.log
mv mult.biomt 1m4x_full.pdb
grep -v "REMARK 350" 1M4X-tmp.pdb >> 1m4x_full.pdb
# assembly #2 icos a.u.
multiplymats 1M4X.cif "(61-88)" >> mult.log
mv mult.biomt 1m4x_icosau.pdb
grep -v "REMARK 350" 1M4X-tmp.pdb >> 1m4x_icosau.pdb
#
# assembly #5  pentasymmetron
# include identity element (chimera doesn't handle mat sets without identity well)
multiplymats 1M4X.cif "(1,(1-5)(63-68))" >> mult.log
mv mult.biomt 1m4x_pentasym.pdb
grep -v "REMARK 350" 1M4X-tmp.pdb >> 1m4x_pentasym.pdb
#
rm -f 1M4X-tmp.pdb
rm -f mult.cif
rm -f mult.log
#
env CHIMERAPDB=`pwd`/1m4x_full.pdb chimera $PTSUITE/chimera/icos_biomt.py  
env CHIMERAPDB=`pwd`/1m4x_icosau.pdb chimera $PTSUITE/chimera/icos_biomt.py  
env CHIMERAPDB=`pwd`/1m4x_pentasym.pdb chimera $PTSUITE/chimera/icos_biomt.py  

display 1m4x_full.pdb.jpg &
display 1m4x_icosau.pdb.jpg &
display 1m4x_pentasym.pdb.jpg &

