#!/bin/csh -f
#  runpt.csh
#  script to process structures with point/helical symmetry
#  CL JAN 2013
###########################################################################
#**********************************************************
#* runpt.csh autoscript
#*
#* all noncrystal entries and MOST crystal entries:
#*
#*     runpt.csh <entry.cif>  <auth-upload-mats> 
#*
#*  <entry.cif> = near-final cif   
#*  <auth-upload-mats> = uploaded matrix file, any format read by importmats.
#*
#**********************************************************
#* complex crystal entries (non crystal frame OR multiple particles per a.u.):
#*
#*     runpt.csh <entry.cif>  <auth-upload-mats> ident <mat1> <mat2> ... 
#*
#*  <entry.cif> = near-final cif   
#*  <auth-upload-mats> = uploaded matrix file, any format read by importmats.
#*  ident keyword <mat#> file(s) = (up to 4) identify multiple independent 
#*    particle positions within a crystal asymmetric unit
#*
#*  ident keyword indicates coordinates in crystal frame for 1st particle position
#*
#* matrix file format is 4x 3 transformation matrix on 3 lines (free text):
#* r11 r12 r13 t1    e.g.    1.0 0.0 0.0 200.0
#* r21 r22 r23 t2            0.0 1.0 0.0 100.0
#* r31 r32 r33 t3            0.0 0.0 1.0 -10.0
#**********************************************************
##########################################################################
#delete old files
if (-e runpt.log)  rm -f runpt.log
if (-e findframe.cif) rm -f findframe.cif
if (-e assembly.cif) rm -f assembly.cif
##########################################################################
    echo "POINTSUITE AUTOSCRIPT runpt.csh "
if ($# < 2) then
    grep "*" $PTSUITE/bin/runpt.csh | grep -v grep
    exit
else if ($# == 2) then
    echo "2 input files: NON-CRYSTAL or SIMPLE CRYSTAL STRUCTURE ASSUMED" 
    echo "2 input files: NON-CRYSTAL or SIMPLE CRYSTAL STRUCTURE ASSUMED"  > runpt.log
    echo " "  >> runpt.log
else if ($# > 6) then
    echo "too many transformations, script and makeassembly program must be modified"
    exit
endif
##########################################################################
set cif_file    = ${1}   # 
set auth_mat    = ${2}   # author matrices in BIOMT format 
set p2xmat1     = ${3}   # crystal transformations
set p2xmat2     = ${4}
set p2xmat3     = ${5} 
set p2xmat4     = ${6} 
##########################################################################
#transform uploaded coordinates
importmats $auth_mat > /dev/null
set icount = `grep -c "BIOMT1" import.biomt`
if (! -e import.biomt|| $icount == 0 ) then
     echo "ERROR: importmats did not run successfully, check author file"
     exit
endif
echo $icount "matrices have been imported"
echo " "
rm -f import.matrix
##########################################################################
#get transformation from deposited to standard icos frame
findframe $cif_file import.biomt >> runpt.log
if (! -e findframe.cif) then
     echo "ERROR: findframe did not run successfully, check matrices"
    exit
endif
##########################################################################
# prepare input file for makeassembly
cp findframe.cif makextrans.cif
makextrans.csh makextrans.cif $p2xmat1 X0
makextrans.csh makextrans.cif $p2xmat2 X1
makextrans.csh makextrans.cif $p2xmat3 X2
makextrans.csh makextrans.cif $p2xmat4 X3
makeassembly $cif_file makextrans.cif >> runpt.log
if (! -e assembly.cif) then
    echo "ERROR: problem creating assembly cif file"
    exit
endif
##########################################################################
# legacymats no longer required for maxit
##########################################################################
# pdb's needed for chimera -- working to remove this pdb dependency
# generate author build pdb
cif2pdb $cif_file build_auth.pdbTMP >/dev/null
cat  import.biomt  > build_auth.pdb
grep -v "REMARK 350 " build_auth.pdbTMP >> build_auth.pdb
rm -f build_auth.pdbTMP
##########################################################################
#generate pointsuite build pdb
cat assembly.biomt   > build_pointsuite.pdb
grep -v "REMARK 350 " build_auth.pdb >> build_pointsuite.pdb
##########################################################################
# generate cif builds for chimera, first remove existing assembly cif if any
clear-assembly.csh $cif_file >> runpt.log
cat edited_file.cif import.cif > build_auth.cif
cat edited_file.cif assembly.cif > build_pointsuite.cif
#successful script
rm -f edited_file.cif
rm -f makextrans.cif
rm -f findframe.cif
rm -f import.biomt
echo 'SUCCESSFUL COMPLETION OF SCRIPT :)  check assembly.cif'
echo ""
echo 'use runchimera.csh to compare author vs. pointsuite builds'
echo ""
echo "RCSBvirusimages:"
echo 'use RCSBvirusimages.csh <final-PDB-file(s)> for single entry'
echo 'for single entry (1 file)  or split entry (multiple files)'
##########################################################################
