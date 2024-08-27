#!/bin/csh -f
#run importmats on each of the test files
foreach mfile (DepUI_output.cif oper.cif 12element biomt_matrix cns-new-2wff.def mtrix ncs.def oldstyle.mtrix pdbset pointsuite.biomt pointsuite.cif )
importmats IMPORT.$mfile
mv import.biomt $mfile.import_biomt
mv import.cif $mfile.import_cif
mv import.matrix $mfile.import_matrix
end
