# PROCESSING VIRUSES USING POINTSUITE  

- - -

## TWO STEPS:   

1.  USE runpt.csh SCRIPT TO CREATE BIOLOGICAL ASSEMBLY CIF
    
2.  USE runchimera.csh TO CHECK RESULT
    

- - -

## RUNPT

CREATE BIOMT/ASSEMBLY.CIF USING runpt.csh

Input files (and what is read from them):  
  
* \<cif\> :  chain id list, crystal symmetry/unit cell info, coordinates  
* \<author matrix file\> : author-provided matrix file in any format that can be read by importmats. The file must define ONE particle only, and must include one identity matrix. If needed can manually convert author-provided format file using "importmats <filename>"  
* \<ident or matfile\> : List of transformations required to place deposited coordinates into all particle positions in crystal ident tag signifies coordinates already in crystal frame (shorthand for identity matrix) matfile contains elements for one non-identity matrix: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/3N7X/3N7X.x0.mat" target="_blank">example</a>. The number of ident  + matfiles = number of particles positioned in the crystal asymmetric unit. If coordinates are deposited in NONCRYSTAL frame, at least one matrix is needed.  
  
  
Command:  
  
MOST CASES: 
```
runpt.csh <cif> <biomt>  
```
  
COMPLEX X-ray: (either out-of-frame coordinates or multiple positions of point symmetry particle in the crystal a.u. )  
```
runpt.csh <cif> <biomt> <matfile1> <matfile2> <matfile3>  
```

  
Output :  
  
On successful run command line will state "successful completion of script:  check assembly.cif"  
  
* runpt.log: log file with processing details.  
* assembly.cif and assembly.biomt: files containing BIOMT information. CIF also describes how to move virus particle to standard icos frame, and how to generate standard icos virus subassemblies.  PTSUITE BIOMT = author matrices adjusted for (1) standard icos. order; (2) exact icosahedral symmetry.  
* build_auth.pdb, build_wwpdb.pdb :  file with author-provided BIOMT + coordinates and  file with final BIOMT + coordinates.  These can be inspected/compared in chimera using pointsuite "runchimera.csh" script, or manually by applying Chimera multiscale model/biomt option.  
* assembly.ncs (X-ray only): these files have pointsuite-calculated NCS.  Useful to compare/check against author-provided NCS, but DO NOT INCLUDE IN PROCESSED ENTRY UNLESS AUTHOR NCS is determined to not be available.  SEE [BELOW](#ncs).  
  
  
## DEMO EXAMPLES

### [2XD8 EM virus](https://www.rcsb.org/structure/2XD8)

UCSF Chimera-created images: author matrices (left) vs PointSuite matrices (right)

<img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/build_auth.cif.jpg"> <img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/build_pointsuite.cif.jpg">

```
runpt.csh 2XD8.cif 2XD8.biomt
runchimera.csh 
```
input files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/2XD8.cif" target="_blank">2XD8.cif</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/2XD8.biomt" target="_blank">2XD8.biomt</a>

output files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/runpt.log" target="_blank">runpt.log</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/assembly.cif" target="_blank">assembly.cif</a>

- - -

### [2W0C X-ray 1 particle/crystal au, coordinates in crystal frame](https://www.rcsb.org/structure/2W0C)
UCSF Chimera-created images: author matrices (left) vs PointSuite matrices (right)

<img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/2W0C/build_auth.cif.jpg"> <img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/2W0C/build_pointsuite.cif.jpg">

```
runpt.csh 2W0C.cif 2W0C.biomt
runchimera.csh 
```

input files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2W0C/2W0C.cif" target="_blank">2W0C.cif</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2W0C/2W0C.biomt" target="_blank">2W0C.biomt</a>

output files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2W0C/runpt.log" target="_blank">runpt.log</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2W0C/assembly.cif" target="_blank">assembly.cif</a>

- - -


### [2VF9 X-ray 2 particles/crystal au, coordinates in crystal frame](https://www.rcsb.org/structure/2VF9)
UCSF Chimera-created images: author matrices (left) vs PointSuite matrices (right)

<img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/build_auth.cif.jpg"><img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/build_pointsuite.cif.jpg">

```
runpt.csh 2VF9.cif 2VF9.biomt ident x1.mat
runchimera.csh 
```

input files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/2VF9.cif" target="_blank">2VF9.cif</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/2VF9.biomt" target="_blank">2VF9.biomt</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/2VF9.x1.mat" target="_blank">2VF9.x1.mat</a>

output files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/runpt.log" target="_blank">runpt.log</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2VF9/assembly.cif" target="_blank">assembly.cif</a>

- - -


### [3N7X X-ray 1/2 particle/crystal au, coordinates NOT in crystal frame](https://www.rcsb.org/structure/3N7X]
(note: example uses version 1 of this entry, in 2023 it was reversioned with coordinates moved to crystal frame)

UCSF Chimera-created images: author matrices (left) vs PointSuite matrices (right)

<img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/3N7X/build_auth.cif.jpg"><img height="200" src="https://github.com/rcsb/PointSuite/blob/master/demo/3N7X/build_pointsuite.cif.jpg">

```
runpt.csh 3N7X.cif 3N7X.biomt x0.mat
runchimera.csh 
```

input files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/3N7X/3N7X.cif" target="_blank">3N7X.cif</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/2XD8/2XD8.biomt" target="_blank">3N7X.biomt</a>

output files: <a href="https://github.com/rcsb/PointSuite/blob/master/demo/3N7X/runpt.log" target="_blank">runpt.log</a> <a href="https://github.com/rcsb/PointSuite/blob/master/demo/3N7X/assembly.cif" target="_blank">assembly.cif</a>

- - -

  

## NCS (MTRIX) RECORDS

FOR ALL X-RAY ENTRIES, NCS (MTRIX) records are handled separately from above.  These matrices should be obtained from the deposited coordinate file and values placed in MTRIX <a href="http://mmcif.pdb.org/dictionaries/mmcif_pdbx.dic/Categories/struct_ncs_oper.html" target="_blank">_struct_ncs_oper</a> records for SF validation. Because pointsuite-generated ncs are based on exact point symmetry operations,  they can differ from author-refined values.   Currently there is no specific pointsuite module for this.   
  
  
  

- - -

  
C. Lawson updated August 2024
