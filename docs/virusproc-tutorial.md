# PROCESSING VIRUSES USING POINTSUITE  

- - -

## TWO STEPS:   

1.  [USE runpt.csh SCRIPT TO CREATE BIOLOGICAL ASSEMBLY CIF](#runpt)  
    
2.  [USE runchimera.csh TO CHECK RESULT](#runchimera)  
    

- - -

## RUNPT

CREATE BIOMT/ASSEMBLY.CIF USING runpt.csh

Input files (and what is read from them):  
  
* <cif> :  chain id list, crystal symmetry/unit cell info, coordinates  
* <author matrix file> : author-provided matrix file in any format that can be read by importmats. The file must define ONE particle only, and must include one identity matrix. If needed can manually convert author-provided format file using "importmats <filename>"  
* <ident or matfile> : List of transformations required to place deposited coordinates into all particle positions in crystal ident tag signifies coordinates already in crystal frame (shorthand for identity matrix) matfile contains elements for one non-identity matrix: [example](../demo/3N7X/3N7X.x0.mat).  #  ident  + matfiles = # of particles positioned in the crystal asymmetric unit. If coordinates are deposited in NONCRYSTAL frame, at least one matrix is needed.  
  
  
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

### 2XD8

 [![2xd8<br>image](../demo/2XD8/build_auth.pdb.jpg)](../demo/2XD8/build_auth.pdb.jpg)[![2xd8 image](../demo/2XD8/build_pointsuite.pdb.jpg)](../demo/2XD8/build_pointsuite.pdb.jpg)  (A) [2XD8](http://www.rcsb.org/pdb/explore/explore.do?structureId=2XD8) -- EM virus  runpt.csh [2xd8.cif](../demo/2XD8/2XD8.cif) [2xd8.biomt](../demo/2XD8/2XD8.biomt)  <br>runchimera.csh  [runpt.log](../demo/2XD8/runpt.log) [assembly.cif](../demo/2XD8/assembly.cif) 

### 2W0C

 [![2w0c<br>image](../demo/2W0C/build_auth.pdb.jpg)](../demo/2W0C/build_auth.pdb.jpg)[![2w0c image](../demo/2W0C/build_pointsuite.pdb.jpg)](../demo/2W0C/build_pointsuite.pdb.jpg)  (B) [2W0C](http://www.rcsb.org/pdb/explore/explore.do?structureId=2W0C) -- X-ray 1 particle/crystal au  <br>coordinates in crystal frame (usually true)  runpt.csh  [2w0c.cif](../demo/2W0C/2W0C.cif) [2w0c.biomt](../demo/2W0C/2W0C.biomt)  <br>runchimera.csh  [runpt.log](../demo/2W0C/runpt.log) [assembly.cif](../demo/2W0C/assembly.cif) 

### 2VF9

 [![2vf9<br>image](../demo/2VF9/build_auth.pdb.jpg)](../demo/2VF9/build_auth.pdb.jpg)[![2vf9 image](../demo/2VF9/build_pointsuite.pdb.jpg)](../demo/2VF9/build_pointsuite.pdb.jpg)  (C) [2VF9](http://www.rcsb.org/pdb/explore/explore.do?structureId=2VF9) -- X-ray [2 particles/crystal au](../html/multiparticle.html)   <br>coordinates in crystal frame  runpt.csh  [2vf9.cif](../demo/2VF9/2VF9.cif) [2vf9.biomt1](../demo/2VF9/2VF9.biomt1)\* ident  [x1.mat](../demo/2VF9/2VF9.x1.mat)\*\*  <br>\*(1st 60 BIOMT matrices in current public file)  <br>\*\*(61st BIOMT matrix in current public file)  <br>runchimera.csh  [runpt.log](../demo/2VF9/runpt.log) [assembly.cif](../demo/2VF9/assembly.cif) 

### 3N7X

 [![3n7x<br>image](../demo/3N7X/build_auth.pdb.jpg)](../demo/3N7X/build_auth.pdb.jpg)[![3n7x image](../demo/3N7X/build_pointsuite.pdb.jpg)](../demo/3N7X/build_pointsuite.pdb.jpg)  (D) [3N7X](http://www.rcsb.org/pdb/explore/explore.do?structureId=3N7X) -- X-ray 1/2 particle/crystal au  <br>coordinates in NONCRYSTAL frame  runpt.csh  [3n7x.cif](../demo/3N7X/3N7X.cif) [3n7x.biomt](../demo/3N7X/3N7X.biomt)  [x0.mat](../demo/3N7X/3N7X.x0.mat)  <br>runchimera.csh  [runpt.log](../demo/3N7X/runpt.log) [assembly.cif](../demo/3N7X/assembly.cif) 

  

- - -

  

## NCS/MTRIX RECORDS

FOR ALL X-RAY ENTRIES, NCS/MTRIX records are handled separately from above.  These matrices should be obtained from the deposited coordinate file and values placed in MTRIX ([\_struct\_ncs\_oper](http://mmcif.pdb.org/dictionaries/mmcif_pdbx.dic/Categories/struct_ncs_oper.html)) records for SF validation. Because pointsuite-generated ncs are based on exact point symmetry operations,  they can differ from author-refined values.   Currently there is no specific pointsuite module for this.   
  
  
  

- - -

  
C. Lawson updated August 2024
