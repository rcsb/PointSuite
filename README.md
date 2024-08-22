2024-August-22 please note we are in the process of updating documentation for this repository in the "docs" folder, and making this repo the location for package download.

# Pointsuite

PointSuite is a set of programs to process macromolecular assemblies described by point and helical symmetry operations, with the goals of uniform annotation, archiving, and viewing.  In order to handle coordinates deposited in any orthogonal Cartesian frame, the relationships between the deposition, standard point and crystal frames are captured as frame transformations.  For example, the transformation required to move icosahedral virus structures from deposited position to the standard frame shown at left is calculated and recorded.  All point symmetries are fully handled; helical entries are handled only for non-crystal cases.

Reference: Lawson CL, Dutta SD, Westbrook JD, Henrick K, Berman HM (2008)   [Representation of viruses in the remediated PDB archive](http://journals.iucr.org/d/issues/2008/08/00/mv5020/index.html), Acta Cryst D, 874-882.Software package for handling matrices of biological assemblies containing point or helical symmetries.

(For now, publicly served at https://iqb.rutgers.edu/pointsuite/)

## Getting Started

Download the latest version here: https://github.com/rcsb/pointsuite/releases
You will also need UCSF Chimera installed and in your path

### Installing
Copy the package in the directory where you want to install the software and type the following commands:
```
tar xvzf pointsuite0.7.tgz
```
To compile:
```
cd pointsuite0.7
make
```
The package is composed of programs written in C along with C-shell-based scripts. It has been extensively tested on linux and mac-intel osx operating systems.

Setup:
```
source setup.csh
```
setup.csh script works for csh/tcsh shells. follow directions to set up your environment more permanently.

To make full use of the package, the graphics program UCSF Chimera should also be installed and in your path.

## Run the demos

```
cd demo
```
run the 1RUG (icosahedral virus) demo:
```
rundemo.csh 1RUG 
```
or, run all of the demos in batch:
```
rundemo.csh all 
```

DEMO LIST: 

-1RUG: Generate archival cif for icosahedral virus crystal structure.  
-1IFD:   Generate archival cif for helical virus fiber diffraction structure.  
-1EI7:   Generate archival cif for D17 symmetry particle.  
-1CGM:  Generate matrix representation for ~900 A length helical TMV-like virus.  
-1M4X:  Generate matrix representations for complex virus particle sub-assemblies.  
-IMPORT:  Generate BIOMT, CIF from typical author-uploaded example input files using importmats.  
-Additional icosahedral virus examples: 2XD8 (EM), 2W0C, 2VF9, 3N7X (X-ray).  
See also the Virus processing tutorial (in the html folder) for more info


## Authors/Contributors

Written/compiled by C. Lawson, with thanks to V.J. Reddy (TSRI) for sharing PDB2VIPER code (findframe);  Tom Goddard (UCSF) for Chimera scripts (runchimera.csh); Huanwang Yang (RCSB PDB) for importmats and cif-handling subroutines.

## Versions

version 0.5.8 (12 June 2007) initial stable release  
  
version 0.6 (20 June 2011) minor updates:  
\*importmats (from H. Yang)  handles additional matrix type (xncsrel) from CNS ncs.def.  
\*update of scripts automating image generation to work with v.1.4 Chimera and higher  
\*when run without arguments, runpt.csh autoscript now prints brief documentation  
\*new virus processing tutorial  
\*additional documentation now provided for utilities: importmats, autoscripts, multiplymats  
\*RCSBvirusimages.csh script to generate set of virus images for web display.  
  
version 0.7 (15 January 2013):  
\*improved cif parsing subroutines added by H. Yang (cifparse.c).  
\*file input reading improvements in importmats, findframe, makeassembly  
\*findframe single input file with matrices and coordinates can now be either PDB or CIF; optional 2nd file in BIOMT format (overrides 1st file matrices)  
\*new program cif2pdb creates simple pdb file (matrices, cryst1 record, coordinates) from cif (H. Yang).  
\*simplified scripts, PDB-dependency removed for runpt.csh  
\*RCSBvirusimages.csh script handles split entry cases (modifications by Ezra Peisach)  
