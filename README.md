# Pointsuite

PointSuite is a set of programs to process macromolecular assemblies described by point and helical symmetry operations, with the goals of uniform annotation, archiving, and viewing.  In order to handle coordinates deposited in any orthogonal Cartesian frame, the relationships between the deposition, standard point and crystal frames are captured as frame transformations.  For example, the transformation required to move icosahedral virus structures from deposited position to the standard frame shown at left is calculated and recorded.  All point symmetries are fully handled; helical entries are handled only for non-crystal cases.

Lawson CL, Dutta SD, Westbrook JD, Henrick K, Berman HM (2008)   [Representation of viruses in the remediated PDB archive](http://journals.iucr.org/d/issues/2008/08/00/mv5020/index.html), Acta Cryst D, 874-882.Software package for handling matrices of biological assemblies containing point or helical symmetries.

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
  *1RUG: Generate archival cif for icosahedral virus crystal structure.⋅⋅
  *1IFD:   Generate archival cif for helical virus fiber diffraction structure.⋅⋅
  *1EI7:   Generate archival cif for D17 symmetry particle.⋅⋅
  *1CGM:  Generate matrix representation for ~900 A length helical TMV-like virus.⋅⋅
  *1M4X:  Generate matrix representations for complex virus particle sub-assemblies.⋅⋅
  *IMPORT:  Generate BIOMT, CIF from typical author-uploaded example input files using importmats.⋅⋅
*Additional icosahedral virus examples: 2XD8 (EM), 2W0C, 2VF9, 3N7X (X-ray). ⋅⋅
See also the Virus processing tutorial (in the html folder) for more info.


## Authors/Contributors

Written/compiled by C. Lawson, with thanks to V.J. Reddy (TSRI) for sharing PDB2VIPER code (findframe);  Tom Goddard (UCSF) for Chimera scripts (runchimera.csh); Huanwang Yang (RCSB PDB) for importmats and cif-handling subroutines.
