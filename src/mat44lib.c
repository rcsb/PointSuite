// pointsuite 0.5.8 
// author Cathy Lawson
// RCSB-PDB
// JUNE 2007

//library of subroutines for pointmats programs
//4x4 matrix representation is used
//for combined rotation/translation transformations 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#define MAXMAT 500
//int stype;
//int zsym;

//get rid of negative signs in front of printed zeros
double checkZero(double x, double precision) {
   double y;
   y = x;
   if (fabs(x)<precision)y=0;
   return y;
}

void multMat(double m1[4][4], double m2[4][4], double ret[4][4]){

  ret[0][0] = m1[0][0]*m2[0][0]+m1[0][1]*m2[1][0]+m1[0][2]*m2[2][0]+m1[0][3]*m2[3][0];
  ret[0][1] = m1[0][0]*m2[0][1]+m1[0][1]*m2[1][1]+m1[0][2]*m2[2][1]+m1[0][3]*m2[3][1];
  ret[0][2] = m1[0][0]*m2[0][2]+m1[0][1]*m2[1][2]+m1[0][2]*m2[2][2]+m1[0][3]*m2[3][2];
  ret[0][3] = m1[0][0]*m2[0][3]+m1[0][1]*m2[1][3]+m1[0][2]*m2[2][3]+m1[0][3]*m2[3][3];

  ret[1][0] = m1[1][0]*m2[0][0]+m1[1][1]*m2[1][0]+m1[1][2]*m2[2][0]+m1[1][3]*m2[3][0];
  ret[1][1] = m1[1][0]*m2[0][1]+m1[1][1]*m2[1][1]+m1[1][2]*m2[2][1]+m1[1][3]*m2[3][1];
  ret[1][2] = m1[1][0]*m2[0][2]+m1[1][1]*m2[1][2]+m1[1][2]*m2[2][2]+m1[1][3]*m2[3][2];
  ret[1][3] = m1[1][0]*m2[0][3]+m1[1][1]*m2[1][3]+m1[1][2]*m2[2][3]+m1[1][3]*m2[3][3];

  ret[2][0] = m1[2][0]*m2[0][0]+m1[2][1]*m2[1][0]+m1[2][2]*m2[2][0]+m1[2][3]*m2[3][0];
  ret[2][1] = m1[2][0]*m2[0][1]+m1[2][1]*m2[1][1]+m1[2][2]*m2[2][1]+m1[2][3]*m2[3][1];
  ret[2][2] = m1[2][0]*m2[0][2]+m1[2][1]*m2[1][2]+m1[2][2]*m2[2][2]+m1[2][3]*m2[3][2];
  ret[2][3] = m1[2][0]*m2[0][3]+m1[2][1]*m2[1][3]+m1[2][2]*m2[2][3]+m1[2][3]*m2[3][3];

  ret[3][0] = m1[3][0]*m2[0][0]+m1[3][1]*m2[1][0]+m1[3][2]*m2[2][0]+m1[3][3]*m2[3][0];
  ret[3][1] = m1[3][0]*m2[0][1]+m1[3][1]*m2[1][1]+m1[3][2]*m2[2][1]+m1[3][3]*m2[3][1];
  ret[3][2] = m1[3][0]*m2[0][2]+m1[3][1]*m2[1][2]+m1[3][2]*m2[2][2]+m1[3][3]*m2[3][2];
  ret[3][3] = m1[3][0]*m2[0][3]+m1[3][1]*m2[1][3]+m1[3][2]*m2[2][3]+m1[3][3]*m2[3][3];

}

double invMat(double m1[4][4], double ret[4][4]) {
//calculate inverse of 4x4 matrix, return determinant
//Laplacian expansion by minors

      double am, det;
      int i,ir,ii,j,jr,jj;
      double c[4][4];
      double x[3][3];

      for (ii=0;ii<4;ii++)
        for(jj=0;jj<4;jj++){
          i=-1;
          for (ir=0;ir<4;ir++){
            if (ir != ii){
              i++;
              j=-1;
              for (jr=0;jr<4;jr++){
                if (jr != jj) {
                j++;
                x[i][j] = m1[ir][jr];  }} }}
          //minor determinant
          am = x[0][0]*x[1][1]*x[2][2] - x[0][0]*x[1][2]*x[2][1] 
             + x[0][1]*x[1][2]*x[2][0] - x[0][1]*x[1][0]*x[2][2] 
             + x[0][2]*x[1][0]*x[2][1] - x[0][2]*x[1][1]*x[2][0];
          //printf ("am: %f\n", am);
          //cofactor 
          c[ii][jj] = am * pow(-1.0,(ii+jj));
        }
//---- calculate determinant
      det = 0.0;
      for (i=0;i<4;i++){
        det = m1[i][0] * c[i][0] + det;
          //printf ("det: %f\n", det);
        }
//---- get inverse matrix
      for (i=0;i<4;i++)
         for (j=0;j<4;j++){
           ret[i][j] = c[j][i]/det;
          //printf ("inverse elem: %f\n", ret[i][j]);
          }
      return det;
}


void defHelixMats(int exp[4], double nmat[4][4][4],double arot, 
                  double arise, int nmxz, int xsym, int pzsym){
//define source matrices for generating helical symmetries
double ar,cr;
int i;

//initialize symmetry generating matrices to identity, exponents to 0
 for (i=0;i<4;i++){
  nmat[i][0][0] = 1.0; nmat[i][0][1] = 0.0; nmat[i][0][2] = 0.0; nmat[i][0][3] = 0.0;
  nmat[i][1][0] = 0.0; nmat[i][1][1] = 1.0; nmat[i][1][2] = 0.0; nmat[i][1][3] = 0.0;
  nmat[i][2][0] = 0.0; nmat[i][2][1] = 0.0; nmat[i][2][2] = 1.0; nmat[i][2][3] = 0.0;
  nmat[i][3][0] = 0.0; nmat[i][3][1] = 0.0; nmat[i][3][2] = 0.0; nmat[i][3][3] = 1.0;
  exp[i]=0;}

//apply C symmetry about z-axis first, then symmetry about x, then rot/trans mat derived from arot, arise

//calculate nmat[0] from zsym
cr = pzsym;
cr = (360.0/cr) * (M_PI/180.0);
nmat[0][0][0] = cos(cr); nmat[0][0][1] = -1.0*sin(cr);    nmat[0][0][2]=0.0;  nmat[0][0][3]=0.0;
nmat[0][1][0] = sin(cr); nmat[0][1][1] = cos(cr);         nmat[0][1][2]=0.0;  nmat[0][1][3]=0.0;
nmat[0][2][0] = 0.0;     nmat[0][2][1] = 0.0;             nmat[0][2][2]=1.0;  nmat[0][2][3]=0.0;
nmat[0][3][0] = 0.0;     nmat[0][3][1] = 0.0;             nmat[0][3][2]=0.0;  nmat[0][3][3]=1.0;
exp[0] = abs(pzsym) - 1;

//nmat[1] = 2fold on x 
if (xsym != 1 && xsym != 2) xsym = 1;
nmat[1][0][0] = 1.0; nmat[1][0][1] = 0.0;  nmat[1][0][2]=0.0;  nmat[1][0][3]=0.0;
nmat[1][1][0] = 0.0; nmat[1][1][1] = -1.0; nmat[1][1][2]=0.0;  nmat[1][1][3]=0.0;
nmat[1][2][0] = 0.0; nmat[1][2][1] = 0.0;  nmat[1][2][2]=-1.0; nmat[1][2][3]=0.0;
nmat[1][3][0] = 0.0; nmat[1][3][1] = 0.0;  nmat[1][3][2]=0.0;  nmat[1][3][3]=1.0;
exp[1] = xsym - 1;

//nmat[2] = from arot, arise 
ar = arot * (M_PI/180.0);
nmat[2][0][0] = cos(ar); nmat[2][0][1] = -1.0*sin(ar);    nmat[2][0][2]=0.0;  nmat[2][0][3]=0.0;
nmat[2][1][0] = sin(ar); nmat[2][1][1] = cos(ar);         nmat[2][1][2]=0.0;  nmat[2][1][3]=0.0;
nmat[2][2][0] = 0.0;     nmat[2][2][1] = 0.0;             nmat[2][2][2]=1.0;  nmat[2][2][3]=arise;
nmat[2][3][0] = 0.0;     nmat[2][3][1] = 0.0;             nmat[2][3][2]=0.0;  nmat[2][3][3]=1.0;

//how many times to do arot/arise multiplication:
//exp[2] = (numbio/xznum) - 1 ;
exp[2] = nmxz ;

}

void defPointMats(int exp[4], double nmat[4][4][4],int *pstype, int *pzsym){
//define source matrices for generating point symmetries
double cr;
int xsym, i;

//for all paths, initialize symmetry generating matrices to identity, exponents to 0
 for (i=0;i<4;i++){
  nmat[i][0][0] = 1.0; nmat[i][0][1] = 0.0; nmat[i][0][2] = 0.0; nmat[i][0][3] = 0.0;
  nmat[i][1][0] = 0.0; nmat[i][1][1] = 1.0; nmat[i][1][2] = 0.0; nmat[i][1][3] = 0.0;
  nmat[i][2][0] = 0.0; nmat[i][2][1] = 0.0; nmat[i][2][2] = 1.0; nmat[i][2][3] = 0.0;
  nmat[i][3][0] = 0.0; nmat[i][3][1] = 0.0; nmat[i][3][2] = 0.0; nmat[i][3][3] = 1.0;
  exp[i]=0;}

//BEGIN HELICAL PATH
//H,C,D symmetry gen 
//if (stype == 'H' || stype == 'C' || stype == 'D') { //disabled H
if (*pstype == 'C' || *pstype == 'D') {

//apply C symmetry about z-axis first, then symmetry about x, then rot/trans mat derived from arot, arise

//calculate nmat[0] from zsym
cr = *pzsym;
cr = (360.0/cr) * (M_PI/180.0);
nmat[0][0][0] = cos(cr); nmat[0][0][1] = -1.0*sin(cr);    nmat[0][0][2]=0.0;  nmat[0][0][3]=0.0;
nmat[0][1][0] = sin(cr); nmat[0][1][1] = cos(cr);         nmat[0][1][2]=0.0;  nmat[0][1][3]=0.0;
nmat[0][2][0] = 0.0;     nmat[0][2][1] = 0.0;             nmat[0][2][2]=1.0;  nmat[0][2][3]=0.0;
nmat[0][3][0] = 0.0;     nmat[0][3][1] = 0.0;             nmat[0][3][2]=0.0;  nmat[0][3][3]=1.0;
exp[0] = abs(*pzsym) - 1;

//nmat[1] = 2fold on x 
if (*pstype == 'D') xsym = 2;
if (*pstype == 'C') xsym = 1;
if (xsym != 1 && xsym != 2) xsym = 1;
nmat[1][0][0] = 1.0; nmat[1][0][1] = 0.0;  nmat[1][0][2]=0.0;  nmat[1][0][3]=0.0;
nmat[1][1][0] = 0.0; nmat[1][1][1] = -1.0; nmat[1][1][2]=0.0;  nmat[1][1][3]=0.0;
nmat[1][2][0] = 0.0; nmat[1][2][1] = 0.0;  nmat[1][2][2]=-1.0; nmat[1][2][3]=0.0;
nmat[1][3][0] = 0.0; nmat[1][3][1] = 0.0;  nmat[1][3][2]=0.0;  nmat[1][3][3]=1.0;
exp[1] = xsym - 1;

/*
//nmat[2] = from arot, arise 
ar = arot * (M_PI/180.0);
nmat[2][0][0] = cos(ar); nmat[2][0][1] = -1.0*sin(ar);    nmat[2][0][2]=0.0;  nmat[2][0][3]=0.0;
nmat[2][1][0] = sin(ar); nmat[2][1][1] = cos(ar);         nmat[2][1][2]=0.0;  nmat[2][1][3]=0.0;
nmat[2][2][0] = 0.0;     nmat[2][2][1] = 0.0;             nmat[2][2][2]=1.0;  nmat[2][2][3]=arise;
nmat[2][3][0] = 0.0;     nmat[2][3][1] = 0.0;             nmat[2][3][2]=0.0;  nmat[2][3][3]=1.0;

//how many times to do arot/arise multiplication:
exp[2] = (numbio/xznum) - 1 ;
*/

} 
//END HELICAL PATH



//BEGIN TETRAHEDRAL PATH: T,O symmetry gen  
if (*pstype == 'T' || *pstype == 'O') {

//define source matrices in order of their application (right to left):

//2-fold on z
  nmat[0][0][0] =-1.0; nmat[0][0][1] = 0.0; nmat[0][0][2] = 0.0; nmat[0][0][3] = 0.0;
  nmat[0][1][0] = 0.0; nmat[0][1][1] =-1.0; nmat[0][1][2] = 0.0; nmat[0][1][3] = 0.0;
  nmat[0][2][0] = 0.0; nmat[0][2][1] = 0.0; nmat[0][2][2] = 1.0; nmat[0][2][3] = 0.0;
  nmat[0][3][0] = 0.0; nmat[0][3][1] = 0.0; nmat[0][3][2] = 0.0; nmat[0][3][3] = 1.0;
  exp[0] = 1;

//2-fold on y
  nmat[1][0][0] =-1.0; nmat[1][0][1] = 0.0; nmat[1][0][2] = 0.0; nmat[1][0][3] = 0.0;
  nmat[1][1][0] = 0.0; nmat[1][1][1] = 1.0; nmat[1][1][2] = 0.0; nmat[1][1][3] = 0.0;
  nmat[1][2][0] = 0.0; nmat[1][2][1] = 0.0; nmat[1][2][2] =-1.0; nmat[1][2][3] = 0.0;
  nmat[1][3][0] = 0.0; nmat[1][3][1] = 0.0; nmat[1][3][2] = 0.0; nmat[1][3][3] = 1.0;
  exp[1] = 1;

//3-fold on body diagonal
  nmat[2][0][0] = 0.0; nmat[2][0][1] = 0.0; nmat[2][0][2] = 1.0; nmat[2][0][3] = 0.0;
  nmat[2][1][0] = 1.0; nmat[2][1][1] = 0.0; nmat[2][1][2] = 0.0; nmat[2][1][3] = 0.0;
  nmat[2][2][0] = 0.0; nmat[2][2][1] = 1.0; nmat[2][2][2] = 0.0; nmat[2][2][3] = 0.0;
  nmat[2][3][0] = 0.0; nmat[2][3][1] = 0.0; nmat[2][3][2] = 0.0; nmat[2][3][3] = 1.0;
  exp[2] = 2;

//2-fold on xy
  nmat[3][0][0] = 0.0; nmat[3][0][1] = 1.0; nmat[3][0][2] = 0.0; nmat[3][0][3] = 0.0;
  nmat[3][1][0] = 1.0; nmat[3][1][1] = 0.0; nmat[3][1][2] = 0.0; nmat[3][1][3] = 0.0;
  nmat[3][2][0] = 0.0; nmat[3][2][1] = 0.0; nmat[3][2][2] =-1.0; nmat[3][2][3] = 0.0;
  nmat[3][3][0] = 0.0; nmat[3][3][1] = 0.0; nmat[3][3][2] = 0.0; nmat[3][3][3] = 1.0;
  if (*pstype == 'O') exp[3] = 1;
  if (*pstype == 'T') exp[3] = 0;   

} 
//END TETRAHEDRAL PATH

//BEGIN ICOSAHEDRAL PATH: I symmetry gen  
if (*pstype == 'I') {

  double hphi = (sqrt(5.0) + 1.0)/4.0;

//5-fold on (0,1,phi)
  nmat[0][0][0] = hphi - 0.5; nmat[0][0][1] = -1.0 * hphi; nmat[0][0][2] = 0.5;        nmat[0][0][3] = 0.0;
  nmat[0][1][0] = hphi;       nmat[0][1][1] = 0.5;         nmat[0][1][2] = hphi - 0.5; nmat[0][1][3] = 0.0;
  nmat[0][2][0] = -0.5;       nmat[0][2][1] = hphi - 0.5;  nmat[0][2][2] = hphi;       nmat[0][2][3] = 0.0;
  nmat[0][3][0] = 0.0;        nmat[0][3][1] = 0.0;         nmat[0][3][2] = 0.0;        nmat[0][3][3] = 1.0;
  exp[0] = 4;

//2-fold on z
  nmat[1][0][0] = -1.0; nmat[1][0][1] = 0.0; nmat[1][0][2] = 0.0; nmat[1][0][3] = 0.0;
  nmat[1][1][0] = 0.0; nmat[1][1][1] = -1.0; nmat[1][1][2] = 0.0; nmat[1][1][3] = 0.0;
  nmat[1][2][0] = 0.0; nmat[1][2][1] = 0.0; nmat[1][2][2] = 1.0; nmat[1][2][3] = 0.0;
  nmat[1][3][0] = 0.0; nmat[1][3][1] = 0.0; nmat[1][3][2] = 0.0; nmat[1][3][3] = 1.0;
  exp[1] = 1;

//2-fold on y
  nmat[2][0][0] = -1.0; nmat[2][0][1] = 0.0; nmat[2][0][2] = 0.0; nmat[2][0][3] = 0.0;
  nmat[2][1][0] = 0.0; nmat[2][1][1] = 1.0; nmat[2][1][2] = 0.0; nmat[2][1][3] = 0.0;
  nmat[2][2][0] = 0.0; nmat[2][2][1] = 0.0; nmat[2][2][2] = -1.0; nmat[2][2][3] = 0.0;
  nmat[2][3][0] = 0.0; nmat[2][3][1] = 0.0; nmat[2][3][2] = 0.0; nmat[2][3][3] = 1.0;
  exp[2] = 1;

//3-fold on body diagonal
  nmat[3][0][0] = 0.0; nmat[3][0][1] = 0.0; nmat[3][0][2] = 1.0; nmat[3][0][3] = 0.0;
  nmat[3][1][0] = 1.0; nmat[3][1][1] = 0.0; nmat[3][1][2] = 0.0; nmat[3][1][3] = 0.0;
  nmat[3][2][0] = 0.0; nmat[3][2][1] = 1.0; nmat[3][2][2] = 0.0; nmat[3][2][3] = 0.0;
  nmat[3][3][0] = 0.0; nmat[3][3][1] = 0.0; nmat[3][3][2] = 0.0; nmat[3][3][3] = 1.0;
  exp[3] = 2;

} 
//END ICOSAHEDRAL PATH
}


void calcIcosMats (double mats[60][4][4]) {
//generate 60 matrices for the icosahedron 
//pentamers occupy positions in std order of a tetrahedron (as in Int Tables s.g. P23)
//
//Let 1 = identity
//    5 = 5-fold on the position (0,1, phi)   (phi = (sqrt(5)+1)/2)
//    2a = 2-fold about the z-axis
//    2b = 2-fold about the y-axis
//    3  = 3-fold about the xyz body diagonal
// then to generate 60 matrices:
//  (1,3,3**2))(1,2b)(1,2a)(1,5,5**2,5**3, 5**4)


  double nmat[4][4][4];
  int i,j,k,n;
  int exp[4], expexp[4];
  int nmats,nlow,nhigh;
  double hphi = (sqrt(5.0) + 1.0)/4.0;

nmats = 4;  //number of source matrices 

//define source matrices in order of their application (right to left):

//5-fold on (0,1,phi)    
  nmat[0][0][0] = hphi - 0.5; nmat[0][0][1] = -1.0 * hphi; nmat[0][0][2] = 0.5;        nmat[0][0][3] = 0.0;
  nmat[0][1][0] = hphi;       nmat[0][1][1] = 0.5;         nmat[0][1][2] = hphi - 0.5; nmat[0][1][3] = 0.0;
  nmat[0][2][0] = -0.5;       nmat[0][2][1] = hphi - 0.5;  nmat[0][2][2] = hphi;       nmat[0][2][3] = 0.0;
  nmat[0][3][0] = 0.0;        nmat[0][3][1] = 0.0;         nmat[0][3][2] = 0.0;        nmat[0][3][3] = 1.0;
//define how many times the source matrix is multiplied by itself. For an n-fold, (n-1); so for a 5-fold, 4
  exp[0] = 4; 

//2-fold on z 
  nmat[1][0][0] = -1.0; nmat[1][0][1] = 0.0; nmat[1][0][2] = 0.0; nmat[1][0][3] = 0.0;
  nmat[1][1][0] = 0.0; nmat[1][1][1] = -1.0; nmat[1][1][2] = 0.0; nmat[1][1][3] = 0.0;
  nmat[1][2][0] = 0.0; nmat[1][2][1] = 0.0; nmat[1][2][2] = 1.0; nmat[1][2][3] = 0.0;
  nmat[1][3][0] = 0.0; nmat[1][3][1] = 0.0; nmat[1][3][2] = 0.0; nmat[1][3][3] = 1.0;
  exp[1] = 1;

//2-fold on y
  nmat[2][0][0] = -1.0; nmat[2][0][1] = 0.0; nmat[2][0][2] = 0.0; nmat[2][0][3] = 0.0;
  nmat[2][1][0] = 0.0; nmat[2][1][1] = 1.0; nmat[2][1][2] = 0.0; nmat[2][1][3] = 0.0;
  nmat[2][2][0] = 0.0; nmat[2][2][1] = 0.0; nmat[2][2][2] = -1.0; nmat[2][2][3] = 0.0;
  nmat[2][3][0] = 0.0; nmat[2][3][1] = 0.0; nmat[2][3][2] = 0.0; nmat[2][3][3] = 1.0;
  exp[2] = 1;

//3-fold on body diagonal 
  nmat[3][0][0] = 0.0; nmat[3][0][1] = 0.0; nmat[3][0][2] = 1.0; nmat[3][0][3] = 0.0;
  nmat[3][1][0] = 1.0; nmat[3][1][1] = 0.0; nmat[3][1][2] = 0.0; nmat[3][1][3] = 0.0;
  nmat[3][2][0] = 0.0; nmat[3][2][1] = 1.0; nmat[3][2][2] = 0.0; nmat[3][2][3] = 0.0;
  nmat[3][3][0] = 0.0; nmat[3][3][1] = 0.0; nmat[3][3][2] = 0.0; nmat[3][3][3] = 1.0;
  exp[3] = 2;


//expexp keeps track of matrices to be multiplied by each source matrix
//for icos, 
//       expexp[0]=1; 
//       expexp[1]=5*1=5; 
//       expexp[2]=2*5=10; 
//       expexp[3]=2*10=20;
//these values are calculated automatically below based on defined exp[i]

expexp[0]= 1;
if (nmats > 1) {for (i=1;i<nmats;i++) {expexp[i] = (exp[i-1] + 1)*expexp[i-1];}}


//return matrices, starting with identity
  mats[0][0][0] = 1.0; mats[0][0][1] = 0.0; mats[0][0][2] = 0.0; mats[0][0][3] = 0.0;
  mats[0][1][0] = 0.0; mats[0][1][1] = 1.0; mats[0][1][2] = 0.0; mats[0][1][3] = 0.0;
  mats[0][2][0] = 0.0; mats[0][2][1] = 0.0; mats[0][2][2] = 1.0; mats[0][2][3] = 0.0;
  mats[0][3][0] = 0.0; mats[0][3][1] = 0.0; mats[0][3][2] = 0.0; mats[0][3][3] = 1.0;

//k is output matrix index
k = 0;
  for (i=0;i<nmats;i++){
   //skip if nmat set is identity only
   if (exp[i]>0){
     for (j=0;j<exp[i];j++){
       //define source mats for multiplication
       nlow = expexp[i] * j;
       nhigh = expexp[i]*(j+1);
       for (n=nlow;n<nhigh;n++){
           k++;
           multMat(nmat[i],mats[n],mats[k]); }}
 }}}


void calcSymmMats (int exp[4], double nmat[4][4][4], double mats[MAXMAT][4][4]) {
//
// nmats are up to four 4x4 source matrices to build up complicated symmetries
//   given in order of application (right to left):
//
//  mats(a,b,c,d) = nmat[3]^(a=0-exp[3]) * nmat[2]^(b=0-exp[2]) 
//                * nmat[1]^(c=0-exp[1]) * nmat[0]^(d=0-exp[0])
//
// where exp[i] is the maximum number of times the source matrix will be applied,
//       and ^0 implies the identity matrix.
// e.g. for a  5-fold rotation, exp[i]=4 ;
// for nmat positions where only the identity matrix is to be applied,
//  set exp[i]=0 . and  nmat[i]  to any rotation
//
//example 1, a helix with C 10 symmetry:
//    nmat[0] = 10-fold about the z-axis           exp[0]=9
//    nmat[1] = 2-fold about the x-axis            exp[1]=0
//    nmat[2] = ang rotation + rise per subunit    exp[2]=#copies of 9 rotation units -1 
//    e.g., exp[2]=9 would result in 100 output matrices
//
//example 2, icosahedron:
//    nmat[0] = 5-fold about (0,1,phi)i            exp[0]=4
//    nmat[1] = 2-fold about the z-axis            exp[1]=1
//    nmat[2] = 2-fold about the y-axis            exp[2]=1
//    nmat[3] = 3-fold about xyz body diagonal     exp[3]=2


  int i,j,k,n;
  int nmats = 4;  //number of source matrices = 1st index of nmat
  int expexp[nmats];
  int nlow,nhigh;



//keep track of # matrices to be multiplied by each source matrix
//for the 1st one it is just identity
expexp[0]= 1;  
//2nd matrix is applied to mat set generated by the 1st matrix
//3rd matrix is applied to mat set generated by the 2nd matrix
//4nd matrix is applied to mat set generated by the 3rd matrix
if (nmats > 1) {for (i=1;i<nmats;i++) {expexp[i] = (exp[i-1] + 1)*expexp[i-1];}}


//return matrices, starting with identity
  mats[0][0][0] = 1.0; mats[0][0][1] = 0.0; mats[0][0][2] = 0.0; mats[0][0][3] = 0.0;
  mats[0][1][0] = 0.0; mats[0][1][1] = 1.0; mats[0][1][2] = 0.0; mats[0][1][3] = 0.0;
  mats[0][2][0] = 0.0; mats[0][2][1] = 0.0; mats[0][2][2] = 1.0; mats[0][2][3] = 0.0;
  mats[0][3][0] = 0.0; mats[0][3][1] = 0.0; mats[0][3][2] = 0.0; mats[0][3][3] = 1.0;

//k is output mats matrix index
k = 0;
  for (i=0;i<nmats;i++){
   //skip if nmat set is identity only
   if (exp[i]>0){
     for (j=0;j<exp[i];j++){
       //define source mats for multiplication
       nlow = expexp[i] * j;
       nhigh = expexp[i]*(j+1);
       for (n=nlow;n<nhigh;n++){
           k++;
           multMat(nmat[i],mats[n],mats[k]); }}
 }}}



double equiMat(double m1[4][4], double m2[4][4]){
//tests for equivalence of 2 matrices
//currently configured to ignore translations
  int i, j;
  double answer;
  answer = 0.0;
   for (i=0; i<4; i++)
     //ignore all translations
     //for (j=0; j<4; j++){
     for (j=0; j<3; j++){
    answer = answer + ((m1[i][j] - m2[i][j]) * (m1[i][j] - m2[i][j]));
    }
    answer = sqrt (answer);
    return answer;
}
    
void applyMat(double mat[4][4], double xyzin[4], double xyzout[4]) {
//applies 4x4 matrix to a coordinate vector
//expect mat as r11 r12 r13 t1 r21 r22 r23 t2 r31 r32 r33 t3 0 0 0 1
//expect coord vector as x y z 1

  xyzout[0] = mat[0][0]*xyzin[0] + mat[0][1]*xyzin[1] + mat[0][2]*xyzin[2] + mat[0][3]*xyzin[3];
  xyzout[1] = mat[1][0]*xyzin[0] + mat[1][1]*xyzin[1] + mat[1][2]*xyzin[2] + mat[1][3]*xyzin[3];
  xyzout[2] = mat[2][0]*xyzin[0] + mat[2][1]*xyzin[1] + mat[2][2]*xyzin[2] + mat[2][3]*xyzin[3];
  xyzout[3] = mat[3][0]*xyzin[0] + mat[3][1]*xyzin[1] + mat[3][2]*xyzin[2] + mat[3][3]*xyzin[3];
   }


double distXYZ (double xyz1[4],double xyz2[4]) {
//returns distance between a pair of coordinate vectors
   double dist;
   int i;
   dist = 0;
   for (i=0;i<3;i++){
     dist += (xyz1[i]-xyz2[i])*(xyz1[i] - xyz2[i]);}  //add squares of delta x, y, z
   dist = sqrt(dist);
   return dist;
   }


void getOrthMat(double uc[6], double orth[4][4]) {
//      a   b(cos(gamma))   c(cos(beta))
//      0   b(sin(gamma))   c(cos(alpha) - cos(beta) cos(gamma)) / sin(gamma)
//      0   0               Vol/(ab sin(gamma))
//Vol = abc(1 - cos**2(alpha) - cos**2(beta) - cos**2(gamma) + 2(cos(alpha) cos(beta) cos(gamma)))**1/2

  double Vol, Volmult;

//move u.c. angles to radians:
  uc[3] = uc[3] * (M_PI/180.0);
  uc[4] = uc[4] * (M_PI/180.0);
  uc[5] = uc[5] * (M_PI/180.0);

//calculate "Vol"
  Volmult =  1.0 - cos(uc[3])*cos(uc[3]) - cos(uc[4])*cos(uc[4]) - cos(uc[5])*cos(uc[5]);
  Volmult =  Volmult + 2*cos(uc[3])* cos(uc[4])*cos(uc[5]);
  Vol = uc[0]*uc[1]*uc[2]*sqrt(Volmult);

//calculate the orthogonalization matrix
  orth[0][0] = uc[0]; orth[0][1]=uc[1]*cos(uc[5]);  orth[0][2]=uc[2]* cos(uc[4]);                                      orth[0][3]=0;
  orth[1][0] = 0;     orth[1][1]=uc[1]*sin(uc[5]);  orth[1][2]=uc[2]*(cos(uc[3]) - cos(uc[4])*cos(uc[5]))/sin(uc[5]);  orth[1][3]=0;
  orth[2][0]=  0;     orth[2][1]=0;                 orth[2][2]=Vol/(uc[0]*uc[1]*sin(uc[5]));                           orth[2][3]=0;
  orth[3][0]=  0;     orth[3][1]=0;                 orth[3][2]=0;                                                      orth[3][3]=1;

//restore u.c. angles to degrees:
  uc[3] = uc[3] * (180.0/M_PI);
  uc[4] = uc[4] * (180.0/M_PI);
  uc[5] = uc[5] * (180.0/M_PI);
}

void getTransVector (double rot[3][3], double trans[3], double tret[3]) {
double  trrottr[4][4];
double  tmp[4][4];
double  vec1[4], vec2[4];
double  det;

trrottr[0][0] = 1 - rot[0][0];
trrottr[0][1] = 0 - rot[0][1];
trrottr[0][2] = 0 - rot[0][2];
trrottr[0][3] = 0;
trrottr[1][0] = 0 - rot[1][0];
trrottr[1][1] = 1 - rot[1][1];
trrottr[1][2] = 0 - rot[1][2];
trrottr[1][3] = 0;
trrottr[2][0] = 0 - rot[2][0];
trrottr[2][1] = 0 - rot[2][1];
trrottr[2][2] = 1 - rot[2][2];
trrottr[2][3] = 0;
trrottr[3][0] = 0;
trrottr[3][1] = 0;
trrottr[3][2] = 0;
trrottr[3][3] = 1;

det = invMat (trrottr, tmp);

vec1[0]=trans[0];
vec1[1]=trans[1];
vec1[2]=trans[2];
vec1[3]=1;

applyMat(tmp,vec1,vec2);

tret[0] = vec2[0];
tret[1] = vec2[1];
tret[2] = vec2[2];

}
