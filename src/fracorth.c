//interactive program to calculate one set of
//orthogonal coords from fractional coords, given x y z and unit cell params
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXMAT 50
#include "mat44lib.c"


 int main () {
 double unitcell[6];
 double orthmat[4][4];
 double fracmat[4][4];
 double frac[3],orth[3];
 double det;
 int i;

 unitcell[0] = 100.0;
 unitcell[1] = 100.0;
 unitcell[2] = 100.0;
 unitcell[3] = 90.;
 unitcell[4] = 90.;
 unitcell[5] = 120.;
 printf ("type in 6 unit cell params a b c alpha beta gamma (Angstroms, degrees) : \n");
 scanf ("%lf", &unitcell[0]); scanf ("%lf", &unitcell[1]); scanf ("%lf", &unitcell[2]);
 scanf ("%lf", &unitcell[3]); scanf ("%lf", &unitcell[4]); scanf ("%lf", &unitcell[5]);

 printf ("type in 3 frac coords for orthogonalization \n");
 scanf ("%lf", &frac[0]); scanf ("%lf", &frac[1]); scanf ("%lf", &frac[2]);

 getOrthMat(unitcell,orthmat);
 det=invMat (orthmat,fracmat);

 printf("unit cell: %f %f %f %f %f %f\n", unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]);
 printf ("\n");
 for (i=0; i<3; i++) {
 printf("orthmat  : %12.6f%12.6f%12.6f    %12.6f\n", orthmat[i][0],orthmat[i][1],orthmat[i][2], orthmat[i][3]);}
 printf ("\n");
 for (i=0; i<3; i++) {
 printf("fracmat  : %12.6f%12.6f%12.6f    %12.6f\n", fracmat[i][0],fracmat[i][1],fracmat[i][2], fracmat[i][3]);}
 printf ("\n"); 

 for (i=0;i<3;i++){
 orth[i]=orthmat[i][0]*frac[0] + orthmat[i][1]*frac[1] + orthmat[i][2] * frac[2];}
 printf("orth coords  : %16.10f %16.10f %16.10f \n", orth[0],orth[1],orth[2]);

 return 1;
}
