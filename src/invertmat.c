//interactive program to calculate the inverse of a 4x4 matrix
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXMAT 50
#include "mat44lib.c"


 int main () {
 double forwardmat[4][4];
 double inversemat[4][4];
 double frac[3],orth[3];
 double det;
 int i;

 printf ("type in 12 matrix elements : \n");
 scanf ("%lf", &forwardmat[0][0]); 
 scanf ("%lf", &forwardmat[0][1]); 
 scanf ("%lf", &forwardmat[0][2]); 
 scanf ("%lf", &forwardmat[0][3]); 
 scanf ("%lf", &forwardmat[1][0]); 
 scanf ("%lf", &forwardmat[1][1]); 
 scanf ("%lf", &forwardmat[1][2]); 
 scanf ("%lf", &forwardmat[1][3]); 
 scanf ("%lf", &forwardmat[2][0]); 
 scanf ("%lf", &forwardmat[2][1]); 
 scanf ("%lf", &forwardmat[2][2]); 
 scanf ("%lf", &forwardmat[2][3]); 

 forwardmat[3][0] = 0.0;
 forwardmat[3][1] = 0.0;
 forwardmat[3][2] = 0.0;
 forwardmat[3][3] = 1.0;


 det=invMat (forwardmat,inversemat);

 for (i=0; i<3; i++) {
 printf("forwardmat  : %12.8f%12.8f%12.8f    %12.8f\n", forwardmat[i][0],forwardmat[i][1],forwardmat[i][2], forwardmat[i][3]);}
 printf ("\n");

 for (i=0; i<3; i++) {
 printf("inversemat  : %12.8f%12.8f%12.8f    %12.8f\n", inversemat[i][0],inversemat[i][1],inversemat[i][2], inversemat[i][3]);}
 printf ("\n"); 

 printf("determinant  : %12.8f\n", det);
 printf ("\n"); 

 return 1;
}
