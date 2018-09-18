//subroutines for findframe (called by findPTV.c)
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define D_PRECISION 0.001

int decon(double t[3][3], double * i, double * j, double * k, double * angle){
  int iz, jz, kz, ijp, ikp, jkp; //boolean flags for zero and pos values
  double mag;
  double d;

  //test for angle = 180 or zero
  if(fabs(t[0][1]-t[1][0])<D_PRECISION && fabs(t[0][2]-t[2][0])<D_PRECISION
      && fabs(t[1][2]-t[2][1])<D_PRECISION){
    //test for identity matrix, if so, just return zero values
    if(fabs(t[0][1]+t[1][0])<D_PRECISION && fabs(t[0][2]+t[2][0])<D_PRECISION
        && fabs(t[1][2]+t[2][1])<D_PRECISION && t[0][0]>(1-D_PRECISION)
        && t[0][0]<(1+D_PRECISION) && t[1][1]>(1-D_PRECISION) && t[1][1]<(1+D_PRECISION)
	&& t[2][2]>(1-D_PRECISION) && t[2][2]<(1+D_PRECISION)){
      *i = 0.0; *j = 0.0; *k = 0.0; *angle = 0.0;
      return 1;
    }
    else{
      //angle is 180
      //printf ("\n");
      //printf ("%8.3f %8.3f %8.3f\n", t[0][0], t[0][1], t[0][2]);
      //printf ("%8.3f %8.3f %8.3f\n", t[1][0], t[1][1], t[1][2]);
      //printf ("%8.3f %8.3f %8.3f\n", t[2][0], t[2][1], t[2][2]);
      *angle = M_PI;
      *i = sqrt((t[0][0]+1)/2);
      *j = sqrt((t[1][1]+1)/2);
      *k = sqrt((t[2][2]+1)/2);
      //printf ("i j k %8.3f %8.3f %8.3f\n", *i, *j, *k);

      //calculate signs of axis
      iz = fabs(*i)<D_PRECISION;
      jz = fabs(*j)<D_PRECISION;
      kz = fabs(*k)<D_PRECISION;
      //printf ("iz  jz  kz  %i %i %i\n", iz,jz,kz);
      ijp = t[0][1] > 0;
      ikp = t[0][2] > 0;
      jkp = t[1][2] > 0;
      //printf ("ijp ikp jkp %i %i %i\n", ijp,ikp,jkp);

      //one axis non-zero:
      if((iz && jz) || (iz && kz) || (jz && kz)) return 1;
      //two axes non-zero:
      if((iz && jkp)|| (jz && ikp) ||(kz && ijp))  return 1;
      if(iz && !jkp) {*j=-*j; return 1;}
      if(jz && !ikp) {*k=-*k; return 1;}
      if(kz && !ijp) {*i=-*i; return 1;}
      //three non-zero axes:
      if (ijp && ikp && jkp) return 1;
      if ((!ijp) && (!ikp) && (!jkp)) return 1;
      if (ijp) {*k=-*k; return 1;}
      if (ikp) {*j=-*j; return 1;}
      if (jkp) {*i=-*i; return 1;}
/*
// this code was buggy...
      if(iz && !jz && !kz) *j = -*j;
      else if (jz && !kz) *k = -*k;
      else if (kz) *i = -*i;
      else if (ijp && ikp && jkp) return;
      else if (jkp) *i = -*i;
      else if (ikp) *j = -*j;
      else if (ijp) *k = -*k;
      *i = -*i; *j = -*j; *k = -*k;
*/
      return 1;
    }
  }
    
  mag = -sqrt((t[2][1]-t[1][2])*(t[2][1]-t[1][2])+
             (t[0][2]-t[2][0])*(t[0][2]-t[2][0])+
            (t[1][0]-t[0][1])*(t[1][0]-t[0][1]));
  d = (t[0][0]+t[1][1]+t[2][2]-1)/2;
  if(d>1)
    d = 1.00;
  if(d<-1)
    d = -1.00;
  *angle = acos(d);
  *i = (t[2][1]-t[1][2])/mag;
  *j = (t[0][2]-t[2][0])/mag;
  *k = (t[1][0]-t[0][1])/mag;
  return 1;
}

int recon(double i, double j, double k, double angle, double t[3][3]){
  double s = sin(angle); double c = cos(angle);
  double a = 1.0-c;
  double x = i; double y = j; double z = k;

  t[0][0] = a*x*x+c; t[0][1] = a*x*y-z*s; t[0][2] = a*x*z+y*s;
  t[1][0] = a*x*y+z*s; t[1][1] = a*y*y+c; t[1][2] = a*y*z-x*s;
  t[2][0] = a*x*z-y*s; t[2][1] = a*y*z+x*s; t[2][2] = a*z*z+c;

  return 1;
}

void crossProduct(double v1[3], double v2[3], double ret[3]){
  ret[0] =     v1[1] * v2[2] - v1[2] * v2[1];
  ret[1] = -1*(v1[0] * v2[2] - v1[2] * v2[0]);
  ret[2] =     v1[0] * v2[1] - v1[1] * v2[0];
}
                                                                                                                             
double dotProduct(double v1[3], double v2[3]){
  double ret;
  ret = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  return ret;
}

/*
void normalize(double v[3]){
  double length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (length == 0.0) length = 0.000000000001;
  v[0] = v[0]/length;
  v[1] = v[1]/length;
  v[2] = v[2]/length;
}
*/

void normalize(double v[3]){
  double length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (length != 0.0){ 
  v[0] = v[0]/length;
  v[1] = v[1]/length;
  v[2] = v[2]/length;}
  else {
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;}
}
                                                                                                                             
void multiplyMatrices(double m1[3][3], double m2[3][3], double ret[3][3]){
  ret[0][0] = m1[0][0]*m2[0][0]+m1[0][1]*m2[1][0]+m1[0][2]*m2[2][0];
  ret[0][1] = m1[0][0]*m2[0][1]+m1[0][1]*m2[1][1]+m1[0][2]*m2[2][1];
  ret[0][2] = m1[0][0]*m2[0][2]+m1[0][1]*m2[1][2]+m1[0][2]*m2[2][2];
  ret[1][0] = m1[1][0]*m2[0][0]+m1[1][1]*m2[1][0]+m1[1][2]*m2[2][0];
  ret[1][1] = m1[1][0]*m2[0][1]+m1[1][1]*m2[1][1]+m1[1][2]*m2[2][1];
  ret[1][2] = m1[1][0]*m2[0][2]+m1[1][1]*m2[1][2]+m1[1][2]*m2[2][2];
  ret[2][0] = m1[2][0]*m2[0][0]+m1[2][1]*m2[1][0]+m1[2][2]*m2[2][0];
  ret[2][1] = m1[2][0]*m2[0][1]+m1[2][1]*m2[1][1]+m1[2][2]*m2[2][1];
  ret[2][2] = m1[2][0]*m2[0][2]+m1[2][1]*m2[1][2]+m1[2][2]*m2[2][2];
}
                                                                                                                             
double radToDeg(double rad){
  return rad*(180.0/M_PI);
}

void multiplyMatVec(double m1[3][3], double v1[3], double v2[3]) {
 v2[0] = m1[0][0] * v1[0] + m1[0][1] * v1[1] + m1[0][2] * v1[2];
 v2[1] = m1[1][0] * v1[0] + m1[1][1] * v1[1] + m1[1][2] * v1[2];
 v2[2] = m1[2][0] * v1[0] + m1[2][1] * v1[1] + m1[2][2] * v1[2];
}

