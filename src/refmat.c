//subroutine for findframe
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int refinePDB2PT(double rot[3][3], 
                 double tr[3], 
                 double pmat[MAXMAT][3][3], 
                 double pt[MAXMAT][3], 
                 int pnum,  
                 double refAtom[3],
                 int *pstype,
                 int *pzsym){
  int i, j, k ;
  double p2i[4][4];
  double matv[4];
  double mat44[4][4], identmat[4][4];
  double allmats[MAXMAT][4][4];
  double xyzfixed[MAXMAT][4];
  double xyzrefin[MAXMAT][4], xyztmp[MAXMAT][4];
  double movcom[4];
  double reftrans[4][4],reftinv[4][4],refrot[4][4];
  double testrot[4][4],testtrans[4][4];
  double det;
  double dist,testdist,dsum1,dsum2;
  int xyzIndex[MAXMAT];
  double umat[3][3],rot33[3][3]; 
  double nmat[4][4][4];
  int    exp[4];

  dsum1=0;
  dsum2=0;
  identmat[0][0]=1; identmat[0][1]=0; identmat[0][2]=0; identmat[0][3]=0;
  identmat[1][0]=0; identmat[1][1]=1; identmat[1][2]=0; identmat[1][3]=0;
  identmat[2][0]=0; identmat[2][1]=0; identmat[2][2]=1; identmat[2][3]=0;
  identmat[3][0]=0; identmat[3][1]=0; identmat[3][2]=0; identmat[3][3]=1;
     
    
//copy center of mass (com) to 1st fixed array vector 
    xyzfixed[0][0]=refAtom[0]; 
    xyzfixed[0][1]=refAtom[1]; 
    xyzfixed[0][2]=refAtom[2]; 
    xyzfixed[0][3]=1;


// move pdb2pt to 4x4 matrix 
  multMat(identmat,identmat,p2i);
  for (i=0; i<3;i++){
    p2i[i][3]=tr[i];
    for (j=0;j<3;j++){
       p2i[i][j] = rot[i][j];}}

// move input mats to 4x4 matrices 
 for(k=0;k<pnum;k++){
  multMat(identmat,identmat,allmats[k]);
  for (i=0; i<3;i++){
    allmats[k][i][3]=pt[k][i];
    for (j=0;j<3;j++){
       allmats[k][i][j] = pmat[k][i][j];}}}

//generate moving constellation in std frame, [p2i][inmats][c.o.m.]
//best fit of moving constellation will be calculated as [RT][p2i][inmats][c.o.m.]

  for (i=0;i<pnum;i++){ 
           multMat(p2i,allmats[i],mat44); 
          applyMat(mat44,xyzfixed[0],xyzrefin[i]);}


//generate fixed constellation in icos frame, [stdmats][p2i][c.o.m.]

    //get the std set of matrices for the point group
    defPointMats(exp,nmat,pstype,pzsym); 
    calcSymmMats(exp,nmat,allmats); 

   //transform c.o.m. by p2i
   applyMat(p2i,xyzfixed[0],matv);  
   for (j=0;j<4;j++) xyzfixed[0][j] = matv[j]; 
   for (i=1;i<pnum;i++) applyMat(allmats[i],xyzfixed[0],xyzfixed[i]); //apply stdmats

//find atom pairs

  //printf("\nCross-Reference Table of Matrix Indices\n  std order  : input file order\n");
   for (i=0;i<pnum;i++){
         testdist=1000000.;
     for (j=0;j<pnum;j++){
         dist = distXYZ(xyzrefin[j],xyzfixed[i]);
         if (dist < testdist) {
             xyzIndex[i]=j;   //index for refine atom that belongs to fixed atom i
             testdist=dist; }}
              }

//double check pairing, must be 1:1 correspondence
    for (i=0;i<pnum;i++)
    for (j=i;j<pnum;j++){
      if(i != j && xyzIndex[i] == xyzIndex[j]) {
        printf("\n\n");
        printf("***                   SEVERE WARNING                               ***\n");
        printf("*** cannot pair up point constellations, no refinement will be done***\n");
        printf("***         PROBLEM WITH INPUT TRANSFORMATIONS                     ***\n");
        printf("***              OR SYMMETRY ASSIGNMENT                            ***\n");
        return 0;}}

//arrange moving constellation in same order as fixed constellation
    for (i=0;i<pnum;i++) for (j=0;j<4;j++) { xyztmp[i][j] = xyzrefin[xyzIndex[i]][j];}


// calculate com of the moving constellation (fixed constellation at origin)
 for(j=0;j<3;j++){ movcom[j]=0; }

 for (i=0;i<pnum;i++) for (j=0;j<3;j++)movcom[j] += xyztmp[i][j]; 

 for(j=0;j<3;j++) movcom[j]=movcom[j]/pnum;

//translation part of fit:
   multMat(identmat,identmat,reftrans);
   for (j=0;j<3;j++) { reftrans[j][3] = 0.0 - movcom[j];}
    //printf("\n  translation change: %18.10f %18.10f %18.10f\n",reftrans[0][3],reftrans[1][3],reftrans[2][3]);

//move moving constellation to origin
 for (i=0;i<pnum;i++) {applyMat(reftrans,xyztmp[i],matv); for (j=0;j<4;j++)xyztmp[i][j]=matv[j];}

//calculate umat needed for qikfit call
    for(i=0;i<3;i++) {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<pnum;j++) {
            switch(i) {
               case 0:
                  umat[i][0] += xyzfixed[j][0] * xyztmp[j][0];
                  umat[i][1] += xyzfixed[j][0] * xyztmp[j][1];
                  umat[i][2] += xyzfixed[j][0] * xyztmp[j][2];
                  break;
               case 1:
                  umat[i][0] += xyzfixed[j][1] * xyztmp[j][0];
                  umat[i][1] += xyzfixed[j][1] * xyztmp[j][1];
                  umat[i][2] += xyzfixed[j][1] * xyztmp[j][2];
                  break;
               case 2:
                  umat[i][0] += xyzfixed[j][2] * xyztmp[j][0];
                  umat[i][1] += xyzfixed[j][2] * xyztmp[j][1];
                  umat[i][2] += xyzfixed[j][2] * xyztmp[j][2];
                  break;
            } } }
   
   qikfit(umat,rot33,1);   //Andrew Martin's coordinate fitting routine (from profit)


    //printf ("\n");
 //for(i=0;i<3;i++){
    //printf ("     rotation change: %18.10f %18.10f %18.10f\n", rot33[i][0], rot33[i][1], rot33[i][2]); }
    printf ("\n");

//move rot33  to 4x4 mat
  multMat(identmat,identmat,refrot);
  for (i=0; i<3;i++){
    for (j=0;j<3;j++){
       refrot[i][j] = rot33[i][j];}}

//translate moving constellation back to original position 
  det=invMat(reftrans,reftinv);
for (i=0;i<pnum;i++){
  applyMat(reftinv,xyztmp[i],matv);   for (j=0;j<4;j++)  xyztmp[i][j]=matv[j];
   }

//initial rmsd
   dsum1=0;
   for (i=0;i<pnum;i++){dist= distXYZ(xyztmp[i],xyzfixed[i]); dsum1+=dist*dist;}
   dsum1= sqrt(dsum1/pnum);

//apply rot/trans  xyz(fitted) = [rot][inv-movcom]*xyz(start)
  for (i=0;i<3;i++) movcom[i] = -1.0 * movcom[i];
for (i=0;i<pnum;i++){
  applyMat(reftrans,xyztmp[i],matv); for (j=0;j<3;j++) xyztmp[i][j]=matv[j];
  applyMat(refrot,  xyztmp[i],matv); for (j=0;j<3;j++) xyztmp[i][j]=matv[j];
          }

//new rmsd
   dsum2=0;
   for (i=0;i<pnum;i++){ dist = distXYZ(xyztmp[i],xyzfixed[i]);
            dsum2+=dist*dist;}
   dsum2= sqrt(dsum2/pnum);

 printf("rmsd (Angstrom) before fit %10.6f  after fit  %10.6f\n", dsum1,dsum2);
 if (dsum2> dsum1) {
      printf("failed to improve fit\n");
      return 0;}
 if (dsum2> 2.0) {
      printf("\n***HIGH RMSD WARNING: input matrices deviate from ideal point geometry***\n");}

//calculate new p2i
  multMat(reftrans,p2i,mat44);
  multMat(refrot,mat44,p2i);

//move back to rot, tr mats
 for (i=0;i<3;i++) tr[i] = p2i[i][3];
 for (i=0;i<3;i++) for (j=0;j<3;j++) rot[i][j] = p2i[i][j];

 return 1;

}
