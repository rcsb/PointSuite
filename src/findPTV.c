//subroutines for findframe
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define VERBOSE 0
#define PRECISION 0.20

double findDist(double v1[3], double v2[3]){
double dist;
dist = sqrt(  (v1[0]-v2[0])*(v1[0]-v2[0]) +
              (v1[1]-v2[1])*(v1[1]-v2[1]) +
              (v1[2]-v2[2])*(v1[2]-v2[2]) );
return dist;
}

double findMag(double v1[3]){
double dist;
dist = sqrt(  (v1[0])*(v1[0]) +
              (v1[1])*(v1[1]) +
              (v1[2])*(v1[2]) );
return dist;
}

int findPTV(int nmat, double ptrans[MAXMAT][3], double prot[MAXMAT][3][3], 
            double refAtom[3], double rotation[3][3], double translation[3],
            int *pstype, int *pzsym){
  int i = 0; int j = 0;
  int imaxfold, itwofold;
  double x,y,z, angle, angle2, angle3, div,divrot, test3fold;
  double mindist, mindistnext;
  double refAtnorm[3];
  double tvector[3], testvec[3];
  double transadjust;
  double t[3][3];
  double t1[3][3];
  double tf[3][3]; //final t matrix
  double * refpt;
  double taxes[3];
  FILE *outmatfile;

  //indices for closest axes
  int p2fold = 0;
  int pxfold = 0; 
  int pxfold2 = 0;
  int pnfold = 0; 

  double pdist[MAXMAT];
  double pvectors[MAXMAT][3];
  double pangles[MAXMAT];
  double pfoldf[MAXMAT];
  double pcoaxtrans[MAXMAT];
  int    pfoldi[MAXMAT];
  int    pindex[MAXMAT];

  double vvectors[4][3];

  double pplane[3];
  double vplane[3];

  double cross[3];
  double cross2[3];
  double dot;

  //find center of translation vectors
  i = 0;
  tvector[0] = 0; tvector[1] = 0; tvector[2] = 0;
  while(i<nmat){
    tvector[0] = tvector[0] + ptrans[i][0];
    tvector[1] = tvector[1] + ptrans[i][1];
    tvector[2] = tvector[2] + ptrans[i][2];
    i++; }
  tvector[0] = -1*tvector[0]/nmat;
  tvector[1] = -1*tvector[1]/nmat;
  tvector[2] = -1*tvector[2]/nmat;

  if(VERBOSE) printf("negated center of translation vectors: %8.3f %8.3f %8.3f\n", 
            tvector[0], tvector[1], tvector[2]);

  //translate refAtom position to origin and normalize
  refAtnorm[0] = refAtom[0] + tvector[0]; 
  refAtnorm[1] = refAtom[1] + tvector[1]; 
  refAtnorm[2] = refAtom[2] + tvector[2]; 
  normalize(refAtnorm);

 
  //process the transformations
  i = 0;
  while(i<nmat){
  
  //get rotation in angle-axis format
    decon(prot[i], &x, &y, &z, &angle);
    pvectors[i][0] = x;
    pvectors[i][1] = y;
    pvectors[i][2] = z;
    pangles[i] = angle;
   
   if (VERBOSE) printf ("pvec:  %f %f %f  angle: %f\n", x,y,z,angle);

  //find translation vector component parallel to rotation axis 
    pcoaxtrans[i]=(dotProduct(ptrans[i],pvectors[i]));

  //Find shortest axis dist from norm reference atom, switch axis direction if needed
      pdist[i] = findDist (pvectors[i], refAtnorm);
              testvec[0] = -1*pvectors[i][0];        
              testvec[1] = -1*pvectors[i][1];        
              testvec[2] = -1*pvectors[i][2];       
      div   =   findDist (testvec, refAtnorm); 
        if (div < pdist[i])  {pdist[i] = div;
             pcoaxtrans[i] = -1*pcoaxtrans[i];
             pvectors[i][0] = testvec[0];
             pvectors[i][1] = testvec[1];
             pvectors[i][2] = testvec[2];} 

    if (pangles[i]>0) pfoldf[i] = 360.0/radToDeg(angle);
    else pfoldf[i]=0;
    pfoldi[i]=0;
    pindex[i]=-1;
   if (VERBOSE) printf ("mat:  %i  pdist: %f\n", i+1, pdist[i]);

    i++;
  }  
//end processing of transformations

// define sort index based on pdist
   j=0;
   mindist = 0.0;
   while (j<nmat){
      mindistnext = 1000.0;
      for (i=0;i<nmat;i++){
              if (pdist[i]==mindist) {pindex[j] =i; j++;}
              else if (pdist[i]>mindist && pdist[i]<mindistnext) mindistnext=pdist[i];}
      mindist=mindistnext;
      }


 //find point symmetry type
  *pstype='U'; 
  i = 0; if (nmat%2 == 0 && nmat>2) 
         while (i<nmat) {if (fabs(pfoldf[i]-(nmat/2)) < PRECISION) *pstype = 'D';i++;} 
  i = 0; while (i<nmat) {if (fabs(pfoldf[i]-nmat)     < PRECISION) *pstype = 'C';i++;} 
  i = 0; if (*pstype=='U' && nmat ==12) 
         while (i<nmat) {if (fabs(pfoldf[i]-3.0)      < PRECISION) *pstype = 'T';i++;} 
  i = 0; if (*pstype=='U' && nmat ==24) 
         while (i<nmat) {if (fabs(pfoldf[i]-4.0)      < PRECISION) *pstype = 'O';i++;} 
  i = 0; if (*pstype=='U' && nmat ==60) 
         while (i<nmat) {if (fabs(pfoldf[i]-5.0)      < PRECISION) *pstype = 'I';i++;} 


//check for helical matrices
//simplistic test (any coax trans in set > 20)
    i=0;
    while (i<nmat) {
      if (fabs(pcoaxtrans[i])> 20.0) {*pstype = 'H'; break; }
      i++;}


//symmetry-type dependent assignments
 switch (*pstype) {
  case 'U': 
    printf ("ERROR: could not determine symmetry, check matrices\n"); 
    return 0; 

  case 'H': 
    //only process helical matrices with helical axis ON z
    // for each transformation
    // pattern must be cos(n) -sin(n)   0  0
    //                 sin(n)  cos(n)   0  0
    //                   0       0      1  trans
    //  (or above mat with 2 fold sym along x applied)
    for (i=0;i<nmat;i++){
      if (fabs(prot[i][0][2]) > PRECISION || fabs(prot[i][1][2]) > PRECISION || 
         fabs(prot[i][2][0]) > PRECISION || fabs(prot[i][2][1]) > PRECISION || 
         (fabs(prot[i][2][2])-1.0) > PRECISION || fabs(ptrans[i][0])> PRECISION ||
         (fabs (prot[i][0][1])-fabs(prot[i][1][0])) > PRECISION ||
         (fabs (prot[i][0][0]) -fabs(prot[i][1][1])) > PRECISION ||
         fabs(ptrans[i][1])> PRECISION ) {
        printf ("ERROR: non-standard helical matrices?\n");
        return 0;}
      }
    printf ("Helical Structure Assumed\n"); 

    //look for exact n-folds, no translations
    imaxfold = 1;
    itwofold = 1;
     for (i=0;i<nmat;i++) {
     if (fabs(pcoaxtrans[i]) > PRECISION) continue;
     for (j=2;j<=nmat;j++) {
      if (fabs(pfoldf[i]-j) < PRECISION) { 
             if (j==2) itwofold=2;
             if (j>imaxfold) imaxfold=j; 
             pfoldi[i] = j;  }}}


    //find "basis matrix"(smallest rotation/translation)
    div = 1000.0; 
    for(i=0;i<nmat;i++) {
      if (fabs(pcoaxtrans[i]) < PRECISION/10.0) continue;
      if (fabs(pcoaxtrans[i]) < div) div = fabs(pcoaxtrans[i]); }

    pnfold = -1;
    divrot =1000.0;
    for(i=0;i<nmat;i++) {
      if (fabs(pcoaxtrans[i]) < PRECISION/10.0) continue;
      if (((fabs(pcoaxtrans[i])- div) < PRECISION/10.0) && (pangles[i] < divrot)) 
            {divrot = pangles[i]; pnfold = i;}}

    if (pnfold == -1){
        printf ("ERROR: cannot determine helix rotation/translation\n");
        return 0;}

    //check hand
    angle = radToDeg(pangles[pnfold]);
    if (prot[pnfold][1][0] < 0) angle = -1.0 * angle; 
    div   = ptrans[pnfold][2];
    if ((angle < 0 && div < 0) || (angle > 0 && div < 0))
       {angle = -1.0 * angle; div = -1.0 * div;}

    //check for 2-fold perpendicular to helix axis 
    if (itwofold == 2) {
      itwofold = 1;
      i=0;
      while (i<nmat) { if (pfoldi[i]==2){
         if (fabs(dotProduct(pvectors[i],pvectors[pnfold])) < PRECISION) {itwofold = 2; break;}}
      i++;}}

    i =  pnfold;
    printf ("\n Basis Matrix: \n");
    for (j=0;j<3;j++) printf (" %10.6f %10.6f %10.6f %10.6f\n", 
           prot[i][j][0],prot[i][j][1], prot[i][j][2], ptrans[i][j]);
    printf ("\n");

    //check rot/trans against basis matrix
    printf("rotation   : %10.3f  degrees\n", angle);
    printf("translation: %10.3f  Angstroms\n", div);
    printf("dyad symm  :      %i\n", itwofold);
    printf("z    symm  :      %i\n\n", imaxfold);
   
/*
    printf ("\nHelix parameter string (#mats rot trans n=1 dyad zsym):\n\n"); 
    printf ("H %i %8.3f %8.3f 1 %i %i\n", nmat,angle,div, itwofold, imaxfold); 
*/

//print out helical parameters to file
  printf("helical symmetry info written to output file *findframe.cif*\n\n");
  outmatfile = fopen("findframe.cif","w");
   //fprintf(outmatfile,  "H %i %8.4f %8.4f 1 %i %i\n", nmat,angle,div, itwofold, imaxfold); 
   fprintf(outmatfile,"_pdbx_helical_symmetry.entry_id                    TEMP\n");
   fprintf(outmatfile,"_pdbx_helical_symmetry.number_of_operations        %i\n",nmat );
   fprintf(outmatfile,"_pdbx_helical_symmetry.rotation_per_n_subunits     %8.3f\n",angle);
   fprintf(outmatfile,"_pdbx_helical_symmetry.rise_per_n_subunits         %8.3f\n",div  );
   fprintf(outmatfile,"_pdbx_helical_symmetry.n_subunits_divisor          %i\n",1    );
   if(itwofold==1)fprintf(outmatfile,"_pdbx_helical_symmetry.dyad_axis                   no\n");
   if(itwofold==2)fprintf(outmatfile,"_pdbx_helical_symmetry.dyad_axis                   yes\n");
   fprintf(outmatfile,"_pdbx_helical_symmetry.circular_symmetry           %i\n",imaxfold);

    fclose(outmatfile);
                                                                                                                  

    return 0; 

  case 'C' : 
    printf ("circular symmetry:  C %i\n", nmat);
    printf ("%i-fold will be aligned to the z-axis\n", nmat);
    printf ("reference atom will be placed on the x-axis\n");
    *pzsym = nmat;

    //define n-folds
    i=0;
    while (i<nmat) {
     for (j=1;j<=nmat/2;j++) {
      div = j; 
      div = nmat/div;
      if (fabs(pfoldf[i]-div)< PRECISION) {pfoldi[i] = nmat; break; }}
      i++;}

    //to move to standard frame:
    //(a) align normal of n-fold-refAtnorm plane to the normal of the Z-X plane
    //(b) rotate about normal to align n-fold to Z

    //std frame vector 0 =   X axis
    vvectors[0][0] = 1.0;
    vvectors[0][1] = 0.0;
    vvectors[0][2] = 0.0;
    normalize(vvectors[0]);

    //std frame vector 1 =  n fold on Z
    vvectors[1][0] = 0.0;
    vvectors[1][1] = 0.0;
    vvectors[1][2] = 1.0;
    normalize(vvectors[1]);

    //find n-fold axis 
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==nmat) { pnfold = pindex[i]; break; } i++; }
    //copy refAtnorm to pvector[nmat] , index nmat is used as the pxfold
    pvectors[nmat][0]=refAtnorm[0];
    pvectors[nmat][1]=refAtnorm[1];
    pvectors[nmat][2]=refAtnorm[2];
    pxfold = nmat;

    //adjust the translation vector so that the refAtom ends up on the x axis
    testvec[0] = refAtom[0] + tvector[0];
    testvec[1] = refAtom[1] + tvector[1];
    testvec[2] = refAtom[2] + tvector[2];
 
    transadjust=(dotProduct(pvectors[pnfold],testvec)); 
    tvector[0] = tvector[0] - transadjust*pvectors[pnfold][0];
    tvector[1] = tvector[1] - transadjust*pvectors[pnfold][1];
    tvector[2] = tvector[2] - transadjust*pvectors[pnfold][2];

    break;

  case 'D' : 
    printf ("dihedral symmetry:  D %i\n", nmat/2);
    *pzsym = nmat/2;
    printf ("%i-fold will be aligned to the z-axis\n", *pzsym);
    printf ("2-fold closest to reference atom will be placed on the x-axis\n");

    //define n-folds
    i=0;
    while (i<nmat) {
     for (j=1;j<=nmat/2;j++) {
      div = j; div = nmat/div;
      if (fabs(pfoldf[i]-div) < PRECISION) {pfoldi[i] = nmat/2; break;}
      if (fabs(pfoldf[i]-2)   < PRECISION) {pfoldi[i] = 2; break;}}
      i++;}

    //strategy:
    //(a) plane alignment: n-fold X (closest 2-fold || to n-fold plane) : Z-X 
    //(b) axis alignment:  n-fold :  Z

    //std frame vector 0 =  2 fold on x
    vvectors[0][0] = 1.0;
    vvectors[0][1] = 0.0;
    vvectors[0][2] = 0.0;
    normalize(vvectors[0]);

    //std frame vector 1 =  n fold on Z
    vvectors[1][0] = 0.0;
    vvectors[1][1] = 0.0;
    vvectors[1][2] = 1.0;
    normalize(vvectors[1]);

    //find n-fold and x-fold axes 
     i = 0;
     while(i<nmat){ if(pfoldi[pindex[i]]==nmat/2) { pnfold = pindex[i]; break; } i++; }
     i = 0;
     while(i<nmat){ if(pfoldi[pindex[i]]==2 && 
           findDist(pvectors[pnfold],pvectors[pindex[i]]) > 1.0) { pxfold = pindex[i]; break; } i++; }

    break;

  case 'T' : 
    printf ("tetrahedral symmetry: T\n");
    printf ("3-fold  closest to reference atom will be aligned to the body diagonal (1,1,1)\n");
    printf ("2-folds closest to reference atom will be aligned to the  x- and z-axes\n");
    //define n-folds
    i=0;
    while (i<nmat) {
      if (fabs(pfoldf[i]-2)   < PRECISION) pfoldi[i] = 2;
      if (fabs(pfoldf[i]-3)   < PRECISION) pfoldi[i] = 3;
      i++;}
    //strategy:
    //(a) plane alignment:  closest-3fold x closest 2-fold : xyz diag x X 
    //(b) axis alignment:  3-fold : xyzdiag 

    //std frame vector 0 =  2 fold on x
    vvectors[0][0] = 1.0;
    vvectors[0][1] = 0.0;
    vvectors[0][2] = 0.0;
    normalize(vvectors[0]);

    //std frame vector 1 =  3 fold on xyz 
    vvectors[1][0] = 1.0;
    vvectors[1][1] = 1.0;
    vvectors[1][2] = 1.0;
    normalize(vvectors[1]);

    //find n-fold and x-fold axes 
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==3) { pnfold = pindex[i]; break; } i++; }
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==2) { pxfold = pindex[i]; break; } i++; }
    while(i<nmat){ 
       if(pfoldi[pindex[i]]==2 
            && findDist(pvectors[pxfold],pvectors[pindex[i]]) > 0.02) 
            { pxfold2 = pindex[i]; break; } i++; }
    printf ("    closest 2-fold index: %3i \n", pxfold);
    printf ("2nd closest 2-fold index: %3i \n", pxfold2);
  
     crossProduct(pvectors[pxfold],pvectors[pxfold2],cross);
     normalize(cross);
     angle3=fabs(acos(dotProduct(cross,pvectors[pnfold])));  
     if (angle3 < 0.5*M_PI) pxfold = pxfold2;  //ensures one 2 fold on z, other on x
    break;

  case 'O' : 
    printf (" octahedral symmetry: O\n");
    printf ("4-fold  closest to reference atom will be aligned to the x-axis\n");
    printf ("3-fold  closest to reference atom will be aligned to the body diagonal (1,1,1)\n");
    //define n-folds
    i=0;
    while (i<nmat) {
      if (fabs(pfoldf[i]-2)   < PRECISION) pfoldi[i] = 2;
      if (fabs(pfoldf[i]-3)   < PRECISION) pfoldi[i] = 3;
      if (fabs(pfoldf[i]-4)   < PRECISION) pfoldi[i] = 4;
      i++;}

    //strategy:
    //(a) plane alignment:  closest-3fold x closest 4-fold : xyz diag x X 
    //(b) axis alignment:  3-fold : xyzdiag 

    //std frame vector 0 =  closest 4 fold on x
    vvectors[0][0] = 1.0;
    vvectors[0][1] = 0.0;
    vvectors[0][2] = 0.0;
    normalize(vvectors[0]);

    //std frame vector 1 =  closest 3 fold on xyz 
    vvectors[1][0] = 1.0;
    vvectors[1][1] = 1.0;
    vvectors[1][2] = 1.0;
    normalize(vvectors[1]);

    //find n-fold and x-fold axes 
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==3) { pnfold = pindex[i]; break; } i++; }
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==4) { pxfold = pindex[i]; break; } i++; }

    break;

  case 'I' : 
    printf ("icosahedral symmetry: I\n");
    //define n-folds
    i=0;
    while (i<nmat) {
      if (fabs(pfoldf[i]-2)   < PRECISION) pfoldi[i] = 2;
      if (fabs(pfoldf[i]-2.5) < PRECISION) pfoldi[i] = 5;
      if (fabs(pfoldf[i]-3)   < PRECISION) pfoldi[i] = 3;
      if (fabs(pfoldf[i]-5)   < PRECISION) pfoldi[i] = 5;
      i++;}
 
    //strategy:
    //(a) plane alignment: closest 5-fold x closest 3-fold* :  viper 5-fold x 3-fold 
    //(b) axis  alignment: 5-fold : 5-fold
    //*closest 3-fold if refAtom is  nearer to 3-fold; 
    // but 3-fold to right of closest 2-fold if refAtom is nearer to 2-fold

    //std frame vector 0 = viper 3fold: x=phi/3, y=0, z=(2*phi+1)/3
    vvectors[0][0] = (1.0 + sqrt(5.0))/6.0;
    vvectors[0][1] = 0.0;
    vvectors[0][2] = (2.0 + sqrt(5.0))/3.0;
    normalize(vvectors[0]);

    //std frame vector 1 = viper5fold: x=0, y = 1, z = phi
    vvectors[1][0] = 0.0;
    vvectors[1][1] = 1.0;
    vvectors[1][2] = (1.0 + sqrt(5.0))/2.0;
    normalize(vvectors[1]);
    
    //find nearest axes 
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==2) { p2fold = pindex[i]; break; } i++; }
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==5) { pnfold = pindex[i]; break; } i++; }
    i = 0;
    while(i<nmat){ if(pfoldi[pindex[i]]==3) { pxfold = pindex[i]; break; } i++; }
    i++;
    while(i<nmat){ 
       if(pfoldi[pindex[i]]==3 
            && findDist(pvectors[pxfold],pvectors[pindex[i]]) > 0.02) 
            { pxfold2 = pindex[i]; break; } i++; }
    if (VERBOSE) {
    printf ("    closest 2-fold index: %3i \n", pxfold);
    printf ("2nd closest 2-fold index: %3i \n", pxfold2);}

    // if icos a.u. is 3-fold-centric, choose the 1st 3 fold
    // if icos a.u. is 2-fold centric, choose the 3fold that is to the "right" of the two-fold
    // test to see if icos a.u. is 2-fold centric and if so, figure out which 3 fold we want
     crossProduct(pvectors[pnfold],pvectors[pxfold],cross);
     normalize(cross);
     angle3=fabs(acos(dotProduct(cross,refAtnorm)) - 0.5*M_PI);  //close to 0 if 3fold centric
     crossProduct(pvectors[pnfold],pvectors[p2fold],cross);
     normalize(cross);
     angle2=fabs(acos(dotProduct(cross,refAtnorm))- 0.5*M_PI);  //close to 0 if 2fold centric
     if (angle2 < angle3){
       printf ("\n Icosahedral asymmetric unit is 2-fold centric");
       printf ("\n closest 5-fold will be aligned to the vector  ( 0, 1, phi)");
       printf ("\n 3-fold to right of 2-fold will be aligned to  ( phi/3, 0, (2*phi+1)/3)\n");
       crossProduct(pvectors[pnfold],pvectors[p2fold],cross);
       normalize(cross);
       angle=acos(dotProduct(cross,pvectors[pxfold]));
       angle2=acos(dotProduct(cross,pvectors[pxfold2]));
       if (angle2 < angle)pxfold=pxfold2;
       }
     else {printf ("\n Icosahedral asymmetric unit is 3-fold centric");
       printf ("\n closest 5-fold will be aligned to the vector  ( 0, 1, phi)");
       printf ("\n closest 3-fold will be aligned to ( phi/3, 0, (2*phi+1)/3)\n");
       } 
    break;
  } //end switch/case

  //print axis info
  if(VERBOSE){
    printf ("\nn fold index: %3i \n", pnfold);
    printf (  "x fold index: %3i \n", pxfold);

      printf("# (normalized rotation axis), rotang, n-fold, coax-trans, dist-to-ref, closest index\n"); 
    for (j=0;j<nmat;j++){
    i = 0;
    while(i<nmat){  if (pindex[i]==j) {
      printf("%3i (%6.3f %6.3f %6.3f) %8.2f %8i  %8.3f -- %10.4f  %8i\n", 
              pindex[i]+1,
              pvectors[pindex[i]][0], 
              pvectors[pindex[i]][1], 
              pvectors[pindex[i]][2],
              radToDeg(pangles[pindex[i]]), 
              pfoldi[pindex[i]], 
              pcoaxtrans[pindex[i]], 
              pdist[pindex[i]],i+1);}
    i++;
    }}
   }
   if (*pstype == 'H') return 0;
 
    //calculate cross product of vectors defining standard frame
    crossProduct(vvectors[1], vvectors[0], vplane);
    normalize(vplane);


  //calculate normal to plane between nfold and xfold
  crossProduct(pvectors[pnfold], pvectors[pxfold], cross);
  normalize(cross);
  pplane[0] = cross[0];
  pplane[1] = cross[1];
  pplane[2] = cross[2];

  if(VERBOSE){
    printf("PPLANE: (%f  %f  %f)\n", pplane[0], pplane[1], pplane[2]);
    printf("VPLANE: (%f  %f  %f)\n", vplane[0], vplane[1], vplane[2]);
  }
  
  //get angle-axis of rotation between two vectors, and turn it into a rotation matrix
  dot = dotProduct(pplane, vplane);
  if (dot<-1) dot = -1;
  if (dot>1) dot = 1;
  angle = acos(dot);
  crossProduct(pplane, vplane, cross);
  normalize(cross);

  //replace zero crossProduct for rare event angle =180 deg (CLL)
  //with arbitrary vector perpendicular to vplane normal
  if (angle == M_PI) {
  div = 1.0 + ((vplane[0]*vplane[0]) / (vplane[1]*vplane[1]));
  cross[0] = sqrt((2.0 - (vplane[0]*vplane[0]) - (vplane[1]*vplane[1]) -(vplane[2]*vplane[2]))/div);
  cross[1] = (vplane[0]*cross[0])/(-1.0 * vplane[1]);
  cross[2] = 0; 
  }
  
  recon(cross[0], cross[1], cross[2], angle, t);

  //apply rotation to place the nfold and xfold axes in the same plane as for standard frame
  refpt = pvectors[pnfold];
  taxes[0] = t[0][0]*refpt[0] + t[0][1]*refpt[1] + t[0][2]*refpt[2];
  taxes[1] = t[1][0]*refpt[0] + t[1][1]*refpt[1] + t[1][2]*refpt[2];
  taxes[2] = t[2][0]*refpt[0] + t[2][1]*refpt[1] + t[2][2]*refpt[2];
  
  
  //calculate rotation angle to place n-fold on standard frame axis
  dot = dotProduct(taxes, vvectors[1]);
  if (dot<-1) dot = -1;
  if (dot>1) dot = 1;
  angle = acos(dot);
  crossProduct(taxes, vvectors[1], cross);
  normalize(cross);

  //replace zero crossProduct for angle =180 deg with normal to vplane
  if (angle == M_PI) {
  cross[0] = vplane[0];
  cross[1] = vplane[1];
  cross[2] = vplane[2]; }
  recon(cross[0], cross[1], cross[2], angle, t1);
  
  //multiply the matrices
  multiplyMatrices(t1, t, rotation);

  translation[0] = rotation[0][0]*tvector[0] + rotation[0][1]*tvector[1] + rotation[0][2]*tvector[2];
  translation[1] = rotation[1][0]*tvector[0] + rotation[1][1]*tvector[1] + rotation[1][2]*tvector[2];
  translation[2] = rotation[2][0]*tvector[0] + rotation[2][1]*tvector[1] + rotation[2][2]*tvector[2];
  return 1;
}


