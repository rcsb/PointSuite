// pointsuite 0.5.8 
// author Cathy Lawson
// RCSB-PDB
// JUNE 2007

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXMAT 3000 //needed for mat44lib.c
#define MAXCHAR 80
#define MAXCHAIN 5
#include "readcif.c"
#include "mat44lib.c"
#define MAKECIF 1
//fabs(rotation elements) < PR8 are reset to zero
#define PR8 0.00000001
//fabs(translation elements) < PR5 are reset to zero
#define PR5 0.00001
//recompile with bigger MAXMAT if needed

int main(int argc, char * argv[]){
  int i, j, k;
  int ierr;
  double cr, ar, arot, arise;
  double copyrot, copyrise;
  int divisor, ident;
  int zsym, xsym, pnum, numbio, xznum, rewind;
  int exp[4]; 
  double nmat[4][4][4];  //4 source matrices (4x4) for symmetry generation
  double pmat[MAXMAT][4][4]; 
  double mmat[MAXMAT][4][4];
  char   idc[MAXMAT][MAXCHAR];
  int    maxmat;
  double skw2pt[4][4]; 
  double pt2skw[4][4]; 
  double  ptemp[4][4],ptemp2[4][4]; 
  double skwdet;
  int stype, mtype;
  char cifstring[80];
  FILE *inmatfile;
  FILE *outmatfile;
  char *matfile;
//following needed for cif out
  char mau[4];
  char mattypestring[30];
  char matframestring[30];
  FILE *outciffile;
  char  pdbid[20],  chain[100][2], entity[100][2], details[4][100][70];
  int chainum;

//defaults
ident = 1;
pnum = 1;
numbio = 1;
divisor = 1;
arot =  0.0; copyrot = 0.0;
arise = 0.0 ;copyrise = 0.0;
xsym = 1;
zsym = 1;
strcpy (mattypestring, "\'point symmetry operation\'");
strcpy (matframestring,"\'transform to point frame\'");
strcpy (mau,"PAU");
mtype = 'P';

  if(argc < 1 || argc > 2) {
    printf("USAGE: pointmats file.cif\n");
    return 0;
  }

strncpy(pdbid,"TEMP\0",20);
ierr=readOneCif(argv[1],"_pdbx_point_symmetry.entry_id",cifstring);
if (!ierr) {strncpy(pdbid,cifstring,20); pdbid[19]='\0';}
ierr=readOneCif(argv[1],"_pdbx_helical_symmetry.entry_id",cifstring);
if (!ierr) {strncpy(pdbid,cifstring,20); pdbid[19]='\0';}
                                                                                                                             

//printf(" Symmetry type (choices H, C, D, T, O, I): \n");
ierr=readOneCif(argv[1],"_pdbx_point_symmetry.Schoenflies_symbol",cifstring);
stype = cifstring[0];
if (ierr==0) { printf("\nSymmetry type: %c\n", stype);}

//assign parameters
 switch (stype) {

  case 'C' : //need only rot
     printf("\n Rotation symmetry along z-axis (integer >=1) :   ");
     // **get following param from cif
     ierr=readOneCif(argv[1],"_pdbx_point_symmetry.circular_symmetry",cifstring);
     zsym = atoi (cifstring);
     //scanf("%d",&zsym);
     if (zsym < 1) zsym = 1;
     arot=0.0;
     arise=0.0;
     xsym=1;
     numbio=zsym;
     xznum = abs(zsym*xsym);    
     printf("\n%i Matrices to be written for point group C%i\n",numbio,zsym);
     break;

  case 'D' : //need only rot
     printf("\n  Rotation symmetry along z-axis (integer >=1) ");
     ierr=readOneCif(argv[1],"_pdbx_point_symmetry.circular_symmetry",cifstring);
     zsym = atoi (cifstring);
     if (zsym < 1) zsym = 1;
     arot=0.0;
     arise=0.0;
     xsym=2;
     numbio=zsym*2;
     xznum = abs(zsym*xsym);    
     printf("\n%i Matrices to be written for point group D%i\n",numbio,zsym);
     break;

 case 'T' : //no input needed
     printf("\n12 Matrices to be written for point group T\n");
     break;

 case 'I' : //no input needed
     printf("\n60 Matrices to be written for point group I\n");
     break;

 case 'O' : //no input needed
     printf("\n24 Matrices to be written for point group O\n");
     break;


  default:
     //check for helical symmetry
     readOneCif(argv[1],"_pdbx_helical_symmetry.rotation_per_n_subunits", cifstring);

     if (strncmp(cifstring,"\0",80)==0) {
     printf("undefined symmetry, writing out identity mat\n"); }

     else {
     stype = 'H';
     printf("\nSymmetry type: %c\n", stype);
     arot = atof (cifstring);
     ierr=readOneCif(argv[1],"_pdbx_helical_symmetry.rise_per_n_subunits", cifstring);
     arise = atof (cifstring);
     ierr=ierr+readOneCif(argv[1],"_pdbx_helical_symmetry.n_subunits_divisor", cifstring);
     divisor = atoi (cifstring);
     ierr=ierr+readOneCif(argv[1],"_pdbx_helical_symmetry.number_of_operations", cifstring);
     numbio = atoi (cifstring);
     ierr=ierr+readOneCif(argv[1],"_pdbx_helical_symmetry.circular_symmetry", cifstring);
     zsym = atoi (cifstring);
     ierr=ierr+readOneCif(argv[1],"_pdbx_helical_symmetry.dyad_axis", cifstring);
     if (ierr > 0) {printf("incomplete helical parameter definition\n"); return 1;}
     if (strncmp(cifstring,"yes",3)==0) xsym = 2;
     if (strncmp(cifstring,"no",2)==0) xsym = 1;

     if (numbio > MAXMAT) numbio=MAXMAT;
     if (divisor == 0) divisor=1;
     copyrot = arot;
     copyrise = arise;
     arot = arot / divisor;
     arise = arise / divisor;
     if (xsym < 1 || xsym > 2) xsym = 1;
     if (zsym < 1) zsym = 1;
     xznum = abs(zsym*xsym);    //how many mats per x/z symmetry
     pnum = (numbio/xznum)*xznum;  //pnum must be multiple of xznum
     mtype = stype;
     strcpy (mau,"HAU");
     strcpy (mattypestring, "\'helical symmetry operation\'");
     strcpy (matframestring,"\'transform to helical frame\'");

     printf ("\n");
     if (arot*arise < 0)printf("%i Matrices requested for a left-handed helix\n",numbio);
     if (arot*arise == 0)printf("%i Matrices requested for an untwisted fiber or disc\n",numbio);
     if (arot*arise > 0)printf("%i Matrices requested for a right-handed helix\n",numbio);
     printf("   %i matrices will be written\n", pnum);
     rewind =   pnum/(2*xznum); //helix will be wound back by rewind arot/arise units 
     ident = rewind*xznum+1;
     printf("   %i =identity matrix index\n", ident);
     printf("   %10.3f degree rotation per subunit around z\n", arot);
     printf("   %10.3f Angstrom rise per subunit along z\n", arise);
     printf("   %i-fold about the x axis\n", xsym);
     printf("   %i-fold about the z axis\n", zsym);
     }

     break;

   }

     printf ("\n");



//for all paths, initialize symmetry generating matrices to identity, exponents to 0
 for (i=0;i<4;i++){
  nmat[i][0][0] = 1.0; nmat[i][0][1] = 0.0; nmat[i][0][2] = 0.0; nmat[i][0][3] = 0.0;
  nmat[i][1][0] = 0.0; nmat[i][1][1] = 1.0; nmat[i][1][2] = 0.0; nmat[i][1][3] = 0.0;
  nmat[i][2][0] = 0.0; nmat[i][2][1] = 0.0; nmat[i][2][2] = 1.0; nmat[i][2][3] = 0.0;
  nmat[i][3][0] = 0.0; nmat[i][3][1] = 0.0; nmat[i][3][2] = 0.0; nmat[i][3][3] = 1.0;
  exp[i]=0;}

//BEGIN HELICAL PATH
//H,C,D symmetry gen 
if (stype == 'H' || stype == 'C' || stype == 'D') {

//apply C symmetry about z-axis first, then symmetry about x, then rot/trans mat derived from arot, arise

//calculate nmat[0] from zsym
cr = zsym;
cr = (360.0/cr) * (M_PI/180.0);
nmat[0][0][0] = cos(cr); nmat[0][0][1] = -1.0*sin(cr);    nmat[0][0][2]=0.0;  nmat[0][0][3]=0.0;
nmat[0][1][0] = sin(cr); nmat[0][1][1] = cos(cr);         nmat[0][1][2]=0.0;  nmat[0][1][3]=0.0;
nmat[0][2][0] = 0.0;     nmat[0][2][1] = 0.0;             nmat[0][2][2]=1.0;  nmat[0][2][3]=0.0;
nmat[0][3][0] = 0.0;     nmat[0][3][1] = 0.0;             nmat[0][3][2]=0.0;  nmat[0][3][3]=1.0;
exp[0] = abs(zsym) - 1;

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
exp[2] = (numbio/xznum) - 1 ;

} 
//END HELICAL PATH



//BEGIN TETRAHEDRAL PATH: T,O symmetry gen  
if (stype == 'T' || stype == 'O') {

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
  if (stype == 'O') exp[3] = 1;
  if (stype == 'T') exp[3] = 0;   

} 
//END TETRAHEDRAL PATH

//BEGIN ICOSAHEDRAL PATH: I symmetry gen  
if (stype == 'I') {

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


//total mats
    pnum = (exp[3]+1) * (exp[2]+1) * (exp[1]+1) * (exp[0]+1);
    //printf ("pnum: %i\n", pnum);

//CALCULATE THE MATRICES 
    calcSymmMats (exp, nmat, pmat); 




//IF HELICAL, put the identity mat in middle of the matset 

  if (stype == 'H') {

   //calculate "rewind" matrix, using  ar and arise as defined in prior HELICAL section 
   rewind = -1.0 * rewind;
   ptemp[0][0] = cos(ar*rewind); ptemp[0][1] = -1.0*sin(ar*rewind);    ptemp[0][2]=0.0;  ptemp[0][3]=0.0;
   ptemp[1][0] = sin(ar*rewind); ptemp[1][1] = cos(ar*rewind);         ptemp[1][2]=0.0;  ptemp[1][3]=0.0;
   ptemp[2][0] = 0.0;            ptemp[2][1] = 0.0;                    ptemp[2][2]=1.0;  ptemp[2][3]=arise*rewind;
   ptemp[3][0] = 0.0;            ptemp[3][1] = 0.0;                    ptemp[3][2]=0.0;  ptemp[3][3]=1.0;

   //rewind the matrix set
   for (i=0;i<pnum;i++) {
       multMat(ptemp,pmat[i],ptemp2); 
          for (j=0;j<4;j++) for (k=0;k<4;k++) {pmat[i][j][k]=ptemp2[j][k]; }}

       } //end rewind HELICAL mats



//apply skew matrix  if present if CIF

//begin skew option

//initialize skw2pt and pt2skw mats to identity in any case
 for (i=0;i<4;i++)
   for (j=0;j<4;j++) { 
     if (i != j) {
     skw2pt[i][j]= 0.0;
     pt2skw[i][j]= 0.0;}
     if (i == j) {
     skw2pt[i][j]= 1.0;
     pt2skw[i][j]= 1.0;} }

//**look for P matrix in input cif
ierr = readAssemblyCif(argv[1], idc, mmat,&maxmat);
for (i=0;i<maxmat;i++)
    {if ((!strncmp("P",idc[i],MAXCHAR)) ||(!strncmp("H",idc[i],MAXCHAR)) ) break;}
if (ierr>0 || i==maxmat) {printf("\nsymmetry operations will be provided in the standard frame\n");}

else {

    for (j=0;j<4;j++)
       for (k=0;k<4;k++){ skw2pt[j][k]=mmat[i][j][k];   }
    skwdet=invMat(skw2pt,pt2skw);


    printf ("symmetry matrices will be transformed by [P-inv][StdMats][P] \n");
    printf ("[P] =  \n");
       for (k=0;k<3;k++) {
           printf("%10.6f%10.6f%10.6f%15.5f\n",
           skw2pt[k][0], skw2pt[k][1], skw2pt[k][2], skw2pt[k][3]); }

  //apply the skew
   // [pt2skw][mat][skw2pt],
    for (i=0;i<pnum; i++) {
       multMat(pt2skw,pmat[i],ptemp);
       multMat(ptemp,skw2pt,pmat[i]);
       }

   }  //end skew option



  //print out BIOMT to biomats
  printf("output file: pointmats.biomt \n");
  outmatfile = fopen("pointmats.biomt","w");
    for (i=0; i<pnum; i++) {
       for (k=0;k<3;k++) {
           fprintf(outmatfile,"REMARK 350   BIOMT%1i%4i%10.6f%10.6f%10.6f%15.5f\n",
           k+1, i+1, pmat[i][k][0], pmat[i][k][1], pmat[i][k][2], pmat[i][k][3]);}}
    fclose(outmatfile);

  //print out cif
  printf("output file: pointmats.cif \n");
  outciffile = fopen("pointmats.cif","w");

//helical symmetry tokens
if (mtype=='H'){   
   fprintf(outciffile,"_pdbx_helical_symmetry.entry_id                    %s\n",pdbid);
   fprintf(outciffile,"_pdbx_helical_symmetry.number_of_operations        %i\n",pnum );
   fprintf(outciffile,"_pdbx_helical_symmetry.rotation_per_n_subunits     %f\n",copyrot);
   fprintf(outciffile,"_pdbx_helical_symmetry.rise_per_n_subunits         %f\n",copyrise);
   fprintf(outciffile,"_pdbx_helical_symmetry.n_subunits_divisor          %i\n",divisor);
   if(xsym==1)fprintf(outciffile,"_pdbx_helical_symmetry.dyad_axis                   no\n"); 
   if(xsym==2)fprintf(outciffile,"_pdbx_helical_symmetry.dyad_axis                   yes\n"); 
   fprintf(outciffile,"_pdbx_helical_symmetry.circular_symmetry           %i\n",zsym);
           }

//point symmetry tokens
if (mtype=='P'){
   fprintf(outciffile,"_pdbx_point_symmetry.entry_id                 %s\n",pdbid);
   fprintf(outciffile,"_pdbx_point_symmetry.Schoenflies_symbol       %c\n",stype);
   if(stype =='C' || stype == 'D') 
  {fprintf(outciffile,"_pdbx_point_symmetry.circular_symmetry        %i\n",zsym); } }


 //matrices
  fprintf(outciffile,"#\nloop_\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.id\n");              
  fprintf(outciffile,"_pdbx_struct_oper_list.type \n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[1][1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[1][2]\n"); 
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[1][3]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.vector[1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[2][1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[2][2]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[2][3] \n");
  fprintf(outciffile,"_pdbx_struct_oper_list.vector[2]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[3][1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[3][2] \n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[3][3]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.vector[3]\n");

  //write transform to helical or point frame matrix
  fprintf(outciffile,"   %c    %30s\n", mtype, matframestring);  
  for (j=0;j<3;j++){
  fprintf(outciffile,"    %14.8f%14.8f%14.8f%20.5f\n", 
            checkZero(skw2pt[j][0], PR8),
            checkZero(skw2pt[j][1], PR8),
            checkZero(skw2pt[j][2], PR8),
            checkZero(skw2pt[j][3], PR5));}

  //write assembly mats  
       for (i=0;i<pnum;i++){
           fprintf(outciffile,"%4i %30s\n",i+1, mattypestring);  
       for (k=0;k<3;k++) {
           fprintf(outciffile,"    %14.8f%14.8f%14.8f%20.5f\n",
              checkZero(pmat[i][k][0], PR8),
              checkZero(pmat[i][k][1], PR8),
              checkZero(pmat[i][k][2], PR8),
              checkZero(pmat[i][k][3], PR5));}}

  fclose(outciffile);
  return 1;
}
