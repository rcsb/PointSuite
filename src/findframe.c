// pointsuite 0.5.8 
// author Cathy Lawson
// RCSB-PDB JUNE 2007
// updated JAN 2013 by Huanwang Yang to accept CIF or PDB as input format
// and by Cathy Lawson to accept an optional 2nd BIOMT file with the matrices
//
//findframe is an RCSB-modified version of pdb2viper
//authors Vijay Reddy and Ian Borelli
// http://viperdb.scripps.edu

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXMAT 500 //used in mat44lib.c, findmat.c, refmat.c

#include "mat44lib.c"
#include "angax.c"
#include "findPTV.c"    
#include "qikfit.c"
#include "refmat.c"
#include "cifparse.c"
#define PR8 0.00000001
#define PR5 0.00001

int main(int argc, char * argv[]){
  int stype, zsym;
  int i = 0, nfiles = 1;
  int idn, iident;
  int row, column,nmat, filetype;
  FILE * file, * outmatfile;
  char  buffer[200];
  char * t;
  char pdb_file[512], pdb_file2[512];
  double pmat[MAXMAT][3][3], pt[MAXMAT][3];
  double rot[3][3], tr[3];
  double finalmat[4][4];
  double refAtom[3]; 
  double x,y,z;
  double xt, yt, zt;
  double test1;

  stype = 'U';
  zsym = 1;

  if((argc != 3 )&& (argc != 2)){
    printf("Usage: %s  <coordfile> <optional:matrixfile>\n\n\n", argv[0]);
    return 0; 
  }

  printf("\n\nProgram findframe\n\n");

//matrix file handling
  if (argc == 3){
    strcpy(pdb_file2,  argv[2]);
    file = fopen(pdb_file2, "r");
  if (file == NULL) {
    printf("ERROR: matrix input file does not exist\n\n"); }
    else { 
    printf("Matrices will be read from 2nd input file\n\n"); 
       nfiles=2;}
    fclose(file); 
     }

//coordinate file handling

  filetype=is_cif(argv[1]);

  if (!filetype){
    printf("Coordinates input is a pdb file\n");
    strcpy(pdb_file,  argv[1]);
  }else{
    printf("Coordinates input is a cif file\n");

    sprintf(pdb_file,"%s.PDB",argv[1]);    
    cif2pdb(argv[1], pdb_file);
  }

  file = fopen(pdb_file, "r");
  if (file == NULL) {
    printf("ERROR: coordinate input file does not exist\n\n");
    return 0;}


  //read in CA or P coords, find center of mass 
  //PDB format 
  i = 0;
  refAtom[0] = 0; refAtom[1] = 0; refAtom[2] = 0;
  while(fgets(buffer, 200, file)){
    if(strncmp(buffer, "ATOM", 4)) continue; 
    t = buffer + 13;
    if(strncmp(t, "CA",2) && strncmp(t,"P ",2)) continue; 
    t = buffer + 30; x = atof(t);
    t = buffer + 38; y = atof(t);
    t = buffer + 46; z = atof(t);
    refAtom[0] += x;
    refAtom[1] += y;
    refAtom[2] += z;
    i++;
  }
  fclose(file);

   if (i==0){ printf("ERROR: no CA or P coordinates read, ending program\n\n"); return 0;}
   else     { printf ("%i CA and/or P coordinates read\n", i);}

  refAtom[0] = refAtom[0]/i;
  refAtom[1] = refAtom[1]/i;
  refAtom[2] = refAtom[2]/i;

    printf ("\ncoordinate center of mass: %f %f %f\n", refAtom[0], refAtom[1], refAtom[2]);

  //Load matrices 
  //try BIOMT first 
if (nfiles==2) { file = fopen(pdb_file2, "r"); } else {file = fopen(pdb_file, "r"); }
  i = 0;
  while(fgets(buffer, 200, file)){
    if (i > MAXMAT*3){
     printf("\nERROR: too many  matrices\n\n");
     return 0;}
    t = buffer + 13;
    if (strncmp(t,"BIOMT",5)!=0)  continue; 
    t = buffer + 24; pmat[i/3][i%3][0] = atof(t);
    t = buffer + 34; pmat[i/3][i%3][1] = atof(t);
    t = buffer + 44; pmat[i/3][i%3][2] = atof(t);
    t = buffer + 58; pt[i/3][i%3] = atof(t);
    i++; }
  fclose(file);
  nmat = i/3; 


  //if # BIOMT < 2, try MTRIX
  if (nmat<2) {
  i = 0;
if (nfiles==2) { file = fopen(pdb_file2, "r"); } else {file = fopen(pdb_file, "r"); }
  while(fgets(buffer, 200, file)){
    t = buffer + 13;
    if(strncmp(buffer,"MTRIX",5)!=0)  continue; 
    t = buffer + 11; pmat[i/3][i%3][0] = atof(t);
    t = buffer + 21; pmat[i/3][i%3][1] = atof(t);
    t = buffer + 31; pmat[i/3][i%3][2] = atof(t);
    t = buffer + 45; pt[i/3][i%3] = atof(t);
    i++; }
  fclose(file);
  nmat = i/3; 
  }

  printf("\n%i matrices read\n",nmat);
   if (nmat < 2){
     printf("\nERROR: must provide at least 2 matrices\n\n");
     return 0;}


  //test for identity element , element limits
  iident=0;
  int ifix=0;
     for (i=0; i<nmat; i++){

         for (row=0;row<3;row++){ 
           for (column=0; column<3;column++) {
               if (pmat[i][row][column] < -1.0) {
                  printf("WARNING: resetting matrix %i element [%i][%i] from %f to -1.0\n", 
                           i+1, row, column, pmat[i][row][column]);
                           pmat[i][row][column] = -1.0; ifix++;}
               if (pmat[i][row][column] >  1.0) {
                  printf("WARNING: resetting matrix %i element [%i][%i] from %f to  1.0\n", 
                           i+1, row, column, pmat[i][row][column]);
                           pmat[i][row][column] =  1.0; ifix++;} }}

         if (pmat[i][0][0]==1 && pmat[i][0][1]==0 && pmat[i][0][2]==0
          && pmat[i][1][0]==0 && pmat[i][1][1]==1 && pmat[i][1][2]==0
          && pmat[i][2][0]==0 && pmat[i][2][1]==0 && pmat[i][2][2]==1
          && pt[i][0]==0 && pt[i][1]==0 && pt[i][2]==0) iident++;}

   if (iident!=1){
     printf("\nERROR: input matrix set must include one identity matrix, %i detected\n\n",iident);
     return 0;}
     else {printf ("\n one identity matrix read in \n");}

   if (ifix>0) printf ("\nWARNING: %i matrix elements outside of -1 to 1 range were reset to limit\n", ifix);

  //Get the rotation matrix(rot), translation vector (tr), symmetry (stype, zsym)
         idn=findPTV(nmat, pt, pmat, refAtom, rot, tr, &stype, &zsym);

  if (idn == 0) {
     printf("\nWARNING: NO FRAME TRANSFORMATION MATRIX CALCULATED\n");
     return 0;}

  
  //display initial result
          printf ("\nFINDFRAME MAT START\n");
          for (row=0; row<3;row++){
               printf("rot: %14.8f%14.8f%14.8f trans: %20.5f\n",
               checkZero(rot[row][0],PR8), 
               checkZero(rot[row][1],PR8), 
               checkZero(rot[row][2],PR8),
               checkZero(tr[row],PR5)); }

       //rotation matrix sanity check
       test1 = rot[0][0]*rot[0][0] + rot[0][1]*rot[0][1] + rot[0][2]*rot[0][2];
       test1 = test1 + rot[1][0]*rot[1][0] + rot[1][1]*rot[1][1] + rot[1][2]*rot[1][2];
       test1 = test1 + rot[2][0]*rot[2][0] + rot[2][1]*rot[2][1] + rot[2][2]*rot[2][2];
       test1 = fabs(test1 - 3.0);
       if (test1 > 0.05){
     printf("\nERROR: PROBLEM WITH MATRIX\n");
     return 0;}

  //refine initial result 
         idn=refinePDB2PT(rot,tr,pmat,pt,nmat,refAtom, &stype, &zsym);

  //display final result
          printf ("\nFINDFRAME MAT FINAL\n");
          for (row=0; row<3;row++){
               printf("rot: %14.8f%14.8f%14.8f trans: %20.5f\n",
               checkZero(rot[row][0],PR8), 
               checkZero(rot[row][1],PR8), 
               checkZero(rot[row][2],PR8),
               checkZero(tr[row],PR5)); }


  //print out final matrix to file 2ptmat

  printf ("\n\n writing assembly cif records to *findframe.cif*\n");

  outmatfile = fopen("findframe.cif","w");
  fprintf(outmatfile,"#\n");
  fprintf(outmatfile,"_pdbx_point_symmetry.entry_id                    TEMP\n");
  fprintf(outmatfile,"_pdbx_point_symmetry.Schoenflies_symbol          %c\n", stype);
  if (stype == 'C' || stype == 'D') {
     fprintf(outmatfile,"_pdbx_point_symmetry.circular_symmetry        %i\n\n",zsym);}
  fprintf(outmatfile,"#\n");
  fprintf(outmatfile,"loop_\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.id\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.type \n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.name \n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[1][1]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[1][2]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[1][3]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.vector[1]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[2][1]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[2][2]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[2][3] \n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.vector[2]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[3][1]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[3][2] \n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.matrix[3][3]\n");
  fprintf(outmatfile,"_pdbx_struct_oper_list.vector[3]\n");
  

  //WRITE OUT pdb2pt mat to cif
           fprintf(outmatfile,"   P   \"transform to point frame\"   . \n" );
           for (row=0;row<3;row++){
           fprintf(outmatfile,"    %14.8f%14.8f%14.8f%20.5f\n",
              checkZero(rot[row][0],PR8),
              checkZero(rot[row][1],PR8),
              checkZero(rot[row][2],PR8), 
              checkZero(tr[row],PR5)); }

  fclose(outmatfile);
  
    
  return 1;
}
