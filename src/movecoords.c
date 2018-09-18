//program reads in matrix file and pdb, outputs transformed coords
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXMAT 50
#include "mat44lib.c"

int main(int argc, char * argv[]){
  int i=0,j,k;
  FILE * file, *file2;
  char  buffer[200];
  char  linestart[32];
  char  lineend[30];
  char * t;
  char * pdb_file;
  char * mat_file;
  double p2i[4][4];
  double coord[4],xyz[4];


  pdb_file = argv[1];
  mat_file = argv[2];

printf("\nprogram movecoords\n");
if (argc != 3) {printf ("use: movecoords pdbfile matfile\n"); return 0;}


    for (j=0;j<4;j++)
       for (k=0;k<4;k++){if (j==k)p2i[j][k]=1.0; else p2i[j][k]=0.0;}

  //read in the transformation matrix 
                                                                                                                             
  file = fopen(mat_file, "r");
  if (file == NULL) { printf("matrix file does not exist\nidentity mat will be used\n"); }
  else {
  i = 0;
    for (j=0;j<4;j++)
       for (k=0;k<4;k++){ if (!feof(file)) fscanf (file,"%lf",&p2i[j][k]);  i++;}
  fclose(file);  
  printf("%i elements read from %s\n", i, mat_file); }
  printf("the following matrix will be applied to xyz coordinates:\n");
  for (i=0;i<4;i++) printf("%12.6f %12.6f %12.6f %12.4f\n", p2i[i][0], p2i[i][1], p2i[i][2], p2i[i][3]);

  file = fopen(pdb_file, "r");
  file2 = fopen("new.pdb", "w");
  if (file ==NULL){printf ("cannot open pdb file\n\n"); return 0;}

  i = 0;
  while(fgets(buffer, 200, file)){
    t = buffer + 0;
    if(!strncmp(t,"TER ",4)) fprintf(file2, "TER\n");
    if(!strncmp(t,"END ",4)) fprintf(file2, "END\n");
    if(strncmp(t,"ATOM ",5)!=0 && strncmp(t,"HETATM",6)!=0)  continue; 

//123456789012345678901234567890123456781234567812345678123456789012345678901234
//ATOM     23  C   ASP A  15      -3.815  49.683  24.980  1.00110.04           C
//30 char start, 3f8.3,26char end

    strncpy(linestart,t,30);
    linestart[30]='\0';
    t = buffer+30; coord[0] = atof(t);
    t = buffer+38; coord[1] = atof(t);
    t = buffer+46; coord[2] = atof(t);
    t = buffer+54; strncpy(lineend,t,29);
    lineend[29]='\0';
    coord[3]=1;
    applyMat(p2i,coord,xyz);
    fprintf(file2,"%30s%8.3f%8.3f%8.3f%-12.12s\n",linestart,xyz[0],xyz[1],xyz[2],lineend);
    i++;
  }
  fclose(file);
  fclose(file2);

}
