// pointsuite 0.5.8 
// author Cathy Lawson
// RCSB-PDB
// JUNE 2007

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void itoa(int i, char *str)
{ sprintf(str, "%d", i); }
                                                                                                                  
                                                                                               
void getString(FILE *file,char string[]) {
    fscanf(file,"%s",string);
    while ((!strncmp(string,"#",1)) && (strncmp(string+strlen(string)-1,"\n",1))) {
    fgets(string,MAXCHAR,file);   //flush comment line
    if (feof(file)) break;
    fscanf(file,"%s", string);}
  }
                                                                                               
void completeString(char string[], FILE *file) {
    int stringend=0;
    char str2[MAXCHAR];
      if (!strncmp(string,"\'",1)) {
          while (strlen(string)==1 || strncmp(string+strlen(string)-1,"\'",1)){
              fscanf(file,"%s", str2);
              strncat(string," ",1);
              strncat(string,str2,(MAXCHAR-1)-strlen(string));}}
      if (!strncmp(string,"\"",1)) {
          while (strlen(string)==1 || strncmp(string+strlen(string)-1,"\"",1)){
              fscanf(file,"%s", str2);
              strncat(string," ",1);
              strncat(string,str2,(MAXCHAR-1)-strlen(string));}}
      if (!strncmp(string,";",1)) {
       //   while (strlen(string)==1 || strncmp(string+strlen(string)-1,";",1)){
          while (strlen(string)==1 || !stringend){
              fscanf(file,"%s", str2);
              if(strlen(string)<MAXCHAR-1) strncat(string,str2,(MAXCHAR-1)-strlen(string));
              if (!strncmp(str2,";",1))stringend=1; }
          string[0]='\"';
          string[strlen(string)-1]='\"';
              }
  }
                                                                                               

int readAssemblyCif(char matfile[80], char idc[MAXMAT][MAXCHAR], 
                    double matall[MAXMAT][4][4], int *maxmat){

  char str[100];
  char tokens[20][40];
  char type[MAXMAT][40];
  char tnam[MAXMAT][40];
  FILE *pFile;
  int crossind[20];
  int checkind[20];
  int i,j,k, maxtoken, iend;


//mmcif items
  strncpy(tokens[0],"_pdbx_struct_oper_list.id\0",40);
  strncpy(tokens[1],"_pdbx_struct_oper_list.type\0",40);
  strncpy(tokens[2],"_pdbx_struct_oper_list.matrix[1][1]\0",40);
  strncpy(tokens[3],"_pdbx_struct_oper_list.matrix[1][2]\0",40);
  strncpy(tokens[4],"_pdbx_struct_oper_list.matrix[1][3]\0",40);
  strncpy(tokens[5],"_pdbx_struct_oper_list.vector[1]\0",40);
  strncpy(tokens[6],"_pdbx_struct_oper_list.matrix[2][1]\0",40);
  strncpy(tokens[7],"_pdbx_struct_oper_list.matrix[2][2]\0",40);
  strncpy(tokens[8],"_pdbx_struct_oper_list.matrix[2][3]\0",40);
  strncpy(tokens[9],"_pdbx_struct_oper_list.vector[2]\0",40);
  strncpy(tokens[10],"_pdbx_struct_oper_list.matrix[3][1]\0",40);
  strncpy(tokens[11],"_pdbx_struct_oper_list.matrix[3][2]\0",40);
  strncpy(tokens[12],"_pdbx_struct_oper_list.matrix[3][3]\0",40);
  strncpy(tokens[13],"_pdbx_struct_oper_list.vector[3]\0",40);
  strncpy(tokens[14],"_pdbx_struct_oper_list.name\0",40);

for (i=0;i<20;i++)checkind[i]=0;
for (i=0;i<20;i++)crossind[i]=-1;

//loop format assumed (even for 1 matrix)
//read once through file, discover token order
  pFile = fopen (matfile,"r");
  if (pFile == NULL) {
    //printf("input matrix set file does not exist\n\n");
    return 1;}
  i = 0;
  while (!feof(pFile)) {
  fscanf (pFile, "%100s", str);
  //printf ("%s \n", str);
 if (strncmp(str,"_pdbx_struct_oper_list",22)==0) {
          crossind[i] = -1;
          for (j=0; j<15;j++){
          if (strncmp(str,tokens[j],35)==0) {crossind[i]=j; checkind[j]=1;}}
   i++;}}
  fclose (pFile);

maxtoken = i;
if (maxtoken <14) return 1;
//first 14 token items must be present, name is optional
for (i=0;i<14;i++) {if (checkind[i]==0) return 1;}

// 2nd readthrough, begin reading data after the tokens 
  pFile = fopen (matfile,"r");
  i = 0;
  while (!feof(pFile)) {
  fscanf (pFile, "%100s", str);
 if (strncmp(str,"_pdbx_struct_oper_list",22)==0) {
  // printf("count: %i, crossind: %i\n", i, crossind[i]);
   i++;
   if(i==maxtoken)break;}}  //beginning of data loop identified
 
  i=0; 
  k=1;
  iend=0;
  while (!feof(pFile)) { 
   for(j=0;j<maxtoken;j++){
      //fscanf(pFile, "%s", str);
      getString(pFile,str);
      if (strncmp(str,"_",1)==0) iend=1;
      if (crossind[j]==0) {strncpy(idc[i],str,strlen(str));
                           idc[i][strlen(str)]='\0';}
      if (crossind[j]==1) { completeString(str,pFile);
                           strncpy(type[i],str,strlen(str)); 
                           type[i][strlen(str)]='\0';
                            }
      if (crossind[j]==14) { completeString(str,pFile);
                           strncpy(type[i],str,strlen(str)); 
                           tnam[i][strlen(str)]='\0';
                            }
      if (crossind[j]==2) matall[i][0][0]=atof(str);
      if (crossind[j]==3) matall[i][0][1]=atof(str);
      if (crossind[j]==4) matall[i][0][2]=atof(str);
      if (crossind[j]==5) matall[i][0][3]=atof(str);
      if (crossind[j]==6) matall[i][1][0]=atof(str);
      if (crossind[j]==7) matall[i][1][1]=atof(str);
      if (crossind[j]==8) matall[i][1][2]=atof(str);
      if (crossind[j]==9) matall[i][1][3]=atof(str);
      if (crossind[j]==10) matall[i][2][0]=atof(str);
      if (crossind[j]==11) matall[i][2][1]=atof(str);
      if (crossind[j]==12) matall[i][2][2]=atof(str);
      if (crossind[j]==13) matall[i][2][3]=atof(str);}
if (feof(pFile)) break;
if (iend==1) break;
i++;
}
  fclose (pFile);

//*maxmat = i-1;
*maxmat = i;
//add 4th line to make 4x4 matrices
for (j=0;j<*maxmat;j++){
 matall[j][3][0]=0; matall[j][3][1]=0; matall[j][3][2]=0; matall[j][3][3]=1;
 }
  return 0;
}

int readOneCif(char matfile[80], char token[80], char retstring[80]) {
  FILE *pFile;
  char str[100];
  
 strncpy(retstring,"\0",80); 
  pFile = fopen (matfile,"r");
  if (pFile == NULL) {return 1;}

while (!feof(pFile)) {
  fscanf (pFile, "%100s", str);
 if (strncmp(str,token,80)==0){ 
   // fscanf (pFile,"%80s", retstring); 
    getString(pFile,retstring);
    completeString(retstring,pFile); 
    fclose(pFile); return 0; }}
  fclose (pFile);
  return 1;
}

int readAsymCif(char ciffile[80], char idc[MAXCHAIN][MAXCHAR], 
                    char ide[MAXCHAIN][MAXCHAR], int *maxch){

  char str[100];
  char tokens[20][40];
  char type[MAXCHAIN][40];
  char tnam[MAXCHAIN][40];
  FILE *pFile;
  int crossind[20];
  int checkind[20];
  int i,j,k, maxtoken, iend, isingle;


//mmcif items
  strncpy(tokens[0],"_struct_asym.id\0",40);
  strncpy(tokens[1],"_struct_asym.entity_id\0",40);
  strncpy(tokens[2],"_struct_asym.pdbx_blank_PDB_chainid_flag\0",40);
  strncpy(tokens[3],"_struct_asym.pdbx_modified\0",40);
  strncpy(tokens[4],"_struct_asym.details\0",40);

for (i=0;i<20;i++)checkind[i]=0;
for (i=0;i<20;i++)crossind[i]=-1;
isingle = 0;

//read once through file, discover token order
  pFile = fopen (ciffile,"r");
  if (pFile == NULL) {
    //printf("input matrix set file does not exist\n\n");
    return 1;}
  i = 0;
  while (!feof(pFile)) {
  fscanf (pFile, "%100s", str);
  //printf ("%s \n", str);
 if (strncmp(str,"_struct_asym.",13)==0) {
          crossind[i] = -1;
          for (j=0; j<5;j++){
          if (strncmp(str,tokens[j],35)==0) {crossind[i]=j; 
          printf ("%s \n", str);
          checkind[j]=1;}}
   i++; }
   else if (i==1) isingle=1;}
  fclose (pFile);

maxtoken = i;
 printf ("# of tokens in _struct_asym: %i \n", maxtoken);
 if (isingle) printf ("non-loop structure detected\n");
if (maxtoken <2) return 1;
//1st 2  token items must be present
for (i=0;i<2;i++) {if (checkind[i]==0) return 1;}

// 2nd readthrough, begin reading data after the tokens 
  pFile = fopen (ciffile,"r");
  if (!isingle) {
  i = 0;
  while (!feof(pFile)) {
  fscanf (pFile, "%100s", str);
  //printf ("%s ", str);
 if (strncmp(str,"_struct_asym.",13)==0) {
   //printf("count: %i, crossind: %i\n", i, crossind[i]);
   i++;
   if(i==maxtoken)break;}}  //beginning of data loop identified
 
  i=0; 
  iend=0;
  while (!feof(pFile)) { 
   for(j=0;j<maxtoken;j++){
      //fscanf(pFile, "%s", str);
      getString(pFile, str);
      if (strncmp(str,"_",1)==0) iend=1;
      if (crossind[j]==0) {strncpy(idc[i],str,strlen(str));
                           idc[i][strlen(str)]='\0';}
      if (crossind[j]==1) {strncpy(ide[i],str,strlen(str)); 
                           ide[i][strlen(str)]='\0';}
     if (crossind[j]==4)  completeString(str,pFile);
      }
if (iend==1) break;
i++;
} }
else { 
//path for singletons
i = 0;
k = 0;
while (!feof(pFile)) {
fscanf (pFile, "%100s", str);
if (k==maxtoken) {i++; break;}
if (strncmp(str,"_struct_asym.",13)==0) {
      fscanf(pFile, "%s", str);
      if (crossind[k]==0) {strncpy(idc[i],str,strlen(str));
                           idc[i][strlen(str)]='\0';}
      if (crossind[k]==1) {strncpy(ide[i],str,strlen(str));
                           ide[i][strlen(str)]='\0';}
     if (crossind[k]==4)  completeString(str,pFile);
     k++; } }
}
  fclose (pFile);

*maxch = i;

  return 0;
}

/*
//simple program to read and report asym id's
int main(int argc, char * argv[]){
  int i, j, istrend, iend, ierror;
  char idc[MAXCHAIN][MAXCHAR];
  char ide[MAXCHAIN][MAXCHAR];
  char *ciffile;
  int maxch,mxmend;
  ierror = 0;

//read in matrices from assembly cif
ciffile = argv[1];
ierror = readAsymCif(ciffile, idc, ide, &maxch);
if (ierror==1) {printf("\nError: Problem reading cif \n"); return 1;}
printf ("number of chains read: %i\n", maxch);

  printf("\nresult:\n\n");
  for (i=0;i<maxch;i++) printf("chain asym id: %s   entity id:  %s\n", idc[i], ide[i]);

  printf("full chain list:   ");
  for (i=0;i<maxch-1;i++) printf("%s,", idc[i]);
  printf("%s\n", idc[maxch-1]);

  return 0;
}
*/
