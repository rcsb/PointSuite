// pointsuite 0.5.8 
// author Cathy Lawson
// RCSB-PDB
// JUNE 2007

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define VERBOSE       0  //more output if not zero
#define MAXSTRINGS 1000  //number of substrings to be interpreted
#define MAXMAT     2000  //items in final matrix multiplication list (1680 needed for string in entry 1M4X) 
#define MAXDELIM     50  //maximum number of any one delimiter (),- in the input string
#define MAXCHAR      80  //maximum size of input string
#define MAXCHAIN      5  //needed for readcif.c but not used
#include "readcif.c"

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


//return positions of char chfind in string str
int strPos (char str[], int chfind, int ipos[]) {
  char * pch;
  int numpos = 0;
  pch=strchr(str,chfind);
  while (pch!=NULL) {
    ipos[numpos] = pch-str+1;
    numpos++;
    pch=strchr(pch+1,chfind);
  }
  return numpos;
}

//parse dash delimited string
int findRange (char str[], char strlist[2][MAXCHAR]){
  int i;
  int ipos[2];
  i=strPos(str,'-',ipos);
  if (i!=1) return 1;
  strncpy(strlist[0],str,ipos[0]-1);
  strlist[0][ipos[0]-1]='\0';
  //printf ("%s\n", strlist[0]);
  strncpy(strlist[1],str+ipos[0],strlen(str)-ipos[0]);
  strlist[1][strlen(str)-ipos[0]]='\0';
  //printf ("%s\n", strlist[1]);
  return 0;
}


//main parsing subroutine, called recursively
//substrings are held in strarray, next available position is stored in istrend
//final list is built in itearray, next available position is stored in iend
int parseString(char str[],char strarray[MAXSTRINGS][MAXCHAR],int *istrend,  
                           int itearray[MAXMAT],int *iend,
                           char idc[MAXMAT][MAXCHAR], 
                           double matall[MAXMAT][4][4], int *maxmat, int *mxmend) {
  char strlims[2][MAXCHAR];
  char strlist[MAXDELIM][MAXCHAR];
  char tmpstr[MAXCHAR];
  int iposl[MAXDELIM],ilevell[MAXDELIM],itopl[MAXDELIM];
  int iposr[MAXDELIM],ilevelr[MAXDELIM],itopr[MAXDELIM];
  int iposc[MAXDELIM],ilevelc[MAXDELIM],itopc[MAXDELIM];
  int iposd[MAXDELIM],ileveld[MAXDELIM],itopd[MAXDELIM];
  int inuml, inumr, inumc, inumd, itopnum, i, j, k;
  int islist,isrange, ilocstart, ilocend;
  int istarta, istartb, ienda, iendb;
  int ierror;

//find delimiter positions
 inuml = strPos(str,'(',iposl);
 inumr = strPos(str,')',iposr);
 inumc = strPos(str,',',iposc);
 inumd = strPos(str,'-',iposd);

//string must have balanced parentheses
if (inuml != inumr) return 1;
 
//scan through string to find delimiting characters at level zero
 i=0;
 for (j=0;j<strlen(str);j++) {
 for (k=0;k<inumd;k++) {if (j == (iposd[k]- 1)) {ileveld[k]=i;}}
 for (k=0;k<inumc;k++) {if (j == (iposc[k]- 1)) {ilevelc[k]=i;}}
 for (k=0;k<inumr;k++) {if (j == (iposr[k]- 1)) {i = i - 1;ilevelr[k]=i;}}
 for (k=0;k<inuml;k++) {if (j == (iposl[k]- 1)) {i = i + 1;ilevell[k]=i-1;}}
  }

//if commas level 0 assume that string is a list
itopc[0]=0;
islist=1;
for (j=0;j<inumc;j++) {
   if (ilevelc[j] == 0) {itopc[islist]=iposc[j]; islist++; }
   }
itopc[islist]=strlen(str)+1;


if (islist>1) { 
                ilocstart=*istrend; 
                i = *istrend;
                for (j=0;j<islist;j++) {
                strncpy(strarray[i],str+(itopc[j]),(itopc[j+1]-itopc[j]-1));
                strarray[i][itopc[j+1]-itopc[j]-1] = '\0';
                if (VERBOSE) printf ("string list item strarray[%i] :  %s \n", i, strarray[i]);
                i++;
                }
                ilocend = i;
                *istrend = i;
                for (i=ilocstart;i<ilocend;i++) {
                if (VERBOSE) printf ("call parse: %s\n", strarray[i]);
                ierror= parseString(strarray[i],strarray,istrend,itearray,iend,idc,matall,maxmat,mxmend); }
                if(ierror) return ierror;
                  return 0;}

//if parentheses at zero level process as multiplication expression
i = 0;
for (j=0;j<inuml;j++) { if (ilevell[j] == 0) {itopl[i]=iposl[j]; i++; } }
itopnum = i;

i = 0;
for (j=0;j<inumr;j++) { if (ilevelr[j] == 0) {itopr[i]=iposr[j]; i++; } }
if (itopnum != i) return 1;

 if (itopnum>0) { 
                ilocstart=*istrend;
                i=*istrend;
                for (j=0;j<itopnum;j++) {
                strncpy(strarray[i],str+(itopl[j]),(itopr[j]-itopl[j]-1));
                strarray[i][itopr[j]-itopl[j]-1] = '\0';
                if (VERBOSE) printf ("string expression strarray[%i] :  %s \n", i, strarray[i]);
                i++;
                   }
                ilocend = i;
                *istrend = i;
                 istarta = *iend;
                 if (VERBOSE) printf ("call parse: %s\n", strarray[ilocend-1]);
                 ierror=parseString(strarray[ilocend-1],strarray,istrend,itearray,iend,idc,matall,maxmat,mxmend);
                 if(ierror) return ierror;
                 ienda = *iend;

                for (i=ilocend-2;i>=ilocstart;i=i-1) {
                 istartb = *iend;
                 if (VERBOSE) printf ("call parse: %s\n", strarray[i]);
                 ierror=parseString(strarray[i],strarray,istrend,itearray,iend,idc,matall,maxmat,mxmend);
                 if(ierror) return ierror;
                 iendb = *iend;

                 //make multiplication strings 
                for (k=istartb;k<iendb;k++){
                  for (j=istarta;j<ienda;j++){
                    strncpy(idc[*mxmend],idc[itearray[k]],strlen(idc[itearray[k]])); 
                    idc[*mxmend][strlen(idc[itearray[k]])]='*';
                    idc[*mxmend][strlen(idc[itearray[k]])+1]='\0';
                    strcat(idc[*mxmend],idc[itearray[j]]);
                    itearray[*iend]=*mxmend;
                    multMat(matall[itearray[k]], matall[itearray[j]], matall[*mxmend] );
                    if(VERBOSE)printf("item assigned to itearray[%i]: %i, %s\n",
                                    *iend,itearray[*iend],idc[itearray[*iend]]);
                    *mxmend = *mxmend + 1;
                    *iend = *iend + 1; }}

                //overwrite the A and B array ranges, reset next avail position
                for (j=0;j<((ienda-istarta)*(iendb-istartb));j++){
                    itearray[istarta+j]=itearray[iendb+j];
                    //strncpy(itearray[istarta+j],itearray[iendb+j],strlen(itearray[iendb+j]));
                    if(VERBOSE)printf("item assigned to itearray[%i]: %i, %s\n",
                                    istarta+j,itearray[istarta+j],idc[itearray[istarta+j]]);
                            }
                    ienda = istarta + (ienda-istarta)*(iendb-istartb);
                    *iend = ienda;

                                                     }
                return 0;
              }


//check for dash delimits
if (inumd > 1) return 2;
if (inumd == 1) {
         findRange(str,strlims);
         istarta = atoi(strlims[0]);
         ienda   = atoi(strlims[1]);
         if (ienda<istarta) return 3;
         if (ienda==0) return 3;
         for (i=istarta;i<=ienda;i++) {
            k=0;
            itoa(i, tmpstr); 
            for (j=0;j<*maxmat;j++){if (!strncmp(tmpstr,idc[j],MAXCHAR)) break;}
            if (j==*maxmat) {printf ("error(1), no transformation labelled %s\n", tmpstr); return 4;}
            itearray[*iend]=j;
            if(VERBOSE) printf("item assigned to itearray[%i]: %i, %s\n", *iend,itearray[*iend],idc[itearray[*iend]]);
            *iend = *iend +1;
                 }
         return 0;
             }

// string is assigned as an "item" in the final list
for (i=0;i<*maxmat;i++){if (!strncmp(str,idc[i],MAXCHAR)) break;}
if (i==*maxmat) {printf ("error(2), no transformation labelled %s\n", str); return 4;}
itearray[*iend]=i;
if(VERBOSE) printf("item assigned to itearray[%i]: %i, %s\n", *iend, itearray[*iend], idc[itearray[*iend]]);
*iend= *iend+1;
//printf ("iend value: %i\n", *iend);

return 0;
}

int main(int argc, char * argv[]){
  int i, j, istrend, iend, ierror;
  char idc[MAXMAT][MAXCHAR];
  double matall[MAXMAT][4][4];
  FILE *outmatfile;
  char *matfile;
  int maxmat,mxmend;
  char strarray[MAXSTRINGS][MAXCHAR];
  int itearray[MAXMAT];
  char str[MAXCHAR];
  istrend = 0;
  iend = 0;
  ierror = 0;

//read in matrices from assembly cif
matfile = argv[1];
ierror = readAssemblyCif(matfile, idc, matall, &maxmat);
if (ierror==1) {printf("\nError: Problem reading cif or cif matrices\n"); return 1;}
printf ("number of matrices read: %i\n", maxmat);
mxmend=maxmat;
/*
//check matrices
for (j=0;j<maxmat;j++){
 printf("%i, %s, \n", j ,idc[j] );
 printf("%f,%f,%f,%f\n", matall[j][0][0],matall[j][0][1],matall[j][0][2],matall[j][0][3]);
 printf("%f,%f,%f,%f\n", matall[j][1][0],matall[j][1][1],matall[j][1][2],matall[j][1][3]);
 printf("%f,%f,%f,%f\n", matall[j][2][0],matall[j][2][1],matall[j][2][2],matall[j][2][3]);
 printf("%f,%f,%f,%f\n\n", matall[j][3][0],matall[j][3][1],matall[j][3][2],matall[j][3][3]); }
*/

//read in string
i = strlen(argv[2]);
if (i>MAXCHAR) {printf("\nError: Multiplication String Too Long\n"); return 1;}
strncpy(str,argv[2],i);
str[i]='\0';
printf ("Matrix multiplication expression to be parsed: %s\n",str);

  if (VERBOSE) printf ("call parse: %s\n",str);
  ierror=parseString (str, strarray, &istrend, itearray, &iend, idc, matall, &maxmat, &mxmend);

  if (ierror==1) {printf("\nError: Unbalanced parentheses in matrix multiplication string\n"); return 1;}
  if (ierror==2) {printf("\nError: More than one dash in matrix multiplication string item range\n"); return 1;}
  if (ierror==3) {printf("\nError: Matrix multiplication item range limits must be low integer-high integer\n"); return 1;}

  printf("\nresult:\n\n");
  for (i=0;i<iend;i++) printf("%s\n", idc[itearray[i]]);




  outmatfile = fopen("mult.cif","w");
  printf ("Writing %i matrices to *mult.cif*\n", iend );

  //WRITE matrix tokens to cif
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


for (j=0;j<iend;j++){
 fprintf(outmatfile, "%i \"general operation\"  %s  \n", j+1 ,idc[itearray[j]] );
 for (i=0;i<3;i++) { fprintf(outmatfile, "%14.8f %14.8f %14.8f %14.5f\n", 
                         matall[itearray[j]][i][0],matall[itearray[j]][i][1],
                         matall[itearray[j]][i][2],matall[itearray[j]][i][3]);}
                           }
   fclose(outmatfile);

  outmatfile = fopen("mult.biomt","w");
  printf ("Writing %i matrices to *mult.biomt*\n", iend );

    for(i=0;i<iend;i++){
      for (j=0;j<3;j++) {
           fprintf(outmatfile,"REMARK 350   BIOMT%1i %4i%10.6f%10.6f%10.6f%15.5f\n",
           j+1, i+1, matall[itearray[i]][j][0], 
                     matall[itearray[i]][j][1], 
                     matall[itearray[i]][j][2], 
                     matall[itearray[i]][j][3]);}
                }
   fclose(outmatfile);

  return 0;
}
